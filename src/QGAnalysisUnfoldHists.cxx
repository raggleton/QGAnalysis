#include "UHH2/QGAnalysis/include/QGAnalysisUnfoldHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisPrinters.h"

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisUnfoldHists::QGAnalysisUnfoldHists(Context & ctx, const string & dirname,
                                             int useNJets, bool doGroomed,
                                             const string & selection,
                                             const string & reco_sel_handle_name, const string & gen_sel_handle_name,
                                             const string & reco_jetlambda_handle_name, const string & gen_jetlambda_handle_name):
  Hists(ctx, dirname),
  useNJets_(useNJets),
  doGroomed_(doGroomed),
  rand_(4357), // random number generator for splitting MC into 2 independent groups for unfold testing
  doMCsplit_(true),
  doJackknifeVariations_(false),
  useBinningValue_(false), // use the centrally determined bining value (e.g. dijet ave), otherwise use jet pT
  N_JACKKNIFE_VARIATIONS(25),
  N_PDF_VARIATIONS(100),
  eventCounter_(-1) // for jackknifing to ensure exclusive samples
  {

  is_mc_ = ctx.get("dataset_type") == "MC";
  doMCsplit_ = is_mc_;
  doJackknifeVariations_ = string2bool(ctx.get("JackknifeVariations", "false"));
  if (doJackknifeVariations_) {
    cout << "Doing jackknife variations in " << dirname << endl;
  } else {
    N_JACKKNIFE_VARIATIONS = 0;
  }
  doPDFvariations_ = (string2bool(ctx.get("PDFvariations", "false")) && is_mc_);
  if (doPDFvariations_) {
    cout << "Doing PDF variations in " << dirname << endl;
    doJackknifeVariations_ = false;
    cout << "Turning off jackknife variations" << endl;
  } else {
    N_PDF_VARIATIONS = 0;
  }

  if (useNJets_ < 0) useNJets_ = 99999; // Do them all
  else if (useNJets_ == 0) throw runtime_error("useNJets should be > 0, or < 0 to use all jets in the event");

  if (selection != "dijet" && selection != "zplusjets") {
    throw runtime_error("selection must be dijet or zplusjets");
  }

  gen_weight_handle = ctx.get_handle<double>("gen_weight");
  pt_binning_reco_handle = ctx.get_handle<double>("pt_binning_reco_value");
  pt_binning_gen_handle = ctx.get_handle<double>("pt_binning_gen_value");

  // setup TUnfold binning for reco & gen
  // ------------------------------------
  // select correct pt binning depending on whether dijet or z+jet
  std::vector<double> pt_bin_edges_reco = Binning::pt_bin_edges_reco;
  int nbins_pt_reco = Binning::nbins_pt_reco;

  std::vector<double> pt_bin_edges_reco_underflow = Binning::pt_bin_edges_reco_underflow;
  int nbins_pt_reco_underflow = Binning::nbins_pt_reco_underflow;

  std::vector<double> pt_bin_edges_gen = Binning::pt_bin_edges_gen;
  int nbins_pt_gen = Binning::nbins_pt_gen;

  std::vector<double> pt_bin_edges_gen_underflow = Binning::pt_bin_edges_gen_underflow;
  int nbins_pt_gen_underflow = Binning::nbins_pt_gen_underflow;

  if (selection == "zplusjets") {
    pt_bin_edges_reco = Binning::pt_bin_edges_zpj_reco;
    nbins_pt_reco = Binning::nbins_pt_zpj_reco;
    pt_bin_edges_reco_underflow = Binning::pt_bin_edges_zpj_reco_underflow;
    nbins_pt_reco_underflow = Binning::nbins_pt_zpj_reco_underflow;
    pt_bin_edges_gen = Binning::pt_bin_edges_zpj_gen;
    nbins_pt_gen = Binning::nbins_pt_zpj_gen;
    pt_bin_edges_gen_underflow = Binning::pt_bin_edges_zpj_gen_underflow;
    nbins_pt_gen_underflow = Binning::nbins_pt_zpj_gen_underflow;
  }

  // remove the bin which fall below reco cut
  auto pt_start = std::find(pt_bin_edges_reco_underflow.begin(), pt_bin_edges_reco_underflow.end(), Cuts::reco_jet_pt_min);
  std::vector<double> pt_bin_edges_reco_underflow_unfold(pt_start, pt_bin_edges_reco_underflow.end());
  cout << "Pt reco bin edges underflow: ";
  for (auto x : pt_bin_edges_reco_underflow_unfold) {
    cout << x << " ";
  }
  cout << endl;
  nbins_pt_reco_underflow = pt_bin_edges_reco_underflow_unfold.size()-1;

  // Remember to update
  std::vector<std::string> descriptions = {
      "!passGen",
      "genjetIndex < 0",
      "genjetIndex >= useNJets",
      "genjet 1 constit"
  };
  h_fake_counter_raw = book<TH1D>("fakes_counter_raw", "Fakes counter unweighted", descriptions.size()+1, -1, descriptions.size());
  h_fake_counter_weighted = book<TH1D>("fakes_counter", "Fakes counter using weights", descriptions.size()+1, -1, descriptions.size());
  h_fake_counter_charged_raw = book<TH1D>("fakes_counter_charged_raw", "Fakes counter for charged lambda unweighted", descriptions.size()+1, -1, descriptions.size());
  h_fake_counter_charged_weighted = book<TH1D>("fakes_countercharged_", "Fakes counter for charged lambda using weights", descriptions.size()+1, -1, descriptions.size());
  for(TAxis * ax : {h_fake_counter_raw->GetXaxis(), h_fake_counter_weighted->GetXaxis(), h_fake_counter_charged_raw->GetXaxis(), h_fake_counter_charged_weighted->GetXaxis()}){
    ax->SetBinLabel(1, "passReco");
    for(size_t i=0; i<descriptions.size(); ++i){
      ax->SetBinLabel(i+2, descriptions[i].c_str());
    }
  }

  // Common flags for TUnfold under/overflow
  bool pt_uf(false), pt_of(true); // don't need underflow, since gen cut is 15, and our binnig starts there
  bool var_uf(false), var_of(true); // don't need underflow bin, TUnfold has global uflow bin for when we have fakes/missing. oflow cos sometimes they go above 1?

  // pT only
  // -------------------------------------
  detector_tu_binning_pt = new TUnfoldBinning("detectorall");
  // add binning scheme for pT underflow regions
  // we handle it ourselves (instead of just having underflow bin in standard detector binning)
  // so that we can also use it to constrain the unfolding
  // (like simultaneously fitting to sideband regions)
  detector_distribution_underflow_pt = detector_tu_binning_pt->AddBinning("detector_underflow");
  detector_distribution_underflow_pt->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false); // reco jets guaranteed to be above lowest bin, so no uflow

  detector_distribution_pt = detector_tu_binning_pt->AddBinning("detector");
  detector_distribution_pt->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of); // no uflow

  generator_tu_binning_pt = new TUnfoldBinning("generator");
  generator_distribution_underflow_pt = generator_tu_binning_pt->AddBinning("signal_underflow");
  generator_distribution_underflow_pt->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false); // do need uflow for gen, as gen cut < first bin

  generator_distribution_pt = generator_tu_binning_pt->AddBinning("signal");
  generator_distribution_pt->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of); // no uflow for signal region (that's why we have our uflow region)

  // Response matrix
  // make tmp copies which we can then copy and use with book<>
  TH2 * h_tu_response_pt_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_pt, detector_tu_binning_pt, "tu_pt_GenReco");
  h_tu_response_pt = copy_book_th2d(h_tu_response_pt_tmp, "tu_pt_GenReco_all");
  h_tu_response_pt_split = copy_book_th2d(h_tu_response_pt_tmp, "tu_pt_GenReco_split");
  delete h_tu_response_pt_tmp;

  // detector histograms
  TH1 * h_tu_reco_pt_tmp = detector_tu_binning_pt->CreateHistogram("hist_pt_reco");
  h_tu_reco_pt = copy_book_th1d(h_tu_reco_pt_tmp, "hist_pt_reco_all");
  h_tu_reco_pt_split = copy_book_th1d(h_tu_reco_pt_tmp, "hist_pt_reco_split");
  // for fakes, detector binning
  h_tu_reco_pt_fake = copy_book_th1d(h_tu_reco_pt_tmp, "hist_pt_reco_fake_all");
  h_tu_reco_pt_fake_split = copy_book_th1d(h_tu_reco_pt_tmp, "hist_pt_reco_fake_split");
  delete h_tu_reco_pt_tmp;

  // truth histograms
  TH1 * h_tu_gen_pt_tmp = generator_tu_binning_pt->CreateHistogram("hist_pt_truth");
  h_tu_gen_pt = copy_book_th1d(h_tu_gen_pt_tmp, "hist_pt_truth_all");
  h_tu_gen_pt_split = copy_book_th1d(h_tu_gen_pt_tmp, "hist_pt_truth_split");
  delete h_tu_gen_pt_tmp;

  // PDF variations
  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_pt_PDF_variations.push_back(copy_book_th1d(h_tu_reco_pt, TString::Format("hist_pt_reco_all_PDF_%d", i).Data()));
      h_tu_gen_pt_PDF_variations.push_back(copy_book_th1d(h_tu_gen_pt, TString::Format("hist_pt_truth_all_PDF_%d", i).Data()));
      h_tu_response_pt_PDF_variations.push_back(copy_book_th2d(h_tu_response_pt, TString::Format("tu_pt_GenReco_all_PDF_%d", i).Data()));
    }
  }

  // LHA
  // -------------------------------------
  // don't give same name as child node, as FindNode() will use this parent node
  detector_tu_binning_LHA = new TUnfoldBinning("detectorall");
  // add binning scheme for pT underflow regions
  // we handle it ourselves (instead of just having underflow bin in standard detector binning)
  // so that we can also use it to constrain the unfolding
  // (like simultaneously fitting to sideband regions)
  detector_distribution_underflow_LHA = detector_tu_binning_LHA->AddBinning("detector_underflow");
  detector_distribution_underflow_LHA->AddAxis("LHA", Binning::nbins_var("LHA", doGroomed_, true), Binning::var_bin_edges("LHA", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_LHA->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_LHA = detector_tu_binning_LHA->AddBinning("detector");
  detector_distribution_LHA->AddAxis("LHA", Binning::nbins_var("LHA", doGroomed_, true), Binning::var_bin_edges("LHA", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_LHA->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_LHA = new TUnfoldBinning("generator");
  generator_distribution_underflow_LHA = generator_tu_binning_LHA->AddBinning("signal_underflow");
  generator_distribution_underflow_LHA->AddAxis("LHA", Binning::nbins_var("LHA", doGroomed_, false), Binning::var_bin_edges("LHA", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_LHA->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_LHA = generator_tu_binning_LHA->AddBinning("signal");
  generator_distribution_LHA->AddAxis("LHA", Binning::nbins_var("LHA", doGroomed_, false), Binning::var_bin_edges("LHA", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_LHA->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  bool doFakeHists = is_mc_;
  LHA_hist_filler.reset(new LambdaHistsFiller(Cuts::lha_args,
                                              detector_tu_binning_LHA,
                                              Cuts::lha_args,
                                              generator_tu_binning_LHA,
                                              "LHA"));
  LHA_hist_filler->setupRecoHists(ctx,
                                  dirname,
                                  doMCsplit_,
                                  doFakeHists,
                                  N_JACKKNIFE_VARIATIONS,
                                  N_PDF_VARIATIONS);
  if (is_mc_) {
    LHA_hist_filler->setupGenHists(ctx,
                                   dirname,
                                   doMCsplit_,
                                   N_JACKKNIFE_VARIATIONS,
                                   N_PDF_VARIATIONS);
    LHA_hist_filler->setupResponseHists(ctx,
                                        dirname,
                                        doMCsplit_,
                                        N_JACKKNIFE_VARIATIONS,
                                        N_PDF_VARIATIONS);
  }

  // Charged LHA
  // -------------------------------------
  detector_tu_binning_LHA_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_LHA_charged = detector_tu_binning_LHA_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_LHA_charged->AddAxis("LHA_charged", Binning::nbins_var("LHA_charged", doGroomed_, true), Binning::var_bin_edges("LHA_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_LHA_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_LHA_charged = detector_tu_binning_LHA_charged->AddBinning("detector");
  detector_distribution_LHA_charged->AddAxis("LHA_charged", Binning::nbins_var("LHA_charged", doGroomed_, true), Binning::var_bin_edges("LHA_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_LHA_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_LHA_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_LHA_charged = generator_tu_binning_LHA_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_LHA_charged->AddAxis("LHA_charged", Binning::nbins_var("LHA_charged", doGroomed_, false), Binning::var_bin_edges("LHA_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_LHA_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_LHA_charged = generator_tu_binning_LHA_charged->AddBinning("signal");
  generator_distribution_LHA_charged->AddAxis("LHA_charged", Binning::nbins_var("LHA_charged", doGroomed_, false), Binning::var_bin_edges("LHA_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_LHA_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  LHA_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::lha_args,
                                                      detector_tu_binning_LHA_charged,
                                                      Cuts::lha_args,
                                                      generator_tu_binning_LHA_charged,
                                                      "LHA_charged"));
  LHA_charged_hist_filler->setupRecoHists(ctx,
                                          dirname,
                                          doMCsplit_,
                                          doFakeHists,
                                          N_JACKKNIFE_VARIATIONS,
                                          N_PDF_VARIATIONS);
  if (is_mc_) {
    LHA_charged_hist_filler->setupGenHists(ctx,
                                           dirname,
                                           doMCsplit_,
                                           N_JACKKNIFE_VARIATIONS,
                                           N_PDF_VARIATIONS);
    LHA_charged_hist_filler->setupResponseHists(ctx,
                                                dirname,
                                                doMCsplit_,
                                                N_JACKKNIFE_VARIATIONS,
                                                N_PDF_VARIATIONS);
  }


  // puppi multiplicity
  // -------------------------------------
  detector_tu_binning_puppiMultiplicity = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_puppiMultiplicity = detector_tu_binning_puppiMultiplicity->AddBinning("detector_underflow");
  detector_distribution_underflow_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_var("puppiMultiplicity", doGroomed_, true), Binning::var_bin_edges("puppiMultiplicity", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_puppiMultiplicity->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_puppiMultiplicity = detector_tu_binning_puppiMultiplicity->AddBinning("detector");
  detector_distribution_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_var("puppiMultiplicity", doGroomed_, true), Binning::var_bin_edges("puppiMultiplicity", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_puppiMultiplicity->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_puppiMultiplicity = new TUnfoldBinning("generator");
  generator_distribution_underflow_puppiMultiplicity = generator_tu_binning_puppiMultiplicity->AddBinning("signal_underflow");
  generator_distribution_underflow_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_var("puppiMultiplicity", doGroomed_, false), Binning::var_bin_edges("puppiMultiplicity", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_puppiMultiplicity->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_puppiMultiplicity = generator_tu_binning_puppiMultiplicity->AddBinning("signal");
  generator_distribution_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_var("puppiMultiplicity", doGroomed_, false), Binning::var_bin_edges("puppiMultiplicity", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_puppiMultiplicity->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  puppiMultiplicity_hist_filler.reset(new LambdaHistsFiller(Cuts::mult_args,
                                                            detector_tu_binning_puppiMultiplicity,
                                                            Cuts::mult_args,
                                                            generator_tu_binning_puppiMultiplicity,
                                                            "puppiMultiplicity"));
  puppiMultiplicity_hist_filler->setupRecoHists(ctx,
                                                dirname,
                                                doMCsplit_,
                                                doFakeHists,
                                                N_JACKKNIFE_VARIATIONS,
                                                N_PDF_VARIATIONS);
  if (is_mc_) {
    puppiMultiplicity_hist_filler->setupGenHists(ctx,
                                                 dirname,
                                                 doMCsplit_,
                                                 N_JACKKNIFE_VARIATIONS,
                                                 N_PDF_VARIATIONS);
    puppiMultiplicity_hist_filler->setupResponseHists(ctx,
                                                      dirname,
                                                      doMCsplit_,
                                                      N_JACKKNIFE_VARIATIONS,
                                                      N_PDF_VARIATIONS);
  }

  // Charged PUPPI multiplicity
  // -------------------------------------
  detector_tu_binning_puppiMultiplicity_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_puppiMultiplicity_charged = detector_tu_binning_puppiMultiplicity_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_var("puppiMultiplicity_charged", doGroomed_, true), Binning::var_bin_edges("puppiMultiplicity_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_puppiMultiplicity_charged = detector_tu_binning_puppiMultiplicity_charged->AddBinning("detector");
  detector_distribution_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_var("puppiMultiplicity_charged", doGroomed_, true), Binning::var_bin_edges("puppiMultiplicity_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_puppiMultiplicity_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_puppiMultiplicity_charged = generator_tu_binning_puppiMultiplicity_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_var("puppiMultiplicity_charged", doGroomed_, false), Binning::var_bin_edges("puppiMultiplicity_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_puppiMultiplicity_charged = generator_tu_binning_puppiMultiplicity_charged->AddBinning("signal");
  generator_distribution_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_var("puppiMultiplicity_charged", doGroomed_, false), Binning::var_bin_edges("puppiMultiplicity_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  puppiMultiplicity_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::mult_args,
                                                                    detector_tu_binning_puppiMultiplicity_charged,
                                                                    Cuts::mult_args,
                                                                    generator_tu_binning_puppiMultiplicity_charged,
                                                                    "puppiMultiplicity_charged"));
  puppiMultiplicity_charged_hist_filler->setupRecoHists(ctx,
                                                        dirname,
                                                        doMCsplit_,
                                                        doFakeHists,
                                                        N_JACKKNIFE_VARIATIONS,
                                                        N_PDF_VARIATIONS);
  if (is_mc_) {
    puppiMultiplicity_charged_hist_filler->setupGenHists(ctx,
                                                         dirname,
                                                         doMCsplit_,
                                                         N_JACKKNIFE_VARIATIONS,
                                                         N_PDF_VARIATIONS);
    puppiMultiplicity_charged_hist_filler->setupResponseHists(ctx,
                                                              dirname,
                                                              doMCsplit_,
                                                              N_JACKKNIFE_VARIATIONS,
                                                              N_PDF_VARIATIONS);
  }

  // pTD
  // -------------------------------------
  detector_tu_binning_pTD = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_pTD = detector_tu_binning_pTD->AddBinning("detector_underflow");
  detector_distribution_underflow_pTD->AddAxis("pTD", Binning::nbins_var("pTD", doGroomed_, true), Binning::var_bin_edges("pTD", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_pTD->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_pTD = detector_tu_binning_pTD->AddBinning("detector");
  detector_distribution_pTD->AddAxis("pTD", Binning::nbins_var("pTD", doGroomed_, true), Binning::var_bin_edges("pTD", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_pTD->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_pTD = new TUnfoldBinning("generator");
  generator_distribution_underflow_pTD = generator_tu_binning_pTD->AddBinning("signal_underflow");
  generator_distribution_underflow_pTD->AddAxis("pTD", Binning::nbins_var("pTD", doGroomed_, false), Binning::var_bin_edges("pTD", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_pTD->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_pTD = generator_tu_binning_pTD->AddBinning("signal");
  generator_distribution_pTD->AddAxis("pTD", Binning::nbins_var("pTD", doGroomed_, false), Binning::var_bin_edges("pTD", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_pTD->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  pTD_hist_filler.reset(new LambdaHistsFiller(Cuts::pTD_args,
                                              detector_tu_binning_pTD,
                                              Cuts::pTD_args,
                                              generator_tu_binning_pTD,
                                              "pTD"));
  pTD_hist_filler->setupRecoHists(ctx,
                                  dirname,
                                  doMCsplit_,
                                  doFakeHists,
                                  N_JACKKNIFE_VARIATIONS,
                                  N_PDF_VARIATIONS);
  if (is_mc_) {
    pTD_hist_filler->setupGenHists(ctx,
                                   dirname,
                                   doMCsplit_,
                                   N_JACKKNIFE_VARIATIONS,
                                   N_PDF_VARIATIONS);
    pTD_hist_filler->setupResponseHists(ctx,
                                        dirname,
                                        doMCsplit_,
                                        N_JACKKNIFE_VARIATIONS,
                                        N_PDF_VARIATIONS);
  }

  // Charged pTD
  // -------------------------------------
  detector_tu_binning_pTD_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_pTD_charged = detector_tu_binning_pTD_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_pTD_charged->AddAxis("pTD_charged", Binning::nbins_var("pTD_charged", doGroomed_, true), Binning::var_bin_edges("pTD_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_pTD_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_pTD_charged = detector_tu_binning_pTD_charged->AddBinning("detector");
  detector_distribution_pTD_charged->AddAxis("pTD_charged", Binning::nbins_var("pTD_charged", doGroomed_, true), Binning::var_bin_edges("pTD_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_pTD_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_pTD_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_pTD_charged = generator_tu_binning_pTD_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_pTD_charged->AddAxis("pTD_charged", Binning::nbins_var("pTD_charged", doGroomed_, false), Binning::var_bin_edges("pTD_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_pTD_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_pTD_charged = generator_tu_binning_pTD_charged->AddBinning("signal");
  generator_distribution_pTD_charged->AddAxis("pTD_charged", Binning::nbins_var("pTD_charged", doGroomed_, false), Binning::var_bin_edges("pTD_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_pTD_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  pTD_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::pTD_args,
                                                      detector_tu_binning_pTD_charged,
                                                      Cuts::pTD_args,
                                                      generator_tu_binning_pTD_charged,
                                                      "pTD_charged"));
  pTD_charged_hist_filler->setupRecoHists(ctx,
                                          dirname,
                                          doMCsplit_,
                                          doFakeHists,
                                          N_JACKKNIFE_VARIATIONS,
                                          N_PDF_VARIATIONS);
  if (is_mc_) {
    pTD_charged_hist_filler->setupGenHists(ctx,
                                           dirname,
                                           doMCsplit_,
                                           N_JACKKNIFE_VARIATIONS,
                                           N_PDF_VARIATIONS);
    pTD_charged_hist_filler->setupResponseHists(ctx,
                                                dirname,
                                                doMCsplit_,
                                                N_JACKKNIFE_VARIATIONS,
                                                N_PDF_VARIATIONS);
  }

  // thrust
  // -------------------------------------
  detector_tu_binning_thrust = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_thrust = detector_tu_binning_thrust->AddBinning("detector_underflow");
  detector_distribution_underflow_thrust->AddAxis("thrust", Binning::nbins_var("thrust", doGroomed_, true), Binning::var_bin_edges("thrust", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_thrust->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_thrust = detector_tu_binning_thrust->AddBinning("detector");
  detector_distribution_thrust->AddAxis("thrust", Binning::nbins_var("thrust", doGroomed_, true), Binning::var_bin_edges("thrust", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_thrust->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_thrust = new TUnfoldBinning("generator");
  generator_distribution_underflow_thrust = generator_tu_binning_thrust->AddBinning("signal_underflow");
  generator_distribution_underflow_thrust->AddAxis("thrust", Binning::nbins_var("thrust", doGroomed_, false), Binning::var_bin_edges("thrust", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_thrust->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_thrust = generator_tu_binning_thrust->AddBinning("signal");
  generator_distribution_thrust->AddAxis("thrust", Binning::nbins_var("thrust", doGroomed_, false), Binning::var_bin_edges("thrust", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_thrust->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  thrust_hist_filler.reset(new LambdaHistsFiller(Cuts::thrust_args,
                                                 detector_tu_binning_thrust,
                                                 Cuts::thrust_args,
                                                 generator_tu_binning_thrust,
                                                 "thrust"));
  thrust_hist_filler->setupRecoHists(ctx,
                                     dirname,
                                     doMCsplit_,
                                     doFakeHists,
                                     N_JACKKNIFE_VARIATIONS,
                                     N_PDF_VARIATIONS);
  if (is_mc_) {
    thrust_hist_filler->setupGenHists(ctx,
                                      dirname,
                                      doMCsplit_,
                                      N_JACKKNIFE_VARIATIONS,
                                      N_PDF_VARIATIONS);
    thrust_hist_filler->setupResponseHists(ctx,
                                           dirname,
                                           doMCsplit_,
                                           N_JACKKNIFE_VARIATIONS,
                                           N_PDF_VARIATIONS);
  }

  // Charged thrust
  // -------------------------------------
  detector_tu_binning_thrust_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_thrust_charged = detector_tu_binning_thrust_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_thrust_charged->AddAxis("thrust_charged", Binning::nbins_var("thrust_charged", doGroomed_, true), Binning::var_bin_edges("thrust_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_thrust_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_thrust_charged = detector_tu_binning_thrust_charged->AddBinning("detector");
  detector_distribution_thrust_charged->AddAxis("thrust_charged", Binning::nbins_var("thrust_charged", doGroomed_, true), Binning::var_bin_edges("thrust_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_thrust_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_thrust_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_thrust_charged = generator_tu_binning_thrust_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_thrust_charged->AddAxis("thrust_charged", Binning::nbins_var("thrust_charged", doGroomed_, false), Binning::var_bin_edges("thrust_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_thrust_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_thrust_charged = generator_tu_binning_thrust_charged->AddBinning("signal");
  generator_distribution_thrust_charged->AddAxis("thrust_charged", Binning::nbins_var("thrust_charged", doGroomed_, false), Binning::var_bin_edges("thrust_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_thrust_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  thrust_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::thrust_args,
                                                         detector_tu_binning_thrust_charged,
                                                         Cuts::thrust_args,
                                                         generator_tu_binning_thrust_charged,
                                                         "thrust_charged"));
  thrust_charged_hist_filler->setupRecoHists(ctx,
                                             dirname,
                                             doMCsplit_,
                                             doFakeHists,
                                             N_JACKKNIFE_VARIATIONS,
                                             N_PDF_VARIATIONS);
  if (is_mc_) {
    thrust_charged_hist_filler->setupGenHists(ctx,
                                              dirname,
                                              doMCsplit_,
                                              N_JACKKNIFE_VARIATIONS,
                                              N_PDF_VARIATIONS);
    thrust_charged_hist_filler->setupResponseHists(ctx,
                                                   dirname,
                                                   doMCsplit_,
                                                   N_JACKKNIFE_VARIATIONS,
                                                   N_PDF_VARIATIONS);
  }

  // width
  // -------------------------------------
  detector_tu_binning_width = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_width = detector_tu_binning_width->AddBinning("detector_underflow");
  detector_distribution_underflow_width->AddAxis("width", Binning::nbins_var("width", doGroomed_, true), Binning::var_bin_edges("width", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_width->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_width = detector_tu_binning_width->AddBinning("detector");
  detector_distribution_width->AddAxis("width", Binning::nbins_var("width", doGroomed_, true), Binning::var_bin_edges("width", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_width->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_width = new TUnfoldBinning("generator");
  generator_distribution_underflow_width = generator_tu_binning_width->AddBinning("signal_underflow");
  generator_distribution_underflow_width->AddAxis("width", Binning::nbins_var("width", doGroomed_, false), Binning::var_bin_edges("width", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_width->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_width = generator_tu_binning_width->AddBinning("signal");
  generator_distribution_width->AddAxis("width", Binning::nbins_var("width", doGroomed_, false), Binning::var_bin_edges("width", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_width->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  width_hist_filler.reset(new LambdaHistsFiller(Cuts::width_args,
                                                detector_tu_binning_width,
                                                Cuts::width_args,
                                                generator_tu_binning_width,
                                                "width"));
  width_hist_filler->setupRecoHists(ctx,
                                    dirname,
                                    doMCsplit_,
                                    doFakeHists,
                                    N_JACKKNIFE_VARIATIONS,
                                    N_PDF_VARIATIONS);
  if (is_mc_) {
    width_hist_filler->setupGenHists(ctx,
                                     dirname,
                                     doMCsplit_,
                                     N_JACKKNIFE_VARIATIONS,
                                     N_PDF_VARIATIONS);
    width_hist_filler->setupResponseHists(ctx,
                                          dirname,
                                          doMCsplit_,
                                          N_JACKKNIFE_VARIATIONS,
                                          N_PDF_VARIATIONS);
  }

  // Charged width
  // -------------------------------------
  detector_tu_binning_width_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_width_charged = detector_tu_binning_width_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_width_charged->AddAxis("width_charged", Binning::nbins_var("width_charged", doGroomed_, true), Binning::var_bin_edges("width_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_underflow_width_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_width_charged = detector_tu_binning_width_charged->AddBinning("detector");
  detector_distribution_width_charged->AddAxis("width_charged", Binning::nbins_var("width_charged", doGroomed_, true), Binning::var_bin_edges("width_charged", doGroomed_, true).data(), var_uf, var_of);
  detector_distribution_width_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_width_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_width_charged = generator_tu_binning_width_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_width_charged->AddAxis("width_charged", Binning::nbins_var("width_charged", doGroomed_, false), Binning::var_bin_edges("width_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_underflow_width_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_width_charged = generator_tu_binning_width_charged->AddBinning("signal");
  generator_distribution_width_charged->AddAxis("width_charged", Binning::nbins_var("width_charged", doGroomed_, false), Binning::var_bin_edges("width_charged", doGroomed_, false).data(), var_uf, var_of);
  generator_distribution_width_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  width_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::width_args,
                                                        detector_tu_binning_width_charged,
                                                        Cuts::width_args,
                                                        generator_tu_binning_width_charged,
                                                        "width_charged"));
  width_charged_hist_filler->setupRecoHists(ctx,
                                            dirname,
                                            doMCsplit_,
                                            doFakeHists,
                                            N_JACKKNIFE_VARIATIONS,
                                            N_PDF_VARIATIONS);
  if (is_mc_) {
    width_charged_hist_filler->setupGenHists(ctx,
                                             dirname,
                                             doMCsplit_,
                                             N_JACKKNIFE_VARIATIONS,
                                             N_PDF_VARIATIONS);
    width_charged_hist_filler->setupResponseHists(ctx,
                                                  dirname,
                                                  doMCsplit_,
                                                  N_JACKKNIFE_VARIATIONS,
                                                  N_PDF_VARIATIONS);
  }


  // Setup other things
  // ---------------------------------------------------------------------------
  if (is_mc_) {
    genJetsLambda_handle = ctx.get_handle< std::vector<GenJetLambdaBundle> > (gen_jetlambda_handle_name);
    pass_gen_handle = ctx.get_handle<bool> (gen_sel_handle_name);
  }

  jetsLambda_handle = ctx.get_handle< std::vector<JetLambdaBundle> > (reco_jetlambda_handle_name);

  pass_reco_handle = ctx.get_handle<bool> (reco_sel_handle_name);

  // TODO: automate setting up of hists - let the LambdaHistsFiller own them?

  chargedPlusNeutralHistFillers = {
    LHA_hist_filler.get(), // use .get() to get raw ptr, since we don't want to own the object
    puppiMultiplicity_hist_filler.get(),
    pTD_hist_filler.get(),
    thrust_hist_filler.get(),
    width_hist_filler.get(),
  };

  chargedOnlyHistFillers = {
    LHA_charged_hist_filler.get(),
    puppiMultiplicity_charged_hist_filler.get(),
    pTD_charged_hist_filler.get(),
    thrust_charged_hist_filler.get(),
    width_charged_hist_filler.get(),
  };

  // combine into one vector of HistFillers
  allHistFillers = chargedPlusNeutralHistFillers;
  allHistFillers.insert(allHistFillers.end(), chargedOnlyHistFillers.begin(), chargedOnlyHistFillers.end());
}

void fill_th1_check(TH1 * h, int value, double weight){
  if (value <= 0) {
    cout << h->GetName() << " : " << value << " : "<< weight << endl;
    throw runtime_error("Filling with <= 0 not allowed");
  }
  h->Fill(value, weight);
}

void fill_th2_check(TH2 * h, int xvalue, int yvalue, double weight){
  if (xvalue <= 0 && yvalue <=0) {
    cout << h->GetName() << " : " << xvalue << " : " << yvalue << " : " << weight << endl;
    throw runtime_error("Filling with <= 0 on x and y not allowed");
  }

  h->Fill(xvalue, yvalue, weight);
}


void QGAnalysisUnfoldHists::fill(const Event & event){
  double weight = event.weight;

  // extract the separate gen & reco weight components, needed for TUnfold
  double gen_weight = event.get(gen_weight_handle);
  double reco_weight = weight / gen_weight;

  // Get selection flags
  bool passReco = event.get(pass_reco_handle);
  bool passGen(false);

  // Get (Gen)Jet-Lambda Bundles
  // ptr since then it can be null
  const std::vector<GenJetLambdaBundle> * genjetLambdas = nullptr;
  if (is_mc_) {
    genjetLambdas = &event.get(genJetsLambda_handle);
    passGen = event.get(pass_gen_handle);
  }

  const std::vector<JetLambdaBundle> & jetLambdas = event.get(jetsLambda_handle);

  // cout << "***" << event.event << endl;

  // Random event flag for MC only filling response matrix
  // This allows us to split the MC into 2 separate samples for testing
  // Make 80% go into response hist so good stats
  bool onlyFillResponse = (rand_.Rndm() > 0.2);

  eventCounter_++;

  // Fill reco jet 1D hists
  // ---------------------------------------------------------------------------
  // At this point, all jet filtering etc should have already been performed
  // The only cut we do is a per-variable cut on # constituents,
  // since each variable can have different constituent cut
  // (the cut itself is hidden away in LambdaHistsFiller::setupReco/Gen)
  // Fill ignoring if there is a genjet match or not
  if (passReco) {
    for (int i = 0; i < useNJets_; i++) {
      const Jet & thisjet = jetLambdas.at(i).jet;
      float jet_pt = useBinningValue_ ? event.get(pt_binning_reco_handle) : thisjet.pt();

      // Update the LambdaHistsFiller objects for this reco jet
      for (auto filler : chargedPlusNeutralHistFillers) {
        filler->setupReco(jet_pt, jetLambdas.at(i).getLambdaCalculator(false, doGroomed_), passReco);
      }
      for (auto filler : chargedOnlyHistFillers) {
        filler->setupReco(jet_pt, jetLambdas.at(i).getLambdaCalculator(true, doGroomed_), passReco);
      }

      // Do all the filling
      for (auto filler : allHistFillers) {
        filler->fillRecoTH1(filler->recoHist(), weight);
        filler->fillRecoTH1GenBinning(filler->recoHistGenBinning(), weight);

        if (!onlyFillResponse && doMCsplit_) {
          filler->fillRecoTH1(filler->recoSplitHist(), weight);
          filler->fillRecoTH1GenBinning(filler->recoSplitHistGenBinning(), weight);
        }

        if (doJackknifeVariations_) {
          // Fill jackknife reco hists
          // we use 100*(1-(1/N_JACKKNIFE_VARIATIONS))% of the sample to fill each gen hist matrix
          // so we fill every hist *except* for the one at index (eventCounter_ % N_JACKKNIFE_VARIATIONS)
          // Then afterwards in the analysis we use RMS * (N_JACKKNIFE_VARIATIONS / (N_JACKKNIFE_VARIATIONS-1))
          // as the uncertainty
          uint index = eventCounter_ % N_JACKKNIFE_VARIATIONS;
          for (uint jk_ind=0; jk_ind<N_JACKKNIFE_VARIATIONS; jk_ind++) {
            if (jk_ind == index) continue;
            filler->fillRecoTH1(filler->recoJackknifeVariations().at(index), weight);
          }
        }

        // Fill PDF variations by varying weight
        if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
          for (int i=0; i<N_PDF_VARIATIONS; i++) {
            double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations, then comes the nominal weight again
            double this_weight = weight * pdf_weight;
            filler->fillRecoTH1(filler->recoPDFVariations().at(i), this_weight);
          }
        }
      } // end loop over hist fillers

      // Manually do pt hist separately
      // -------------------------------------
      int recBinPt(0); // default to 0 for underflow
      bool isUnderflow = (jet_pt < Binning::pt_bin_edges_reco[0]);
      // here we get bin number based on if underflow or not
      if (isUnderflow) {
        recBinPt = detector_distribution_underflow_pt->GetGlobalBinNumber(jet_pt);
      } else {
        recBinPt = detector_distribution_pt->GetGlobalBinNumber(jet_pt);
      }
      fill_th1_check(h_tu_reco_pt, recBinPt, weight);

      if (!onlyFillResponse && doMCsplit_) {
        h_tu_reco_pt_split->Fill(recBinPt, weight);
      }

      // Fill PDF variations by varying weight
      if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
        for (int i=0; i<N_PDF_VARIATIONS; i++) {
          double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations, then comes the nominal weight again
          double this_weight = weight * pdf_weight;
          h_tu_reco_pt_PDF_variations.at(i)->Fill(recBinPt, this_weight);
        }
      }

      // Fill fakes hists
      // ----------------
      if (is_mc_) {
        // Fake = not pass Gen, or pass Gen but matching genjet not one of the
        // selection ones (i.e. doesn't count as a match)
        // Or the genjet fails #constits > 1
        if ((thisjet.genjet_index() >= 0) && (thisjet.genjet_index() < (int)genjetLambdas->size())) {
          // if there is a matching genjet, set up the LambdaHistsFiller with it so we can test the # constits
          for (auto filler : chargedPlusNeutralHistFillers) {
            filler->setupGen(genjetLambdas->at(thisjet.genjet_index()).jet.pt(),
                             genjetLambdas->at(thisjet.genjet_index()).getLambdaCalculator(false, doGroomed_),
                             passGen);
          }
          for (auto filler : chargedOnlyHistFillers) {
            filler->setupGen(genjetLambdas->at(thisjet.genjet_index()).jet.pt(),
                             genjetLambdas->at(thisjet.genjet_index()).getLambdaCalculator(true, doGroomed_),
                             passGen);
          }
        }

        for (auto filler: allHistFillers) {
          if (!passGen || (thisjet.genjet_index() < 0 || thisjet.genjet_index() >= useNJets_) || !filler->passGen()) {
            filler->fillRecoTH1(filler->recoHistFakes(), weight);
            filler->fillRecoTH1GenBinning(filler->recoHistFakesGenBinning(), weight);
            if (!onlyFillResponse && doMCsplit_) { filler->fillRecoTH1(filler->recoSplitHistFakes(), weight); }
          }
        }

        // fill fakes cutflow & pt fakes
        if (!passGen || (thisjet.genjet_index() < 0 || thisjet.genjet_index() >= useNJets_)) {
          if (!passGen) {
            int ind=0;
            h_fake_counter_raw->Fill(ind);
            h_fake_counter_weighted->Fill(ind, weight);
          }
          else if (thisjet.genjet_index() < 0) {
            int ind=1;
            h_fake_counter_raw->Fill(ind);
            h_fake_counter_weighted->Fill(ind, weight);
          }
          else if (thisjet.genjet_index() >= useNJets_) {
            int ind=2;
            h_fake_counter_raw->Fill(ind);
            h_fake_counter_weighted->Fill(ind, weight);
          }
          // else if (!genConstit) {
          //   int ind=3;
          //   h_fake_counter_raw->Fill(ind);
          //   h_fake_counter_weighted->Fill(ind, weight);
          // }

          fill_th1_check(h_tu_reco_pt_fake, recBinPt, weight);
          if (!onlyFillResponse && doMCsplit_) {
            h_tu_reco_pt_fake_split->Fill(recBinPt, weight);
          }
        } // end of test if failed gen

      } // end if mc / fakes
    } // end recojet loop
  } // end if passReco

  // Fill gen jet 1D hists & response matrices
  // ---------------------------------------------------------------------------
  if (is_mc_ && passGen) {
    for (int i = 0; i < useNJets_; i++) {
      const GenJet & thisjet = genjetLambdas->at(i).jet;
      float genjet_pt = useBinningValue_ ? event.get(pt_binning_gen_handle) : thisjet.pt();

      // Update the LambdaHistsFiller objects for this gen jet
      for (auto filler : chargedPlusNeutralHistFillers) {
        filler->setupGen(genjet_pt, genjetLambdas->at(i).getLambdaCalculator(false, doGroomed_), passGen);
      }
      for (auto filler : chargedOnlyHistFillers) {
        filler->setupGen(genjet_pt, genjetLambdas->at(i).getLambdaCalculator(true, doGroomed_), passGen);
      }

      // Do all the filling
      for (auto filler : allHistFillers) {
        filler->fillGenTH1(filler->genHist(), gen_weight);

        if (!onlyFillResponse && doMCsplit_) {
          filler->fillGenTH1(filler->genSplitHist(), gen_weight);
        }

        // Fill jackknife hists
        if (doJackknifeVariations_) {
          uint index = eventCounter_ % N_JACKKNIFE_VARIATIONS;
          for (uint jk_ind=0; jk_ind<N_JACKKNIFE_VARIATIONS; jk_ind++) {
            if (jk_ind == index) continue;
            filler->fillGenTH1(filler->genJackknifeVariations().at(index), gen_weight);
          }
        }

        // Fill PDF variations by varying weight
        if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
          for (int i=0; i<N_PDF_VARIATIONS; i++) {
            double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations, then comes the nominal weight again
            double this_weight = gen_weight * pdf_weight;
            filler->fillGenTH1(filler->genPDFVariations().at(i), this_weight);
          }
        }
      } // end loop over hist fillers

      // Fill 1D gen pt hists manually
      // -----------------------------------------------------------------------
      int genBinPt(0);
      bool isUnderflowGen = (genjet_pt < Binning::pt_bin_edges_gen[0]);
      if (isUnderflowGen) {
        genBinPt = generator_distribution_underflow_pt->GetGlobalBinNumber(genjet_pt);
      } else {
        genBinPt = generator_distribution_pt->GetGlobalBinNumber(genjet_pt);
      }

      h_tu_gen_pt->Fill(genBinPt, gen_weight);

      if (!onlyFillResponse && doMCsplit_) {
        h_tu_gen_pt_split->Fill(genBinPt, gen_weight);
      }

      // Fill PDF variations by varying weight
      if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
        for (int i=0; i<N_PDF_VARIATIONS; i++) {
          double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations, then comes the nominal weight again
          double this_weight = pdf_weight * gen_weight;
          h_tu_gen_pt_PDF_variations.at(i)->Fill(genBinPt, this_weight);
        }
      }

      // Now fill response matrices
      // -----------------------------------------------------------------------
      // Loop through genjets, since if there isn't a reco jet then it's a miss-reco (= bin 0),
      // whereas a recojet & no genjet is a fake, which we don't want in our migration matrix
      int recBinPt(0);
      int recoInd = -1;
      if (passReco) { // only valid match if reco selection also passed
        // First find matching recojet - annoyingly matches are held in the Jet class, not GenJet
        // default to -1 for underflow
        for (int j=0; j<useNJets_; j++) {
          if (jetLambdas.at(j).jet.genjet_index() == i) { recoInd = j; break; }
        }

        if (recoInd < 0) {
          // out of interest, try to figure out which is the match
          for (uint j=0; j<jetLambdas.size(); j++) {
            if (jetLambdas.at(j).jet.genjet_index() == i) {
              cout << "Match for genjet " << i << " actually found with reco jet " << j << endl;
              // printJets(*event.jets);
              // printGenJets(*event.genjets);
            }
          }
        }

        if (recoInd >= 0 && recoInd < useNJets_) {
          const Jet & thisjet = jetLambdas.at(recoInd).jet;
          float jet_pt = useBinningValue_ ? event.get(pt_binning_reco_handle) : thisjet.pt();

          for (auto filler : chargedPlusNeutralHistFillers) {
            // check if the reco jet already in this filler is the same one (saves time)
            // and if not, update it
            if (filler->recoJetPt() != jet_pt) {
              filler->setupReco(jet_pt, jetLambdas.at(recoInd).getLambdaCalculator(false, doGroomed_), passReco);
            }
          }
          for (auto filler : chargedOnlyHistFillers) {
            // check if the reco jet already in this filler is the same one (Saves time)
            // and if not, update it
            if (filler->recoJetPt() != jet_pt) {
              filler->setupReco(jet_pt, jetLambdas.at(recoInd).getLambdaCalculator(true, doGroomed_), passReco);
            }
          }

          bool isUnderflow = (jet_pt < Binning::pt_bin_edges_reco[0]);
          // here we get bin number based on if underflow or not
          if (isUnderflow) {
            recBinPt = detector_distribution_underflow_pt->GetGlobalBinNumber(jet_pt);
          } else {
            recBinPt = detector_distribution_pt->GetGlobalBinNumber(jet_pt);
          }
        } else {
          // mark as failed, since we didn't get a match
          for (auto filler : allHistFillers) {
            filler->setPassReco(false);
          }
        } // end if recoInd >= 0
      } else {
        // Setup all fillers to mark that they failed reco
        for (auto filler : allHistFillers) {
          filler->setPassReco(false);
        }
      }// end if passReco

      // Fill TUnfold 2D response maps
      // We fill it so long as we have a passGen, whether or not we have passReco
      // This ensures we account for the efficiency correctly
      // -----------------------------------------------------------------------
      fill_th2_check(h_tu_response_pt, genBinPt, recBinPt, weight);

      for (auto filler : allHistFillers) {
        filler->fillResponseTH2(filler->responseHist(), reco_weight, gen_weight);

        if (onlyFillResponse && doMCsplit_) filler->fillResponseTH2(filler->responseSplitHist(), reco_weight, gen_weight);

        if (doJackknifeVariations_) {
          uint index = eventCounter_ % N_JACKKNIFE_VARIATIONS;
          for (uint jk_ind=0; jk_ind<N_JACKKNIFE_VARIATIONS; jk_ind++) {
            if (jk_ind == index) continue;
            filler->fillResponseTH2(filler->responseJackknifeVariations().at(index), reco_weight, gen_weight);
          }
        }

        if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
          for (int i=0; i<N_PDF_VARIATIONS; i++) {
            double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations + nominal again
            double this_gen_weight = gen_weight * pdf_weight;
            filler->fillResponseTH2(filler->responsePDFVariations().at(i), reco_weight, this_gen_weight);
          }
        }
      }

      // we need to add in an extra part, such that the 1D projection on the gen axis
      // agrees with the 1D gen histogram
      // i.e. account for the difference between the reco_weight (goes into event.weight) and gen_weight
      // we use the underflow bin for this
      // (0 as bin edges are tunfold bin numbers, not physical values, so the first bin is 0-1,
      // not to be confused with ROOT binning, starting at 1 :s)
      // double corr_weight = gen_weight * (1 - reco_weight);
      double corr_weight = gen_weight - weight;
      int underflow_bin = 0;
      fill_th2_check(h_tu_response_pt, genBinPt, underflow_bin, corr_weight);

      // do split bit
      if (onlyFillResponse && doMCsplit_) {
        h_tu_response_pt_split->Fill(genBinPt, recBinPt, weight);
        // extra part for correct weighting
        h_tu_response_pt_split->Fill(genBinPt, underflow_bin, corr_weight);
      } // end of doMCsplit_

      // do PDF variations
      if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
        for (int i=0; i<N_PDF_VARIATIONS; i++) {
          double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations + nominal again
          double this_gen_weight = gen_weight * pdf_weight;
          double this_weight = this_gen_weight * reco_weight;
          double this_corr_weight = this_gen_weight * (1 - reco_weight);
          h_tu_response_pt_PDF_variations.at(i)->Fill(genBinPt, recBinPt, this_weight);
          h_tu_response_pt_PDF_variations.at(i)->Fill(genBinPt, recBinPt, this_corr_weight);
        }
      } // end of PDF variations

    } // end of for loop over genjets
  } // end if is_mc_ && passGen

}


TH1D * QGAnalysisUnfoldHists::copy_book_th1d(TH1 * h, const std::string & newName) {
  // DO NOT USE h->GetXaxis()->GetXbins()->GetArray() to clone bins, it just doesn't work
  return book<TH1D>(newName.c_str(),
                    h->GetTitle(),
                    h->GetNbinsX(),
                    h->GetXaxis()->GetXmin(),
                    h->GetXaxis()->GetXmax());
}


TH2D * QGAnalysisUnfoldHists::copy_book_th2d(TH2 * h, const std::string & newName) {
  return book<TH2D>(newName.c_str(),
                    h->GetTitle(),
                    h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                    h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

}


QGAnalysisUnfoldHists::~QGAnalysisUnfoldHists(){}

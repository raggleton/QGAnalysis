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
  }
  doPDFvariations_ = (string2bool(ctx.get("PDFvariations", "false")) && is_mc_);
  if (doPDFvariations_) {
    cout << "Doing PDF variations in " << dirname << endl;
    doJackknifeVariations_ = false;
    cout << "Turning off jackknife variations" << endl;
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

  // remove the first 2 bins (15, 22) since we have no entries there in reco
  std::vector<double> pt_bin_edges_reco_underflow_unfold(pt_bin_edges_reco_underflow.begin()+2, pt_bin_edges_reco_underflow.end());
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
  detector_distribution_underflow_LHA->AddAxis("LHA", Binning::nbins_lha_reco, Binning::lha_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_LHA->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_LHA = detector_tu_binning_LHA->AddBinning("detector");
  detector_distribution_LHA->AddAxis("LHA", Binning::nbins_lha_reco, Binning::lha_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_LHA->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_LHA = new TUnfoldBinning("generator");
  generator_distribution_underflow_LHA = generator_tu_binning_LHA->AddBinning("signal_underflow");
  generator_distribution_underflow_LHA->AddAxis("LHA", Binning::nbins_lha_gen, Binning::lha_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_LHA->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_LHA = generator_tu_binning_LHA->AddBinning("signal");
  generator_distribution_LHA->AddAxis("LHA", Binning::nbins_lha_gen, Binning::lha_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_LHA->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  // Response matrix
  // make tmp copies which we can then copy and use with book<>
  TH2 * h_tu_response_LHA_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_LHA, detector_tu_binning_LHA, "tu_LHA_GenReco");
  h_tu_response_LHA = copy_book_th2d(h_tu_response_LHA_tmp, "tu_LHA_GenReco_all");
  h_tu_response_LHA_split = copy_book_th2d(h_tu_response_LHA_tmp, "tu_LHA_GenReco_split");
  delete h_tu_response_LHA_tmp;

  // detector histograms
  TH1 * h_tu_reco_LHA_tmp = detector_tu_binning_LHA->CreateHistogram("hist_LHA_reco");
  h_tu_reco_LHA = copy_book_th1d(h_tu_reco_LHA_tmp, "hist_LHA_reco_all");
  h_tu_reco_LHA_split = copy_book_th1d(h_tu_reco_LHA_tmp, "hist_LHA_reco_split");
  // for fakes, detector binning
  h_tu_reco_LHA_fake = copy_book_th1d(h_tu_reco_LHA_tmp, "hist_LHA_reco_fake_all");
  h_tu_reco_LHA_fake_split = copy_book_th1d(h_tu_reco_LHA_tmp, "hist_LHA_reco_fake_split");
  delete h_tu_reco_LHA_tmp;

  // truth histograms
  TH1 * h_tu_gen_LHA_tmp = generator_tu_binning_LHA->CreateHistogram("hist_LHA_truth");
  h_tu_gen_LHA = copy_book_th1d(h_tu_gen_LHA_tmp, "hist_LHA_truth_all");
  h_tu_gen_LHA_split = copy_book_th1d(h_tu_gen_LHA_tmp, "hist_LHA_truth_split");
  delete h_tu_gen_LHA_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_response_LHA_jackknife_variations.push_back(copy_book_th2d(h_tu_response_LHA, TString::Format("tu_LHA_GenReco_all_jackknife_%d", i).Data()));
      h_tu_reco_LHA_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_LHA, TString::Format("hist_LHA_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_LHA_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_LHA, TString::Format("hist_LHA_truth_all_jackknife_%d", i).Data()));
    }
  }

  // PDF variations
  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_LHA_PDF_variations.push_back(copy_book_th1d(h_tu_reco_LHA, TString::Format("hist_LHA_reco_all_PDF_%d", i).Data()));
      h_tu_gen_LHA_PDF_variations.push_back(copy_book_th1d(h_tu_gen_LHA, TString::Format("hist_LHA_truth_all_PDF_%d", i).Data()));
      h_tu_response_LHA_PDF_variations.push_back(copy_book_th2d(h_tu_response_LHA, TString::Format("tu_LHA_GenReco_all_PDF_%d", i).Data()));
    }
  }

  // detector variable, but using gen binning for comparison later
  h_tu_reco_LHA_gen_binning = copy_book_th1d(h_tu_gen_LHA, "hist_LHA_reco_gen_binning");
  h_tu_reco_LHA_gen_binning_split = copy_book_th1d(h_tu_gen_LHA, "hist_LHA_reco_gen_binning_split");
  // for fakes, gen binning
  h_tu_reco_LHA_fake_gen_binning = copy_book_th1d(h_tu_gen_LHA, "hist_LHA_reco_fake_gen_binning");


  // Charged LHA
  // -------------------------------------
  detector_tu_binning_LHA_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_LHA_charged = detector_tu_binning_LHA_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_LHA_charged->AddAxis("LHA_charged", Binning::nbins_lha_charged_reco, Binning::lha_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_LHA_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_LHA_charged = detector_tu_binning_LHA_charged->AddBinning("detector");
  detector_distribution_LHA_charged->AddAxis("LHA_charged", Binning::nbins_lha_charged_reco, Binning::lha_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_LHA_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_LHA_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_LHA_charged = generator_tu_binning_LHA_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_LHA_charged->AddAxis("LHA_charged", Binning::nbins_lha_charged_gen, Binning::lha_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_LHA_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_LHA_charged = generator_tu_binning_LHA_charged->AddBinning("signal");
  generator_distribution_LHA_charged->AddAxis("LHA_charged", Binning::nbins_lha_charged_gen, Binning::lha_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_LHA_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_LHA_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_LHA_charged, detector_tu_binning_LHA_charged, "tu_LHA_charged_GenReco");
  h_tu_response_LHA_charged = copy_book_th2d(h_tu_response_LHA_charged_tmp, "tu_LHA_charged_GenReco_all");
  h_tu_response_LHA_charged_split = copy_book_th2d(h_tu_response_LHA_charged_tmp, "tu_LHA_charged_GenReco_split");
  delete h_tu_response_LHA_charged_tmp;

  TH1 * h_tu_reco_LHA_charged_tmp = detector_tu_binning_LHA_charged->CreateHistogram("hist_LHA_charged_reco");
  h_tu_reco_LHA_charged = copy_book_th1d(h_tu_reco_LHA_charged_tmp, "hist_LHA_charged_reco_all");
  h_tu_reco_LHA_charged_split = copy_book_th1d(h_tu_reco_LHA_charged_tmp, "hist_LHA_charged_reco_split");
  h_tu_reco_LHA_charged_fake = copy_book_th1d(h_tu_reco_LHA_charged_tmp, "hist_LHA_charged_reco_fake_all");
  h_tu_reco_LHA_charged_fake_split = copy_book_th1d(h_tu_reco_LHA_charged_tmp, "hist_LHA_charged_reco_fake_split");
  delete h_tu_reco_LHA_charged_tmp;

  TH1 * h_tu_gen_LHA_charged_tmp = generator_tu_binning_LHA_charged->CreateHistogram("hist_LHA_charged_truth");
  h_tu_gen_LHA_charged = copy_book_th1d(h_tu_gen_LHA_charged_tmp, "hist_LHA_charged_truth_all");
  h_tu_gen_LHA_charged_split = copy_book_th1d(h_tu_gen_LHA_charged_tmp, "hist_LHA_charged_truth_split");
  delete h_tu_gen_LHA_charged_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_LHA_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_LHA_charged, TString::Format("hist_LHA_charged_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_LHA_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_LHA_charged, TString::Format("hist_LHA_charged_truth_all_jackknife_%d", i).Data()));
      h_tu_response_LHA_charged_jackknife_variations.push_back(copy_book_th2d(h_tu_response_LHA_charged, TString::Format("tu_LHA_charged_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_LHA_charged_PDF_variations.push_back(copy_book_th1d(h_tu_reco_LHA_charged, TString::Format("hist_LHA_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_LHA_charged_PDF_variations.push_back(copy_book_th1d(h_tu_gen_LHA_charged, TString::Format("hist_LHA_charged_truth_all_PDF_%d", i).Data()));
      h_tu_response_LHA_charged_PDF_variations.push_back(copy_book_th2d(h_tu_response_LHA_charged, TString::Format("tu_LHA_charged_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_LHA_charged_gen_binning = copy_book_th1d(h_tu_gen_LHA_charged, "hist_LHA_charged_reco_gen_binning");
  h_tu_reco_LHA_charged_gen_binning_split = copy_book_th1d(h_tu_gen_LHA_charged, "hist_LHA_charged_reco_gen_binning_split");
  h_tu_reco_LHA_charged_fake_gen_binning = copy_book_th1d(h_tu_gen_LHA_charged, "hist_LHA_charged_reco_fake_gen_binning");

  // puppi multiplicity
  // -------------------------------------
  detector_tu_binning_puppiMultiplicity = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_puppiMultiplicity = detector_tu_binning_puppiMultiplicity->AddBinning("detector_underflow");
  detector_distribution_underflow_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_puppiMultiplicity_reco, Binning::puppiMultiplicity_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_puppiMultiplicity->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_puppiMultiplicity = detector_tu_binning_puppiMultiplicity->AddBinning("detector");
  detector_distribution_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_puppiMultiplicity_reco, Binning::puppiMultiplicity_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_puppiMultiplicity->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_puppiMultiplicity = new TUnfoldBinning("generator");
  generator_distribution_underflow_puppiMultiplicity = generator_tu_binning_puppiMultiplicity->AddBinning("signal_underflow");
  generator_distribution_underflow_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_puppiMultiplicity_gen, Binning::puppiMultiplicity_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_puppiMultiplicity->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_puppiMultiplicity = generator_tu_binning_puppiMultiplicity->AddBinning("signal");
  generator_distribution_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_puppiMultiplicity_gen, Binning::puppiMultiplicity_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_puppiMultiplicity->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_puppiMultiplicity_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_puppiMultiplicity, detector_tu_binning_puppiMultiplicity, "tu_puppiMultiplicity_GenReco");
  h_tu_response_puppiMultiplicity = copy_book_th2d(h_tu_response_puppiMultiplicity_tmp, "tu_puppiMultiplicity_GenReco_all");
  h_tu_response_puppiMultiplicity_split = copy_book_th2d(h_tu_response_puppiMultiplicity_tmp, "tu_puppiMultiplicity_GenReco_split");
  delete h_tu_response_puppiMultiplicity_tmp;

  TH1 * h_tu_reco_puppiMultiplicity_tmp = detector_tu_binning_puppiMultiplicity->CreateHistogram("hist_puppiMultiplicity_reco");
  h_tu_reco_puppiMultiplicity = copy_book_th1d(h_tu_reco_puppiMultiplicity_tmp, "hist_puppiMultiplicity_reco_all");
  h_tu_reco_puppiMultiplicity_split = copy_book_th1d(h_tu_reco_puppiMultiplicity_tmp, "hist_puppiMultiplicity_reco_split");
  h_tu_reco_puppiMultiplicity_fake = copy_book_th1d(h_tu_reco_puppiMultiplicity_tmp, "hist_puppiMultiplicity_reco_fake_all");
  h_tu_reco_puppiMultiplicity_fake_split = copy_book_th1d(h_tu_reco_puppiMultiplicity_tmp, "hist_puppiMultiplicity_reco_fake_split");
  delete h_tu_reco_puppiMultiplicity_tmp;

  TH1 * h_tu_gen_puppiMultiplicity_tmp = generator_tu_binning_puppiMultiplicity->CreateHistogram("hist_puppiMultiplicity_truth");
  h_tu_gen_puppiMultiplicity = copy_book_th1d(h_tu_gen_puppiMultiplicity_tmp, "hist_puppiMultiplicity_truth_all");
  h_tu_gen_puppiMultiplicity_split = copy_book_th1d(h_tu_gen_puppiMultiplicity_tmp, "hist_puppiMultiplicity_truth_split");
  delete h_tu_gen_puppiMultiplicity_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_puppiMultiplicity_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_puppiMultiplicity, TString::Format("hist_puppiMultiplicity_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_puppiMultiplicity_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_puppiMultiplicity, TString::Format("hist_puppiMultiplicity_truth_all_jackknife_%d", i).Data()));
      h_tu_response_puppiMultiplicity_jackknife_variations.push_back(copy_book_th2d(h_tu_response_puppiMultiplicity, TString::Format("tu_puppiMultiplicity_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_puppiMultiplicity_PDF_variations.push_back(copy_book_th1d(h_tu_reco_puppiMultiplicity, TString::Format("hist_puppiMultiplicity_reco_all_PDF_%d", i).Data()));
      h_tu_gen_puppiMultiplicity_PDF_variations.push_back(copy_book_th1d(h_tu_gen_puppiMultiplicity, TString::Format("hist_puppiMultiplicity_truth_all_PDF_%d", i).Data()));
      h_tu_response_puppiMultiplicity_PDF_variations.push_back(copy_book_th2d(h_tu_response_puppiMultiplicity, TString::Format("tu_puppiMultiplicity_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_puppiMultiplicity_gen_binning = copy_book_th1d(h_tu_gen_puppiMultiplicity, "hist_puppiMultiplicity_reco_gen_binning");
  h_tu_reco_puppiMultiplicity_gen_binning_split = copy_book_th1d(h_tu_gen_puppiMultiplicity, "hist_puppiMultiplicity_reco_gen_binning_split");

  h_tu_reco_puppiMultiplicity_fake_gen_binning = copy_book_th1d(h_tu_gen_puppiMultiplicity, "hist_puppiMultiplicity_reco_fake_gen_binning");

  // Charged PUPPI multiplicity
  // -------------------------------------
  detector_tu_binning_puppiMultiplicity_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_puppiMultiplicity_charged = detector_tu_binning_puppiMultiplicity_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_puppiMultiplicity_charged_reco, Binning::puppiMultiplicity_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_puppiMultiplicity_charged = detector_tu_binning_puppiMultiplicity_charged->AddBinning("detector");
  detector_distribution_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_puppiMultiplicity_charged_reco, Binning::puppiMultiplicity_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_puppiMultiplicity_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_puppiMultiplicity_charged = generator_tu_binning_puppiMultiplicity_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_puppiMultiplicity_charged_gen, Binning::puppiMultiplicity_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_puppiMultiplicity_charged = generator_tu_binning_puppiMultiplicity_charged->AddBinning("signal");
  generator_distribution_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_puppiMultiplicity_charged_gen, Binning::puppiMultiplicity_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_puppiMultiplicity_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_puppiMultiplicity_charged, detector_tu_binning_puppiMultiplicity_charged, "tu_puppiMultiplicity_charged_GenReco");
  h_tu_response_puppiMultiplicity_charged = copy_book_th2d(h_tu_response_puppiMultiplicity_charged_tmp, "tu_puppiMultiplicity_charged_GenReco_all");
  h_tu_response_puppiMultiplicity_charged_split = copy_book_th2d(h_tu_response_puppiMultiplicity_charged_tmp, "tu_puppiMultiplicity_charged_GenReco_split");
  delete h_tu_response_puppiMultiplicity_charged_tmp;

  TH1 * h_tu_reco_puppiMultiplicity_charged_tmp = detector_tu_binning_puppiMultiplicity_charged->CreateHistogram("hist_puppiMultiplicity_charged_reco");
  h_tu_reco_puppiMultiplicity_charged = copy_book_th1d(h_tu_reco_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_reco_all");
  h_tu_reco_puppiMultiplicity_charged_split = copy_book_th1d(h_tu_reco_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_reco_split");
  h_tu_reco_puppiMultiplicity_charged_fake = copy_book_th1d(h_tu_reco_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_reco_fake_all");
  h_tu_reco_puppiMultiplicity_charged_fake_split = copy_book_th1d(h_tu_reco_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_reco_fake_split");
  delete h_tu_reco_puppiMultiplicity_charged_tmp;

  TH1 * h_tu_gen_puppiMultiplicity_charged_tmp = generator_tu_binning_puppiMultiplicity_charged->CreateHistogram("hist_puppiMultiplicity_charged_truth");
  h_tu_gen_puppiMultiplicity_charged = copy_book_th1d(h_tu_gen_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_truth_all");
  h_tu_gen_puppiMultiplicity_charged_split = copy_book_th1d(h_tu_gen_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_truth_split");
  delete h_tu_gen_puppiMultiplicity_charged_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_puppiMultiplicity_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_puppiMultiplicity_charged, TString::Format("hist_puppiMultiplicity_charged_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_puppiMultiplicity_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_puppiMultiplicity_charged, TString::Format("hist_puppiMultiplicity_charged_truth_all_jackknife_%d", i).Data()));
      h_tu_response_puppiMultiplicity_charged_jackknife_variations.push_back(copy_book_th2d(h_tu_response_puppiMultiplicity_charged, TString::Format("tu_puppiMultiplicity_charged_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_puppiMultiplicity_charged_PDF_variations.push_back(copy_book_th1d(h_tu_reco_puppiMultiplicity_charged, TString::Format("hist_puppiMultiplicity_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_puppiMultiplicity_charged_PDF_variations.push_back(copy_book_th1d(h_tu_gen_puppiMultiplicity_charged, TString::Format("hist_puppiMultiplicity_charged_truth_all_PDF_%d", i).Data()));
      h_tu_response_puppiMultiplicity_charged_PDF_variations.push_back(copy_book_th2d(h_tu_response_puppiMultiplicity_charged, TString::Format("tu_puppiMultiplicity_charged_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_puppiMultiplicity_charged_gen_binning = copy_book_th1d(h_tu_gen_puppiMultiplicity_charged, "hist_puppiMultiplicity_charged_reco_gen_binning");
  h_tu_reco_puppiMultiplicity_charged_gen_binning_split = copy_book_th1d(h_tu_gen_puppiMultiplicity_charged, "hist_puppiMultiplicity_charged_reco_gen_binning_split");

  h_tu_reco_puppiMultiplicity_charged_fake_gen_binning = copy_book_th1d(h_tu_gen_puppiMultiplicity_charged, "hist_puppiMultiplicity_charged_reco_fake_gen_binning");

  // pTD
  // -------------------------------------
  detector_tu_binning_pTD = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_pTD = detector_tu_binning_pTD->AddBinning("detector_underflow");
  detector_distribution_underflow_pTD->AddAxis("pTD", Binning::nbins_pTD_reco, Binning::pTD_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_pTD->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_pTD = detector_tu_binning_pTD->AddBinning("detector");
  detector_distribution_pTD->AddAxis("pTD", Binning::nbins_pTD_reco, Binning::pTD_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_pTD->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_pTD = new TUnfoldBinning("generator");
  generator_distribution_underflow_pTD = generator_tu_binning_pTD->AddBinning("signal_underflow");
  generator_distribution_underflow_pTD->AddAxis("pTD", Binning::nbins_pTD_gen, Binning::pTD_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_pTD->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_pTD = generator_tu_binning_pTD->AddBinning("signal");
  generator_distribution_pTD->AddAxis("pTD", Binning::nbins_pTD_gen, Binning::pTD_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_pTD->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_pTD_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_pTD, detector_tu_binning_pTD, "tu_pTD_GenReco");
  h_tu_response_pTD = copy_book_th2d(h_tu_response_pTD_tmp, "tu_pTD_GenReco_all");
  h_tu_response_pTD_split = copy_book_th2d(h_tu_response_pTD_tmp, "tu_pTD_GenReco_split");
  delete h_tu_response_pTD_tmp;

  TH1 * h_tu_reco_pTD_tmp = detector_tu_binning_pTD->CreateHistogram("hist_pTD_reco");
  h_tu_reco_pTD = copy_book_th1d(h_tu_reco_pTD_tmp, "hist_pTD_reco_all");
  h_tu_reco_pTD_split = copy_book_th1d(h_tu_reco_pTD_tmp, "hist_pTD_reco_split");
  h_tu_reco_pTD_fake = copy_book_th1d(h_tu_reco_pTD_tmp, "hist_pTD_reco_fake_all");
  h_tu_reco_pTD_fake_split = copy_book_th1d(h_tu_reco_pTD_tmp, "hist_pTD_reco_fake_split");
  delete h_tu_reco_pTD_tmp;

  TH1 * h_tu_gen_pTD_tmp = generator_tu_binning_pTD->CreateHistogram("hist_pTD_truth");
  h_tu_gen_pTD = copy_book_th1d(h_tu_gen_pTD_tmp, "hist_pTD_truth_all");
  h_tu_gen_pTD_split = copy_book_th1d(h_tu_gen_pTD_tmp, "hist_pTD_truth_split");
  delete h_tu_gen_pTD_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_pTD_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_pTD, TString::Format("hist_pTD_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_pTD_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_pTD, TString::Format("hist_pTD_truth_all_jackknife_%d", i).Data()));
      h_tu_response_pTD_jackknife_variations.push_back(copy_book_th2d(h_tu_response_pTD, TString::Format("tu_pTD_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_pTD_PDF_variations.push_back(copy_book_th1d(h_tu_reco_pTD, TString::Format("hist_pTD_reco_all_PDF_%d", i).Data()));
      h_tu_gen_pTD_PDF_variations.push_back(copy_book_th1d(h_tu_gen_pTD, TString::Format("hist_pTD_truth_all_PDF_%d", i).Data()));
      h_tu_response_pTD_PDF_variations.push_back(copy_book_th2d(h_tu_response_pTD, TString::Format("tu_pTD_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_pTD_gen_binning = copy_book_th1d(h_tu_gen_pTD, "hist_pTD_reco_gen_binning");
  h_tu_reco_pTD_gen_binning_split = copy_book_th1d(h_tu_gen_pTD, "hist_pTD_reco_gen_binning_split");

  h_tu_reco_pTD_fake_gen_binning = copy_book_th1d(h_tu_gen_pTD, "hist_pTD_reco_fake_gen_binning");

  // Charged pTD
  // -------------------------------------
  detector_tu_binning_pTD_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_pTD_charged = detector_tu_binning_pTD_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_pTD_charged->AddAxis("pTD_charged", Binning::nbins_pTD_charged_reco, Binning::pTD_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_pTD_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_pTD_charged = detector_tu_binning_pTD_charged->AddBinning("detector");
  detector_distribution_pTD_charged->AddAxis("pTD_charged", Binning::nbins_pTD_charged_reco, Binning::pTD_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_pTD_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_pTD_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_pTD_charged = generator_tu_binning_pTD_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_pTD_charged->AddAxis("pTD_charged", Binning::nbins_pTD_charged_gen, Binning::pTD_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_pTD_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_pTD_charged = generator_tu_binning_pTD_charged->AddBinning("signal");
  generator_distribution_pTD_charged->AddAxis("pTD_charged", Binning::nbins_pTD_charged_gen, Binning::pTD_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_pTD_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_pTD_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_pTD_charged, detector_tu_binning_pTD_charged, "tu_pTD_charged_GenReco");
  h_tu_response_pTD_charged = copy_book_th2d(h_tu_response_pTD_charged_tmp, "tu_pTD_charged_GenReco_all");
  h_tu_response_pTD_charged_split = copy_book_th2d(h_tu_response_pTD_charged_tmp, "tu_pTD_charged_GenReco_split");
  delete h_tu_response_pTD_charged_tmp;

  TH1 * h_tu_reco_pTD_charged_tmp = detector_tu_binning_pTD_charged->CreateHistogram("hist_pTD_charged_reco");
  h_tu_reco_pTD_charged = copy_book_th1d(h_tu_reco_pTD_charged_tmp, "hist_pTD_charged_reco_all");
  h_tu_reco_pTD_charged_split = copy_book_th1d(h_tu_reco_pTD_charged_tmp, "hist_pTD_charged_reco_split");
  h_tu_reco_pTD_charged_fake = copy_book_th1d(h_tu_reco_pTD_charged_tmp, "hist_pTD_charged_reco_fake_all");
  h_tu_reco_pTD_charged_fake_split = copy_book_th1d(h_tu_reco_pTD_charged_tmp, "hist_pTD_charged_reco_fake_split");
  delete h_tu_reco_pTD_charged_tmp;

  TH1 * h_tu_gen_pTD_charged_tmp = generator_tu_binning_pTD_charged->CreateHistogram("hist_pTD_charged_truth");
  h_tu_gen_pTD_charged = copy_book_th1d(h_tu_gen_pTD_charged_tmp, "hist_pTD_charged_truth_all");
  h_tu_gen_pTD_charged_split = copy_book_th1d(h_tu_gen_pTD_charged_tmp, "hist_pTD_charged_truth_split");
  delete h_tu_gen_pTD_charged_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_pTD_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_pTD_charged, TString::Format("hist_pTD_charged_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_pTD_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_pTD_charged, TString::Format("hist_pTD_charged_truth_all_jackknife_%d", i).Data()));
      h_tu_response_pTD_charged_jackknife_variations.push_back(copy_book_th2d(h_tu_response_pTD_charged, TString::Format("tu_pTD_charged_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_pTD_charged_PDF_variations.push_back(copy_book_th1d(h_tu_reco_pTD_charged, TString::Format("hist_pTD_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_pTD_charged_PDF_variations.push_back(copy_book_th1d(h_tu_gen_pTD_charged, TString::Format("hist_pTD_charged_truth_all_PDF_%d", i).Data()));
      h_tu_response_pTD_charged_PDF_variations.push_back(copy_book_th2d(h_tu_response_pTD_charged, TString::Format("tu_pTD_charged_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_pTD_charged_gen_binning = copy_book_th1d(h_tu_gen_pTD_charged, "hist_pTD_charged_reco_gen_binning");
  h_tu_reco_pTD_charged_gen_binning_split = copy_book_th1d(h_tu_gen_pTD_charged, "hist_pTD_charged_reco_gen_binning_split");

  h_tu_reco_pTD_charged_fake_gen_binning = copy_book_th1d(h_tu_gen_pTD_charged, "hist_pTD_charged_reco_fake_gen_binning");

  // thrust
  // -------------------------------------
  detector_tu_binning_thrust = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_thrust = detector_tu_binning_thrust->AddBinning("detector_underflow");
  detector_distribution_underflow_thrust->AddAxis("thrust", Binning::nbins_thrust_reco, Binning::thrust_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_thrust->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_thrust = detector_tu_binning_thrust->AddBinning("detector");
  detector_distribution_thrust->AddAxis("thrust", Binning::nbins_thrust_reco, Binning::thrust_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_thrust->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_thrust = new TUnfoldBinning("generator");
  generator_distribution_underflow_thrust = generator_tu_binning_thrust->AddBinning("signal_underflow");
  generator_distribution_underflow_thrust->AddAxis("thrust", Binning::nbins_thrust_gen, Binning::thrust_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_thrust->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_thrust = generator_tu_binning_thrust->AddBinning("signal");
  generator_distribution_thrust->AddAxis("thrust", Binning::nbins_thrust_gen, Binning::thrust_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_thrust->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_thrust_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_thrust, detector_tu_binning_thrust, "tu_thrust_GenReco");
  h_tu_response_thrust = copy_book_th2d(h_tu_response_thrust_tmp, "tu_thrust_GenReco_all");
  h_tu_response_thrust_split = copy_book_th2d(h_tu_response_thrust_tmp, "tu_thrust_GenReco_split");
  delete h_tu_response_thrust_tmp;

  TH1 * h_tu_reco_thrust_tmp = detector_tu_binning_thrust->CreateHistogram("hist_thrust_reco");
  h_tu_reco_thrust = copy_book_th1d(h_tu_reco_thrust_tmp, "hist_thrust_reco_all");
  h_tu_reco_thrust_split = copy_book_th1d(h_tu_reco_thrust_tmp, "hist_thrust_reco_split");
  h_tu_reco_thrust_fake = copy_book_th1d(h_tu_reco_thrust_tmp, "hist_thrust_reco_fake_all");
  h_tu_reco_thrust_fake_split = copy_book_th1d(h_tu_reco_thrust_tmp, "hist_thrust_reco_fake_split");
  delete h_tu_reco_thrust_tmp;

  TH1 * h_tu_gen_thrust_tmp = generator_tu_binning_thrust->CreateHistogram("hist_thrust_truth");
  h_tu_gen_thrust = copy_book_th1d(h_tu_gen_thrust_tmp, "hist_thrust_truth_all");
  h_tu_gen_thrust_split = copy_book_th1d(h_tu_gen_thrust_tmp, "hist_thrust_truth_split");
  delete h_tu_gen_thrust_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_thrust_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_thrust, TString::Format("hist_thrust_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_thrust_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_thrust, TString::Format("hist_thrust_truth_all_jackknife_%d", i).Data()));
      h_tu_response_thrust_jackknife_variations.push_back(copy_book_th2d(h_tu_response_thrust, TString::Format("tu_thrust_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_thrust_PDF_variations.push_back(copy_book_th1d(h_tu_reco_thrust, TString::Format("hist_thrust_reco_all_PDF_%d", i).Data()));
      h_tu_gen_thrust_PDF_variations.push_back(copy_book_th1d(h_tu_gen_thrust, TString::Format("hist_thrust_truth_all_PDF_%d", i).Data()));
      h_tu_response_thrust_PDF_variations.push_back(copy_book_th2d(h_tu_response_thrust, TString::Format("tu_thrust_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_thrust_gen_binning = copy_book_th1d(h_tu_gen_thrust, "hist_thrust_reco_gen_binning");
  h_tu_reco_thrust_gen_binning_split = copy_book_th1d(h_tu_gen_thrust, "hist_thrust_reco_gen_binning_split");

  h_tu_reco_thrust_fake_gen_binning = copy_book_th1d(h_tu_gen_thrust, "hist_thrust_reco_fake_gen_binning");

  // Charged thrust
  // -------------------------------------
  detector_tu_binning_thrust_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_thrust_charged = detector_tu_binning_thrust_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_thrust_charged->AddAxis("thrust_charged", Binning::nbins_thrust_charged_reco, Binning::thrust_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_thrust_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_thrust_charged = detector_tu_binning_thrust_charged->AddBinning("detector");
  detector_distribution_thrust_charged->AddAxis("thrust_charged", Binning::nbins_thrust_charged_reco, Binning::thrust_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_thrust_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_thrust_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_thrust_charged = generator_tu_binning_thrust_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_thrust_charged->AddAxis("thrust_charged", Binning::nbins_thrust_charged_gen, Binning::thrust_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_thrust_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_thrust_charged = generator_tu_binning_thrust_charged->AddBinning("signal");
  generator_distribution_thrust_charged->AddAxis("thrust_charged", Binning::nbins_thrust_charged_gen, Binning::thrust_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_thrust_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_thrust_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_thrust_charged, detector_tu_binning_thrust_charged, "tu_thrust_charged_GenReco");
  h_tu_response_thrust_charged = copy_book_th2d(h_tu_response_thrust_charged_tmp, "tu_thrust_charged_GenReco_all");
  h_tu_response_thrust_charged_split = copy_book_th2d(h_tu_response_thrust_charged_tmp, "tu_thrust_charged_GenReco_split");
  delete h_tu_response_thrust_charged_tmp;

  TH1 * h_tu_reco_thrust_charged_tmp = detector_tu_binning_thrust_charged->CreateHistogram("hist_thrust_charged_reco");
  h_tu_reco_thrust_charged = copy_book_th1d(h_tu_reco_thrust_charged_tmp, "hist_thrust_charged_reco_all");
  h_tu_reco_thrust_charged_split = copy_book_th1d(h_tu_reco_thrust_charged_tmp, "hist_thrust_charged_reco_split");
  h_tu_reco_thrust_charged_fake = copy_book_th1d(h_tu_reco_thrust_charged_tmp, "hist_thrust_charged_reco_fake_all");
  h_tu_reco_thrust_charged_fake_split = copy_book_th1d(h_tu_reco_thrust_charged_tmp, "hist_thrust_charged_reco_fake_split");
  delete h_tu_reco_thrust_charged_tmp;

  TH1 * h_tu_gen_thrust_charged_tmp = generator_tu_binning_thrust_charged->CreateHistogram("hist_thrust_charged_truth");
  h_tu_gen_thrust_charged = copy_book_th1d(h_tu_gen_thrust_charged_tmp, "hist_thrust_charged_truth_all");
  h_tu_gen_thrust_charged_split = copy_book_th1d(h_tu_gen_thrust_charged_tmp, "hist_thrust_charged_truth_split");
  delete h_tu_gen_thrust_charged_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_thrust_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_thrust_charged, TString::Format("hist_thrust_charged_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_thrust_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_thrust_charged, TString::Format("hist_thrust_charged_truth_all_jackknife_%d", i).Data()));
      h_tu_response_thrust_charged_jackknife_variations.push_back(copy_book_th2d(h_tu_response_thrust_charged, TString::Format("tu_thrust_charged_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_thrust_charged_PDF_variations.push_back(copy_book_th1d(h_tu_reco_thrust_charged, TString::Format("hist_thrust_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_thrust_charged_PDF_variations.push_back(copy_book_th1d(h_tu_gen_thrust_charged, TString::Format("hist_thrust_charged_truth_all_PDF_%d", i).Data()));
      h_tu_response_thrust_charged_PDF_variations.push_back(copy_book_th2d(h_tu_response_thrust_charged, TString::Format("tu_thrust_charged_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_thrust_charged_gen_binning = copy_book_th1d(h_tu_gen_thrust_charged, "hist_thrust_charged_reco_gen_binning");
  h_tu_reco_thrust_charged_gen_binning_split = copy_book_th1d(h_tu_gen_thrust_charged, "hist_thrust_charged_reco_gen_binning_split");

  h_tu_reco_thrust_charged_fake_gen_binning = copy_book_th1d(h_tu_gen_thrust_charged, "hist_thrust_charged_reco_fake_gen_binning");

  // width
  // -------------------------------------
  detector_tu_binning_width = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_width = detector_tu_binning_width->AddBinning("detector_underflow");
  detector_distribution_underflow_width->AddAxis("width", Binning::nbins_width_reco, Binning::width_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_width->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_width = detector_tu_binning_width->AddBinning("detector");
  detector_distribution_width->AddAxis("width", Binning::nbins_width_reco, Binning::width_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_width->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_width = new TUnfoldBinning("generator");
  generator_distribution_underflow_width = generator_tu_binning_width->AddBinning("signal_underflow");
  generator_distribution_underflow_width->AddAxis("width", Binning::nbins_width_gen, Binning::width_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_width->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_width = generator_tu_binning_width->AddBinning("signal");
  generator_distribution_width->AddAxis("width", Binning::nbins_width_gen, Binning::width_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_width->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_width_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_width, detector_tu_binning_width, "tu_width_GenReco");
  h_tu_response_width = copy_book_th2d(h_tu_response_width_tmp, "tu_width_GenReco_all");
  h_tu_response_width_split = copy_book_th2d(h_tu_response_width_tmp, "tu_width_GenReco_split");
  delete h_tu_response_width_tmp;

  TH1 * h_tu_reco_width_tmp = detector_tu_binning_width->CreateHistogram("hist_width_reco");
  h_tu_reco_width = copy_book_th1d(h_tu_reco_width_tmp, "hist_width_reco_all");
  h_tu_reco_width_split = copy_book_th1d(h_tu_reco_width_tmp, "hist_width_reco_split");
  h_tu_reco_width_fake = copy_book_th1d(h_tu_reco_width_tmp, "hist_width_reco_fake_all");
  h_tu_reco_width_fake_split = copy_book_th1d(h_tu_reco_width_tmp, "hist_width_reco_fake_split");
  delete h_tu_reco_width_tmp;

  TH1 * h_tu_gen_width_tmp = generator_tu_binning_width->CreateHistogram("hist_width_truth");
  h_tu_gen_width = copy_book_th1d(h_tu_gen_width_tmp, "hist_width_truth_all");
  h_tu_gen_width_split = copy_book_th1d(h_tu_gen_width_tmp, "hist_width_truth_split");
  delete h_tu_gen_width_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_width_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_width, TString::Format("hist_width_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_width_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_width, TString::Format("hist_width_truth_all_jackknife_%d", i).Data()));
      h_tu_response_width_jackknife_variations.push_back(copy_book_th2d(h_tu_response_width, TString::Format("tu_width_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_width_PDF_variations.push_back(copy_book_th1d(h_tu_reco_width, TString::Format("hist_width_reco_all_PDF_%d", i).Data()));
      h_tu_gen_width_PDF_variations.push_back(copy_book_th1d(h_tu_gen_width, TString::Format("hist_width_truth_all_PDF_%d", i).Data()));
      h_tu_response_width_PDF_variations.push_back(copy_book_th2d(h_tu_response_width, TString::Format("tu_width_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_width_gen_binning = copy_book_th1d(h_tu_gen_width, "hist_width_reco_gen_binning");
  h_tu_reco_width_gen_binning_split = copy_book_th1d(h_tu_gen_width, "hist_width_reco_gen_binning_split");

  h_tu_reco_width_fake_gen_binning = copy_book_th1d(h_tu_gen_width, "hist_width_reco_fake_gen_binning");

  // Charged width
  // -------------------------------------
  detector_tu_binning_width_charged = new TUnfoldBinning("detectorall");
  detector_distribution_underflow_width_charged = detector_tu_binning_width_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_width_charged->AddAxis("width_charged", Binning::nbins_width_charged_reco, Binning::width_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_width_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow_unfold.data(), false, false);

  detector_distribution_width_charged = detector_tu_binning_width_charged->AddBinning("detector");
  detector_distribution_width_charged->AddAxis("width_charged", Binning::nbins_width_charged_reco, Binning::width_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_width_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), false, pt_of);

  generator_tu_binning_width_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_width_charged = generator_tu_binning_width_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_width_charged->AddAxis("width_charged", Binning::nbins_width_charged_gen, Binning::width_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_width_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, false);

  generator_distribution_width_charged = generator_tu_binning_width_charged->AddBinning("signal");
  generator_distribution_width_charged->AddAxis("width_charged", Binning::nbins_width_charged_gen, Binning::width_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_width_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), false, pt_of);

  TH2 * h_tu_response_width_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_width_charged, detector_tu_binning_width_charged, "tu_width_charged_GenReco");
  h_tu_response_width_charged = copy_book_th2d(h_tu_response_width_charged_tmp, "tu_width_charged_GenReco_all");
  h_tu_response_width_charged_split = copy_book_th2d(h_tu_response_width_charged_tmp, "tu_width_charged_GenReco_split");
  delete h_tu_response_width_charged_tmp;

  TH1 * h_tu_reco_width_charged_tmp = detector_tu_binning_width_charged->CreateHistogram("hist_width_charged_reco");
  h_tu_reco_width_charged = copy_book_th1d(h_tu_reco_width_charged_tmp, "hist_width_charged_reco_all");
  h_tu_reco_width_charged_split = copy_book_th1d(h_tu_reco_width_charged_tmp, "hist_width_charged_reco_split");
  h_tu_reco_width_charged_fake = copy_book_th1d(h_tu_reco_width_charged_tmp, "hist_width_charged_reco_fake_all");
  h_tu_reco_width_charged_fake_split = copy_book_th1d(h_tu_reco_width_charged_tmp, "hist_width_charged_reco_fake_split");
  delete h_tu_reco_width_charged_tmp;

  TH1 * h_tu_gen_width_charged_tmp = generator_tu_binning_width_charged->CreateHistogram("hist_width_charged_truth");
  h_tu_gen_width_charged = copy_book_th1d(h_tu_gen_width_charged_tmp, "hist_width_charged_truth_all");
  h_tu_gen_width_charged_split = copy_book_th1d(h_tu_gen_width_charged_tmp, "hist_width_charged_truth_split");
  delete h_tu_gen_width_charged_tmp;

  if (doJackknifeVariations_) {
    for (uint i=0; i < N_JACKKNIFE_VARIATIONS; i++) {
      h_tu_reco_width_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_reco_width_charged, TString::Format("hist_width_charged_reco_all_jackknife_%d", i).Data()));
      h_tu_gen_width_charged_jackknife_variations.push_back(copy_book_th1d(h_tu_gen_width_charged, TString::Format("hist_width_charged_truth_all_jackknife_%d", i).Data()));
      h_tu_response_width_charged_jackknife_variations.push_back(copy_book_th2d(h_tu_response_width_charged, TString::Format("tu_width_charged_GenReco_all_jackknife_%d", i).Data()));
    }
  }

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_width_charged_PDF_variations.push_back(copy_book_th1d(h_tu_reco_width_charged, TString::Format("hist_width_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_width_charged_PDF_variations.push_back(copy_book_th1d(h_tu_gen_width_charged, TString::Format("hist_width_charged_truth_all_PDF_%d", i).Data()));
      h_tu_response_width_charged_PDF_variations.push_back(copy_book_th2d(h_tu_response_width_charged, TString::Format("tu_width_charged_GenReco_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_width_charged_gen_binning = copy_book_th1d(h_tu_gen_width_charged, "hist_width_charged_reco_gen_binning");
  h_tu_reco_width_charged_gen_binning_split = copy_book_th1d(h_tu_gen_width_charged, "hist_width_charged_reco_gen_binning_split");

  h_tu_reco_width_charged_fake_gen_binning = copy_book_th1d(h_tu_gen_width_charged, "hist_width_charged_reco_fake_gen_binning");

  if (is_mc_) {
    genJetsLambda_handle = ctx.get_handle< std::vector<GenJetLambdaBundle> > (gen_jetlambda_handle_name);
    pass_gen_handle = ctx.get_handle<bool> (gen_sel_handle_name);
  }

  jetsLambda_handle = ctx.get_handle< std::vector<JetLambdaBundle> > (reco_jetlambda_handle_name);

  pass_reco_handle = ctx.get_handle<bool> (reco_sel_handle_name);

  // TODO: automate setting up of hists - let the LambdaHistsFiller own them?
  LHA_hist_filler.reset(new LambdaHistsFiller(Cuts::lha_pf_args,
                                              detector_tu_binning_LHA,
                                              Cuts::lha_gen_args,
                                              generator_tu_binning_LHA));
  LHA_hist_filler->assignRecoHists(h_tu_reco_LHA,
                                   h_tu_reco_LHA_split,
                                   h_tu_reco_LHA_gen_binning,
                                   h_tu_reco_LHA_gen_binning_split,
                                   h_tu_reco_LHA_fake,
                                   h_tu_reco_LHA_fake_split,
                                   h_tu_reco_LHA_fake_gen_binning,
                                   &h_tu_reco_LHA_jackknife_variations,
                                   &h_tu_reco_LHA_PDF_variations);
  LHA_hist_filler->assignGenHists(h_tu_gen_LHA,
                                  h_tu_gen_LHA_split,
                                  &h_tu_gen_LHA_jackknife_variations,
                                  &h_tu_gen_LHA_PDF_variations);
  LHA_hist_filler->assignResponseHists(h_tu_response_LHA,
                                       h_tu_response_LHA_split,
                                       &h_tu_response_LHA_jackknife_variations,
                                       &h_tu_response_LHA_PDF_variations);

  puppiMultiplicity_hist_filler.reset(new LambdaHistsFiller(Cuts::mult_pf_args,
                                                            detector_tu_binning_puppiMultiplicity,
                                                            Cuts::mult_gen_args,
                                                            generator_tu_binning_puppiMultiplicity));
  puppiMultiplicity_hist_filler->assignRecoHists(h_tu_reco_puppiMultiplicity,
                                                 h_tu_reco_puppiMultiplicity_split,
                                                 h_tu_reco_puppiMultiplicity_gen_binning,
                                                 h_tu_reco_puppiMultiplicity_gen_binning_split,
                                                 h_tu_reco_puppiMultiplicity_fake,
                                                 h_tu_reco_puppiMultiplicity_fake_split,
                                                 h_tu_reco_puppiMultiplicity_fake_gen_binning,
                                                 &h_tu_reco_puppiMultiplicity_jackknife_variations,
                                                 &h_tu_reco_puppiMultiplicity_PDF_variations);
  puppiMultiplicity_hist_filler->assignGenHists(h_tu_gen_puppiMultiplicity,
                                                h_tu_gen_puppiMultiplicity_split,
                                                &h_tu_gen_puppiMultiplicity_jackknife_variations,
                                                &h_tu_gen_puppiMultiplicity_PDF_variations);
  puppiMultiplicity_hist_filler->assignResponseHists(h_tu_response_puppiMultiplicity,
                                                     h_tu_response_puppiMultiplicity_split,
                                                     &h_tu_response_puppiMultiplicity_jackknife_variations,
                                                     &h_tu_response_puppiMultiplicity_PDF_variations);

  pTD_hist_filler.reset(new LambdaHistsFiller(Cuts::pTD_pf_args,
                                              detector_tu_binning_pTD,
                                              Cuts::pTD_gen_args,
                                              generator_tu_binning_pTD));
  pTD_hist_filler->assignRecoHists(h_tu_reco_pTD,
                                  h_tu_reco_pTD_split,
                                  h_tu_reco_pTD_gen_binning,
                                  h_tu_reco_pTD_gen_binning_split,
                                  h_tu_reco_pTD_fake,
                                  h_tu_reco_pTD_fake_split,
                                  h_tu_reco_pTD_fake_gen_binning,
                                  &h_tu_reco_pTD_jackknife_variations,
                                  &h_tu_reco_pTD_PDF_variations);
  pTD_hist_filler->assignGenHists(h_tu_gen_pTD,
                                  h_tu_gen_pTD_split,
                                  &h_tu_gen_pTD_jackknife_variations,
                                  &h_tu_gen_pTD_PDF_variations);
  pTD_hist_filler->assignResponseHists(h_tu_response_pTD,
                                       h_tu_response_pTD_split,
                                       &h_tu_response_pTD_jackknife_variations,
                                       &h_tu_response_pTD_PDF_variations);

  thrust_hist_filler.reset(new LambdaHistsFiller(Cuts::thrust_pf_args,
                                                 detector_tu_binning_thrust,
                                                 Cuts::thrust_gen_args,
                                                 generator_tu_binning_thrust));
  thrust_hist_filler->assignRecoHists(h_tu_reco_thrust,
                                     h_tu_reco_thrust_split,
                                     h_tu_reco_thrust_gen_binning,
                                     h_tu_reco_thrust_gen_binning_split,
                                     h_tu_reco_thrust_fake,
                                     h_tu_reco_thrust_fake_split,
                                     h_tu_reco_thrust_fake_gen_binning,
                                     &h_tu_reco_thrust_jackknife_variations,
                                     &h_tu_reco_thrust_PDF_variations);
  thrust_hist_filler->assignGenHists(h_tu_gen_thrust,
                                     h_tu_gen_thrust_split,
                                     &h_tu_gen_thrust_jackknife_variations,
                                     &h_tu_gen_thrust_PDF_variations);
  thrust_hist_filler->assignResponseHists(h_tu_response_thrust,
                                          h_tu_response_thrust_split,
                                          &h_tu_response_thrust_jackknife_variations,
                                          &h_tu_response_thrust_PDF_variations);

  width_hist_filler.reset(new LambdaHistsFiller(Cuts::width_pf_args,
                                                detector_tu_binning_width,
                                                Cuts::width_gen_args,
                                                generator_tu_binning_width));
  width_hist_filler->assignRecoHists(h_tu_reco_width,
                                    h_tu_reco_width_split,
                                    h_tu_reco_width_gen_binning,
                                    h_tu_reco_width_gen_binning_split,
                                    h_tu_reco_width_fake,
                                    h_tu_reco_width_fake_split,
                                    h_tu_reco_width_fake_gen_binning,
                                    &h_tu_reco_width_jackknife_variations,
                                    &h_tu_reco_width_PDF_variations);
  width_hist_filler->assignGenHists(h_tu_gen_width,
                                    h_tu_gen_width_split,
                                    &h_tu_gen_width_jackknife_variations,
                                    &h_tu_gen_width_PDF_variations);
  width_hist_filler->assignResponseHists(h_tu_response_width,
                                         h_tu_response_width_split,
                                         &h_tu_response_width_jackknife_variations,
                                         &h_tu_response_width_PDF_variations);

  chargedPlusNeutralHistFillers = {
    LHA_hist_filler.get(), // use .get() to get raw ptr, since we don't want to own the object
    puppiMultiplicity_hist_filler.get(),
    pTD_hist_filler.get(),
    thrust_hist_filler.get(),
    width_hist_filler.get(),
  };

  LHA_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::lha_pf_args,
                                                      detector_tu_binning_LHA_charged,
                                                      Cuts::lha_gen_args,
                                                      generator_tu_binning_LHA_charged));
  LHA_charged_hist_filler->assignRecoHists(h_tu_reco_LHA_charged,
                                           h_tu_reco_LHA_charged_split,
                                           h_tu_reco_LHA_charged_gen_binning,
                                           h_tu_reco_LHA_charged_gen_binning_split,
                                           h_tu_reco_LHA_charged_fake,
                                           h_tu_reco_LHA_charged_fake_split,
                                           h_tu_reco_LHA_charged_fake_gen_binning,
                                           &h_tu_reco_LHA_charged_jackknife_variations,
                                           &h_tu_reco_LHA_charged_PDF_variations);
  LHA_charged_hist_filler->assignResponseHists(h_tu_response_LHA_charged,
                                               h_tu_response_LHA_charged_split,
                                               &h_tu_response_LHA_charged_jackknife_variations,
                                               &h_tu_response_LHA_charged_PDF_variations);

  LHA_charged_hist_filler->assignGenHists(h_tu_gen_LHA_charged,
                                          h_tu_gen_LHA_charged_split,
                                          &h_tu_gen_LHA_charged_jackknife_variations,
                                          &h_tu_gen_LHA_charged_PDF_variations);

  puppiMultiplicity_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::mult_pf_args,
                                                                    detector_tu_binning_puppiMultiplicity_charged,
                                                                    Cuts::mult_gen_args,
                                                                    generator_tu_binning_puppiMultiplicity_charged));
  puppiMultiplicity_charged_hist_filler->assignRecoHists(h_tu_reco_puppiMultiplicity_charged,
                                                         h_tu_reco_puppiMultiplicity_charged_split,
                                                         h_tu_reco_puppiMultiplicity_charged_gen_binning,
                                                         h_tu_reco_puppiMultiplicity_charged_gen_binning_split,
                                                         h_tu_reco_puppiMultiplicity_charged_fake,
                                                         h_tu_reco_puppiMultiplicity_charged_fake_split,
                                                         h_tu_reco_puppiMultiplicity_charged_fake_gen_binning,
                                                         &h_tu_reco_puppiMultiplicity_charged_jackknife_variations,
                                                         &h_tu_reco_puppiMultiplicity_charged_PDF_variations);
  puppiMultiplicity_charged_hist_filler->assignGenHists(h_tu_gen_puppiMultiplicity_charged,
                                                        h_tu_gen_puppiMultiplicity_charged_split,
                                                        &h_tu_gen_puppiMultiplicity_charged_jackknife_variations,
                                                        &h_tu_gen_puppiMultiplicity_charged_PDF_variations);
  puppiMultiplicity_charged_hist_filler->assignResponseHists(h_tu_response_puppiMultiplicity_charged,
                                                             h_tu_response_puppiMultiplicity_charged_split,
                                                             &h_tu_response_puppiMultiplicity_charged_jackknife_variations,
                                                             &h_tu_response_puppiMultiplicity_charged_PDF_variations);

  pTD_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::pTD_pf_args,
                                                      detector_tu_binning_pTD_charged,
                                                      Cuts::pTD_gen_args,
                                                      generator_tu_binning_pTD_charged));
  pTD_charged_hist_filler->assignRecoHists(h_tu_reco_pTD_charged,
                                           h_tu_reco_pTD_charged_split,
                                           h_tu_reco_pTD_charged_gen_binning,
                                           h_tu_reco_pTD_charged_gen_binning_split,
                                           h_tu_reco_pTD_charged_fake,
                                           h_tu_reco_pTD_charged_fake_split,
                                           h_tu_reco_pTD_charged_fake_gen_binning,
                                           &h_tu_reco_pTD_charged_jackknife_variations,
                                           &h_tu_reco_pTD_charged_PDF_variations);
  pTD_charged_hist_filler->assignGenHists(h_tu_gen_pTD_charged,
                                          h_tu_gen_pTD_charged_split,
                                          &h_tu_gen_pTD_charged_jackknife_variations,
                                          &h_tu_gen_pTD_charged_PDF_variations);
  pTD_charged_hist_filler->assignResponseHists(h_tu_response_pTD_charged,
                                               h_tu_response_pTD_charged_split,
                                               &h_tu_response_pTD_charged_jackknife_variations,
                                               &h_tu_response_pTD_charged_PDF_variations);

  thrust_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::thrust_pf_args,
                                                         detector_tu_binning_thrust_charged,
                                                         Cuts::thrust_gen_args,
                                                         generator_tu_binning_thrust_charged));
  thrust_charged_hist_filler->assignRecoHists(h_tu_reco_thrust_charged,
                                              h_tu_reco_thrust_charged_split,
                                              h_tu_reco_thrust_charged_gen_binning,
                                              h_tu_reco_thrust_charged_gen_binning_split,
                                              h_tu_reco_thrust_charged_fake,
                                              h_tu_reco_thrust_charged_fake_split,
                                              h_tu_reco_thrust_charged_fake_gen_binning,
                                              &h_tu_reco_thrust_charged_jackknife_variations,
                                              &h_tu_reco_thrust_charged_PDF_variations);
  thrust_charged_hist_filler->assignGenHists(h_tu_gen_thrust_charged,
                                             h_tu_gen_thrust_charged_split,
                                             &h_tu_gen_thrust_charged_jackknife_variations,
                                             &h_tu_gen_thrust_charged_PDF_variations);
  thrust_charged_hist_filler->assignResponseHists(h_tu_response_thrust_charged,
                                                  h_tu_response_thrust_charged_split,
                                                  &h_tu_response_thrust_charged_jackknife_variations,
                                                  &h_tu_response_thrust_charged_PDF_variations);

  width_charged_hist_filler.reset(new LambdaHistsFiller(Cuts::width_pf_args,
                                                        detector_tu_binning_width_charged,
                                                        Cuts::width_gen_args,
                                                        generator_tu_binning_width_charged));
  width_charged_hist_filler->assignRecoHists(h_tu_reco_width_charged,
                                             h_tu_reco_width_charged_split,
                                             h_tu_reco_width_charged_gen_binning,
                                             h_tu_reco_width_charged_gen_binning_split,
                                             h_tu_reco_width_charged_fake,
                                             h_tu_reco_width_charged_fake_split,
                                             h_tu_reco_width_charged_fake_gen_binning,
                                             &h_tu_reco_width_charged_jackknife_variations,
                                             &h_tu_reco_width_charged_PDF_variations);
  width_charged_hist_filler->assignGenHists(h_tu_gen_width_charged,
                                            h_tu_gen_width_charged_split,
                                            &h_tu_gen_width_charged_jackknife_variations,
                                            &h_tu_gen_width_charged_PDF_variations);
  width_charged_hist_filler->assignResponseHists(h_tu_response_width_charged,
                                                 h_tu_response_width_charged_split,
                                                 &h_tu_response_width_charged_jackknife_variations,
                                                 &h_tu_response_width_charged_PDF_variations);
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
            filler->fillRecoTH1(filler->recoJackknifeVariations()->at(index), weight);
          }
        }

        // Fill PDF variations by varying weight
        if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
          for (int i=0; i<N_PDF_VARIATIONS; i++) {
            double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations, then comes the nominal weight again
            double this_weight = weight * pdf_weight;
            filler->fillRecoTH1(filler->recoPDFVariations()->at(i), this_weight);
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
            filler->fillGenTH1(filler->genJackknifeVariations()->at(index), gen_weight);
          }
        }

        // Fill PDF variations by varying weight
        if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
          for (int i=0; i<N_PDF_VARIATIONS; i++) {
            double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations, then comes the nominal weight again
            double this_weight = gen_weight * pdf_weight;
            filler->fillGenTH1(filler->genPDFVariations()->at(i), this_weight);
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
      float lha, jet_pt;
      bool thisPassReco = passReco;
      bool thisPassRecoCharged = passReco;
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

        if (recoInd >= 0) {
          const Jet & thisjet = jetLambdas.at(recoInd).jet;
          float jet_pt = useBinningValue_ ? event.get(pt_binning_reco_handle) : thisjet.pt();

          for (auto filler : chargedPlusNeutralHistFillers) {
            // check if the reco jet already in this filler is the same one
            // and if not, update it
            if (filler->recoJetPt() != jet_pt) {
              filler->setupReco(jet_pt, jetLambdas.at(recoInd).getLambdaCalculator(false, doGroomed_), passReco);
            }
          }
          for (auto filler : chargedOnlyHistFillers) {
            // check if the reco jet already in this filler is the same one
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
            filler->fillResponseTH2(filler->responseJackknifeVariations()->at(index), reco_weight, gen_weight);
          }
        }

        if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
          for (int i=0; i<N_PDF_VARIATIONS; i++) {
            double pdf_weight = event.genInfo->systweights().at(i+10) / event.genInfo->systweights().at(9); // +10 as first 9 are scale variations + nominal again
            double this_gen_weight = gen_weight * pdf_weight;
            filler->fillResponseTH2(filler->responsePDFVariations()->at(i), reco_weight, this_gen_weight);
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

LambdaHistsFiller::LambdaHistsFiller(const PFLambdaArgs & pfLambdaArgs,
                                     TUnfoldBinning * recoBinning,
                                     const GenLambdaArgs & genLambdaArgs,
                                     TUnfoldBinning * genBinning):
pfLambdaArgs_(pfLambdaArgs),
recoBinning_(recoBinning),
passReco_(false),
genLambdaArgs_(genLambdaArgs),
genBinning_(genBinning),
passGen_(false),
recoBin_(0),
recoBinGenBinning_(0),
genBin_(0),
recoJetPt_(0),
genJetPt_(0)
{}

void LambdaHistsFiller::setPassReco(bool passReco) {
  passReco_ = passReco;
  if (!passReco_) {
    // if fail, then reset bins
    recoBin_ = 0;
    recoBinGenBinning_ = 0;
  }
}

void LambdaHistsFiller::setupReco(float recoJetPt,
                                  const LambdaCalculator<PFParticle> & recoJetCalc,
                                  bool passReco) {
  // reset bin numbers, bools
  recoBin_ = 0;
  recoBinGenBinning_ = 0;
  setPassReco(passReco);
  recoJetPt_ = recoJetPt;

  // calculate reco lambda variable, if those selections passed
  // if value is < 0 then indicates not enough constituents, so fail that selection
  // Get the correct TUnfoldBinning node (ie if underflow or not)
  // also calculate global bin numbers for TUnfold hists
  if (passReco_) {
    double recoLambda = recoJetCalc.getLambda(pfLambdaArgs_.kappa, pfLambdaArgs_.beta, pfLambdaArgs_.id);
    if (recoLambda < 0) {
      setPassReco(false);
    } else {
      bool isUnderflow = recoJetPt_ < Binning::pt_bin_edges_reco[0];
      if (recoBinning_ == nullptr) throw std::runtime_error("LambdaHistsFiller::recoBinning is nullptr");
      std::string nodeName = isUnderflow ? "detector_underflow" : "detector";
      const TUnfoldBinning * thisRecoBinning = recoBinning_->FindNode(nodeName.c_str());
      if (thisRecoBinning == nullptr) throw std::runtime_error("Cannot find TUnfoldBinning node with name " + nodeName);
      recoBin_ = thisRecoBinning->GetGlobalBinNumber(recoLambda, recoJetPt_);

      nodeName = isUnderflow ? "signal_underflow" : "signal";
      if (genBinning_ == nullptr) throw std::runtime_error("LambdaHistsFiller::genBinning is nullptr");
      const TUnfoldBinning * thisGenBinning = genBinning_->FindNode(nodeName.c_str());
      if (thisGenBinning == nullptr) throw std::runtime_error("Cannot find TUnfoldBinning node with name " + nodeName);
      recoBinGenBinning_ = thisGenBinning->GetGlobalBinNumber(recoLambda, recoJetPt_);
    }
  }
}

void LambdaHistsFiller::setPassGen(bool passGen) {
  passGen_ = passGen;
  if (!passGen_) {
    // if fail, then reset bin
    genBin_ = 0;
  }
}

void LambdaHistsFiller::setupGen(float genJetPt,
                                 const LambdaCalculator<GenParticle> & genJetCalc,
                                 bool passGen) {
  // reset bin numbers, bools
  genBin_ = 0;
  setPassGen(passGen);
  genJetPt_ = genJetPt;

  // calculate gen lambda variable, if those selections passed
  // if value is < 0 then indicates not enough constituents, so fail that selection
  // Get the correct TUnfoldBinning node (ie if underflow or not)
  // also calculate global bin numbers for TUnfold hists
  if (passGen_) {
    double genLambda = genJetCalc.getLambda(genLambdaArgs_.kappa, genLambdaArgs_.beta, genLambdaArgs_.id);
    if (genLambda < 0) {
      setPassGen(false);
    } else {
      bool isUnderflow = genJetPt_ < Binning::pt_bin_edges_gen[0];
      std::string nodeName = isUnderflow ? "signal_underflow" : "signal";
      const TUnfoldBinning * thisGenBinning = genBinning_->FindNode(nodeName.c_str());
      if (thisGenBinning == nullptr) throw std::runtime_error("Cannot find TUnfoldBinning node with name " + nodeName);
      genBin_ = thisGenBinning->GetGlobalBinNumber(genLambda, genJetPt_);
    }
  }
}

void LambdaHistsFiller::assignRecoHists(TH1D * allReco,
                                        TH1D * splitReco,
                                        TH1D * allRecoGenBinning,
                                        TH1D * splitRecoGenBinning,
                                        TH1D * allRecoFakes,
                                        TH1D * splitRecoFakes,
                                        TH1D * allRecoFakesGenBinning,
                                        std::vector<TH1D*> * jackknifeVariations,
                                        std::vector<TH1D*> * PDFVariations) {
  recoHist_ = allReco;
  recoSplitHist_ = splitReco;
  recoHistGenBinning_ = allRecoGenBinning;
  recoSplitHistGenBinning_ = splitRecoGenBinning;
  recoHistFakes_ = allRecoFakes;
  recoSplitHistFakes_ = splitRecoFakes;
  recoHistFakesGenBinning_ = allRecoFakesGenBinning;
  recoJackknifeVariations_ = jackknifeVariations;
  recoPDFVariations_ = PDFVariations;
}

void LambdaHistsFiller::assignGenHists(TH1D * allGen,
                                       TH1D * splitGen,
                                       std::vector<TH1D*> * jackknifeVariations,
                                       std::vector<TH1D*> * PDFVariations) {
  genHist_ = allGen;
  genSplitHist_ = splitGen;
  genJackknifeVariations_ = jackknifeVariations;
  genPDFVariations_ = PDFVariations;
}

void LambdaHistsFiller::assignResponseHists(TH2D * allResponse,
                                            TH2D * splitResponse,
                                            std::vector<TH2D*> * jackknifeVariations,
                                            std::vector<TH2D*> * PDFVariations) {
  responseHist_ = allResponse;
  responseSplitHist_ = splitResponse;
  responseJackknifeVariations_ = jackknifeVariations;
  responsePDFVariations_ = PDFVariations;
}

void LambdaHistsFiller::fillRecoTH1(TH1 * h, double weight) {
  if (h != nullptr && passReco()) {
    h->Fill(recoBin_, weight);
  }
}

void LambdaHistsFiller::fillRecoTH1GenBinning(TH1 * h, double weight) {
  if (h != nullptr && passReco()) {
    h->Fill(recoBinGenBinning_, weight);
  }
}

void LambdaHistsFiller::fillGenTH1(TH1 * h, double weight) {
  if (h != nullptr && passGen()) {
    h->Fill(genBin_, weight);
  }
}

void LambdaHistsFiller::fillResponseTH2(TH2 * h, double recoWeight, double genWeight) {
  // How to handle the fact that the matching reco jet may not be the one stored?
  // Assume the user has re-called setupReco with the matched reco jet
  if (h != nullptr && passGen()) {
    // fill if we pass gen requirements, irrespective of reco requirements
    // (will be bin 0 if fail)
    // fill twice: once normally with total event weight,
    // then again but in global underflow with (gen weight - total event weight).
    // This ensures that the projection on the gen axis matches the 1D gen distribution
    // (which would only have the gen weight)
    double totalWeight = recoWeight * genWeight;
    h->Fill(genBin_, recoBin_, totalWeight);
    double corrWeight = genWeight - totalWeight;
    int underflowBin = 0;
    // this is the underflow bin, since the histogram's bin edges are the global bin numbers,
    // not physical values
    h->Fill(genBin_, underflowBin, corrWeight);
  }
}

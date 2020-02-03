#include "UHH2/QGAnalysis/include/QGAnalysisUnfoldHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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
  useBinningValue_(false), // use the centrally determined bining value (e.g. dijet ave), otherwise use jet pT
  N_PDF_VARIATIONS(100)
  {

  is_mc_ = ctx.get("dataset_type") == "MC";
  doMCsplit_ = is_mc_;
  doPDFvariations_ = (string2bool(ctx.get("PDFvariations", "false")) && is_mc_);
  if (doPDFvariations_) {
    cout << "Doing PDF variations" << endl;
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
  h_fake_counter_raw = book<TH1F>("fakes_counter_raw", "Fakes counter unweighted", descriptions.size()+1, -1, descriptions.size());
  h_fake_counter_weighted = book<TH1F>("fakes_counter", "Fakes counter using weights", descriptions.size()+1, -1, descriptions.size());
  h_fake_counter_charged_raw = book<TH1F>("fakes_counter_charged_raw", "Fakes counter for charged lambda unweighted", descriptions.size()+1, -1, descriptions.size());
  h_fake_counter_charged_weighted = book<TH1F>("fakes_countercharged_", "Fakes counter for charged lambda using weights", descriptions.size()+1, -1, descriptions.size());
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
  detector_tu_binning_pt = new TUnfoldBinning("detector");
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
  h_tu_response_pt = copy_book_th2f(h_tu_response_pt_tmp, "tu_pt_GenReco_all");
  h_tu_response_pt_split = copy_book_th2f(h_tu_response_pt_tmp, "tu_pt_GenReco_split");
  delete h_tu_response_pt_tmp;

  // detector histograms
  TH1 * h_tu_reco_pt_tmp = detector_tu_binning_pt->CreateHistogram("hist_pt_reco");
  h_tu_reco_pt = copy_book_th1f(h_tu_reco_pt_tmp, "hist_pt_reco_all");
  h_tu_reco_pt_split = copy_book_th1f(h_tu_reco_pt_tmp, "hist_pt_reco_split");
  // for fakes, detector binning
  h_tu_reco_pt_fake = copy_book_th1f(h_tu_reco_pt_tmp, "hist_pt_reco_fake_all");
  h_tu_reco_pt_fake_split = copy_book_th1f(h_tu_reco_pt_tmp, "hist_pt_reco_fake_split");
  delete h_tu_reco_pt_tmp;

  // truth histograms
  TH1 * h_tu_gen_pt_tmp = generator_tu_binning_pt->CreateHistogram("hist_pt_truth");
  h_tu_gen_pt = copy_book_th1f(h_tu_gen_pt_tmp, "hist_pt_truth_all");
  h_tu_gen_pt_split = copy_book_th1f(h_tu_gen_pt_tmp, "hist_pt_truth_split");
  delete h_tu_gen_pt_tmp;

  // LHA
  // -------------------------------------
  detector_tu_binning_LHA = new TUnfoldBinning("detector");
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
  h_tu_response_LHA = copy_book_th2f(h_tu_response_LHA_tmp, "tu_LHA_GenReco_all");
  h_tu_response_LHA_split = copy_book_th2f(h_tu_response_LHA_tmp, "tu_LHA_GenReco_split");
  delete h_tu_response_LHA_tmp;

  // detector histograms
  TH1 * h_tu_reco_LHA_tmp = detector_tu_binning_LHA->CreateHistogram("hist_LHA_reco");
  h_tu_reco_LHA = copy_book_th1f(h_tu_reco_LHA_tmp, "hist_LHA_reco_all");
  h_tu_reco_LHA_split = copy_book_th1f(h_tu_reco_LHA_tmp, "hist_LHA_reco_split");
  // for fakes, detector binning
  h_tu_reco_LHA_fake = copy_book_th1f(h_tu_reco_LHA_tmp, "hist_LHA_reco_fake_all");
  h_tu_reco_LHA_fake_split = copy_book_th1f(h_tu_reco_LHA_tmp, "hist_LHA_reco_fake_split");
  delete h_tu_reco_LHA_tmp;

  // truth histograms
  TH1 * h_tu_gen_LHA_tmp = generator_tu_binning_LHA->CreateHistogram("hist_LHA_truth");
  h_tu_gen_LHA = copy_book_th1f(h_tu_gen_LHA_tmp, "hist_LHA_truth_all");
  h_tu_gen_LHA_split = copy_book_th1f(h_tu_gen_LHA_tmp, "hist_LHA_truth_split");
  delete h_tu_gen_LHA_tmp;

  // PDF variations
  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_LHA_PDF_variations.push_back(copy_book_th1f(h_tu_reco_LHA, TString::Format("hist_LHA_reco_all_PDF_%d", i).Data()));
      h_tu_gen_LHA_PDF_variations.push_back(copy_book_th1f(h_tu_gen_LHA, TString::Format("hist_LHA_gen_all_PDF_%d", i).Data()));
    }
  }

  // detector variable, but using gen binning for comparison later
  h_tu_reco_LHA_gen_binning = copy_book_th1f(h_tu_gen_LHA, "hist_LHA_reco_gen_binning");
  h_tu_reco_LHA_gen_binning_split = copy_book_th1f(h_tu_gen_LHA, "hist_LHA_reco_gen_binning_split");
  // for fakes, gen binning
  h_tu_reco_LHA_fake_gen_binning = copy_book_th1f(h_tu_gen_LHA, "hist_LHA_reco_fake_gen_binning");
  h_tu_reco_LHA_fake_gen_binning_split = copy_book_th1f(h_tu_gen_LHA, "hist_LHA_reco_fake_gen_binning_split");


  // Charged LHA
  // -------------------------------------
  detector_tu_binning_LHA_charged = new TUnfoldBinning("detector");
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
  h_tu_response_LHA_charged = copy_book_th2f(h_tu_response_LHA_charged_tmp, "tu_LHA_charged_GenReco_all");
  h_tu_response_LHA_charged_split = copy_book_th2f(h_tu_response_LHA_charged_tmp, "tu_LHA_charged_GenReco_split");
  delete h_tu_response_LHA_charged_tmp;

  TH1 * h_tu_reco_LHA_charged_tmp = detector_tu_binning_LHA_charged->CreateHistogram("hist_LHA_charged_reco");
  h_tu_reco_LHA_charged = copy_book_th1f(h_tu_reco_LHA_charged_tmp, "hist_LHA_charged_reco_all");
  h_tu_reco_LHA_charged_split = copy_book_th1f(h_tu_reco_LHA_charged_tmp, "hist_LHA_charged_reco_split");
  h_tu_reco_LHA_charged_fake = copy_book_th1f(h_tu_reco_LHA_charged_tmp, "hist_LHA_charged_reco_fake_all");
  h_tu_reco_LHA_charged_fake_split = copy_book_th1f(h_tu_reco_LHA_charged_tmp, "hist_LHA_charged_reco_fake_split");
  delete h_tu_reco_LHA_charged_tmp;

  TH1 * h_tu_gen_LHA_charged_tmp = generator_tu_binning_LHA_charged->CreateHistogram("hist_LHA_charged_truth");
  h_tu_gen_LHA_charged = copy_book_th1f(h_tu_gen_LHA_charged_tmp, "hist_LHA_charged_truth_all");
  h_tu_gen_LHA_charged_split = copy_book_th1f(h_tu_gen_LHA_charged_tmp, "hist_LHA_charged_truth_split");
  delete h_tu_gen_LHA_charged_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_LHA_charged_PDF_variations.push_back(copy_book_th1f(h_tu_reco_LHA_charged, TString::Format("hist_LHA_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_LHA_charged_PDF_variations.push_back(copy_book_th1f(h_tu_gen_LHA_charged, TString::Format("hist_LHA_charged_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_LHA_charged_gen_binning = copy_book_th1f(h_tu_gen_LHA_charged, "hist_LHA_charged_reco_gen_binning");
  h_tu_reco_LHA_charged_gen_binning_split = copy_book_th1f(h_tu_gen_LHA_charged, "hist_LHA_charged_reco_gen_binning_split");
  h_tu_reco_LHA_charged_fake_gen_binning = copy_book_th1f(h_tu_gen_LHA_charged, "hist_LHA_charged_reco_fake_gen_binning");
  h_tu_reco_LHA_charged_fake_gen_binning_split = copy_book_th1f(h_tu_gen_LHA_charged, "hist_LHA_charged_reco_fake_gen_binning_split");

  // puppi multiplicity
  // -------------------------------------
  detector_tu_binning_puppiMultiplicity = new TUnfoldBinning("detector");
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
  h_tu_response_puppiMultiplicity = copy_book_th2f(h_tu_response_puppiMultiplicity_tmp, "tu_puppiMultiplicity_GenReco_all");
  h_tu_response_puppiMultiplicity_split = copy_book_th2f(h_tu_response_puppiMultiplicity_tmp, "tu_puppiMultiplicity_GenReco_split");
  delete h_tu_response_puppiMultiplicity_tmp;

  TH1 * h_tu_reco_puppiMultiplicity_tmp = detector_tu_binning_puppiMultiplicity->CreateHistogram("hist_puppiMultiplicity_reco");
  h_tu_reco_puppiMultiplicity = copy_book_th1f(h_tu_reco_puppiMultiplicity_tmp, "hist_puppiMultiplicity_reco_all");
  h_tu_reco_puppiMultiplicity_split = copy_book_th1f(h_tu_reco_puppiMultiplicity_tmp, "hist_puppiMultiplicity_reco_split");
  h_tu_reco_puppiMultiplicity_fake = copy_book_th1f(h_tu_reco_puppiMultiplicity_tmp, "hist_puppiMultiplicity_reco_fake_all");
  h_tu_reco_puppiMultiplicity_fake_split = copy_book_th1f(h_tu_reco_puppiMultiplicity_tmp, "hist_puppiMultiplicity_reco_fake_split");
  delete h_tu_reco_puppiMultiplicity_tmp;

  TH1 * h_tu_gen_puppiMultiplicity_tmp = generator_tu_binning_puppiMultiplicity->CreateHistogram("hist_puppiMultiplicity_truth");
  h_tu_gen_puppiMultiplicity = copy_book_th1f(h_tu_gen_puppiMultiplicity_tmp, "hist_puppiMultiplicity_truth_all");
  h_tu_gen_puppiMultiplicity_split = copy_book_th1f(h_tu_gen_puppiMultiplicity_tmp, "hist_puppiMultiplicity_truth_split");
  delete h_tu_gen_puppiMultiplicity_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_puppiMultiplicity_PDF_variations.push_back(copy_book_th1f(h_tu_reco_puppiMultiplicity, TString::Format("hist_puppiMultiplicity_reco_all_PDF_%d", i).Data()));
      h_tu_gen_puppiMultiplicity_PDF_variations.push_back(copy_book_th1f(h_tu_gen_puppiMultiplicity, TString::Format("hist_puppiMultiplicity_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_puppiMultiplicity_gen_binning = copy_book_th1f(h_tu_gen_puppiMultiplicity, "hist_puppiMultiplicity_reco_gen_binning");
  h_tu_reco_puppiMultiplicity_gen_binning_split = copy_book_th1f(h_tu_gen_puppiMultiplicity, "hist_puppiMultiplicity_reco_gen_binning_split");

  h_tu_reco_puppiMultiplicity_fake_gen_binning = copy_book_th1f(h_tu_gen_puppiMultiplicity, "hist_puppiMultiplicity_reco_fake_gen_binning");
  h_tu_reco_puppiMultiplicity_fake_gen_binning_split = copy_book_th1f(h_tu_gen_puppiMultiplicity, "hist_puppiMultiplicity_reco_fake_gen_binning_split");

  // Charged PUPPI multiplicity
  // -------------------------------------
  detector_tu_binning_puppiMultiplicity_charged = new TUnfoldBinning("detector");
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
  h_tu_response_puppiMultiplicity_charged = copy_book_th2f(h_tu_response_puppiMultiplicity_charged_tmp, "tu_puppiMultiplicity_charged_GenReco_all");
  h_tu_response_puppiMultiplicity_charged_split = copy_book_th2f(h_tu_response_puppiMultiplicity_charged_tmp, "tu_puppiMultiplicity_charged_GenReco_split");
  delete h_tu_response_puppiMultiplicity_charged_tmp;

  TH1 * h_tu_reco_puppiMultiplicity_charged_tmp = detector_tu_binning_puppiMultiplicity_charged->CreateHistogram("hist_puppiMultiplicity_charged_reco");
  h_tu_reco_puppiMultiplicity_charged = copy_book_th1f(h_tu_reco_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_reco_all");
  h_tu_reco_puppiMultiplicity_charged_split = copy_book_th1f(h_tu_reco_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_reco_split");
  h_tu_reco_puppiMultiplicity_charged_fake = copy_book_th1f(h_tu_reco_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_reco_fake_all");
  h_tu_reco_puppiMultiplicity_charged_fake_split = copy_book_th1f(h_tu_reco_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_reco_fake_split");
  delete h_tu_reco_puppiMultiplicity_charged_tmp;

  TH1 * h_tu_gen_puppiMultiplicity_charged_tmp = generator_tu_binning_puppiMultiplicity_charged->CreateHistogram("hist_puppiMultiplicity_charged_truth");
  h_tu_gen_puppiMultiplicity_charged = copy_book_th1f(h_tu_gen_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_truth_all");
  h_tu_gen_puppiMultiplicity_charged_split = copy_book_th1f(h_tu_gen_puppiMultiplicity_charged_tmp, "hist_puppiMultiplicity_charged_truth_split");
  delete h_tu_gen_puppiMultiplicity_charged_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_puppiMultiplicity_charged_PDF_variations.push_back(copy_book_th1f(h_tu_reco_puppiMultiplicity_charged, TString::Format("hist_puppiMultiplicity_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_puppiMultiplicity_charged_PDF_variations.push_back(copy_book_th1f(h_tu_gen_puppiMultiplicity_charged, TString::Format("hist_puppiMultiplicity_charged_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_puppiMultiplicity_charged_gen_binning = copy_book_th1f(h_tu_gen_puppiMultiplicity_charged, "hist_puppiMultiplicity_charged_reco_gen_binning");
  h_tu_reco_puppiMultiplicity_charged_gen_binning_split = copy_book_th1f(h_tu_gen_puppiMultiplicity_charged, "hist_puppiMultiplicity_charged_reco_gen_binning_split");

  h_tu_reco_puppiMultiplicity_charged_fake_gen_binning = copy_book_th1f(h_tu_gen_puppiMultiplicity_charged, "hist_puppiMultiplicity_charged_reco_fake_gen_binning");
  h_tu_reco_puppiMultiplicity_charged_fake_gen_binning_split = copy_book_th1f(h_tu_gen_puppiMultiplicity_charged, "hist_puppiMultiplicity_charged_reco_fake_gen_binning_split");

  // pTD
  // -------------------------------------
  detector_tu_binning_pTD = new TUnfoldBinning("detector");
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
  h_tu_response_pTD = copy_book_th2f(h_tu_response_pTD_tmp, "tu_pTD_GenReco_all");
  h_tu_response_pTD_split = copy_book_th2f(h_tu_response_pTD_tmp, "tu_pTD_GenReco_split");
  delete h_tu_response_pTD_tmp;

  TH1 * h_tu_reco_pTD_tmp = detector_tu_binning_pTD->CreateHistogram("hist_pTD_reco");
  h_tu_reco_pTD = copy_book_th1f(h_tu_reco_pTD_tmp, "hist_pTD_reco_all");
  h_tu_reco_pTD_split = copy_book_th1f(h_tu_reco_pTD_tmp, "hist_pTD_reco_split");
  h_tu_reco_pTD_fake = copy_book_th1f(h_tu_reco_pTD_tmp, "hist_pTD_reco_fake_all");
  h_tu_reco_pTD_fake_split = copy_book_th1f(h_tu_reco_pTD_tmp, "hist_pTD_reco_fake_split");
  delete h_tu_reco_pTD_tmp;

  TH1 * h_tu_gen_pTD_tmp = generator_tu_binning_pTD->CreateHistogram("hist_pTD_truth");
  h_tu_gen_pTD = copy_book_th1f(h_tu_gen_pTD_tmp, "hist_pTD_truth_all");
  h_tu_gen_pTD_split = copy_book_th1f(h_tu_gen_pTD_tmp, "hist_pTD_truth_split");
  delete h_tu_gen_pTD_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_pTD_PDF_variations.push_back(copy_book_th1f(h_tu_reco_pTD, TString::Format("hist_pTD_reco_all_PDF_%d", i).Data()));
      h_tu_gen_pTD_PDF_variations.push_back(copy_book_th1f(h_tu_gen_pTD, TString::Format("hist_pTD_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_pTD_gen_binning = copy_book_th1f(h_tu_gen_pTD, "hist_pTD_reco_gen_binning");
  h_tu_reco_pTD_gen_binning_split = copy_book_th1f(h_tu_gen_pTD, "hist_pTD_reco_gen_binning_split");

  h_tu_reco_pTD_fake_gen_binning = copy_book_th1f(h_tu_gen_pTD, "hist_pTD_reco_fake_gen_binning");
  h_tu_reco_pTD_fake_gen_binning_split = copy_book_th1f(h_tu_gen_pTD, "hist_pTD_reco_fake_gen_binning_split");

  // Charged pTD
  // -------------------------------------
  detector_tu_binning_pTD_charged = new TUnfoldBinning("detector");
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
  h_tu_response_pTD_charged = copy_book_th2f(h_tu_response_pTD_charged_tmp, "tu_pTD_charged_GenReco_all");
  h_tu_response_pTD_charged_split = copy_book_th2f(h_tu_response_pTD_charged_tmp, "tu_pTD_charged_GenReco_split");
  delete h_tu_response_pTD_charged_tmp;

  TH1 * h_tu_reco_pTD_charged_tmp = detector_tu_binning_pTD_charged->CreateHistogram("hist_pTD_charged_reco");
  h_tu_reco_pTD_charged = copy_book_th1f(h_tu_reco_pTD_charged_tmp, "hist_pTD_charged_reco_all");
  h_tu_reco_pTD_charged_split = copy_book_th1f(h_tu_reco_pTD_charged_tmp, "hist_pTD_charged_reco_split");
  h_tu_reco_pTD_charged_fake = copy_book_th1f(h_tu_reco_pTD_charged_tmp, "hist_pTD_charged_reco_fake_all");
  h_tu_reco_pTD_charged_fake_split = copy_book_th1f(h_tu_reco_pTD_charged_tmp, "hist_pTD_charged_reco_fake_split");
  delete h_tu_reco_pTD_charged_tmp;

  TH1 * h_tu_gen_pTD_charged_tmp = generator_tu_binning_pTD_charged->CreateHistogram("hist_pTD_charged_truth");
  h_tu_gen_pTD_charged = copy_book_th1f(h_tu_gen_pTD_charged_tmp, "hist_pTD_charged_truth_all");
  h_tu_gen_pTD_charged_split = copy_book_th1f(h_tu_gen_pTD_charged_tmp, "hist_pTD_charged_truth_split");
  delete h_tu_gen_pTD_charged_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_pTD_charged_PDF_variations.push_back(copy_book_th1f(h_tu_reco_pTD_charged, TString::Format("hist_pTD_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_pTD_charged_PDF_variations.push_back(copy_book_th1f(h_tu_gen_pTD_charged, TString::Format("hist_pTD_charged_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_pTD_charged_gen_binning = copy_book_th1f(h_tu_gen_pTD_charged, "hist_pTD_charged_reco_gen_binning");
  h_tu_reco_pTD_charged_gen_binning_split = copy_book_th1f(h_tu_gen_pTD_charged, "hist_pTD_charged_reco_gen_binning_split");

  h_tu_reco_pTD_charged_fake_gen_binning = copy_book_th1f(h_tu_gen_pTD_charged, "hist_pTD_charged_reco_fake_gen_binning");
  h_tu_reco_pTD_charged_fake_gen_binning_split = copy_book_th1f(h_tu_gen_pTD_charged, "hist_pTD_charged_reco_fake_gen_binning_split");

  // thrust
  // -------------------------------------
  detector_tu_binning_thrust = new TUnfoldBinning("detector");
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
  h_tu_response_thrust = copy_book_th2f(h_tu_response_thrust_tmp, "tu_thrust_GenReco_all");
  h_tu_response_thrust_split = copy_book_th2f(h_tu_response_thrust_tmp, "tu_thrust_GenReco_split");
  delete h_tu_response_thrust_tmp;

  TH1 * h_tu_reco_thrust_tmp = detector_tu_binning_thrust->CreateHistogram("hist_thrust_reco");
  h_tu_reco_thrust = copy_book_th1f(h_tu_reco_thrust_tmp, "hist_thrust_reco_all");
  h_tu_reco_thrust_split = copy_book_th1f(h_tu_reco_thrust_tmp, "hist_thrust_reco_split");
  h_tu_reco_thrust_fake = copy_book_th1f(h_tu_reco_thrust_tmp, "hist_thrust_reco_fake_all");
  h_tu_reco_thrust_fake_split = copy_book_th1f(h_tu_reco_thrust_tmp, "hist_thrust_reco_fake_split");
  delete h_tu_reco_thrust_tmp;

  TH1 * h_tu_gen_thrust_tmp = generator_tu_binning_thrust->CreateHistogram("hist_thrust_truth");
  h_tu_gen_thrust = copy_book_th1f(h_tu_gen_thrust_tmp, "hist_thrust_truth_all");
  h_tu_gen_thrust_split = copy_book_th1f(h_tu_gen_thrust_tmp, "hist_thrust_truth_split");
  delete h_tu_gen_thrust_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_thrust_PDF_variations.push_back(copy_book_th1f(h_tu_reco_thrust, TString::Format("hist_thrust_reco_all_PDF_%d", i).Data()));
      h_tu_gen_thrust_PDF_variations.push_back(copy_book_th1f(h_tu_gen_thrust, TString::Format("hist_thrust_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_thrust_gen_binning = copy_book_th1f(h_tu_gen_thrust, "hist_thrust_reco_gen_binning");
  h_tu_reco_thrust_gen_binning_split = copy_book_th1f(h_tu_gen_thrust, "hist_thrust_reco_gen_binning_split");

  h_tu_reco_thrust_fake_gen_binning = copy_book_th1f(h_tu_gen_thrust, "hist_thrust_reco_fake_gen_binning");
  h_tu_reco_thrust_fake_gen_binning_split = copy_book_th1f(h_tu_gen_thrust, "hist_thrust_reco_fake_gen_binning_split");

  // Charged thrust
  // -------------------------------------
  detector_tu_binning_thrust_charged = new TUnfoldBinning("detector");
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
  h_tu_response_thrust_charged = copy_book_th2f(h_tu_response_thrust_charged_tmp, "tu_thrust_charged_GenReco_all");
  h_tu_response_thrust_charged_split = copy_book_th2f(h_tu_response_thrust_charged_tmp, "tu_thrust_charged_GenReco_split");
  delete h_tu_response_thrust_charged_tmp;

  TH1 * h_tu_reco_thrust_charged_tmp = detector_tu_binning_thrust_charged->CreateHistogram("hist_thrust_charged_reco");
  h_tu_reco_thrust_charged = copy_book_th1f(h_tu_reco_thrust_charged_tmp, "hist_thrust_charged_reco_all");
  h_tu_reco_thrust_charged_split = copy_book_th1f(h_tu_reco_thrust_charged_tmp, "hist_thrust_charged_reco_split");
  h_tu_reco_thrust_charged_fake = copy_book_th1f(h_tu_reco_thrust_charged_tmp, "hist_thrust_charged_reco_fake_all");
  h_tu_reco_thrust_charged_fake_split = copy_book_th1f(h_tu_reco_thrust_charged_tmp, "hist_thrust_charged_reco_fake_split");
  delete h_tu_reco_thrust_charged_tmp;

  TH1 * h_tu_gen_thrust_charged_tmp = generator_tu_binning_thrust_charged->CreateHistogram("hist_thrust_charged_truth");
  h_tu_gen_thrust_charged = copy_book_th1f(h_tu_gen_thrust_charged_tmp, "hist_thrust_charged_truth_all");
  h_tu_gen_thrust_charged_split = copy_book_th1f(h_tu_gen_thrust_charged_tmp, "hist_thrust_charged_truth_split");
  delete h_tu_gen_thrust_charged_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_thrust_charged_PDF_variations.push_back(copy_book_th1f(h_tu_reco_thrust_charged, TString::Format("hist_thrust_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_thrust_charged_PDF_variations.push_back(copy_book_th1f(h_tu_gen_thrust_charged, TString::Format("hist_thrust_charged_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_thrust_charged_gen_binning = copy_book_th1f(h_tu_gen_thrust_charged, "hist_thrust_charged_reco_gen_binning");
  h_tu_reco_thrust_charged_gen_binning_split = copy_book_th1f(h_tu_gen_thrust_charged, "hist_thrust_charged_reco_gen_binning_split");

  h_tu_reco_thrust_charged_fake_gen_binning = copy_book_th1f(h_tu_gen_thrust_charged, "hist_thrust_charged_reco_fake_gen_binning");
  h_tu_reco_thrust_charged_fake_gen_binning_split = copy_book_th1f(h_tu_gen_thrust_charged, "hist_thrust_charged_reco_fake_gen_binning_split");

  // width
  // -------------------------------------
  detector_tu_binning_width = new TUnfoldBinning("detector");
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
  h_tu_response_width = copy_book_th2f(h_tu_response_width_tmp, "tu_width_GenReco_all");
  h_tu_response_width_split = copy_book_th2f(h_tu_response_width_tmp, "tu_width_GenReco_split");
  delete h_tu_response_width_tmp;

  TH1 * h_tu_reco_width_tmp = detector_tu_binning_width->CreateHistogram("hist_width_reco");
  h_tu_reco_width = copy_book_th1f(h_tu_reco_width_tmp, "hist_width_reco_all");
  h_tu_reco_width_split = copy_book_th1f(h_tu_reco_width_tmp, "hist_width_reco_split");
  h_tu_reco_width_fake = copy_book_th1f(h_tu_reco_width_tmp, "hist_width_reco_fake_all");
  h_tu_reco_width_fake_split = copy_book_th1f(h_tu_reco_width_tmp, "hist_width_reco_fake_split");
  delete h_tu_reco_width_tmp;

  TH1 * h_tu_gen_width_tmp = generator_tu_binning_width->CreateHistogram("hist_width_truth");
  h_tu_gen_width = copy_book_th1f(h_tu_gen_width_tmp, "hist_width_truth_all");
  h_tu_gen_width_split = copy_book_th1f(h_tu_gen_width_tmp, "hist_width_truth_split");
  delete h_tu_gen_width_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_width_PDF_variations.push_back(copy_book_th1f(h_tu_reco_width, TString::Format("hist_width_reco_all_PDF_%d", i).Data()));
      h_tu_gen_width_PDF_variations.push_back(copy_book_th1f(h_tu_gen_width, TString::Format("hist_width_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_width_gen_binning = copy_book_th1f(h_tu_gen_width, "hist_width_reco_gen_binning");
  h_tu_reco_width_gen_binning_split = copy_book_th1f(h_tu_gen_width, "hist_width_reco_gen_binning_split");

  h_tu_reco_width_fake_gen_binning = copy_book_th1f(h_tu_gen_width, "hist_width_reco_fake_gen_binning");
  h_tu_reco_width_fake_gen_binning_split = copy_book_th1f(h_tu_gen_width, "hist_width_reco_fake_gen_binning_split");

  // Charged width
  // -------------------------------------
  detector_tu_binning_width_charged = new TUnfoldBinning("detector");
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
  h_tu_response_width_charged = copy_book_th2f(h_tu_response_width_charged_tmp, "tu_width_charged_GenReco_all");
  h_tu_response_width_charged_split = copy_book_th2f(h_tu_response_width_charged_tmp, "tu_width_charged_GenReco_split");
  delete h_tu_response_width_charged_tmp;

  TH1 * h_tu_reco_width_charged_tmp = detector_tu_binning_width_charged->CreateHistogram("hist_width_charged_reco");
  h_tu_reco_width_charged = copy_book_th1f(h_tu_reco_width_charged_tmp, "hist_width_charged_reco_all");
  h_tu_reco_width_charged_split = copy_book_th1f(h_tu_reco_width_charged_tmp, "hist_width_charged_reco_split");
  h_tu_reco_width_charged_fake = copy_book_th1f(h_tu_reco_width_charged_tmp, "hist_width_charged_reco_fake_all");
  h_tu_reco_width_charged_fake_split = copy_book_th1f(h_tu_reco_width_charged_tmp, "hist_width_charged_reco_fake_split");
  delete h_tu_reco_width_charged_tmp;

  TH1 * h_tu_gen_width_charged_tmp = generator_tu_binning_width_charged->CreateHistogram("hist_width_charged_truth");
  h_tu_gen_width_charged = copy_book_th1f(h_tu_gen_width_charged_tmp, "hist_width_charged_truth_all");
  h_tu_gen_width_charged_split = copy_book_th1f(h_tu_gen_width_charged_tmp, "hist_width_charged_truth_split");
  delete h_tu_gen_width_charged_tmp;

  if (doPDFvariations_) {
    for (int i=0; i < N_PDF_VARIATIONS; i++) {
      h_tu_reco_width_charged_PDF_variations.push_back(copy_book_th1f(h_tu_reco_width_charged, TString::Format("hist_width_charged_reco_all_PDF_%d", i).Data()));
      h_tu_gen_width_charged_PDF_variations.push_back(copy_book_th1f(h_tu_gen_width_charged, TString::Format("hist_width_charged_gen_all_PDF_%d", i).Data()));
    }
  }

  h_tu_reco_width_charged_gen_binning = copy_book_th1f(h_tu_gen_width_charged, "hist_width_charged_reco_gen_binning");
  h_tu_reco_width_charged_gen_binning_split = copy_book_th1f(h_tu_gen_width_charged, "hist_width_charged_reco_gen_binning_split");

  h_tu_reco_width_charged_fake_gen_binning = copy_book_th1f(h_tu_gen_width_charged, "hist_width_charged_reco_fake_gen_binning");
  h_tu_reco_width_charged_fake_gen_binning_split = copy_book_th1f(h_tu_gen_width_charged, "hist_width_charged_reco_fake_gen_binning_split");

  if (is_mc_) {
    genJetsLambda_handle = ctx.get_handle< std::vector<GenJetLambdaBundle> > (gen_jetlambda_handle_name);
    pass_gen_handle = ctx.get_handle<bool> (gen_sel_handle_name);
  }

  jetsLambda_handle = ctx.get_handle< std::vector<JetLambdaBundle> > (reco_jetlambda_handle_name);

  pass_reco_handle = ctx.get_handle<bool> (reco_sel_handle_name);
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

  const std::vector<JetLambdaBundle> jetLambdas = event.get(jetsLambda_handle);

  // cout << "***" << event.event << endl;

  // Random event flag for MC only filling response matrix
  // This allows us to split the MC into 2 separate samples for testing
  // Make 80% go into response hist so good stats
  bool onlyFillResponse = (rand_.Rndm() > 0.2);

  // Fill reco jet 1D hists
  // ---------------------------------------------------------------------------
  // At this point, all jet filtering etc should have already been performed
  // Fill ignoring if there is a genjet match or not
  if (passReco) {
    for (int i = 0; i < useNJets_; i++) {
      const Jet & thisjet = jetLambdas.at(i).jet;

      LambdaCalculator<PFParticle> recoJetCalc = jetLambdas.at(i).getLambdaCalculator(false, doGroomed_); // can't be const as getLambda() modifies it

      // To account for Lambda Caclulators with 1 constituent from charged-only or grooming,
      // (which isn't really a jet) we test per jet, and treat it otherwise as a fail
      // TODO: should this be done outside?
      bool thisPassReco = (passReco && recoJetCalc.constits().size() > 1);

      // FIXME check this corresponds to same jet as normal lambdas?
      LambdaCalculator<PFParticle> recoJetCalcCharged = jetLambdas.at(i).getLambdaCalculator(true, doGroomed_);
      bool thisPassRecoCharged = (passReco && recoJetCalcCharged.constits().size() > 1);
      if (!thisPassReco && !thisPassRecoCharged) continue;

      float lha(0.), puppiMult(0.), ptd(0.), width(0.), thrust(0.);
      if (thisPassReco) {
        lha = recoJetCalc.getLambda(1, 0.5);
        puppiMult = recoJetCalc.getLambda(0, 0);
        ptd = recoJetCalc.getLambda(2, 0);
        width = recoJetCalc.getLambda(1, 1);
        thrust = recoJetCalc.getLambda(1, 2);
      }

      float lha_charged(0), puppiMult_charged(0), ptd_charged(0), width_charged(0), thrust_charged(0);
      if (thisPassRecoCharged) {
        lha_charged = recoJetCalcCharged.getLambda(1, 0.5);
        puppiMult_charged = recoJetCalcCharged.getLambda(0, 0);
        ptd_charged = recoJetCalcCharged.getLambda(2, 0);
        width_charged = recoJetCalcCharged.getLambda(1, 1);
        thrust_charged = recoJetCalcCharged.getLambda(1, 2);
      }

      float jet_pt = useBinningValue_ ? event.get(pt_binning_reco_handle) : thisjet.pt();

      // default to 0 for underflow
      int recBinPt(0);
      int recBinLHA(0), recBinPuppiMult(0), recBinpTD(0), recBinThrust(0), recBinWidth(0);
      int recBinLHACharged(0), recBinPuppiMultCharged(0), recBinpTDCharged(0), recBinThrustCharged(0), recBinWidthCharged(0);
      // generator-bin for reco-level quantities (for displaying later)
      int genBinLHA(0), genBinPuppiMult(0), genBinpTD(0), genBinThrust(0), genBinWidth(0);
      int genBinLHACharged(0), genBinPuppiMultCharged(0), genBinpTDCharged(0), genBinThrustCharged(0), genBinWidthCharged(0);
      bool isUnderflow = (jet_pt < Binning::pt_bin_edges_reco[0]);
      // here we get bin number based on if underflow or not
      if (isUnderflow) {
        recBinPt = detector_distribution_underflow_pt->GetGlobalBinNumber(jet_pt);

        if (thisPassReco) {
          recBinLHA = detector_distribution_underflow_LHA->GetGlobalBinNumber(lha, jet_pt);
          recBinPuppiMult = detector_distribution_underflow_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
          recBinpTD = detector_distribution_underflow_pTD->GetGlobalBinNumber(ptd, jet_pt);
          recBinThrust = detector_distribution_underflow_thrust->GetGlobalBinNumber(thrust, jet_pt);
          recBinWidth = detector_distribution_underflow_width->GetGlobalBinNumber(width, jet_pt);

          genBinLHA = generator_distribution_underflow_LHA->GetGlobalBinNumber(lha, jet_pt);
          genBinPuppiMult = generator_distribution_underflow_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
          genBinpTD = generator_distribution_underflow_pTD->GetGlobalBinNumber(ptd, jet_pt);
          genBinThrust = generator_distribution_underflow_thrust->GetGlobalBinNumber(thrust, jet_pt);
          genBinWidth = generator_distribution_underflow_width->GetGlobalBinNumber(width, jet_pt);
        }

        if (thisPassRecoCharged) {
          recBinLHACharged = detector_distribution_underflow_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
          recBinPuppiMultCharged = detector_distribution_underflow_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
          recBinpTDCharged = detector_distribution_underflow_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
          recBinThrustCharged = detector_distribution_underflow_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
          recBinWidthCharged = detector_distribution_underflow_width_charged->GetGlobalBinNumber(width_charged, jet_pt);

          genBinLHACharged = generator_distribution_underflow_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
          genBinPuppiMultCharged = generator_distribution_underflow_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
          genBinpTDCharged = generator_distribution_underflow_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
          genBinThrustCharged = generator_distribution_underflow_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
          genBinWidthCharged = generator_distribution_underflow_width_charged->GetGlobalBinNumber(width_charged, jet_pt);
        }
      } else {
        recBinPt = detector_distribution_pt->GetGlobalBinNumber(jet_pt);

        if (thisPassReco) {
          recBinLHA = detector_distribution_LHA->GetGlobalBinNumber(lha, jet_pt);
          recBinPuppiMult = detector_distribution_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
          recBinpTD = detector_distribution_pTD->GetGlobalBinNumber(ptd, jet_pt);
          recBinThrust = detector_distribution_thrust->GetGlobalBinNumber(thrust, jet_pt);
          recBinWidth = detector_distribution_width->GetGlobalBinNumber(width, jet_pt);

          genBinLHA = generator_distribution_LHA->GetGlobalBinNumber(lha, jet_pt);
          genBinPuppiMult = generator_distribution_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
          genBinpTD = generator_distribution_pTD->GetGlobalBinNumber(ptd, jet_pt);
          genBinThrust = generator_distribution_thrust->GetGlobalBinNumber(thrust, jet_pt);
          genBinWidth = generator_distribution_width->GetGlobalBinNumber(width, jet_pt);
        }

        if (thisPassRecoCharged) {
          recBinLHACharged = detector_distribution_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
          recBinPuppiMultCharged = detector_distribution_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
          recBinpTDCharged = detector_distribution_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
          recBinThrustCharged = detector_distribution_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
          recBinWidthCharged = detector_distribution_width_charged->GetGlobalBinNumber(width_charged, jet_pt);

          genBinLHACharged = generator_distribution_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
          genBinPuppiMultCharged = generator_distribution_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
          genBinpTDCharged = generator_distribution_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
          genBinThrustCharged = generator_distribution_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
          genBinWidthCharged = generator_distribution_width_charged->GetGlobalBinNumber(width_charged, jet_pt);
        }
      }

      h_tu_reco_pt->Fill(recBinPt, weight);

      if (thisPassReco) {
        h_tu_reco_LHA->Fill(recBinLHA, weight);
        h_tu_reco_puppiMultiplicity->Fill(recBinPuppiMult, weight);
        h_tu_reco_pTD->Fill(recBinpTD, weight);
        h_tu_reco_thrust->Fill(recBinThrust, weight);
        h_tu_reco_width->Fill(recBinWidth, weight);

        h_tu_reco_LHA_gen_binning->Fill(genBinLHA, weight);
        h_tu_reco_puppiMultiplicity_gen_binning->Fill(genBinPuppiMult, weight);
        h_tu_reco_pTD_gen_binning->Fill(genBinpTD, weight);
        h_tu_reco_thrust_gen_binning->Fill(genBinThrust, weight);
        h_tu_reco_width_gen_binning->Fill(genBinWidth, weight);
      }

      if (thisPassRecoCharged) {
        h_tu_reco_LHA_charged->Fill(recBinLHACharged, weight);
        h_tu_reco_puppiMultiplicity_charged->Fill(recBinPuppiMultCharged, weight);
        h_tu_reco_pTD_charged->Fill(recBinpTDCharged, weight);
        h_tu_reco_thrust_charged->Fill(recBinThrustCharged, weight);
        h_tu_reco_width_charged->Fill(recBinWidthCharged, weight);

        h_tu_reco_LHA_charged_gen_binning->Fill(genBinLHACharged, weight);
        h_tu_reco_puppiMultiplicity_charged_gen_binning->Fill(genBinPuppiMultCharged, weight);
        h_tu_reco_pTD_charged_gen_binning->Fill(genBinpTDCharged, weight);
        h_tu_reco_thrust_charged_gen_binning->Fill(genBinThrustCharged, weight);
        h_tu_reco_width_charged_gen_binning->Fill(genBinWidthCharged, weight);
      }

      if (!onlyFillResponse && doMCsplit_) {
        h_tu_reco_pt_split->Fill(recBinPt, weight);

        if (thisPassReco) {
          h_tu_reco_LHA_split->Fill(recBinLHA, weight);
          h_tu_reco_puppiMultiplicity_split->Fill(recBinPuppiMult, weight);
          h_tu_reco_pTD_split->Fill(recBinpTD, weight);
          h_tu_reco_thrust_split->Fill(recBinThrust, weight);
          h_tu_reco_width_split->Fill(recBinWidth, weight);

          h_tu_reco_LHA_gen_binning_split->Fill(genBinLHA, weight);
          h_tu_reco_puppiMultiplicity_gen_binning_split->Fill(genBinPuppiMult, weight);
          h_tu_reco_pTD_gen_binning_split->Fill(genBinpTD, weight);
          h_tu_reco_thrust_gen_binning_split->Fill(genBinThrust, weight);
          h_tu_reco_width_gen_binning_split->Fill(genBinWidth, weight);
        }

        if (thisPassRecoCharged) {
          h_tu_reco_LHA_charged_split->Fill(recBinLHACharged, weight);
          h_tu_reco_puppiMultiplicity_charged_split->Fill(recBinPuppiMultCharged, weight);
          h_tu_reco_pTD_charged_split->Fill(recBinpTDCharged, weight);
          h_tu_reco_thrust_charged_split->Fill(recBinThrustCharged, weight);
          h_tu_reco_width_charged_split->Fill(recBinWidthCharged, weight);

          h_tu_reco_LHA_charged_gen_binning_split->Fill(genBinLHACharged, weight);
          h_tu_reco_puppiMultiplicity_charged_gen_binning_split->Fill(genBinPuppiMultCharged, weight);
          h_tu_reco_pTD_charged_gen_binning_split->Fill(genBinpTDCharged, weight);
          h_tu_reco_thrust_charged_gen_binning_split->Fill(genBinThrustCharged, weight);
          h_tu_reco_width_charged_gen_binning_split->Fill(genBinWidthCharged, weight);
        }
      }

      // Fill PDF variations by varying weight
      if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
        for (int i=0; i<N_PDF_VARIATIONS; i++) {
          double pdf_weight = event.genInfo->systweights().at(i+9) / event.genInfo->originalXWGTUP(); // +9 as first 9 are scale variations
          double this_weight = weight * pdf_weight;
          if (thisPassReco) {
            h_tu_reco_LHA_PDF_variations.at(i)->Fill(recBinLHA, this_weight);
            h_tu_reco_puppiMultiplicity_PDF_variations.at(i)->Fill(recBinPuppiMult, this_weight);
            h_tu_reco_pTD_PDF_variations.at(i)->Fill(recBinpTD, this_weight);
            h_tu_reco_thrust_PDF_variations.at(i)->Fill(recBinThrust, this_weight);
            h_tu_reco_width_PDF_variations.at(i)->Fill(recBinWidth, this_weight);
          }
          if (thisPassRecoCharged) {
            h_tu_reco_LHA_charged_PDF_variations.at(i)->Fill(recBinLHACharged, this_weight);
            h_tu_reco_puppiMultiplicity_charged_PDF_variations.at(i)->Fill(recBinPuppiMultCharged, this_weight);
            h_tu_reco_pTD_charged_PDF_variations.at(i)->Fill(recBinpTDCharged, this_weight);
            h_tu_reco_thrust_charged_PDF_variations.at(i)->Fill(recBinThrustCharged, this_weight);
            h_tu_reco_width_charged_PDF_variations.at(i)->Fill(recBinWidthCharged, this_weight);
          }
        }
      }

      if (is_mc_) {
        // Fill fakes, both in detector & truth binnings
        // Fake = not pass Gen, or pass Gen but matching genjet not one of the
        // selection ones (i.e. doesn't count as a match)
        // Or the genjet fails #constits > 1
        // -----------------------------------------------------------------------
        if (thisPassReco) {
          h_fake_counter_raw->Fill(-1);
          h_fake_counter_weighted->Fill(-1, weight);

          bool genConstit = false;
          if ((thisjet.genjet_index() >= 0) && (thisjet.genjet_index() < genjetLambdas->size())) {
            genConstit = (genjetLambdas->at(thisjet.genjet_index()).getLambdaCalculator(false, doGroomed_).constits().size() > 1);
          }

          if (!passGen || (thisjet.genjet_index() < 0 || thisjet.genjet_index() >= useNJets_) || !genConstit) {
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
            else if (!genConstit) {
              int ind=3;
              h_fake_counter_raw->Fill(ind);
              h_fake_counter_weighted->Fill(ind, weight);
            }

            h_tu_reco_pt_fake->Fill(recBinPt, weight);

            h_tu_reco_LHA_fake->Fill(recBinLHA, weight);
            h_tu_reco_puppiMultiplicity_fake->Fill(recBinPuppiMult, weight);
            h_tu_reco_pTD_fake->Fill(recBinpTD, weight);
            h_tu_reco_thrust_fake->Fill(recBinThrust, weight);
            h_tu_reco_width_fake->Fill(recBinWidth, weight);

            h_tu_reco_LHA_fake_gen_binning->Fill(genBinLHA, weight);
            h_tu_reco_puppiMultiplicity_fake_gen_binning->Fill(genBinPuppiMult, weight);
            h_tu_reco_pTD_fake_gen_binning->Fill(genBinpTD, weight);
            h_tu_reco_thrust_fake_gen_binning->Fill(genBinThrust, weight);
            h_tu_reco_width_fake_gen_binning->Fill(genBinWidth, weight);

            if (!onlyFillResponse && doMCsplit_) {
              h_tu_reco_pt_fake_split->Fill(recBinPt, weight);

              h_tu_reco_LHA_fake_split->Fill(recBinLHA, weight);
              h_tu_reco_puppiMultiplicity_fake_split->Fill(recBinPuppiMult, weight);
              h_tu_reco_pTD_fake_split->Fill(recBinpTD, weight);
              h_tu_reco_thrust_fake_split->Fill(recBinThrust, weight);
              h_tu_reco_width_fake_split->Fill(recBinWidth, weight);

              h_tu_reco_LHA_fake_gen_binning_split->Fill(genBinLHA, weight);
              h_tu_reco_puppiMultiplicity_fake_gen_binning_split->Fill(genBinPuppiMult, weight);
              h_tu_reco_pTD_fake_gen_binning_split->Fill(genBinpTD, weight);
              h_tu_reco_thrust_fake_gen_binning_split->Fill(genBinThrust, weight);
              h_tu_reco_width_fake_gen_binning_split->Fill(genBinWidth, weight);
            }
          } // end of test if failed gen
        } // end of if thisPassReco

        if (thisPassRecoCharged) {
          h_fake_counter_charged_raw->Fill(-1);
          h_fake_counter_charged_weighted->Fill(-1, weight);

          bool genConstit = false;
          if ((thisjet.genjet_index() >= 0) && (thisjet.genjet_index() < genjetLambdas->size())) {
            genConstit = (genjetLambdas->at(thisjet.genjet_index()).getLambdaCalculator(true, doGroomed_).constits().size() > 1);
          }
          if (!passGen || (thisjet.genjet_index() < 0 || thisjet.genjet_index() >= useNJets_) || !genConstit) {
            if (!passGen) {
              int ind=0;
              h_fake_counter_charged_raw->Fill(ind);
              h_fake_counter_charged_weighted->Fill(ind, weight);
            }
            else if (thisjet.genjet_index() < 0) {
              int ind=1;
              h_fake_counter_charged_raw->Fill(ind);
              h_fake_counter_charged_weighted->Fill(ind, weight);
            }
            else if (thisjet.genjet_index() >= useNJets_) {
              int ind=2;
              h_fake_counter_charged_raw->Fill(ind);
              h_fake_counter_charged_weighted->Fill(ind, weight);
            }
            else if (!genConstit) {
              int ind=3;
              h_fake_counter_charged_raw->Fill(ind);
              h_fake_counter_charged_weighted->Fill(ind, weight);
            }

            h_tu_reco_LHA_charged_fake->Fill(recBinLHACharged, weight);
            h_tu_reco_puppiMultiplicity_charged_fake->Fill(recBinPuppiMultCharged, weight);
            h_tu_reco_pTD_charged_fake->Fill(recBinpTDCharged, weight);
            h_tu_reco_thrust_charged_fake->Fill(recBinThrustCharged, weight);
            h_tu_reco_width_charged_fake->Fill(recBinWidthCharged, weight);

            h_tu_reco_LHA_charged_fake_gen_binning->Fill(genBinLHACharged, weight);
            h_tu_reco_puppiMultiplicity_charged_fake_gen_binning->Fill(genBinPuppiMultCharged, weight);
            h_tu_reco_pTD_charged_fake_gen_binning->Fill(genBinpTDCharged, weight);
            h_tu_reco_thrust_charged_fake_gen_binning->Fill(genBinThrustCharged, weight);
            h_tu_reco_width_charged_fake_gen_binning->Fill(genBinWidthCharged, weight);

            if (!onlyFillResponse && doMCsplit_) {
              h_tu_reco_LHA_charged_fake_split->Fill(recBinLHACharged, weight);
              h_tu_reco_puppiMultiplicity_charged_fake_split->Fill(recBinPuppiMultCharged, weight);
              h_tu_reco_pTD_charged_fake_split->Fill(recBinpTDCharged, weight);
              h_tu_reco_thrust_charged_fake_split->Fill(recBinThrustCharged, weight);
              h_tu_reco_width_charged_fake_split->Fill(recBinWidthCharged, weight);

              h_tu_reco_LHA_charged_fake_gen_binning_split->Fill(genBinLHACharged, weight);
              h_tu_reco_puppiMultiplicity_charged_fake_gen_binning_split->Fill(genBinPuppiMultCharged, weight);
              h_tu_reco_pTD_charged_fake_gen_binning_split->Fill(genBinpTDCharged, weight);
              h_tu_reco_thrust_charged_fake_gen_binning_split->Fill(genBinThrustCharged, weight);
              h_tu_reco_width_charged_fake_gen_binning_split->Fill(genBinWidthCharged, weight);
            }
          } // end of test if failed gen charged
        } // end if thisPassRecoCharged
      } // end if mc
    } // end recojet loop
  } // end if passReco

  // Fill gen jet 1D hists & response matrices
  // ---------------------------------------------------------------------------
  if (is_mc_ && passGen) {
    for (int i = 0; i < useNJets_; i++) {
      const GenJetWithParts & thisjet = genjetLambdas->at(i).jet;
      LambdaCalculator<GenParticle> genJetCalc = genjetLambdas->at(i).getLambdaCalculator(false, doGroomed_); // can't be const as getLambda() modifies it
      bool thisPassGen = (passGen && genJetCalc.constits().size() > 1);
      // FIXME check this corresponds to same jet as normal lambdas?
      LambdaCalculator<GenParticle> genJetCalcCharged = genjetLambdas->at(i).getLambdaCalculator(true, doGroomed_);
      bool thisPassGenCharged = (passGen && genJetCalcCharged.constits().size() > 1);

      if (!thisPassGen && !thisPassGenCharged) continue;

      float gen_lha(0.), gen_mult(0.), gen_ptd(0.), gen_width(0.), gen_thrust(0.);
      if (thisPassGen) {
        gen_lha = genJetCalc.getLambda(1, 0.5);
        gen_mult = genJetCalc.getLambda(0, 0);
        gen_ptd = genJetCalc.getLambda(2, 0);
        gen_width = genJetCalc.getLambda(1, 1);
        gen_thrust = genJetCalc.getLambda(1, 2);
      }

      float gen_lha_charged(0.), gen_mult_charged(0.), gen_ptd_charged(0.), gen_width_charged(0.), gen_thrust_charged(0.);
      if (thisPassGenCharged) {
        gen_lha_charged = genJetCalcCharged.getLambda(1, 0.5);
        gen_mult_charged = genJetCalcCharged.getLambda(0, 0);
        gen_ptd_charged = genJetCalcCharged.getLambda(2, 0);
        gen_width_charged = genJetCalcCharged.getLambda(1, 1);
        gen_thrust_charged = genJetCalcCharged.getLambda(1, 2);
      }

      float genjet_pt = useBinningValue_ ? event.get(pt_binning_gen_handle) : thisjet.pt();

      // Get TUnfold gen bins
      // default to 0 for underflow
      int genBinPt(0);
      int genBinLHA(0), genBinPuppiMult(0), genBinpTD(0), genBinThrust(0), genBinWidth(0);
      int genBinLHACharged(0), genBinPuppiMultCharged(0), genBinpTDCharged(0), genBinThrustCharged(0), genBinWidthCharged(0);
      bool isUnderflowGen = (genjet_pt < Binning::pt_bin_edges_gen[0]);
      // here we get bin number based on if underflow or not
      if (isUnderflowGen) {
        genBinPt = generator_distribution_underflow_pt->GetGlobalBinNumber(genjet_pt);
        if (thisPassGen) {
          genBinLHA = generator_distribution_underflow_LHA->GetGlobalBinNumber(gen_lha, genjet_pt);
          genBinPuppiMult = generator_distribution_underflow_puppiMultiplicity->GetGlobalBinNumber(gen_mult, genjet_pt);
          genBinpTD = generator_distribution_underflow_pTD->GetGlobalBinNumber(gen_ptd, genjet_pt);
          genBinThrust = generator_distribution_underflow_thrust->GetGlobalBinNumber(gen_thrust, genjet_pt);
          genBinWidth = generator_distribution_underflow_width->GetGlobalBinNumber(gen_width, genjet_pt);
        }
        if (thisPassGenCharged) {
          genBinLHACharged = generator_distribution_underflow_LHA_charged->GetGlobalBinNumber(gen_lha_charged, genjet_pt);
          genBinPuppiMultCharged = generator_distribution_underflow_puppiMultiplicity_charged->GetGlobalBinNumber(gen_mult_charged, genjet_pt);
          genBinpTDCharged = generator_distribution_underflow_pTD_charged->GetGlobalBinNumber(gen_ptd_charged, genjet_pt);
          genBinThrustCharged = generator_distribution_underflow_thrust_charged->GetGlobalBinNumber(gen_thrust_charged, genjet_pt);
          genBinWidthCharged = generator_distribution_underflow_width_charged->GetGlobalBinNumber(gen_width_charged, genjet_pt);
        }
      } else {
        genBinPt = generator_distribution_pt->GetGlobalBinNumber(genjet_pt);
        if (thisPassGen) {
          genBinLHA = generator_distribution_LHA->GetGlobalBinNumber(gen_lha, genjet_pt);
          genBinPuppiMult = generator_distribution_puppiMultiplicity->GetGlobalBinNumber(gen_mult, genjet_pt);
          genBinpTD = generator_distribution_pTD->GetGlobalBinNumber(gen_ptd, genjet_pt);
          genBinThrust = generator_distribution_thrust->GetGlobalBinNumber(gen_thrust, genjet_pt);
          genBinWidth = generator_distribution_width->GetGlobalBinNumber(gen_width, genjet_pt);
        }
        if (thisPassGenCharged) {
          genBinLHACharged = generator_distribution_LHA_charged->GetGlobalBinNumber(gen_lha_charged, genjet_pt);
          genBinPuppiMultCharged = generator_distribution_puppiMultiplicity_charged->GetGlobalBinNumber(gen_mult_charged, genjet_pt);
          genBinpTDCharged = generator_distribution_pTD_charged->GetGlobalBinNumber(gen_ptd_charged, genjet_pt);
          genBinThrustCharged = generator_distribution_thrust_charged->GetGlobalBinNumber(gen_thrust_charged, genjet_pt);
          genBinWidthCharged = generator_distribution_width_charged->GetGlobalBinNumber(gen_width_charged, genjet_pt);
        }
      }

      // Fill 1D gen hist
      // -----------------------------------------------------------------------
      h_tu_gen_pt->Fill(genBinPt, gen_weight);

      if (thisPassGen) {
        h_tu_gen_puppiMultiplicity->Fill(genBinPuppiMult, gen_weight);
        h_tu_gen_LHA->Fill(genBinLHA, gen_weight);
        h_tu_gen_pTD->Fill(genBinpTD, gen_weight);
        h_tu_gen_width->Fill(genBinWidth, gen_weight);
        h_tu_gen_thrust->Fill(genBinThrust, gen_weight);
      }

      if (thisPassGenCharged) {
        h_tu_gen_puppiMultiplicity_charged->Fill(genBinPuppiMultCharged, gen_weight);
        h_tu_gen_LHA_charged->Fill(genBinLHACharged, gen_weight);
        h_tu_gen_pTD_charged->Fill(genBinpTDCharged, gen_weight);
        h_tu_gen_width_charged->Fill(genBinWidthCharged, gen_weight);
        h_tu_gen_thrust_charged->Fill(genBinThrustCharged, gen_weight);
      }

      if (!onlyFillResponse && doMCsplit_) {
        // Fill 1D gen hist
        h_tu_gen_pt_split->Fill(genBinPt, gen_weight);
        if (thisPassGen) {
          h_tu_gen_puppiMultiplicity_split->Fill(genBinPuppiMult, gen_weight);
          h_tu_gen_LHA_split->Fill(genBinLHA, gen_weight);
          h_tu_gen_pTD_split->Fill(genBinpTD, gen_weight);
          h_tu_gen_width_split->Fill(genBinWidth, gen_weight);
          h_tu_gen_thrust_split->Fill(genBinThrust, gen_weight);
        }
        if (thisPassGenCharged) {
          h_tu_gen_puppiMultiplicity_charged_split->Fill(genBinPuppiMultCharged, gen_weight);
          h_tu_gen_LHA_charged_split->Fill(genBinLHACharged, gen_weight);
          h_tu_gen_pTD_charged_split->Fill(genBinpTDCharged, gen_weight);
          h_tu_gen_width_charged_split->Fill(genBinWidthCharged, gen_weight);
          h_tu_gen_thrust_charged_split->Fill(genBinThrustCharged, gen_weight);
        }
      }

      // Fill PDF variations by varying weight
      if (doPDFvariations_ && event.genInfo->systweights().size() > 0) {
        for (int i=0; i<N_PDF_VARIATIONS; i++) {
          double pdf_weight = event.genInfo->systweights().at(i+9) / event.genInfo->originalXWGTUP(); // +9 as first 9 are scale variations
          double this_weight = gen_weight * pdf_weight;
          if (thisPassGen) {
            h_tu_gen_LHA_PDF_variations.at(i)->Fill(genBinLHA, this_weight);
            h_tu_gen_puppiMultiplicity_PDF_variations.at(i)->Fill(genBinPuppiMult, this_weight);
            h_tu_gen_pTD_PDF_variations.at(i)->Fill(genBinpTD, this_weight);
            h_tu_gen_thrust_PDF_variations.at(i)->Fill(genBinThrust, this_weight);
            h_tu_gen_width_PDF_variations.at(i)->Fill(genBinWidth, this_weight);
          }
          if (thisPassGenCharged) {
            h_tu_gen_LHA_charged_PDF_variations.at(i)->Fill(genBinLHACharged, this_weight);
            h_tu_gen_puppiMultiplicity_charged_PDF_variations.at(i)->Fill(genBinPuppiMultCharged, this_weight);
            h_tu_gen_pTD_charged_PDF_variations.at(i)->Fill(genBinpTDCharged, this_weight);
            h_tu_gen_thrust_charged_PDF_variations.at(i)->Fill(genBinThrustCharged, this_weight);
            h_tu_gen_width_charged_PDF_variations.at(i)->Fill(genBinWidthCharged, this_weight);
          }
        }
      }

      // Now fill response matrices
      // -----------------------------------------------------------------------
      // Loop through genjets, since if there isn't a reco jet then it's a miss-reco,
      // whereas a recojet & no genjet is a fake, which we don't want in our migration matrix
      int recBinPt(0);
      int recBinLHA(0), recBinPuppiMult(0), recBinpTD(0), recBinThrust(0), recBinWidth(0);
      int recBinLHACharged(0), recBinPuppiMultCharged(0), recBinpTDCharged(0), recBinThrustCharged(0), recBinWidthCharged(0);
      int recoInd = -1;
      if (passReco) { // only valid match if reco selection also passed
        // First find matching recojet - annoyingly matches are held in the Jet class, not GenJet
        // default to -1 for underflow
        for (int j=0; j<useNJets_; j++) {
          if (jetLambdas.at(j).jet.genjet_index() == i) { recoInd = j; break; }
        }
        if (recoInd < 0) {
          for (uint j=0; j<jetLambdas.size(); j++) {
            if (jetLambdas.at(j).jet.genjet_index() == i) { cout << "Match for genjet " << i << " actually found with reco jet " << j << endl;}
          }
        }
        if (recoInd >= 0) {
          const Jet & thisjet = jetLambdas.at(recoInd).jet;
          LambdaCalculator<PFParticle> recoJetCalc = jetLambdas.at(recoInd).getLambdaCalculator(false, doGroomed_); // can't be const as getLambda() modifies it
          // FIXME check this corresponds to same jet as normal lambdas?
          LambdaCalculator<PFParticle> recoJetCalcCharged = jetLambdas.at(recoInd).getLambdaCalculator(true, doGroomed_);

          thisPassReco = (passReco && recoJetCalc.constits().size() > 1);
          thisPassRecoCharged = (passReco && recoJetCalcCharged.constits().size() > 1);
          if (!thisPassReco && !thisPassRecoCharged) continue;

          float lha(0.), puppiMult(0.), ptd(0.), width(0.), thrust(0.);
          if (thisPassReco) {
            lha = recoJetCalc.getLambda(1, 0.5);
            puppiMult = recoJetCalc.getLambda(0, 0);
            ptd = recoJetCalc.getLambda(2, 0);
            width = recoJetCalc.getLambda(1, 1);
            thrust = recoJetCalc.getLambda(1, 2);
          }

          float lha_charged(0), puppiMult_charged(0), ptd_charged(0), width_charged(0), thrust_charged(0);
          if (thisPassRecoCharged) {
            lha_charged = recoJetCalcCharged.getLambda(1, 0.5);
            puppiMult_charged = recoJetCalcCharged.getLambda(0, 0);
            ptd_charged = recoJetCalcCharged.getLambda(2, 0);
            width_charged = recoJetCalcCharged.getLambda(1, 1);
            thrust_charged = recoJetCalcCharged.getLambda(1, 2);
          }

          float jet_pt = useBinningValue_ ? event.get(pt_binning_reco_handle) : thisjet.pt();

          bool isUnderflow = (jet_pt < Binning::pt_bin_edges_reco[0]);
          // here we get bin number based on if underflow or not
          if (isUnderflow) {
            recBinPt = detector_distribution_underflow_pt->GetGlobalBinNumber(jet_pt);
            if (thisPassReco) {
              recBinLHA = detector_distribution_underflow_LHA->GetGlobalBinNumber(lha, jet_pt);
              recBinPuppiMult = detector_distribution_underflow_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
              recBinpTD = detector_distribution_underflow_pTD->GetGlobalBinNumber(ptd, jet_pt);
              recBinThrust = detector_distribution_underflow_thrust->GetGlobalBinNumber(thrust, jet_pt);
              recBinWidth = detector_distribution_underflow_width->GetGlobalBinNumber(width, jet_pt);
            }
            if (thisPassRecoCharged) {
              recBinLHACharged = detector_distribution_underflow_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
              recBinPuppiMultCharged = detector_distribution_underflow_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
              recBinpTDCharged = detector_distribution_underflow_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
              recBinThrustCharged = detector_distribution_underflow_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
              recBinWidthCharged = detector_distribution_underflow_width_charged->GetGlobalBinNumber(width_charged, jet_pt);
            }
          } else {
            recBinPt = detector_distribution_pt->GetGlobalBinNumber(jet_pt);
            if (thisPassReco) {
              recBinLHA = detector_distribution_LHA->GetGlobalBinNumber(lha, jet_pt);
              recBinPuppiMult = detector_distribution_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
              recBinpTD = detector_distribution_pTD->GetGlobalBinNumber(ptd, jet_pt);
              recBinThrust = detector_distribution_thrust->GetGlobalBinNumber(thrust, jet_pt);
              recBinWidth = detector_distribution_width->GetGlobalBinNumber(width, jet_pt);
            }
            if (thisPassRecoCharged) {
              recBinLHACharged = detector_distribution_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
              recBinPuppiMultCharged = detector_distribution_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
              recBinpTDCharged = detector_distribution_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
              recBinThrustCharged = detector_distribution_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
              recBinWidthCharged = detector_distribution_width_charged->GetGlobalBinNumber(width_charged, jet_pt);
            }
          }
        }
      } // end if passReco

      // Fill TUnfold 2D response maps
      // ---------------------------------------------------------------
      h_tu_response_pt->Fill(genBinPt, recBinPt, weight);
      if (thisPassGen) {
        h_tu_response_puppiMultiplicity->Fill(genBinPuppiMult, recBinPuppiMult, weight);
        h_tu_response_LHA->Fill(genBinLHA, recBinLHA, weight);
        if (genBinLHA == 0 && recBinLHA == 0) {
          cout << "BADBIN event " << event.event << " genBinLHA == 0 && recBinLHA = 0: " << gen_lha << " : " << genjet_pt << " reco: "<< lha << " : " << jet_pt << endl;
        }
        h_tu_response_pTD->Fill(genBinpTD, recBinpTD, weight);
        h_tu_response_width->Fill(genBinWidth, recBinWidth, weight);
        h_tu_response_thrust->Fill(genBinThrust, recBinThrust, weight);
      }

      if (thisPassGenCharged) {
        h_tu_response_puppiMultiplicity_charged->Fill(genBinPuppiMultCharged, recBinPuppiMultCharged, weight);
        h_tu_response_LHA_charged->Fill(genBinLHACharged, recBinLHACharged, weight);
        h_tu_response_pTD_charged->Fill(genBinpTDCharged, recBinpTDCharged, weight);
        h_tu_response_width_charged->Fill(genBinWidthCharged, recBinWidthCharged, weight);
        h_tu_response_thrust_charged->Fill(genBinThrustCharged, recBinThrustCharged, weight);
      }
      // we need to add in an extra part, such that the 1D projection on the gen axis
      // agrees with the 1D gen histogram
      // i.e. account for the difference between the reco_weight (goes into event.weight) and gen_weight
      // we use the underflow bin for this
      // (0 as bin edges are tunfold bin numbers, not physical values, so the first bin is 0-1,
      // not to be confused with ROOT binning, starting at 1 :s)
      double corr_weight = gen_weight * (1 - reco_weight);
      int underflow_bin = 0;
      h_tu_response_pt->Fill(genBinPt, underflow_bin, corr_weight);
      if (thisPassGen) {
        h_tu_response_puppiMultiplicity->Fill(genBinPuppiMult, underflow_bin, corr_weight);
        h_tu_response_LHA->Fill(genBinLHA, underflow_bin, corr_weight);
        h_tu_response_pTD->Fill(genBinpTD, underflow_bin, corr_weight);
        h_tu_response_width->Fill(genBinWidth, underflow_bin, corr_weight);
        h_tu_response_thrust->Fill(genBinThrust, underflow_bin, corr_weight);
      }
      if (thisPassGenCharged) {
        h_tu_response_puppiMultiplicity_charged->Fill(genBinPuppiMultCharged, underflow_bin, corr_weight);
        h_tu_response_LHA_charged->Fill(genBinLHACharged, underflow_bin, corr_weight);
        h_tu_response_pTD_charged->Fill(genBinpTDCharged, underflow_bin, corr_weight);
        h_tu_response_width_charged->Fill(genBinWidthCharged, underflow_bin, corr_weight);
        h_tu_response_thrust_charged->Fill(genBinThrustCharged, underflow_bin, corr_weight);
      }

      if (onlyFillResponse && doMCsplit_) {
        h_tu_response_pt_split->Fill(genBinPt, recBinPt, weight);
        // extra part for correct weighting
        h_tu_response_pt_split->Fill(genBinPt, underflow_bin, corr_weight);

        if (thisPassGen) {
          h_tu_response_puppiMultiplicity_split->Fill(genBinPuppiMult, recBinPuppiMult, weight);
          h_tu_response_LHA_split->Fill(genBinLHA, recBinLHA, weight);
          h_tu_response_pTD_split->Fill(genBinpTD, recBinpTD, weight);
          h_tu_response_width_split->Fill(genBinWidth, recBinWidth, weight);
          h_tu_response_thrust_split->Fill(genBinThrust, recBinThrust, weight);
          // correct weighting part
          h_tu_response_puppiMultiplicity_split->Fill(genBinPuppiMult, underflow_bin, corr_weight);
          h_tu_response_LHA_split->Fill(genBinLHA, underflow_bin, corr_weight);
          h_tu_response_pTD_split->Fill(genBinpTD, underflow_bin, corr_weight);
          h_tu_response_width_split->Fill(genBinWidth, underflow_bin, corr_weight);
          h_tu_response_thrust_split->Fill(genBinThrust, underflow_bin, corr_weight);
        }
        if (thisPassGenCharged) {
          h_tu_response_puppiMultiplicity_charged_split->Fill(genBinPuppiMultCharged, recBinPuppiMultCharged, weight);
          h_tu_response_LHA_charged_split->Fill(genBinLHACharged, recBinLHACharged, weight);
          h_tu_response_pTD_charged_split->Fill(genBinpTDCharged, recBinpTDCharged, weight);
          h_tu_response_width_charged_split->Fill(genBinWidthCharged, recBinWidthCharged, weight);
          h_tu_response_thrust_charged_split->Fill(genBinThrustCharged, recBinThrustCharged, weight);
          // correct weighting part
          h_tu_response_puppiMultiplicity_charged_split->Fill(genBinPuppiMultCharged, underflow_bin, corr_weight);
          h_tu_response_LHA_charged_split->Fill(genBinLHACharged, underflow_bin, corr_weight);
          h_tu_response_pTD_charged_split->Fill(genBinpTDCharged, underflow_bin, corr_weight);
          h_tu_response_width_charged_split->Fill(genBinWidthCharged, underflow_bin, corr_weight);
          h_tu_response_thrust_charged_split->Fill(genBinThrustCharged, underflow_bin, corr_weight);
        }

      }
    } // end of for loop over genjets
  } // end if passGen

}


TH1F * QGAnalysisUnfoldHists::copy_book_th1f(TH1 * h, const std::string & newName) {
  // DO NOT USE h->GetXaxis()->GetXbins()->GetArray() to clone bins, it just doesn't work
  return book<TH1F>(newName.c_str(),
                    h->GetTitle(),
                    h->GetNbinsX(),
                    h->GetXaxis()->GetXmin(),
                    h->GetXaxis()->GetXmax());
}


TH2F * QGAnalysisUnfoldHists::copy_book_th2f(TH2 * h, const std::string & newName) {
  return book<TH2F>(newName.c_str(),
                    h->GetTitle(),
                    h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                    h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

}


QGAnalysisUnfoldHists::~QGAnalysisUnfoldHists(){}

#include "UHH2/QGAnalysis/include/QGAnalysisUnfoldHists.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisUnfoldHists::QGAnalysisUnfoldHists(Context & ctx, const string & dirname, int useNJets, const string & selection, const string & reco_sel_handle_name, const string & gen_sel_handle_name):
  Hists(ctx, dirname),
  useNJets_(useNJets),
  rand_(4357), // random number generator for splitting MC into 2 independent groups for unfold testing
  doMCsplit_(true)
  {

  is_mc_ = ctx.get("dataset_type") == "MC";
  doMCsplit_ = is_mc_;

  if (useNJets_ < 0) useNJets_ = 99999; // Do them all

  if (useNJets_ == 0) throw runtime_error("useNJets should be > 0, or < 0 to use all jets in the event");

  if (selection != "dijet" && selection != "zplusjets") {
    throw runtime_error("selection must be dijet or zplusjets");
  }

  doHerwigReweighting = ctx.get("herwig_reweight_file", "") != "";
  if (doHerwigReweighting) {
    TFile f_weight(ctx.get("herwig_reweight_file", "").c_str());
    if (selection == "dijet")
      reweightHist = (TH1F*) f_weight.Get("dijet_reco");
    else if (selection == "zplusjets")
      reweightHist = (TH1F*) f_weight.Get("zpj_reco");

    if (reweightHist == nullptr) {
      doHerwigReweighting = false;
      cout << "WARNING: could not find reweight hist - not reweighting AnalysisHists!" << endl;
    } else {
      reweightHist->SetDirectory(0);
    }
  }

  gen_weight_handle = ctx.get_handle<double>("gen_weight");

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

  // Common flags for TUnfold under/overflow
  bool pt_uf(false), pt_of(false); // handle under/overflow ourselves
  bool var_uf(false), var_of(false); // don't need underflow bin, TUnfold has global uflow bin for when we have fakes/missing

  // LHA
  // -------------------------------------
  detector_tu_binning_LHA = new TUnfoldBinning("detector");
  // add binning scheme for pT underflow regions
  // we handle it ourselves (instead of just having underflow bin in standard detector binning)
  // so that we can also use it to constrain the unfolding
  // (like simultaneously fitting to sideband regions)
  detector_distribution_underflow_LHA = detector_tu_binning_LHA->AddBinning("detector_underflow");
  detector_distribution_underflow_LHA->AddAxis("LHA", Binning::nbins_lha_reco, Binning::lha_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_LHA->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_LHA = detector_tu_binning_LHA->AddBinning("detector");
  detector_distribution_LHA->AddAxis("LHA", Binning::nbins_lha_reco, Binning::lha_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_LHA->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_LHA = new TUnfoldBinning("generator");
  generator_distribution_underflow_LHA = generator_tu_binning_LHA->AddBinning("signal_underflow");
  generator_distribution_underflow_LHA->AddAxis("LHA", Binning::nbins_lha_gen, Binning::lha_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_LHA->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_LHA = generator_tu_binning_LHA->AddBinning("signal");
  generator_distribution_LHA->AddAxis("LHA", Binning::nbins_lha_gen, Binning::lha_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_LHA->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  // Response matrix
  // make tmp copies which we can then copy and use with book<>
  TH2 * h_tu_response_LHA_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_LHA, detector_tu_binning_LHA, "tu_LHA_GenReco");
  h_tu_response_LHA = copy_book_th2f(h_tu_response_LHA_tmp, "_all");
  h_tu_response_LHA_half = copy_book_th2f(h_tu_response_LHA_tmp, "_half");
  delete h_tu_response_LHA_tmp;

  // detector histograms
  TH1 * h_tu_reco_LHA_tmp = detector_tu_binning_LHA->CreateHistogram("hist_LHA_reco");
  h_tu_reco_LHA = copy_book_th1f(h_tu_reco_LHA_tmp, "_all");
  h_tu_reco_LHA_half = copy_book_th1f(h_tu_reco_LHA_tmp, "_half");
  // for fakes, detector binning
  h_tu_reco_LHA_fake = copy_book_th1f(h_tu_reco_LHA_tmp, "_fake_all");
  h_tu_reco_LHA_fake_half = copy_book_th1f(h_tu_reco_LHA_tmp, "_fake_half");
  delete h_tu_reco_LHA_tmp;

  // truth histograms
  TH1 * h_tu_gen_LHA_tmp = generator_tu_binning_LHA->CreateHistogram("hist_LHA_truth");
  h_tu_gen_LHA = copy_book_th1f(h_tu_gen_LHA_tmp, "_all");
  h_tu_gen_LHA_half = copy_book_th1f(h_tu_gen_LHA_tmp, "_half");
  delete h_tu_gen_LHA_tmp;

  // detector variable, but using gen binning for comparison later
  h_tu_reco_LHA_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_LHA->Clone("hist_LHA_reco_gen_binning"));
  h_tu_reco_LHA_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_LHA->Clone("hist_LHA_reco_gen_binning_half"));
  // for fakes, gen binning
  h_tu_reco_LHA_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_LHA->Clone("hist_LHA_reco_fake_gen_binning"));
  h_tu_reco_LHA_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_LHA->Clone("hist_LHA_reco_fake_gen_binning_half"));


  // Charged LHA
  // -------------------------------------
  detector_tu_binning_LHA_charged = new TUnfoldBinning("detector");
  detector_distribution_underflow_LHA_charged = detector_tu_binning_LHA_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_LHA_charged->AddAxis("LHA_charged", Binning::nbins_lha_charged_reco, Binning::lha_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_LHA_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_LHA_charged = detector_tu_binning_LHA_charged->AddBinning("detector");
  detector_distribution_LHA_charged->AddAxis("LHA_charged", Binning::nbins_lha_charged_reco, Binning::lha_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_LHA_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);


  generator_tu_binning_LHA_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_LHA_charged = generator_tu_binning_LHA_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_LHA_charged->AddAxis("LHA_charged", Binning::nbins_lha_charged_gen, Binning::lha_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_LHA_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_LHA_charged = generator_tu_binning_LHA_charged->AddBinning("signal");
  generator_distribution_LHA_charged->AddAxis("LHA_charged", Binning::nbins_lha_charged_gen, Binning::lha_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_LHA_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);


  TH2 * h_tu_response_LHA_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_LHA_charged, detector_tu_binning_LHA_charged, "tu_LHA_charged_GenReco");
  h_tu_response_LHA_charged = copy_book_th2f(h_tu_response_LHA_charged_tmp, "_all");
  h_tu_response_LHA_charged_half = copy_book_th2f(h_tu_response_LHA_charged_tmp, "_half");
  delete h_tu_response_LHA_charged_tmp;

  TH1 * h_tu_reco_LHA_charged_tmp = detector_tu_binning_LHA_charged->CreateHistogram("hist_LHA_charged_reco");
  h_tu_reco_LHA_charged = copy_book_th1f(h_tu_reco_LHA_charged_tmp, "_all");
  h_tu_reco_LHA_charged_half = copy_book_th1f(h_tu_reco_LHA_charged_tmp, "_half");
  h_tu_reco_LHA_charged_fake = copy_book_th1f(h_tu_reco_LHA_charged_tmp, "_fake_all");
  h_tu_reco_LHA_charged_fake_half = copy_book_th1f(h_tu_reco_LHA_charged_tmp, "_fake_half");
  delete h_tu_reco_LHA_charged_tmp;

  TH1 * h_tu_gen_LHA_charged_tmp = generator_tu_binning_LHA_charged->CreateHistogram("hist_LHA_charged_truth");
  h_tu_gen_LHA_charged = copy_book_th1f(h_tu_gen_LHA_charged_tmp, "_all");
  h_tu_gen_LHA_charged_half = copy_book_th1f(h_tu_gen_LHA_charged_tmp, "_half");
  delete h_tu_gen_LHA_charged_tmp;

  h_tu_reco_LHA_charged_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_LHA_charged->Clone("hist_LHA_charged_reco_gen_binning"));
  h_tu_reco_LHA_charged_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_LHA_charged->Clone("hist_LHA_charged_reco_gen_binning_half"));
  h_tu_reco_LHA_charged_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_LHA_charged->Clone("hist_LHA_charged_reco_fake_gen_binning"));
  h_tu_reco_LHA_charged_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_LHA_charged->Clone("hist_LHA_charged_reco_fake_gen_binning_half"));

  // puppi multiplicity
  // -------------------------------------
  detector_tu_binning_puppiMultiplicity = new TUnfoldBinning("detector");

  detector_distribution_underflow_puppiMultiplicity = detector_tu_binning_puppiMultiplicity->AddBinning("detector_underflow");
  detector_distribution_underflow_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_puppiMultiplicity_reco, Binning::puppiMultiplicity_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_puppiMultiplicity->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_puppiMultiplicity = detector_tu_binning_puppiMultiplicity->AddBinning("detector");
  detector_distribution_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_puppiMultiplicity_reco, Binning::puppiMultiplicity_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_puppiMultiplicity->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_puppiMultiplicity = new TUnfoldBinning("generator");

  generator_distribution_underflow_puppiMultiplicity = generator_tu_binning_puppiMultiplicity->AddBinning("signal_underflow");
  generator_distribution_underflow_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_puppiMultiplicity_gen, Binning::puppiMultiplicity_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_puppiMultiplicity->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_puppiMultiplicity = generator_tu_binning_puppiMultiplicity->AddBinning("signal");
  generator_distribution_puppiMultiplicity->AddAxis("puppiMultiplicity", Binning::nbins_puppiMultiplicity_gen, Binning::puppiMultiplicity_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_puppiMultiplicity->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  TH2 * h_tu_response_puppiMultiplicity_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_puppiMultiplicity, detector_tu_binning_puppiMultiplicity, "tu_puppiMultiplicity_GenReco");
  h_tu_response_puppiMultiplicity = copy_book_th2f(h_tu_response_puppiMultiplicity_tmp, "_all");
  h_tu_response_puppiMultiplicity_half = copy_book_th2f(h_tu_response_puppiMultiplicity_tmp, "_half");
  delete h_tu_response_puppiMultiplicity_tmp;

  TH1 * h_tu_reco_puppiMultiplicity_tmp = detector_tu_binning_puppiMultiplicity->CreateHistogram("hist_puppiMultiplicity_reco");
  h_tu_reco_puppiMultiplicity = copy_book_th1f(h_tu_reco_puppiMultiplicity_tmp, "_all");
  h_tu_reco_puppiMultiplicity_half = copy_book_th1f(h_tu_reco_puppiMultiplicity_tmp, "_half");
  h_tu_reco_puppiMultiplicity_fake = copy_book_th1f(h_tu_reco_puppiMultiplicity_tmp, "_fake_all");
  h_tu_reco_puppiMultiplicity_fake_half = copy_book_th1f(h_tu_reco_puppiMultiplicity_tmp, "_fake_half");
  delete h_tu_reco_puppiMultiplicity_tmp;

  TH1 * h_tu_gen_puppiMultiplicity_tmp = generator_tu_binning_puppiMultiplicity->CreateHistogram("hist_puppiMultiplicity_truth");
  h_tu_gen_puppiMultiplicity = copy_book_th1f(h_tu_gen_puppiMultiplicity_tmp, "_all");
  h_tu_gen_puppiMultiplicity_half = copy_book_th1f(h_tu_gen_puppiMultiplicity_tmp, "_half");
  delete h_tu_gen_puppiMultiplicity_tmp;

  h_tu_reco_puppiMultiplicity_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_puppiMultiplicity->Clone("hist_puppiMultiplicity_reco_gen_binning"));
  h_tu_reco_puppiMultiplicity_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_puppiMultiplicity->Clone("hist_puppiMultiplicity_reco_gen_binning_half"));

  h_tu_reco_puppiMultiplicity_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_puppiMultiplicity->Clone("hist_puppiMultiplicity_reco_fake_gen_binning"));
  h_tu_reco_puppiMultiplicity_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_puppiMultiplicity->Clone("hist_puppiMultiplicity_reco_fake_gen_binning_half"));

  // Charged PUPPI multiplicity
  // -------------------------------------
  detector_tu_binning_puppiMultiplicity_charged = new TUnfoldBinning("detector");

  detector_distribution_underflow_puppiMultiplicity_charged = detector_tu_binning_puppiMultiplicity_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_puppiMultiplicity_charged_reco, Binning::puppiMultiplicity_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_puppiMultiplicity_charged = detector_tu_binning_puppiMultiplicity_charged->AddBinning("detector");
  detector_distribution_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_puppiMultiplicity_charged_reco, Binning::puppiMultiplicity_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_puppiMultiplicity_charged = new TUnfoldBinning("generator");

  generator_distribution_underflow_puppiMultiplicity_charged = generator_tu_binning_puppiMultiplicity_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_puppiMultiplicity_charged_gen, Binning::puppiMultiplicity_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_puppiMultiplicity_charged = generator_tu_binning_puppiMultiplicity_charged->AddBinning("signal");
  generator_distribution_puppiMultiplicity_charged->AddAxis("puppiMultiplicity_charged", Binning::nbins_puppiMultiplicity_charged_gen, Binning::puppiMultiplicity_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_puppiMultiplicity_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  TH2 * h_tu_response_puppiMultiplicity_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_puppiMultiplicity_charged, detector_tu_binning_puppiMultiplicity_charged, "tu_puppiMultiplicity_charged_GenReco");
  h_tu_response_puppiMultiplicity_charged = copy_book_th2f(h_tu_response_puppiMultiplicity_charged_tmp, "_all");
  h_tu_response_puppiMultiplicity_charged_half = copy_book_th2f(h_tu_response_puppiMultiplicity_charged_tmp, "_half");
  delete h_tu_response_puppiMultiplicity_charged_tmp;

  TH1 * h_tu_reco_puppiMultiplicity_charged_tmp = detector_tu_binning_puppiMultiplicity_charged->CreateHistogram("hist_puppiMultiplicity_charged_reco");
  h_tu_reco_puppiMultiplicity_charged = copy_book_th1f(h_tu_reco_puppiMultiplicity_charged_tmp, "_all");
  h_tu_reco_puppiMultiplicity_charged_half = copy_book_th1f(h_tu_reco_puppiMultiplicity_charged_tmp, "_half");
  h_tu_reco_puppiMultiplicity_charged_fake = copy_book_th1f(h_tu_reco_puppiMultiplicity_charged_tmp, "_fake_all");
  h_tu_reco_puppiMultiplicity_charged_fake_half = copy_book_th1f(h_tu_reco_puppiMultiplicity_charged_tmp, "_fake_half");
  delete h_tu_reco_puppiMultiplicity_charged_tmp;

  TH1 * h_tu_gen_puppiMultiplicity_charged_tmp = generator_tu_binning_puppiMultiplicity_charged->CreateHistogram("hist_puppiMultiplicity_charged_truth");
  h_tu_gen_puppiMultiplicity_charged = copy_book_th1f(h_tu_gen_puppiMultiplicity_charged_tmp, "_all");
  h_tu_gen_puppiMultiplicity_charged_half = copy_book_th1f(h_tu_gen_puppiMultiplicity_charged_tmp, "_half");
  delete h_tu_gen_puppiMultiplicity_charged_tmp;

  h_tu_reco_puppiMultiplicity_charged_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_puppiMultiplicity_charged->Clone("hist_puppiMultiplicity_charged_reco_gen_binning"));
  h_tu_reco_puppiMultiplicity_charged_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_puppiMultiplicity_charged->Clone("hist_puppiMultiplicity_charged_reco_gen_binning_half"));

  h_tu_reco_puppiMultiplicity_charged_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_puppiMultiplicity_charged->Clone("hist_puppiMultiplicity_charged_reco_fake_gen_binning"));
  h_tu_reco_puppiMultiplicity_charged_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_puppiMultiplicity_charged->Clone("hist_puppiMultiplicity_charged_reco_fake_gen_binning_half"));

  // pTD
  // -------------------------------------
  detector_tu_binning_pTD = new TUnfoldBinning("detector");
  detector_distribution_underflow_pTD = detector_tu_binning_pTD->AddBinning("detector_underflow");
  detector_distribution_underflow_pTD->AddAxis("pTD", Binning::nbins_pTD_reco, Binning::pTD_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_pTD->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_pTD = detector_tu_binning_pTD->AddBinning("detector");
  detector_distribution_pTD->AddAxis("pTD", Binning::nbins_pTD_reco, Binning::pTD_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_pTD->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_pTD = new TUnfoldBinning("generator");
  generator_distribution_underflow_pTD = generator_tu_binning_pTD->AddBinning("signal_underflow");
  generator_distribution_underflow_pTD->AddAxis("pTD", Binning::nbins_pTD_gen, Binning::pTD_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_pTD->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_pTD = generator_tu_binning_pTD->AddBinning("signal");
  generator_distribution_pTD->AddAxis("pTD", Binning::nbins_pTD_gen, Binning::pTD_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_pTD->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  TH2 * h_tu_response_pTD_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_pTD, detector_tu_binning_pTD, "tu_pTD_GenReco");
  h_tu_response_pTD = copy_book_th2f(h_tu_response_pTD_tmp, "_all");
  h_tu_response_pTD_half = copy_book_th2f(h_tu_response_pTD_tmp, "_half");
  delete h_tu_response_pTD_tmp;

  TH1 * h_tu_reco_pTD_tmp = detector_tu_binning_pTD->CreateHistogram("hist_pTD_reco");
  h_tu_reco_pTD = copy_book_th1f(h_tu_reco_pTD_tmp, "_all");
  h_tu_reco_pTD_half = copy_book_th1f(h_tu_reco_pTD_tmp, "_half");
  h_tu_reco_pTD_fake = copy_book_th1f(h_tu_reco_pTD_tmp, "_fake_all");
  h_tu_reco_pTD_fake_half = copy_book_th1f(h_tu_reco_pTD_tmp, "_fake_half");
  delete h_tu_reco_pTD_tmp;

  TH1 * h_tu_gen_pTD_tmp = generator_tu_binning_pTD->CreateHistogram("hist_pTD_truth");
  h_tu_gen_pTD = copy_book_th1f(h_tu_gen_pTD_tmp, "_all");
  h_tu_gen_pTD_half = copy_book_th1f(h_tu_gen_pTD_tmp, "_half");
  delete h_tu_gen_pTD_tmp;

  h_tu_reco_pTD_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_pTD->Clone("hist_pTD_reco_gen_binning"));
  h_tu_reco_pTD_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_pTD->Clone("hist_pTD_reco_gen_binning_half"));

  h_tu_reco_pTD_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_pTD->Clone("hist_pTD_reco_fake_gen_binning"));
  h_tu_reco_pTD_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_pTD->Clone("hist_pTD_reco_fake_gen_binning_half"));

  // Charged pTD
  // -------------------------------------
  detector_tu_binning_pTD_charged = new TUnfoldBinning("detector");
  detector_distribution_underflow_pTD_charged = detector_tu_binning_pTD_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_pTD_charged->AddAxis("pTD_charged", Binning::nbins_pTD_charged_reco, Binning::pTD_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_pTD_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_pTD_charged = detector_tu_binning_pTD_charged->AddBinning("detector");
  detector_distribution_pTD_charged->AddAxis("pTD_charged", Binning::nbins_pTD_charged_reco, Binning::pTD_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_pTD_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_pTD_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_pTD_charged = generator_tu_binning_pTD_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_pTD_charged->AddAxis("pTD_charged", Binning::nbins_pTD_charged_gen, Binning::pTD_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_pTD_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_pTD_charged = generator_tu_binning_pTD_charged->AddBinning("signal");
  generator_distribution_pTD_charged->AddAxis("pTD_charged", Binning::nbins_pTD_charged_gen, Binning::pTD_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_pTD_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  TH2 * h_tu_response_pTD_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_pTD_charged, detector_tu_binning_pTD_charged, "tu_pTD_charged_GenReco");
  h_tu_response_pTD_charged = copy_book_th2f(h_tu_response_pTD_charged_tmp, "_all");
  h_tu_response_pTD_charged_half = copy_book_th2f(h_tu_response_pTD_charged_tmp, "_half");
  delete h_tu_response_pTD_charged_tmp;

  TH1 * h_tu_reco_pTD_charged_tmp = detector_tu_binning_pTD_charged->CreateHistogram("hist_pTD_charged_reco");
  h_tu_reco_pTD_charged = copy_book_th1f(h_tu_reco_pTD_charged_tmp, "_all");
  h_tu_reco_pTD_charged_half = copy_book_th1f(h_tu_reco_pTD_charged_tmp, "_half");
  h_tu_reco_pTD_charged_fake = copy_book_th1f(h_tu_reco_pTD_charged_tmp, "_fake_all");
  h_tu_reco_pTD_charged_fake_half = copy_book_th1f(h_tu_reco_pTD_charged_tmp, "_fake_half");
  delete h_tu_reco_pTD_charged_tmp;

  TH1 * h_tu_gen_pTD_charged_tmp = generator_tu_binning_pTD_charged->CreateHistogram("hist_pTD_charged_truth");
  h_tu_gen_pTD_charged = copy_book_th1f(h_tu_gen_pTD_charged_tmp, "_all");
  h_tu_gen_pTD_charged_half = copy_book_th1f(h_tu_gen_pTD_charged_tmp, "_half");
  delete h_tu_gen_pTD_charged_tmp;

  h_tu_reco_pTD_charged_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_pTD_charged->Clone("hist_pTD_charged_reco_gen_binning"));
  h_tu_reco_pTD_charged_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_pTD_charged->Clone("hist_pTD_charged_reco_gen_binning_half"));

  h_tu_reco_pTD_charged_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_pTD_charged->Clone("hist_pTD_charged_reco_fake_gen_binning"));
  h_tu_reco_pTD_charged_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_pTD_charged->Clone("hist_pTD_charged_reco_fake_gen_binning_half"));

  // thrust
  // -------------------------------------
  detector_tu_binning_thrust = new TUnfoldBinning("detector");
  detector_distribution_underflow_thrust = detector_tu_binning_thrust->AddBinning("detector_underflow");
  detector_distribution_underflow_thrust->AddAxis("thrust", Binning::nbins_thrust_reco, Binning::thrust_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_thrust->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_thrust = detector_tu_binning_thrust->AddBinning("detector");
  detector_distribution_thrust->AddAxis("thrust", Binning::nbins_thrust_reco, Binning::thrust_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_thrust->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_thrust = new TUnfoldBinning("generator");
  generator_distribution_underflow_thrust = generator_tu_binning_thrust->AddBinning("signal_underflow");
  generator_distribution_underflow_thrust->AddAxis("thrust", Binning::nbins_thrust_gen, Binning::thrust_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_thrust->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_thrust = generator_tu_binning_thrust->AddBinning("signal");
  generator_distribution_thrust->AddAxis("thrust", Binning::nbins_thrust_gen, Binning::thrust_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_thrust->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  TH2 * h_tu_response_thrust_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_thrust, detector_tu_binning_thrust, "tu_thrust_GenReco");
  h_tu_response_thrust = copy_book_th2f(h_tu_response_thrust_tmp, "_all");
  h_tu_response_thrust_half = copy_book_th2f(h_tu_response_thrust_tmp, "_half");
  delete h_tu_response_thrust_tmp;

  TH1 * h_tu_reco_thrust_tmp = detector_tu_binning_thrust->CreateHistogram("hist_thrust_reco");
  h_tu_reco_thrust = copy_book_th1f(h_tu_reco_thrust_tmp, "_all");
  h_tu_reco_thrust_half = copy_book_th1f(h_tu_reco_thrust_tmp, "_half");
  h_tu_reco_thrust_fake = copy_book_th1f(h_tu_reco_thrust_tmp, "_fake_all");
  h_tu_reco_thrust_fake_half = copy_book_th1f(h_tu_reco_thrust_tmp, "_fake_half");
  delete h_tu_reco_thrust_tmp;

  TH1 * h_tu_gen_thrust_tmp = generator_tu_binning_thrust->CreateHistogram("hist_thrust_truth");
  h_tu_gen_thrust = copy_book_th1f(h_tu_gen_thrust_tmp, "_all");
  h_tu_gen_thrust_half = copy_book_th1f(h_tu_gen_thrust_tmp, "_half");
  delete h_tu_gen_thrust_tmp;

  h_tu_reco_thrust_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_thrust->Clone("hist_thrust_reco_gen_binning"));
  h_tu_reco_thrust_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_thrust->Clone("hist_thrust_reco_gen_binning_half"));

  h_tu_reco_thrust_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_thrust->Clone("hist_thrust_reco_fake_gen_binning"));
  h_tu_reco_thrust_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_thrust->Clone("hist_thrust_reco_fake_gen_binning_half"));

  // Charged thrust
  // -------------------------------------
  detector_tu_binning_thrust_charged = new TUnfoldBinning("detector");
  detector_distribution_underflow_thrust_charged = detector_tu_binning_thrust_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_thrust_charged->AddAxis("thrust_charged", Binning::nbins_thrust_charged_reco, Binning::thrust_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_thrust_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_thrust_charged = detector_tu_binning_thrust_charged->AddBinning("detector");
  detector_distribution_thrust_charged->AddAxis("thrust_charged", Binning::nbins_thrust_charged_reco, Binning::thrust_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_thrust_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_thrust_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_thrust_charged = generator_tu_binning_thrust_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_thrust_charged->AddAxis("thrust_charged", Binning::nbins_thrust_charged_gen, Binning::thrust_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_thrust_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_thrust_charged = generator_tu_binning_thrust_charged->AddBinning("signal");
  generator_distribution_thrust_charged->AddAxis("thrust_charged", Binning::nbins_thrust_charged_gen, Binning::thrust_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_thrust_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  TH2 * h_tu_response_thrust_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_thrust_charged, detector_tu_binning_thrust_charged, "tu_thrust_charged_GenReco");
  h_tu_response_thrust_charged = copy_book_th2f(h_tu_response_thrust_charged_tmp, "_all");
  h_tu_response_thrust_charged_half = copy_book_th2f(h_tu_response_thrust_charged_tmp, "_half");
  delete h_tu_response_thrust_charged_tmp;

  TH1 * h_tu_reco_thrust_charged_tmp = detector_tu_binning_thrust_charged->CreateHistogram("hist_thrust_charged_reco");
  h_tu_reco_thrust_charged = copy_book_th1f(h_tu_reco_thrust_charged_tmp, "_all");
  h_tu_reco_thrust_charged_half = copy_book_th1f(h_tu_reco_thrust_charged_tmp, "_half");
  h_tu_reco_thrust_charged_fake = copy_book_th1f(h_tu_reco_thrust_charged_tmp, "_fake_all");
  h_tu_reco_thrust_charged_fake_half = copy_book_th1f(h_tu_reco_thrust_charged_tmp, "_fake_half");
  delete h_tu_reco_thrust_charged_tmp;

  TH1 * h_tu_gen_thrust_charged_tmp = generator_tu_binning_thrust_charged->CreateHistogram("hist_thrust_charged_truth");
  h_tu_gen_thrust_charged = copy_book_th1f(h_tu_gen_thrust_charged_tmp, "_all");
  h_tu_gen_thrust_charged_half = copy_book_th1f(h_tu_gen_thrust_charged_tmp, "_half");
  delete h_tu_gen_thrust_charged_tmp;

  h_tu_reco_thrust_charged_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_thrust_charged->Clone("hist_thrust_charged_reco_gen_binning"));
  h_tu_reco_thrust_charged_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_thrust_charged->Clone("hist_thrust_charged_reco_gen_binning_half"));

  h_tu_reco_thrust_charged_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_thrust_charged->Clone("hist_thrust_charged_reco_fake_gen_binning"));
  h_tu_reco_thrust_charged_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_thrust_charged->Clone("hist_thrust_charged_reco_fake_gen_binning_half"));

  // width
  // -------------------------------------
  detector_tu_binning_width = new TUnfoldBinning("detector");
  detector_distribution_underflow_width = detector_tu_binning_width->AddBinning("detector_underflow");
  detector_distribution_underflow_width->AddAxis("width", Binning::nbins_width_reco, Binning::width_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_width->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_width = detector_tu_binning_width->AddBinning("detector");
  detector_distribution_width->AddAxis("width", Binning::nbins_width_reco, Binning::width_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_width->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_width = new TUnfoldBinning("generator");
  generator_distribution_underflow_width = generator_tu_binning_width->AddBinning("signal_underflow");
  generator_distribution_underflow_width->AddAxis("width", Binning::nbins_width_gen, Binning::width_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_width->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_width = generator_tu_binning_width->AddBinning("signal");
  generator_distribution_width->AddAxis("width", Binning::nbins_width_gen, Binning::width_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_width->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  TH2 * h_tu_response_width_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_width, detector_tu_binning_width, "tu_width_GenReco");
  h_tu_response_width = copy_book_th2f(h_tu_response_width_tmp, "_all");
  h_tu_response_width_half = copy_book_th2f(h_tu_response_width_tmp, "_half");
  delete h_tu_response_width_tmp;

  TH1 * h_tu_reco_width_tmp = detector_tu_binning_width->CreateHistogram("hist_width_reco");
  h_tu_reco_width = copy_book_th1f(h_tu_reco_width_tmp, "_all");
  h_tu_reco_width_half = copy_book_th1f(h_tu_reco_width_tmp, "_half");
  h_tu_reco_width_fake = copy_book_th1f(h_tu_reco_width_tmp, "_fake_all");
  h_tu_reco_width_fake_half = copy_book_th1f(h_tu_reco_width_tmp, "_fake_half");
  delete h_tu_reco_width_tmp;

  TH1 * h_tu_gen_width_tmp = generator_tu_binning_width->CreateHistogram("hist_width_truth");
  h_tu_gen_width = copy_book_th1f(h_tu_gen_width_tmp, "_all");
  h_tu_gen_width_half = copy_book_th1f(h_tu_gen_width_tmp, "_half");
  delete h_tu_gen_width_tmp;

  h_tu_reco_width_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_width->Clone("hist_width_reco_gen_binning"));
  h_tu_reco_width_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_width->Clone("hist_width_reco_gen_binning_half"));

  h_tu_reco_width_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_width->Clone("hist_width_reco_fake_gen_binning"));
  h_tu_reco_width_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_width->Clone("hist_width_reco_fake_gen_binning_half"));

  // Charged width
  // -------------------------------------
  detector_tu_binning_width_charged = new TUnfoldBinning("detector");
  detector_distribution_underflow_width_charged = detector_tu_binning_width_charged->AddBinning("detector_underflow");
  detector_distribution_underflow_width_charged->AddAxis("width_charged", Binning::nbins_width_charged_reco, Binning::width_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_underflow_width_charged->AddAxis("pt", nbins_pt_reco_underflow, pt_bin_edges_reco_underflow.data(), pt_uf, pt_of);

  detector_distribution_width_charged = detector_tu_binning_width_charged->AddBinning("detector");
  detector_distribution_width_charged->AddAxis("width_charged", Binning::nbins_width_charged_reco, Binning::width_charged_bin_edges_reco.data(), var_uf, var_of);
  detector_distribution_width_charged->AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco.data(), pt_uf, pt_of);

  generator_tu_binning_width_charged = new TUnfoldBinning("generator");
  generator_distribution_underflow_width_charged = generator_tu_binning_width_charged->AddBinning("signal_underflow");
  generator_distribution_underflow_width_charged->AddAxis("width_charged", Binning::nbins_width_charged_gen, Binning::width_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_underflow_width_charged->AddAxis("pt", nbins_pt_gen_underflow, pt_bin_edges_gen_underflow.data(), pt_uf, pt_of);

  generator_distribution_width_charged = generator_tu_binning_width_charged->AddBinning("signal");
  generator_distribution_width_charged->AddAxis("width_charged", Binning::nbins_width_charged_gen, Binning::width_charged_bin_edges_gen.data(), var_uf, var_of);
  generator_distribution_width_charged->AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen.data(), pt_uf, pt_of);

  TH2 * h_tu_response_width_charged_tmp = TUnfoldBinning::CreateHistogramOfMigrations(generator_tu_binning_width_charged, detector_tu_binning_width_charged, "tu_width_charged_GenReco");
  h_tu_response_width_charged = copy_book_th2f(h_tu_response_width_charged_tmp, "_all");
  h_tu_response_width_charged_half = copy_book_th2f(h_tu_response_width_charged_tmp, "_half");
  delete h_tu_response_width_charged_tmp;

  TH1 * h_tu_reco_width_charged_tmp = detector_tu_binning_width_charged->CreateHistogram("hist_width_charged_reco");
  h_tu_reco_width_charged = copy_book_th1f(h_tu_reco_width_charged_tmp, "_all");
  h_tu_reco_width_charged_half = copy_book_th1f(h_tu_reco_width_charged_tmp, "_half");
  h_tu_reco_width_charged_fake = copy_book_th1f(h_tu_reco_width_charged_tmp, "_fake_all");
  h_tu_reco_width_charged_fake_half = copy_book_th1f(h_tu_reco_width_charged_tmp, "_fake_half");
  delete h_tu_reco_width_charged_tmp;

  TH1 * h_tu_gen_width_charged_tmp = generator_tu_binning_width_charged->CreateHistogram("hist_width_charged_truth");
  h_tu_gen_width_charged = copy_book_th1f(h_tu_gen_width_charged_tmp, "_all");
  h_tu_gen_width_charged_half = copy_book_th1f(h_tu_gen_width_charged_tmp, "_half");
  delete h_tu_gen_width_charged_tmp;

  h_tu_reco_width_charged_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_width_charged->Clone("hist_width_charged_reco_gen_binning"));
  h_tu_reco_width_charged_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_width_charged->Clone("hist_width_charged_reco_gen_binning_half"));

  h_tu_reco_width_charged_fake_gen_binning = copy_book_th1f((TH1F*) h_tu_gen_width_charged->Clone("hist_width_charged_reco_fake_gen_binning"));
  h_tu_reco_width_charged_fake_gen_binning_half = copy_book_th1f((TH1F*) h_tu_gen_width_charged->Clone("hist_width_charged_reco_fake_gen_binning_half"));

  if (is_mc_) {
    genJetsLambda_handle = ctx.get_handle< std::vector<GenJetLambdaBundle> > ("GoodGenJetLambdas");
    genJetsChargedLambda_handle = ctx.get_handle< std::vector<GenJetLambdaBundle> > ("GoodGenJetChargedLambdas");
    pass_gen_handle = ctx.get_handle<bool> (gen_sel_handle_name);
  }

  jetsLambda_handle = ctx.get_handle< std::vector<JetLambdaBundle> > ("JetLambdas");
  jetsChargedLambda_handle = ctx.get_handle< std::vector<JetLambdaBundle> > ("JetChargedLambdas");

  pass_reco_handle = ctx.get_handle<bool> (reco_sel_handle_name);
}


void QGAnalysisUnfoldHists::fill(const Event & event){
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();

  // Optionally apply weight to Herwig to ensure spectrum matches Pythia spectrum
  float herwig_weight = 1.;
  if (doHerwigReweighting && Njets >= 1) {
    float pt = jets->at(0).pt();
    if (pt >= reweightHist->GetXaxis()->GetXmax()) {
      pt = reweightHist->GetXaxis()->GetXmax() - 0.1;
    }
    int bin_num = reweightHist->GetXaxis()->FindBin(pt);
    herwig_weight = reweightHist->GetBinContent(bin_num);
  }

  double weight = event.weight * herwig_weight;

  // extract the separate gen & reco weight components, needed for TUnfold
  double gen_weight = event.get(gen_weight_handle);
  gen_weight *= herwig_weight;
  double reco_weight = weight / gen_weight;

  // Get selection flags
  bool passReco = event.get(pass_reco_handle);
  bool passGen(false);

  // Get (Gen)Jet-Lambda Bundles
  // ptr since then it can be null
  const std::vector<GenJetLambdaBundle> * genjetLambdas = nullptr;
  const std::vector<GenJetLambdaBundle> * genjetChargedLambdas = nullptr;
  if (is_mc_) {
    genjetLambdas = &event.get(genJetsLambda_handle);
    genjetChargedLambdas = &event.get(genJetsChargedLambda_handle);
    passGen = event.get(pass_gen_handle);
  }

  const std::vector<JetLambdaBundle> jetLambdas = event.get(jetsLambda_handle);
  const std::vector<JetLambdaBundle> jetChargedLambdas = event.get(jetsChargedLambda_handle);

  // cout << "***" << event.event << endl;

  // Random event flag for MC only filling response matrix
  // This allows us to split the MC into 2 separate samples for testing
  bool onlyResponseHalf = (rand_.Rndm() > 0.5);

  // figure out ave pt across N jets
  // float avePtReco = 0.;
  // float avePtGen = 0.;
  // int nGenJets = 0;
  // for (int i = 0; i < useNJets_; i++) {
  //   const Jet & thisjet = jets->at(i);
  //   avePtReco += thisjet.pt();
  //   if (is_mc_) {
  //     if (thisjet.genjet_index() > -1) {
  //       nGenJets++;
  //       const GenJetWithParts & genjet = genjets->at(thisjet.genjet_index());
  //       avePtGen += genjet.pt();
  //     }
  //   }
  // }

  // avePtReco /= useNJets_;
  // avePtGen /= nGenJets;
  std::cout << " Filling unfold reco jet hists" << std::endl;

  // Fill reco jet 1D hists
  // ---------------------------------------------------------------------------
  // At this point, all jet filtering etc should have already been performed
  // Fill ignoring if there is a genjet match or not
  if (passReco) {
    for (int i = 0; i < useNJets_; i++) {
      const Jet & thisjet = jetLambdas.at(i).jet;
      LambdaCalculator<PFParticle> recoJetCalc = jetLambdas.at(i).lambda; // can't be const as getLambda() modifies it
      // FIXME check this corresponds to same jet as normal lambdas?
      LambdaCalculator<PFParticle> recoJetCalcCharged = jetChargedLambdas.at(i).lambda;

      float lha = recoJetCalc.getLambda(1, 0.5);
      float puppiMult = recoJetCalc.getLambda(0, 0);
      float ptd = recoJetCalc.getLambda(2, 0);
      float width = recoJetCalc.getLambda(1, 1);
      float thrust = recoJetCalc.getLambda(1, 2);

      float lha_charged = recoJetCalcCharged.getLambda(1, 0.5);
      float puppiMult_charged = recoJetCalcCharged.getLambda(0, 0);
      float ptd_charged = recoJetCalcCharged.getLambda(2, 0);
      float width_charged = recoJetCalcCharged.getLambda(1, 1);
      float thrust_charged = recoJetCalcCharged.getLambda(1, 2);

      float jet_pt = thisjet.pt();

      // default to 0 for underflow
      int recBinLHA(0), recBinPuppiMult(0), recBinpTD(0), recBinThrust(0), recBinWidth(0);
      int recBinLHACharged(0), recBinPuppiMultCharged(0), recBinpTDCharged(0), recBinThrustCharged(0), recBinWidthCharged(0);
      // generator-bin for reco-level quantities (for displaying later)
      int genBinLHA(0), genBinPuppiMult(0), genBinpTD(0), genBinThrust(0), genBinWidth(0);
      int genBinLHACharged(0), genBinPuppiMultCharged(0), genBinpTDCharged(0), genBinThrustCharged(0), genBinWidthCharged(0);
      bool isUnderflow = (jet_pt < Binning::pt_bin_edges_reco[0]);
      // here we get bin number based on if underflow or not
      if (isUnderflow) {
        recBinLHA = detector_distribution_underflow_LHA->GetGlobalBinNumber(lha, jet_pt);
        recBinPuppiMult = detector_distribution_underflow_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
        recBinpTD = detector_distribution_underflow_pTD->GetGlobalBinNumber(ptd, jet_pt);
        recBinThrust = detector_distribution_underflow_thrust->GetGlobalBinNumber(thrust, jet_pt);
        recBinWidth = detector_distribution_underflow_width->GetGlobalBinNumber(width, jet_pt);

        recBinLHACharged = detector_distribution_underflow_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
        recBinPuppiMultCharged = detector_distribution_underflow_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
        recBinpTDCharged = detector_distribution_underflow_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
        recBinThrustCharged = detector_distribution_underflow_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
        recBinWidthCharged = detector_distribution_underflow_width_charged->GetGlobalBinNumber(width_charged, jet_pt);

        genBinLHA = generator_distribution_underflow_LHA->GetGlobalBinNumber(lha, jet_pt);
        genBinPuppiMult = generator_distribution_underflow_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
        genBinpTD = generator_distribution_underflow_pTD->GetGlobalBinNumber(ptd, jet_pt);
        genBinThrust = generator_distribution_underflow_thrust->GetGlobalBinNumber(thrust, jet_pt);
        genBinWidth = generator_distribution_underflow_width->GetGlobalBinNumber(width, jet_pt);

        genBinLHACharged = generator_distribution_underflow_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
        genBinPuppiMultCharged = generator_distribution_underflow_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
        genBinpTDCharged = generator_distribution_underflow_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
        genBinThrustCharged = generator_distribution_underflow_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
        genBinWidthCharged = generator_distribution_underflow_width_charged->GetGlobalBinNumber(width_charged, jet_pt);
      } else {
        recBinLHA = detector_distribution_LHA->GetGlobalBinNumber(lha, jet_pt);
        recBinPuppiMult = detector_distribution_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
        recBinpTD = detector_distribution_pTD->GetGlobalBinNumber(ptd, jet_pt);
        recBinThrust = detector_distribution_thrust->GetGlobalBinNumber(thrust, jet_pt);
        recBinWidth = detector_distribution_width->GetGlobalBinNumber(width, jet_pt);

        recBinLHACharged = detector_distribution_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
        recBinPuppiMultCharged = detector_distribution_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
        recBinpTDCharged = detector_distribution_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
        recBinThrustCharged = detector_distribution_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
        recBinWidthCharged = detector_distribution_width_charged->GetGlobalBinNumber(width_charged, jet_pt);

        genBinLHA = generator_distribution_LHA->GetGlobalBinNumber(lha, jet_pt);
        genBinPuppiMult = generator_distribution_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
        genBinpTD = generator_distribution_pTD->GetGlobalBinNumber(ptd, jet_pt);
        genBinThrust = generator_distribution_thrust->GetGlobalBinNumber(thrust, jet_pt);
        genBinWidth = generator_distribution_width->GetGlobalBinNumber(width, jet_pt);

        genBinLHACharged = generator_distribution_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
        genBinPuppiMultCharged = generator_distribution_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
        genBinpTDCharged = generator_distribution_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
        genBinThrustCharged = generator_distribution_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
        genBinWidthCharged = generator_distribution_width_charged->GetGlobalBinNumber(width_charged, jet_pt);
      }

      h_tu_reco_LHA->Fill(recBinLHA, weight);
      h_tu_reco_puppiMultiplicity->Fill(recBinPuppiMult, weight);
      h_tu_reco_pTD->Fill(recBinpTD, weight);
      h_tu_reco_thrust->Fill(recBinThrust, weight);
      h_tu_reco_width->Fill(recBinWidth, weight);

      h_tu_reco_LHA_charged->Fill(recBinLHACharged, weight);
      h_tu_reco_puppiMultiplicity_charged->Fill(recBinPuppiMultCharged, weight);
      h_tu_reco_pTD_charged->Fill(recBinpTDCharged, weight);
      h_tu_reco_thrust_charged->Fill(recBinThrustCharged, weight);
      h_tu_reco_width_charged->Fill(recBinWidthCharged, weight);

      h_tu_reco_LHA_gen_binning->Fill(genBinLHA, weight);
      h_tu_reco_puppiMultiplicity_gen_binning->Fill(genBinPuppiMult, weight);
      h_tu_reco_pTD_gen_binning->Fill(genBinpTD, weight);
      h_tu_reco_thrust_gen_binning->Fill(genBinThrust, weight);
      h_tu_reco_width_gen_binning->Fill(genBinWidth, weight);

      h_tu_reco_LHA_charged_gen_binning->Fill(genBinLHACharged, weight);
      h_tu_reco_puppiMultiplicity_charged_gen_binning->Fill(genBinPuppiMultCharged, weight);
      h_tu_reco_pTD_charged_gen_binning->Fill(genBinpTDCharged, weight);
      h_tu_reco_thrust_charged_gen_binning->Fill(genBinThrustCharged, weight);
      h_tu_reco_width_charged_gen_binning->Fill(genBinWidthCharged, weight);

      if (!onlyResponseHalf && doMCsplit_) {
        h_tu_reco_LHA_half->Fill(recBinLHA, weight);
        h_tu_reco_puppiMultiplicity_half->Fill(recBinPuppiMult, weight);
        h_tu_reco_pTD_half->Fill(recBinpTD, weight);
        h_tu_reco_thrust_half->Fill(recBinThrust, weight);
        h_tu_reco_width_half->Fill(recBinWidth, weight);

        h_tu_reco_LHA_charged_half->Fill(recBinLHACharged, weight);
        h_tu_reco_puppiMultiplicity_charged_half->Fill(recBinPuppiMultCharged, weight);
        h_tu_reco_pTD_charged_half->Fill(recBinpTDCharged, weight);
        h_tu_reco_thrust_charged_half->Fill(recBinThrustCharged, weight);
        h_tu_reco_width_charged_half->Fill(recBinWidthCharged, weight);

        h_tu_reco_LHA_gen_binning_half->Fill(genBinLHA, weight);
        h_tu_reco_puppiMultiplicity_gen_binning_half->Fill(genBinPuppiMult, weight);
        h_tu_reco_pTD_gen_binning_half->Fill(genBinpTD, weight);
        h_tu_reco_thrust_gen_binning_half->Fill(genBinThrust, weight);
        h_tu_reco_width_gen_binning_half->Fill(genBinWidth, weight);

        h_tu_reco_LHA_charged_gen_binning_half->Fill(genBinLHACharged, weight);
        h_tu_reco_puppiMultiplicity_charged_gen_binning_half->Fill(genBinPuppiMultCharged, weight);
        h_tu_reco_pTD_charged_gen_binning_half->Fill(genBinpTDCharged, weight);
        h_tu_reco_thrust_charged_gen_binning_half->Fill(genBinThrustCharged, weight);
        h_tu_reco_width_charged_gen_binning_half->Fill(genBinWidthCharged, weight);
      }

      // Fill fakes, both in detector & truth binnings
      // Fake = not pass Gen, or pass Gen but matching genjet not one of the
      // selection ones (i.e. doesn't count as a match)
      // -----------------------------------------------------------------------
      if (!passGen || (thisjet.genjet_index() < 0 || thisjet.genjet_index() >= useNJets_)) {
        h_tu_reco_LHA_fake->Fill(recBinLHA, weight);
        h_tu_reco_puppiMultiplicity_fake->Fill(recBinPuppiMult, weight);
        h_tu_reco_pTD_fake->Fill(recBinpTD, weight);
        h_tu_reco_thrust_fake->Fill(recBinThrust, weight);
        h_tu_reco_width_fake->Fill(recBinWidth, weight);

        h_tu_reco_LHA_charged_fake->Fill(recBinLHACharged, weight);
        h_tu_reco_puppiMultiplicity_charged_fake->Fill(recBinPuppiMultCharged, weight);
        h_tu_reco_pTD_charged_fake->Fill(recBinpTDCharged, weight);
        h_tu_reco_thrust_charged_fake->Fill(recBinThrustCharged, weight);
        h_tu_reco_width_charged_fake->Fill(recBinWidthCharged, weight);

        h_tu_reco_LHA_fake_gen_binning->Fill(genBinLHA, weight);
        h_tu_reco_puppiMultiplicity_fake_gen_binning->Fill(genBinPuppiMult, weight);
        h_tu_reco_pTD_fake_gen_binning->Fill(genBinpTD, weight);
        h_tu_reco_thrust_fake_gen_binning->Fill(genBinThrust, weight);
        h_tu_reco_width_fake_gen_binning->Fill(genBinWidth, weight);

        h_tu_reco_LHA_charged_fake_gen_binning->Fill(genBinLHACharged, weight);
        h_tu_reco_puppiMultiplicity_charged_fake_gen_binning->Fill(genBinPuppiMultCharged, weight);
        h_tu_reco_pTD_charged_fake_gen_binning->Fill(genBinpTDCharged, weight);
        h_tu_reco_thrust_charged_fake_gen_binning->Fill(genBinThrustCharged, weight);
        h_tu_reco_width_charged_fake_gen_binning->Fill(genBinWidthCharged, weight);

        if (!onlyResponseHalf && doMCsplit_) {
          h_tu_reco_LHA_fake_half->Fill(recBinLHA, weight);
          h_tu_reco_puppiMultiplicity_fake_half->Fill(recBinPuppiMult, weight);
          h_tu_reco_pTD_fake_half->Fill(recBinpTD, weight);
          h_tu_reco_thrust_fake_half->Fill(recBinThrust, weight);
          h_tu_reco_width_fake_half->Fill(recBinWidth, weight);

          h_tu_reco_LHA_charged_fake_half->Fill(recBinLHACharged, weight);
          h_tu_reco_puppiMultiplicity_charged_fake_half->Fill(recBinPuppiMultCharged, weight);
          h_tu_reco_pTD_charged_fake_half->Fill(recBinpTDCharged, weight);
          h_tu_reco_thrust_charged_fake_half->Fill(recBinThrustCharged, weight);
          h_tu_reco_width_charged_fake_half->Fill(recBinWidthCharged, weight);

          h_tu_reco_LHA_fake_gen_binning_half->Fill(genBinLHA, weight);
          h_tu_reco_puppiMultiplicity_fake_gen_binning_half->Fill(genBinPuppiMult, weight);
          h_tu_reco_pTD_fake_gen_binning_half->Fill(genBinpTD, weight);
          h_tu_reco_thrust_fake_gen_binning_half->Fill(genBinThrust, weight);
          h_tu_reco_width_fake_gen_binning_half->Fill(genBinWidth, weight);

          h_tu_reco_LHA_charged_fake_gen_binning_half->Fill(genBinLHACharged, weight);
          h_tu_reco_puppiMultiplicity_charged_fake_gen_binning_half->Fill(genBinPuppiMultCharged, weight);
          h_tu_reco_pTD_charged_fake_gen_binning_half->Fill(genBinpTDCharged, weight);
          h_tu_reco_thrust_charged_fake_gen_binning_half->Fill(genBinThrustCharged, weight);
          h_tu_reco_width_charged_fake_gen_binning_half->Fill(genBinWidthCharged, weight);
        }
      }
    }
  }

  // Fill gen jet 1D hists
  // ---------------------------------------------------------------------------
  if (is_mc_ && passGen) {
    for (int i = 0; i < useNJets_; i++) {
      const GenJetWithParts & thisjet = genjetLambdas->at(i).jet;
      LambdaCalculator<GenParticle> genJetCalc = genjetLambdas->at(i).lambda; // can't be const as getLambda() modifies it
      // FIXME check this corresponds to same jet as normal lambdas?
      LambdaCalculator<GenParticle> genJetCalcCharged = genjetChargedLambdas->at(i).lambda;

      float gen_lha = genJetCalc.getLambda(1, 0.5);
      float gen_mult = genJetCalc.getLambda(0, 0);
      float gen_ptd = genJetCalc.getLambda(2, 0);
      float gen_width = genJetCalc.getLambda(1, 1);
      float gen_thrust = genJetCalc.getLambda(1, 2);

      float gen_lha_charged = genJetCalcCharged.getLambda(1, 0.5);
      float gen_mult_charged = genJetCalcCharged.getLambda(0, 0);
      float gen_ptd_charged = genJetCalcCharged.getLambda(2, 0);
      float gen_width_charged = genJetCalcCharged.getLambda(1, 1);
      float gen_thrust_charged = genJetCalcCharged.getLambda(1, 2);

      float genjet_pt = thisjet.pt();

      // Get TUnfold gen bins
      // default to 0 for underflow
      int genBinLHA(0), genBinPuppiMult(0), genBinpTD(0), genBinThrust(0), genBinWidth(0);
      int genBinLHACharged(0), genBinPuppiMultCharged(0), genBinpTDCharged(0), genBinThrustCharged(0), genBinWidthCharged(0);
      bool isUnderflowGen = (genjet_pt < Binning::pt_bin_edges_gen[0]);
      // here we get bin number based on if underflow or not
      if (isUnderflowGen) {
        genBinLHA = generator_distribution_underflow_LHA->GetGlobalBinNumber(gen_lha, genjet_pt);
        genBinPuppiMult = generator_distribution_underflow_puppiMultiplicity->GetGlobalBinNumber(gen_mult, genjet_pt);
        genBinpTD = generator_distribution_underflow_pTD->GetGlobalBinNumber(gen_ptd, genjet_pt);
        genBinThrust = generator_distribution_underflow_thrust->GetGlobalBinNumber(gen_thrust, genjet_pt);
        genBinWidth = generator_distribution_underflow_width->GetGlobalBinNumber(gen_width, genjet_pt);

        genBinLHACharged = generator_distribution_underflow_LHA_charged->GetGlobalBinNumber(gen_lha_charged, genjet_pt);
        genBinPuppiMultCharged = generator_distribution_underflow_puppiMultiplicity_charged->GetGlobalBinNumber(gen_mult_charged, genjet_pt);
        genBinpTDCharged = generator_distribution_underflow_pTD_charged->GetGlobalBinNumber(gen_ptd_charged, genjet_pt);
        genBinThrustCharged = generator_distribution_underflow_thrust_charged->GetGlobalBinNumber(gen_thrust_charged, genjet_pt);
        genBinWidthCharged = generator_distribution_underflow_width_charged->GetGlobalBinNumber(gen_width_charged, genjet_pt);
      } else {
        genBinLHA = generator_distribution_LHA->GetGlobalBinNumber(gen_lha, genjet_pt);
        genBinPuppiMult = generator_distribution_puppiMultiplicity->GetGlobalBinNumber(gen_mult, genjet_pt);
        genBinpTD = generator_distribution_pTD->GetGlobalBinNumber(gen_ptd, genjet_pt);
        genBinThrust = generator_distribution_thrust->GetGlobalBinNumber(gen_thrust, genjet_pt);
        genBinWidth = generator_distribution_width->GetGlobalBinNumber(gen_width, genjet_pt);

        genBinLHACharged = generator_distribution_LHA_charged->GetGlobalBinNumber(gen_lha_charged, genjet_pt);
        genBinPuppiMultCharged = generator_distribution_puppiMultiplicity_charged->GetGlobalBinNumber(gen_mult_charged, genjet_pt);
        genBinpTDCharged = generator_distribution_pTD_charged->GetGlobalBinNumber(gen_ptd_charged, genjet_pt);
        genBinThrustCharged = generator_distribution_thrust_charged->GetGlobalBinNumber(gen_thrust_charged, genjet_pt);
        genBinWidthCharged = generator_distribution_width_charged->GetGlobalBinNumber(gen_width_charged, genjet_pt);
      }

      // Fill 1D gen hist
      // -----------------------------------------------------------------------
      h_tu_gen_puppiMultiplicity->Fill(genBinPuppiMult, gen_weight);
      h_tu_gen_LHA->Fill(genBinLHA, gen_weight);
      h_tu_gen_pTD->Fill(genBinpTD, gen_weight);
      h_tu_gen_width->Fill(genBinWidth, gen_weight);
      h_tu_gen_thrust->Fill(genBinThrust, gen_weight);

      h_tu_gen_puppiMultiplicity_charged->Fill(genBinPuppiMultCharged, gen_weight);
      h_tu_gen_LHA_charged->Fill(genBinLHACharged, gen_weight);
      h_tu_gen_pTD_charged->Fill(genBinpTDCharged, gen_weight);
      h_tu_gen_width_charged->Fill(genBinWidthCharged, gen_weight);
      h_tu_gen_thrust_charged->Fill(genBinThrustCharged, gen_weight);

      if (!onlyResponseHalf && doMCsplit_) {
        // Fill 1D gen hist
        h_tu_gen_puppiMultiplicity_half->Fill(genBinPuppiMult, gen_weight);
        h_tu_gen_LHA_half->Fill(genBinLHA, gen_weight);
        h_tu_gen_pTD_half->Fill(genBinpTD, gen_weight);
        h_tu_gen_width_half->Fill(genBinWidth, gen_weight);
        h_tu_gen_thrust_half->Fill(genBinThrust, gen_weight);

        h_tu_gen_puppiMultiplicity_charged_half->Fill(genBinPuppiMultCharged, gen_weight);
        h_tu_gen_LHA_charged_half->Fill(genBinLHACharged, gen_weight);
        h_tu_gen_pTD_charged_half->Fill(genBinpTDCharged, gen_weight);
        h_tu_gen_width_charged_half->Fill(genBinWidthCharged, gen_weight);
        h_tu_gen_thrust_charged_half->Fill(genBinThrustCharged, gen_weight);
      }

      // Now fill response matrices
      // -----------------------------------------------------------------------
      // Loop through genjets, since if there isn't a reco jet then it's a miss-reco,
      // whereas a recojet & no genjet is a fake, which we don't want in our migration matrix
      int recBinLHA(0), recBinPuppiMult(0), recBinpTD(0), recBinThrust(0), recBinWidth(0);
      int recBinLHACharged(0), recBinPuppiMultCharged(0), recBinpTDCharged(0), recBinThrustCharged(0), recBinWidthCharged(0);
      int recoInd = -1;
      if (passReco) { // only valid match if reco selection also passed
        // First find matching recojet - annoyingly matches are held in the Jet class, not GenJet
        // default to 0 for underflow
        for (int j=0; j<useNJets_; j++) {
          if (jetLambdas.at(j).jet.genjet_index() == i) { recoInd = j; break; }
        }
        if (recoInd >= 0) {
          const Jet & thisjet = jetLambdas.at(recoInd).jet;
          LambdaCalculator<PFParticle> recoJetCalc = jetLambdas.at(recoInd).lambda; // can't be const as getLambda() modifies it
          // FIXME check this corresponds to same jet as normal lambdas?
          LambdaCalculator<PFParticle> recoJetCalcCharged = jetChargedLambdas.at(recoInd).lambda;

          float lha = recoJetCalc.getLambda(1, 0.5);
          float puppiMult = recoJetCalc.getLambda(0, 0);
          float ptd = recoJetCalc.getLambda(2, 0);
          float width = recoJetCalc.getLambda(1, 1);
          float thrust = recoJetCalc.getLambda(1, 2);

          float lha_charged = recoJetCalcCharged.getLambda(1, 0.5);
          float puppiMult_charged = recoJetCalcCharged.getLambda(0, 0);
          float ptd_charged = recoJetCalcCharged.getLambda(2, 0);
          float width_charged = recoJetCalcCharged.getLambda(1, 1);
          float thrust_charged = recoJetCalcCharged.getLambda(1, 2);

          float jet_pt = thisjet.pt();

          bool isUnderflow = (jet_pt < Binning::pt_bin_edges_reco[0]);
          // here we get bin number based on if underflow or not
          if (isUnderflow) {
            recBinLHA = detector_distribution_underflow_LHA->GetGlobalBinNumber(lha, jet_pt);
            recBinPuppiMult = detector_distribution_underflow_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
            recBinpTD = detector_distribution_underflow_pTD->GetGlobalBinNumber(ptd, jet_pt);
            recBinThrust = detector_distribution_underflow_thrust->GetGlobalBinNumber(thrust, jet_pt);
            recBinWidth = detector_distribution_underflow_width->GetGlobalBinNumber(width, jet_pt);

            recBinLHACharged = detector_distribution_underflow_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
            recBinPuppiMultCharged = detector_distribution_underflow_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
            recBinpTDCharged = detector_distribution_underflow_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
            recBinThrustCharged = detector_distribution_underflow_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
            recBinWidthCharged = detector_distribution_underflow_width_charged->GetGlobalBinNumber(width_charged, jet_pt);
          } else {
            recBinLHA = detector_distribution_LHA->GetGlobalBinNumber(lha, jet_pt);
            recBinPuppiMult = detector_distribution_puppiMultiplicity->GetGlobalBinNumber(puppiMult, jet_pt);
            recBinpTD = detector_distribution_pTD->GetGlobalBinNumber(ptd, jet_pt);
            recBinThrust = detector_distribution_thrust->GetGlobalBinNumber(thrust, jet_pt);
            recBinWidth = detector_distribution_width->GetGlobalBinNumber(width, jet_pt);

            recBinLHACharged = detector_distribution_LHA_charged->GetGlobalBinNumber(lha_charged, jet_pt);
            recBinPuppiMultCharged = detector_distribution_puppiMultiplicity_charged->GetGlobalBinNumber(puppiMult_charged, jet_pt);
            recBinpTDCharged = detector_distribution_pTD_charged->GetGlobalBinNumber(ptd_charged, jet_pt);
            recBinThrustCharged = detector_distribution_thrust_charged->GetGlobalBinNumber(thrust_charged, jet_pt);
            recBinWidthCharged = detector_distribution_width_charged->GetGlobalBinNumber(width_charged, jet_pt);
          }
        }
      } // end if passReco

      // Fill TUnfold 2D response maps
      // ---------------------------------------------------------------
      h_tu_response_puppiMultiplicity->Fill(genBinPuppiMult, recBinPuppiMult, weight);
      h_tu_response_LHA->Fill(genBinLHA, recBinLHA, weight);
      h_tu_response_pTD->Fill(genBinpTD, recBinpTD, weight);
      h_tu_response_width->Fill(genBinWidth, recBinWidth, weight);
      h_tu_response_thrust->Fill(genBinThrust, recBinThrust, weight);

      h_tu_response_puppiMultiplicity_charged->Fill(genBinPuppiMultCharged, recBinPuppiMultCharged, weight);
      h_tu_response_LHA_charged->Fill(genBinLHACharged, recBinLHACharged, weight);
      h_tu_response_pTD_charged->Fill(genBinpTDCharged, recBinpTDCharged, weight);
      h_tu_response_width_charged->Fill(genBinWidthCharged, recBinWidthCharged, weight);
      h_tu_response_thrust_charged->Fill(genBinThrustCharged, recBinThrustCharged, weight);

      // we need to add in an extra part, such that the 1D projection on the gen axis
      // agrees with the 1D gen histogram
      // i.e. account for the difference between the reco_weight (goes into event.weight) and gen_weight
      // we use the underflow bin for this
      // (0 as bin edges are tunfold bin numbers, not physical values, so the first bin is 0-1,
      // not to be confused with ROOT binning, starting at 1 :s)
      double corr_weight = gen_weight * (1 - reco_weight);
      int underflow_bin = 0;
      h_tu_response_puppiMultiplicity->Fill(genBinPuppiMult, underflow_bin, corr_weight);
      h_tu_response_LHA->Fill(genBinLHA, underflow_bin, corr_weight);
      h_tu_response_pTD->Fill(genBinpTD, underflow_bin, corr_weight);
      h_tu_response_width->Fill(genBinWidth, underflow_bin, corr_weight);
      h_tu_response_thrust->Fill(genBinThrust, underflow_bin, corr_weight);

      h_tu_response_puppiMultiplicity_charged->Fill(genBinPuppiMultCharged, underflow_bin, corr_weight);
      h_tu_response_LHA_charged->Fill(genBinLHACharged, underflow_bin, corr_weight);
      h_tu_response_pTD_charged->Fill(genBinpTDCharged, underflow_bin, corr_weight);
      h_tu_response_width_charged->Fill(genBinWidthCharged, underflow_bin, corr_weight);
      h_tu_response_thrust_charged->Fill(genBinThrustCharged, underflow_bin, corr_weight);

      if (onlyResponseHalf && doMCsplit_) {
        h_tu_response_puppiMultiplicity_half->Fill(genBinPuppiMult, recBinPuppiMult, weight);
        h_tu_response_LHA_half->Fill(genBinLHA, recBinLHA, weight);
        h_tu_response_pTD_half->Fill(genBinpTD, recBinpTD, weight);
        h_tu_response_width_half->Fill(genBinWidth, recBinWidth, weight);
        h_tu_response_thrust_half->Fill(genBinThrust, recBinThrust, weight);

        h_tu_response_puppiMultiplicity_charged_half->Fill(genBinPuppiMultCharged, recBinPuppiMultCharged, weight);
        h_tu_response_LHA_charged_half->Fill(genBinLHACharged, recBinLHACharged, weight);
        h_tu_response_pTD_charged_half->Fill(genBinpTDCharged, recBinpTDCharged, weight);
        h_tu_response_width_charged_half->Fill(genBinWidthCharged, recBinWidthCharged, weight);
        h_tu_response_thrust_charged_half->Fill(genBinThrustCharged, recBinThrustCharged, weight);

        // extra part for correct weighting
        h_tu_response_puppiMultiplicity_half->Fill(genBinPuppiMult, underflow_bin, corr_weight);
        h_tu_response_LHA_half->Fill(genBinLHA, underflow_bin, corr_weight);
        h_tu_response_pTD_half->Fill(genBinpTD, underflow_bin, corr_weight);
        h_tu_response_width_half->Fill(genBinWidth, underflow_bin, corr_weight);
        h_tu_response_thrust_half->Fill(genBinThrust, underflow_bin, corr_weight);

        h_tu_response_puppiMultiplicity_charged_half->Fill(genBinPuppiMultCharged, underflow_bin, corr_weight);
        h_tu_response_LHA_charged_half->Fill(genBinLHACharged, underflow_bin, corr_weight);
        h_tu_response_pTD_charged_half->Fill(genBinpTDCharged, underflow_bin, corr_weight);
        h_tu_response_width_charged_half->Fill(genBinWidthCharged, underflow_bin, corr_weight);
        h_tu_response_thrust_charged_half->Fill(genBinThrustCharged, underflow_bin, corr_weight);
      }
    } // end of for loop over genjets
  }

}


TH1F * QGAnalysisUnfoldHists::copy_book_th1f(TH1 * h, const std::string & append) {
  // DO NOT USE h->GetXaxis()->GetXbins()->GetArray() to clone bins, it just doesn't work
  return book<TH1F>((std::string(h->GetName())+append).c_str(),
                                 h->GetTitle(),
                                 h->GetNbinsX(),
                                 h->GetXaxis()->GetXmin(),
                                 h->GetXaxis()->GetXmax());
}


TH2F * QGAnalysisUnfoldHists::copy_book_th2f(TH2 * h, const std::string & append) {
  return book<TH2F>((std::string(h->GetName())+append).c_str(),
                                 h->GetTitle(),
                                 h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                                 h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

}


QGAnalysisUnfoldHists::~QGAnalysisUnfoldHists(){}

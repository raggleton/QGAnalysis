#include "UHH2/QGAnalysis/include/QGAnalysisHists.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisHists::QGAnalysisHists(Context & ctx, const string & dirname, int useNJets, const string & selection):
  Hists(ctx, dirname),
  useNJets_(useNJets),
  bins_pt_response(get_pt_bin_edges(3500, 1.3)),
  nbins_pt_response(bins_pt_response.size() - 1),
  neutral_pf_hadron_shift_(0.)
  {

  is_mc_ = ctx.get("dataset_type") == "MC";
  useGenPartonFlav_ = (is_mc_ && ctx.get("useGenPartonFlav") == "true");
  doPuppi_ = (ctx.get("PURemoval") == "PUPPI");

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

  // book all histograms here
  int nWeightBins = 11;
  double weightBins [nWeightBins+1] = {1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 10};

  h_weights = book<TH1F>("weights", ";weight;", nWeightBins, weightBins);
  h_weights_vs_pt = book<TH2F>("weights_vs_pt", ";weight;", nWeightBins, weightBins, 200, 0., 2000.);
  h_pthat_vs_weight = book<TH2F>("pthat_vs_weight", ";weight;ptHat", nWeightBins, weightBins, 200, 0, 2000);
  h_pthat_vs_jet_pt = book<TH2F>("pthat_vs_jet_pt", ";p_{T}^{leading jet};ptHat", 200, 0, 2000, 200, 0, 2000);

  // RECO jet hists
  // --------------
  // For all jets
  int nPtBins = 2000;
  float ptMin(0.), ptMax(2000.);
  h_jet_pt = book<TH1F>("jet_pt", ";p_{T}^{j} [GeV];", 2*nPtBins, ptMin, ptMax); // use as a catch-all to see we haven't missed highPT jets
  h_jet_pt_unweighted = book<TH1F>("jet_pt_unweighted", ";p_{T}^{j} [GeV];", 2*nPtBins, ptMin, ptMax); // use as a catch-all to see we haven't missed highPT jets

  int nEtaBins = 50;
  float etaMin(-5), etaMax(5);
  h_jet_eta = book<TH1F>("jet_eta", ";#eta^{j};", nEtaBins, etaMin, etaMax);

  h_jet_flavour = book<TH1F>("jet_flavour", "jet flavour;PDGID;", 23, -0.5, 22.5);
  h_jet_genParton_flavour = book<TH1F>("jet_genParton_flavour", "jet flavour (genParton);PDGID;", 23, -0.5, 22.5);

  int nMultBins = 150;
  h_jet_multiplicity = book<TH1F>("jet_multiplicity", ";# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity = book<TH1F>("jet_puppiMultiplicity", ";# of PUPPI-weighted constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);

  int nBins = 100;
  h_jet_LHA = book<TH1F>("jet_LHA", ";LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_jet_pTD = book<TH1F>("jet_pTD", ";p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_jet_width = book<TH1F>("jet_width", ";Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_jet_thrust = book<TH1F>("jet_thrust", ";Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // gen-reco response hists
  float rsp_max = 5.;
  int nBinsNormRsp = 500;
  // hist names must always end in _repsone or _rel_response
  h_jet_multiplicity_response = book<TH2F>("jet_multiplicity_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_charged_response = book<TH2F>("jet_multiplicity_charged_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO), charged only", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_rel_response = book<TH2F>("jet_multiplicity_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_multiplicity_charged_rel_response = book<TH2F>("jet_multiplicity_charged_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  // 2D hist for each pT bin
  // multiplicity_response_binned.resize(nbins_pt_response, std::vector<TH2F*>(nbins_pt_response));
  // for (int i=0; i < nbins_pt_response; i++) {
  //   for (int j=0; j < nbins_pt_response; j++) {
  //     multiplicity_response_binned[i][j] = book<TH2F>(TString::Format("jet_multiplicity_response_bin_%d_%d", i, j),
  //                                                     TString::Format("%g < p_{T}^{Gen} < %g GeV, %g < p_{T}^{Reco} < %g GeV;# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO)",
  //                                                                     bins_pt_response[i], bins_pt_response[i+1], bins_pt_response[j], bins_pt_response[j+1]),
  //                                                     nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  //   }
  // }

  h_jet_puppiMultiplicity_response = book<TH2F>("jet_puppiMultiplicity_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_charged_response = book<TH2F>("jet_puppiMultiplicity_charged_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO, charged only)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_rel_response = book<TH2F>("jet_puppiMultiplicity_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_puppiMultiplicity_charged_rel_response = book<TH2F>("jet_puppiMultiplicity_charged_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  // puppiMultiplicity_response_binned.resize(nbins_pt_response, std::vector<TH2F*>(nbins_pt_response));
  // for (int i=0; i < nbins_pt_response; i++) {
  //   for (int j=0; j < nbins_pt_response; j++) {
  //     puppiMultiplicity_response_binned[i][j] = book<TH2F>(TString::Format("jet_puppiMultiplicity_response_bin_%d_%d", i, j),
  //                                                     TString::Format("%g < p_{T}^{Gen} < %g GeV, %g < p_{T}^{Reco} < %g GeV;# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO)",
  //                                                                     bins_pt_response[i], bins_pt_response[i+1], bins_pt_response[j], bins_pt_response[j+1]),
  //                                                     nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  //   }
  // }

  h_jet_LHA_response = book<TH2F>("jet_LHA_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_charged_response = book<TH2F>("jet_LHA_charged_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_rel_response = book<TH2F>("jet_LHA_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_LHA_charged_rel_response = book<TH2F>("jet_LHA_charged_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  // LHA_response_binned.resize(nbins_pt_response, std::vector<TH2F*>(nbins_pt_response));
  // for (int i=0; i < nbins_pt_response; i++) {
  //   for (int j=0; j < nbins_pt_response; j++) {
  //     LHA_response_binned[i][j] = book<TH2F>(TString::Format("jet_lha_response_bin_%d_%d", i, j),
  //                                                     TString::Format("%g < p_{T}^{Gen} < %g GeV, %g < p_{T}^{Reco} < %g GeV;LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO)",
  //                                                                     bins_pt_response[i], bins_pt_response[i+1], bins_pt_response[j], bins_pt_response[j+1]),
  //                                                     nBins, 0, 1, nBins, 0, 1);
  //   }
  // }


  h_jet_pTD_response = book<TH2F>("jet_pTD_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_charged_response = book<TH2F>("jet_pTD_charged_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN, charged);p_{T}^{D} (#lambda_{0}^{2}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_rel_response = book<TH2F>("jet_pTD_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_pTD_charged_rel_response = book<TH2F>("jet_pTD_charged_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  // pTD_response_binned.resize(nbins_pt_response, std::vector<TH2F*>(nbins_pt_response));
  // for (int i=0; i < nbins_pt_response; i++) {
  //   for (int j=0; j < nbins_pt_response; j++) {
  //     pTD_response_binned[i][j] = book<TH2F>(TString::Format("jet_pTD_response_bin_%d_%d", i, j),
  //                                                     TString::Format("%g < p_{T}^{Gen} < %g GeV, %g < p_{T}^{Reco} < %g GeV;p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO)",
  //                                                                     bins_pt_response[i], bins_pt_response[i+1], bins_pt_response[j], bins_pt_response[j+1]),
  //                                                     nBins, 0, 1, nBins, 0, 1);
  //   }
  // }

  h_jet_width_response = book<TH2F>("jet_width_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_width_charged_response = book<TH2F>("jet_width_charged_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_width_rel_response = book<TH2F>("jet_width_rel_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_width_charged_rel_response = book<TH2F>("jet_width_charged_rel_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  // width_response_binned.resize(nbins_pt_response, std::vector<TH2F*>(nbins_pt_response));
  // for (int i=0; i < nbins_pt_response; i++) {
  //   for (int j=0; j < nbins_pt_response; j++) {
  //     width_response_binned[i][j] = book<TH2F>(TString::Format("jet_width_response_bin_%d_%d", i, j),
  //                                                     TString::Format("%g < p_{T}^{Gen} < %g GeV, %g < p_{T}^{Reco} < %g GeV;Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO)",
  //                                                                     bins_pt_response[i], bins_pt_response[i+1], bins_pt_response[j], bins_pt_response[j+1]),
  //                                                     nBins, 0, 1, nBins, 0, 1);
  //   }
  // }

  h_jet_thrust_response = book<TH2F>("jet_thrust_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_thrust_charged_response = book<TH2F>("jet_thrust_charged_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_thrust_rel_response = book<TH2F>("jet_thrust_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_thrust_charged_rel_response = book<TH2F>("jet_thrust_charged_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  // thrust_response_binned.resize(nbins_pt_response, std::vector<TH2F*>(nbins_pt_response));
  // for (int i=0; i < nbins_pt_response; i++) {
  //   for (int j=0; j < nbins_pt_response; j++) {
  //     thrust_response_binned[i][j] = book<TH2F>(TString::Format("jet_thrust_response_bin_%d_%d", i, j),
  //                                                     TString::Format("%g < p_{T}^{Gen} < %g GeV, %g < p_{T}^{Reco} < %g GeV;Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO)",
  //                                                                     bins_pt_response[i], bins_pt_response[i+1], bins_pt_response[j], bins_pt_response[j+1]),
  //                                                     nBins, 0, 1, nBins, 0, 1);
  //   }
  // }

  // q jet only
  h_qjet_multiplicity = book<TH1F>("qjet_multiplicity", "q-flavour;# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_qjet_LHA = book<TH1F>("qjet_LHA", "q-flavour;LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_qjet_pTD = book<TH1F>("qjet_pTD", "q-flavour;p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_qjet_width = book<TH1F>("qjet_width", "q-flavour;Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_qjet_thrust = book<TH1F>("qjet_thrust", "q-flavour jet;Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // g jet only
  h_gjet_multiplicity = book<TH1F>("gjet_multiplicity", "g-flavour;# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_gjet_LHA = book<TH1F>("gjet_LHA", "g-flavour;LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_gjet_pTD = book<TH1F>("gjet_pTD", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_gjet_width = book<TH1F>("gjet_width", "g-flavour;Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_gjet_thrust = book<TH1F>("gjet_thrust", "g-flavour jet;Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // 2D versions vs PT
  // -----------------
  // All flavs
  h_jet_multiplicity_vs_pt = book<TH2F>("jet_multiplicity_vs_pt", ";# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_jet_LHA_vs_pt = book<TH2F>("jet_LHA_vs_pt", ";LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_pTD_vs_pt = book<TH2F>("jet_pTD_vs_pt", ";p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_width_vs_pt = book<TH2F>("jet_width_vs_pt", ";Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_thrust_vs_pt = book<TH2F>("jet_thrust_vs_pt", ";Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  int rspNbins = 500;
  float rspMax = 5;
  h_jet_response_vs_genjet_pt = book<TH2F>("jet_response_vs_genjet_pt", ";response;GenJet p_{T} [GeV]", rspNbins, 0, rspMax, nPtBins, ptMin, ptMax);

  h_jet_flavour_vs_pt = book<TH2F>("jet_flavour_vs_pt", "jet flavour;PDGID;Jet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet_genParton_flavour_vs_pt = book<TH2F>("jet_genParton_flavour_vs_pt", "jet flavour;PDGID;Jet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet1_flavour_vs_pt = book<TH2F>("jet1_flavour_vs_pt", "jet1 flavour;PDGID;Jet1 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet1_genParton_flavour_vs_pt = book<TH2F>("jet1_genParton_flavour_vs_pt", "jet1 flavour;PDGID;Jet1 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet2_flavour_vs_pt = book<TH2F>("jet2_flavour_vs_pt", "jet2 flavour;PDGID;Jet2 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet2_genParton_flavour_vs_pt = book<TH2F>("jet2_genParton_flavour_vs_pt", "jet2 flavour;PDGID;Jet2 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);

  h_jet_flavour_vs_eta = book<TH2F>("jet_flavour_vs_eta", "jet flavour;PDGID;Jet #eta", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);
  h_jet_genParton_flavour_vs_eta = book<TH2F>("jet_genParton_flavour_vs_eta", "jet flavour;PDGID;Jet #eta", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);

  // q jet only
  h_qjet_multiplicity_vs_pt = book<TH2F>("qjet_multiplicity_vs_pt", "q-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_qjet_LHA_vs_pt = book<TH2F>("qjet_LHA_vs_pt", "q-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_pTD_vs_pt = book<TH2F>("qjet_pTD_vs_pt", "q-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_width_vs_pt = book<TH2F>("qjet_width_vs_pt", "q-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_thrust_vs_pt = book<TH2F>("qjet_thrust_vs_pt", "q-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_response_vs_genjet_pt = book<TH2F>("qjet_response_vs_genjet_pt", ";response;GenJet p_{T} [GeV]", rspNbins, 0, rspMax, nPtBins, ptMin, ptMax);

  // g jet only
  h_gjet_multiplicity_vs_pt = book<TH2F>("gjet_multiplicity_vs_pt", "g-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_gjet_LHA_vs_pt = book<TH2F>("gjet_LHA_vs_pt", "g-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_pTD_vs_pt = book<TH2F>("gjet_pTD_vs_pt", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_width_vs_pt = book<TH2F>("gjet_width_vs_pt", "g-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_thrust_vs_pt = book<TH2F>("gjet_thrust_vs_pt", "g-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_response_vs_genjet_pt = book<TH2F>("gjet_response_vs_genjet_pt", ";response;GenJet p_{T} [GeV]", rspNbins, 0, rspMax, nPtBins, ptMin, ptMax);

  // puppi ones
  h_jet_puppiMultiplicity_vs_pt = book<TH2F>("jet_puppiMultiplicity_vs_pt", ";# of constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_qjet_puppiMultiplicity_vs_pt = book<TH2F>("qjet_puppiMultiplicity_vs_pt", "q-flavour;# of constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_gjet_puppiMultiplicity_vs_pt = book<TH2F>("gjet_puppiMultiplicity_vs_pt", "g-flavour;# of constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);

  // Variable correlation hists
  // --------------------------
  // All flavs
  // h_jet_multiplicity_vs_LHA = book<TH2F>("jet_multiplicity_vs_LHA", ";multiplicity (#lambda_{0}^{0});LHA (#lambda_{0.5}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_jet_multiplicity_vs_pTD = book<TH2F>("jet_multiplicity_vs_pTD", ";multiplicity (#lambda_{0}^{0});{p_{T}^{D}}^2 (#lambda_{0}^{2})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_jet_multiplicity_vs_width = book<TH2F>("jet_multiplicity_vs_width", ";multiplicity (#lambda_{0}^{0});width (#lambda_{1}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_jet_multiplicity_vs_thrust = book<TH2F>("jet_multiplicity_vs_thrust", ";multiplicity (#lambda_{0}^{0});thrust (#lambda_{2}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_jet_LHA_vs_pTD = book<TH2F>("jet_LHA_vs_pTD", ";LHA (#lambda_{0.5}^{1});{p_{T}^{D}}^2 (#lambda_{0}^{2})", nBins, 0, 1, nBins, 0, 1);
  // h_jet_LHA_vs_width = book<TH2F>("jet_LHA_vs_width", ";LHA (#lambda_{0.5}^{1});width (#lambda_{1}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_jet_LHA_vs_thrust = book<TH2F>("jet_LHA_vs_thrust", ";LHA (#lambda_{0.5}^{1});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_jet_pTD_vs_width = book<TH2F>("jet_pTD_vs_width", ";{p_{T}^{D}}^2 (#lambda_{0}^{2});width (#lambda_{1}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_jet_pTD_vs_thrust = book<TH2F>("jet_pTD_vs_thrust", ";{p_{T}^{D}}^2 (#lambda_{0}^{2});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_jet_width_vs_thrust = book<TH2F>("jet_width_vs_thrust", ";width (#lambda_{1}^{1});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);

  // q jet only
  // h_qjet_multiplicity_vs_LHA = book<TH2F>("qjet_multiplicity_vs_LHA", "q-flavour;multiplicity (#lambda_{0}^{0});LHA (#lambda_{0.5}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_qjet_multiplicity_vs_pTD = book<TH2F>("qjet_multiplicity_vs_pTD", "q-flavour;multiplicity (#lambda_{0}^{0});{p_{T}^{D}}^2 (#lambda_{0}^{2})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_qjet_multiplicity_vs_width = book<TH2F>("qjet_multiplicity_vs_width", "q-flavour;multiplicity (#lambda_{0}^{0});width (#lambda_{1}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_qjet_multiplicity_vs_thrust = book<TH2F>("qjet_multiplicity_vs_thrust", "q-flavour;multiplicity (#lambda_{0}^{0});thrust (#lambda_{2}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_qjet_LHA_vs_pTD = book<TH2F>("qjet_LHA_vs_pTD", "q-flavour;LHA (#lambda_{0.5}^{1});{p_{T}^{D}}^2 (#lambda_{0}^{2})", nBins, 0, 1, nBins, 0, 1);
  // h_qjet_LHA_vs_width = book<TH2F>("qjet_LHA_vs_width", "q-flavour;LHA (#lambda_{0.5}^{1});width (#lambda_{1}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_qjet_LHA_vs_thrust = book<TH2F>("qjet_LHA_vs_thrust", "q-flavour;LHA (#lambda_{0.5}^{1});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_qjet_pTD_vs_width = book<TH2F>("qjet_pTD_vs_width", "q-flavour;{p_{T}^{D}}^2 (#lambda_{0}^{2});width (#lambda_{1}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_qjet_pTD_vs_thrust = book<TH2F>("qjet_pTD_vs_thrust", "q-flavour;{p_{T}^{D}}^2 (#lambda_{0}^{2});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_qjet_width_vs_thrust = book<TH2F>("qjet_width_vs_thrust", "q-flavour;width (#lambda_{1}^{1});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);

  // g jet only
  // h_gjet_multiplicity_vs_LHA = book<TH2F>("gjet_multiplicity_vs_LHA", "g-flavour;multiplicity (#lambda_{0}^{0});LHA (#lambda_{0.5}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_gjet_multiplicity_vs_pTD = book<TH2F>("gjet_multiplicity_vs_pTD", "g-flavour;multiplicity (#lambda_{0}^{0});{p_{T}^{D}}^2 (#lambda_{0}^{2})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_gjet_multiplicity_vs_width = book<TH2F>("gjet_multiplicity_vs_width", "g-flavour;multiplicity (#lambda_{0}^{0});width (#lambda_{1}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_gjet_multiplicity_vs_thrust = book<TH2F>("gjet_multiplicity_vs_thrust", "g-flavour;multiplicity (#lambda_{0}^{0});thrust (#lambda_{2}^{1})", nMultBins, 0, nMultBins, nBins, 0, 1);
  // h_gjet_LHA_vs_pTD = book<TH2F>("gjet_LHA_vs_pTD", "g-flavour;LHA (#lambda_{0.5}^{1});{p_{T}^{D}}^2 (#lambda_{0}^{2})", nBins, 0, 1, nBins, 0, 1);
  // h_gjet_LHA_vs_width = book<TH2F>("gjet_LHA_vs_width", "g-flavour;LHA (#lambda_{0.5}^{1});width (#lambda_{1}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_gjet_LHA_vs_thrust = book<TH2F>("gjet_LHA_vs_thrust", "g-flavour;LHA (#lambda_{0.5}^{1});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_gjet_pTD_vs_width = book<TH2F>("gjet_pTD_vs_width", "g-flavour;{p_{T}^{D}}^2 (#lambda_{0}^{2});width (#lambda_{1}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_gjet_pTD_vs_thrust = book<TH2F>("gjet_pTD_vs_thrust", "g-flavour;{p_{T}^{D}}^2 (#lambda_{0}^{2});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);
  // h_gjet_width_vs_thrust = book<TH2F>("gjet_width_vs_thrust", "g-flavour;width (#lambda_{1}^{1});thrust (#lambda_{2}^{1})", nBins, 0, 1, nBins, 0, 1);


  // GENJET hists
  // ------------
  // For all jets
  h_genjet_pt = book<TH1F>("genjet_pt", ";p_{T}^{j} [GeV];", 2*nPtBins, ptMin, ptMax);
  h_genjet_eta = book<TH1F>("genjet_eta", ";#eta^{j};", nEtaBins, etaMin, etaMax);

  h_genjet_multiplicity = book<TH1F>("genjet_multiplicity", ";# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_genjet_LHA = book<TH1F>("genjet_LHA", ";LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_genjet_pTD = book<TH1F>("genjet_pTD", ";p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_genjet_width = book<TH1F>("genjet_width", ";Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_genjet_thrust = book<TH1F>("genjet_thrust", ";Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // Flavour-tagged
  h_qgenjet_multiplicity = book<TH1F>("qgenjet_multiplicity", "q-flavour;# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_qgenjet_LHA = book<TH1F>("qgenjet_LHA", "q-flavour;LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_qgenjet_pTD = book<TH1F>("qgenjet_pTD", "q-flavour;p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_qgenjet_width = book<TH1F>("qgenjet_width", "q-flavour;Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_qgenjet_thrust = book<TH1F>("qgenjet_thrust", "q-flavour jet;Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  h_ggenjet_multiplicity = book<TH1F>("ggenjet_multiplicity", "g-flavour;# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_ggenjet_LHA = book<TH1F>("ggenjet_LHA", "g-flavour;LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_ggenjet_pTD = book<TH1F>("ggenjet_pTD", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_ggenjet_width = book<TH1F>("ggenjet_width", "g-flavour;Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_ggenjet_thrust = book<TH1F>("ggenjet_thrust", "g-flavour jet;Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // 2D versions vs PT
  h_genjet_multiplicity_vs_pt = book<TH2F>("genjet_multiplicity_vs_pt", ";# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_genjet_LHA_vs_pt = book<TH2F>("genjet_LHA_vs_pt", ";LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_genjet_pTD_vs_pt = book<TH2F>("genjet_pTD_vs_pt", ";p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_genjet_width_vs_pt = book<TH2F>("genjet_width_vs_pt", ";Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_genjet_thrust_vs_pt = book<TH2F>("genjet_thrust_vs_pt", ";Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  // Some normalisations that weren't done in the initial calculation
  string jet_cone = ctx.get("JetCone", "AK4");
  if (jet_cone.find("AK4") != string::npos)
    jetRadius = 0.4;
  else if (jet_cone.find("AK8") != string::npos)
    jetRadius = 0.8;
  else if (jet_cone.find("ca15") != string::npos)
    jetRadius = 1.5;
  else
    throw runtime_error("Cannot determine jetRadius in QGAnalysisHists");

  LHA_rescale = pow(jetRadius, 0.5);
  width_rescale = jetRadius;
  thrust_rescale = pow(jetRadius, 2.0);

  if (is_mc_) genJets_handle = ctx.get_handle< std::vector<GenJetWithParts> > ("GoodGenJets");

  std::string neutralPfShift = ctx.get("neutralHadronShift", "nominal");
  if (neutralPfShift == "nominal") {
    neutral_pf_hadron_shift_ = 0.;
  } else if (neutralPfShift == "up") {
    neutral_pf_hadron_shift_ = 1.;
  } else if (neutralPfShift == "down") {
    neutral_pf_hadron_shift_ = -1;
  } else {
    throw runtime_error("neutralPfShift must be nominal, up, or down");
  }
}


void QGAnalysisHists::fill(const Event & event){
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

  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  // Don't forget to always use the weight when filling.
  double weight = event.weight * herwig_weight;

  h_weights->Fill(weight);

  const std::vector<GenJetWithParts> * genjets = nullptr;
  if (is_mc_) {

    if (event.genInfo->binningValues().size() > 0)
      h_pthat_vs_weight->Fill(weight, event.genInfo->binningValues()[0]);

    genjets = &event.get(genJets_handle);
  }

  // Fill reco jet hists
  for (int i = 0; i < useNJets_; i++) {
    const Jet & thisjet = jets->at(i);

    std::vector<PFParticle*> orig_daughters = get_jet_pfparticles(thisjet, event.pfparticles);

    // Do uncertainty shifting of neutral hadron components
    if (neutral_pf_hadron_shift_ != 0) {
      std::vector<PFParticle*> daughters_copy = create_copy(orig_daughters);
      shift_neutral_hadron_pfparticles(daughters_copy, neutral_pf_hadron_shift_, 0.1);
      orig_daughters = daughters_copy;
    }

    std::vector<PFParticle*> daughters;
    float puppiMult = 0;
    for (auto dau : orig_daughters) {
      if (dau->pt() > 1.) {
        daughters.push_back(dau);
        puppiMult += dau->puppiWeight();
      }
    }

    float mult = daughters.size();

    LambdaCalculator<PFParticle> recoJetCalc(daughters, jetRadius, thisjet.v4(), doPuppi_);
    float lha = recoJetCalc.getLambda(1, 0.5);
    float ptd = recoJetCalc.getLambda(2, 0);
    float width = recoJetCalc.getLambda(1, 1);
    float thrust = recoJetCalc.getLambda(1, 2);

    float jet_pt = thisjet.pt();
    h_weights_vs_pt->Fill(weight, jet_pt);
    // if (i == 0) {
    h_weights_vs_pt->Fill(weight, thisjet.pt());
    if (is_mc_) {
      if (event.genInfo->binningValues().size() > 0)
        h_pthat_vs_jet_pt->Fill(thisjet.pt(), event.genInfo->binningValues()[0]);
    }
    // }
    h_jet_pt_unweighted->Fill(jet_pt);
    h_jet_pt->Fill(jet_pt, weight);
    h_jet_eta->Fill(thisjet.eta(), weight);

    h_jet_multiplicity->Fill(mult, weight);
    h_jet_LHA->Fill(lha, weight);
    h_jet_pTD->Fill(ptd, weight);
    h_jet_width->Fill(width, weight);
    h_jet_thrust->Fill(thrust, weight);

    h_jet_multiplicity_vs_pt->Fill(mult, jet_pt, weight);
    h_jet_LHA_vs_pt->Fill(lha, jet_pt, weight);
    h_jet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
    h_jet_width_vs_pt->Fill(width, jet_pt, weight);
    h_jet_thrust_vs_pt->Fill(thrust, jet_pt, weight);

    // float puppiMult = thisjet.get_tag(Jet::puppiMultiplicity);
    h_jet_puppiMultiplicity->Fill(puppiMult, weight);
    h_jet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);

    if (is_mc_) {
      // Store variables for matched GenJet
      bool matchedJet = false;
      float genjet_pt = -1.;
      float response = -1.;
      if (thisjet.genjet_index() > -1) {
        matchedJet = true;
        const GenJetWithParts & genjet = genjets->at(thisjet.genjet_index());
        genjet_pt = genjet.pt();
        response = jet_pt/genjet_pt;
        h_jet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);

        std::vector<GenParticle*> orig_matchedDaughters = get_genjet_genparticles(genjet, event.genparticles);

        std::vector<GenParticle*> matchedDaughters;
        for (auto dau : orig_matchedDaughters) {
          if (dau->pt() > 1) {
            matchedDaughters.push_back(dau);
          }
        }

        // Response hists over all pt
        float gen_mult = matchedDaughters.size();
        h_jet_multiplicity_response->Fill(gen_mult, mult, weight);
        h_jet_puppiMultiplicity_response->Fill(gen_mult, puppiMult, weight);
        h_jet_multiplicity_rel_response->Fill(gen_mult, mult/gen_mult, weight);
        h_jet_puppiMultiplicity_rel_response->Fill(gen_mult, puppiMult/gen_mult, weight);

        LambdaCalculator<GenParticle> matchedGenJetCalc(matchedDaughters, jetRadius, genjet.v4(), false);
        float gen_lha = matchedGenJetCalc.getLambda(1, 0.5);
        h_jet_LHA_response->Fill(gen_lha, lha, weight);
        h_jet_LHA_rel_response->Fill(gen_lha, lha/gen_lha, weight);
        float gen_ptd = matchedGenJetCalc.getLambda(2, 0);
        h_jet_pTD_response->Fill(gen_ptd, ptd, weight);
        h_jet_pTD_rel_response->Fill(gen_ptd, ptd/gen_ptd, weight);
        float gen_width = matchedGenJetCalc.getLambda(1, 1);
        h_jet_width_response->Fill(gen_width, width, weight);
        h_jet_width_rel_response->Fill(gen_width, width/gen_width, weight);
        float gen_thrust = matchedGenJetCalc.getLambda(1, 2);
        h_jet_thrust_response->Fill(gen_thrust, thrust, weight);
        h_jet_thrust_rel_response->Fill(gen_thrust, thrust/gen_thrust, weight);

        // Charged-only daughters version
        std::vector<PFParticle*> chargedDaughters;
        float puppiMult_charged = 0;
        for (auto dau : daughters) {
          if (dau->charge() != 0) {
            chargedDaughters.push_back(dau);
            puppiMult_charged += dau->puppiWeight();
          }
        }
        LambdaCalculator<PFParticle> recoJetCalcCharged(chargedDaughters, jetRadius, thisjet.v4(), doPuppi_);

        std::vector<GenParticle*> matchedChargedDaughters;
        for (auto dau : matchedDaughters) {
          if (dau->charge() != 0) {
            matchedChargedDaughters.push_back(dau);
          }
        }
        float gen_mult_charged = matchedChargedDaughters.size();
        float mult_charged = chargedDaughters.size();
        h_jet_multiplicity_charged_response->Fill(gen_mult_charged, mult_charged, weight);
        h_jet_multiplicity_charged_rel_response->Fill(gen_mult_charged, mult_charged/gen_mult_charged, weight);

        h_jet_puppiMultiplicity_charged_response->Fill(gen_mult_charged, puppiMult_charged, weight);
        h_jet_puppiMultiplicity_charged_rel_response->Fill(gen_mult_charged, puppiMult_charged/gen_mult_charged, weight);

        LambdaCalculator<GenParticle> matchedGenJetCalcCharged(matchedChargedDaughters, jetRadius, genjet.v4(), false);
        float lha_charged = recoJetCalcCharged.getLambda(1, 0.5);
        float gen_lha_charged = matchedGenJetCalcCharged.getLambda(1, 0.5);
        h_jet_LHA_charged_response->Fill(gen_lha_charged, lha_charged, weight);
        h_jet_LHA_charged_rel_response->Fill(gen_lha_charged, lha_charged/gen_lha_charged, weight);

        float ptd_charged = recoJetCalcCharged.getLambda(2, 0);
        float gen_ptd_charged = matchedGenJetCalcCharged.getLambda(2, 0);
        h_jet_pTD_charged_response->Fill(gen_ptd_charged, ptd_charged, weight);
        h_jet_pTD_charged_rel_response->Fill(gen_ptd_charged, ptd_charged/gen_ptd_charged, weight);

        float width_charged = recoJetCalcCharged.getLambda(1, 1);
        float gen_width_charged = matchedGenJetCalcCharged.getLambda(1, 1);
        h_jet_width_charged_response->Fill(gen_width_charged, width_charged, weight);
        h_jet_width_charged_rel_response->Fill(gen_width_charged, width_charged/gen_width_charged, weight);

        float thrust_charged = recoJetCalcCharged.getLambda(1, 2);
        float gen_thrust_charged = matchedGenJetCalcCharged.getLambda(1, 2);
        h_jet_thrust_charged_response->Fill(gen_thrust_charged, thrust_charged, weight);
        h_jet_thrust_charged_rel_response->Fill(gen_thrust_charged, thrust_charged/gen_thrust_charged, weight);

        // Response hists per pt bin
        // auto genit = std::lower_bound(bins_pt_response.begin(), bins_pt_response.end(), genjet_pt);
        // int genind = genit - bins_pt_response.begin() - 1; // need the -1 to get the offset correct
        // if (genind < 0) {
        //   throw std::runtime_error("Invalid genind pos - bins start at too high a value");
        // } else if (genind >= (int) multiplicity_response_binned.size()) {
        //   throw std::runtime_error("Invalid genind pos - bins don't extend far enough");
        // }
        // auto recoit = std::lower_bound(bins_pt_response.begin(), bins_pt_response.end(), jet_pt);
        // int recoind = recoit - bins_pt_response.begin() - 1;
        // if (recoind < 0) {
        //   throw std::runtime_error("Invalid recoind pos - bins start at too high a value");
        // } else if (recoind >= (int) multiplicity_response_binned[genind].size()) {
        //   throw std::runtime_error("Invalid recoind pos - bins don't extend far enough");
        // }
        // multiplicity_response_binned[genind][recoind]->Fill(gen_mult, mult, weight);
        // puppiMultiplicity_response_binned[genind][recoind]->Fill(gen_mult, puppiMult, weight);
        // LHA_response_binned[genind][recoind]->Fill(gen_lha, lha, weight);
        // pTD_response_binned[genind][recoind]->Fill(gen_ptd, ptd, weight);
        // width_response_binned[genind][recoind]->Fill(gen_width, width, weight);
        // thrust_response_binned[genind][recoind]->Fill(gen_thrust, thrust, weight);
      }

      // int jet_flav = get_jet_flavour(thisjet, event.genparticles);
      // int jet_flav = abs(thisjet.genPartonFlavor());
      int jet_flav = useGenPartonFlav_ ? abs(thisjet.genPartonFlavor()) : abs(thisjet.flavor());

      // Split by actual jet flavour - these only make sense for MC
      if (jet_flav == 21) { // gluon jets
        h_gjet_multiplicity->Fill(mult, weight);
        h_gjet_LHA->Fill(lha, weight);
        h_gjet_pTD->Fill(ptd, weight);
        h_gjet_width->Fill(width, weight);
        h_gjet_thrust->Fill(thrust, weight);

        h_gjet_multiplicity_vs_pt->Fill(mult, jet_pt, weight);
        h_gjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
        h_gjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
        h_gjet_width_vs_pt->Fill(width, jet_pt, weight);
        h_gjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);
        if (matchedJet) h_gjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
        h_gjet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);
      } else if ((jet_flav <= 3) && (jet_flav > 0)){ // uds jets
        h_qjet_multiplicity->Fill(mult, weight);
        h_qjet_LHA->Fill(lha, weight);
        h_qjet_pTD->Fill(ptd, weight);
        h_qjet_width->Fill(width, weight);
        h_qjet_thrust->Fill(thrust, weight);

        h_qjet_multiplicity_vs_pt->Fill(mult, jet_pt, weight);
        h_qjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
        h_qjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
        h_qjet_width_vs_pt->Fill(width, jet_pt, weight);
        h_qjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);
        if (matchedJet) h_qjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
        h_qjet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);
      }

      h_jet_flavour->Fill(jet_flav, weight);
      h_jet_genParton_flavour->Fill(abs(thisjet.genPartonFlavor()), weight);
      h_jet_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
      h_jet_genParton_flavour_vs_pt->Fill(abs(thisjet.genPartonFlavor()), jet_pt, weight);
      if (i == 0) {
        h_jet1_genParton_flavour_vs_pt->Fill(abs(thisjet.genPartonFlavor()), jet_pt, weight);
        h_jet1_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
      }
      else if (i == 1) {
        h_jet2_genParton_flavour_vs_pt->Fill(abs(thisjet.genPartonFlavor()), jet_pt, weight);
        h_jet2_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
      }

      h_jet_flavour_vs_eta->Fill(jet_flav, thisjet.eta(), weight);
      h_jet_genParton_flavour_vs_eta->Fill(abs(thisjet.genPartonFlavor()), thisjet.eta(), weight);
    }
    //  Do lambda correlation hists
    // if (jet_pt > 100 && jet_pt < 200) {
    //   h_jet_multiplicity_vs_LHA->Fill(mult, lha, weight);
    //   h_jet_multiplicity_vs_pTD->Fill(mult, ptd, weight);
    //   h_jet_multiplicity_vs_width->Fill(mult, width, weight);
    //   h_jet_multiplicity_vs_thrust->Fill(mult, thrust, weight);
    //   h_jet_LHA_vs_pTD->Fill(lha, ptd, weight);
    //   h_jet_LHA_vs_width->Fill(lha, width, weight);
    //   h_jet_LHA_vs_thrust->Fill(lha, thrust, weight);
    //   h_jet_pTD_vs_width->Fill(ptd, width, weight);
    //   h_jet_pTD_vs_thrust->Fill(ptd, thrust, weight);
    //   h_jet_width_vs_thrust->Fill(width, thrust, weight);

    //   if (abs(jet_flav) == 21) { // gluon jets
    //     h_gjet_multiplicity_vs_LHA->Fill(mult, lha, weight);
    //     h_gjet_multiplicity_vs_pTD->Fill(mult, ptd, weight);
    //     h_gjet_multiplicity_vs_width->Fill(mult, width, weight);
    //     h_gjet_multiplicity_vs_thrust->Fill(mult, thrust, weight);
    //     h_gjet_LHA_vs_pTD->Fill(lha, ptd, weight);
    //     h_gjet_LHA_vs_width->Fill(lha, width, weight);
    //     h_gjet_LHA_vs_thrust->Fill(lha, thrust, weight);
    //     h_gjet_pTD_vs_width->Fill(ptd, width, weight);
    //     h_gjet_pTD_vs_thrust->Fill(ptd, thrust, weight);
    //     h_gjet_width_vs_thrust->Fill(width, thrust, weight);
    //   } else if ((abs(jet_flav) <= 3) && (abs(jet_flav) > 0)){ // uds jets
    //     h_qjet_multiplicity_vs_LHA->Fill(mult, lha, weight);
    //     h_qjet_multiplicity_vs_pTD->Fill(mult, ptd, weight);
    //     h_qjet_multiplicity_vs_width->Fill(mult, width, weight);
    //     h_qjet_multiplicity_vs_thrust->Fill(mult, thrust, weight);
    //     h_qjet_LHA_vs_pTD->Fill(lha, ptd, weight);
    //     h_qjet_LHA_vs_width->Fill(lha, width, weight);
    //     h_qjet_LHA_vs_thrust->Fill(lha, thrust, weight);
    //     h_qjet_pTD_vs_width->Fill(ptd, width, weight);
    //     h_qjet_pTD_vs_thrust->Fill(ptd, thrust, weight);
    //     h_qjet_width_vs_thrust->Fill(width, thrust, weight);
    //   }
    // }
  }

  if (is_mc_) {
    // Fill GenJet hists
    // std::vector<GenJetWithParts>* genjets = event.genjets;
    std::vector<GenParticle>* genparticles = event.genparticles;

    // for (int i = 0; i < useNJets_; i++) {
    //   const GenJetWithParts & thisjet = genjets->at(i);
    int counter = 0;
    for (const auto & thisjet : *genjets) {
      if (!(thisjet.pt() > 30 && fabs(thisjet.eta()) < 2.4))
        continue;

      counter++;
      if (counter > useNJets_)
        break;

      std::vector<GenParticle*> daughters = get_genjet_genparticles(thisjet, genparticles);

      h_genjet_pt->Fill(thisjet.pt(), weight);
      h_genjet_eta->Fill(thisjet.eta(), weight);

      // do special vars according to 1704.03878
      LambdaCalculator<GenParticle> genJetCalc(daughters, jetRadius, thisjet.v4(), false);
      float lha = genJetCalc.getLambda(1, 0.5);
      float ptd = genJetCalc.getLambda(2, 0);
      float width = genJetCalc.getLambda(1, 1);
      float thrust = genJetCalc.getLambda(1, 2);
      uint mult = daughters.size();

      h_genjet_multiplicity->Fill(mult, weight);
      h_genjet_LHA->Fill(lha, weight);
      h_genjet_pTD->Fill(ptd, weight);
      h_genjet_width->Fill(width, weight);
      h_genjet_thrust->Fill(thrust, weight);

      float jet_pt = thisjet.pt();
      h_genjet_multiplicity_vs_pt->Fill(mult, jet_pt, weight);
      h_genjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
      h_genjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
      h_genjet_width_vs_pt->Fill(width, jet_pt, weight);
      h_genjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);

      // if (thisjet.flavor() == 21) { // gluon jets
      //   h_ggenjet_multiplicity->Fill(mult, weight);
      //   h_ggenjet_LHA->Fill(lha / LHA_rescale, weight);
      //   h_ggenjet_pTD->Fill(ptd, weight);
      //   h_ggenjet_width->Fill(width / width_rescale, weight);
      //   h_ggenjet_thrust->Fill(thrust / thrust_rescale, weight);
      // } else if ((abs(thisjet.flavor()) <= 3) && (abs(thisjet.flavor()) > 0)){ // uds jets
      //   h_qgenjet_multiplicity->Fill(mult, weight);
      //   h_qgenjet_LHA->Fill(lha / LHA_rescale, weight);
      //   h_qgenjet_pTD->Fill(ptd, weight);
      //   h_qgenjet_width->Fill(width / width_rescale, weight);
      //   h_qgenjet_thrust->Fill(thrust / thrust_rescale, weight);
      // }
    }
  }
}


/**
 * Get the collection of GenParticle*s for a given GenJet
 */
std::vector<GenParticle*> QGAnalysisHists::get_genjet_genparticles(const GenJetWithParts & jet, std::vector<GenParticle>* genparticles) {
  std::vector<GenParticle*> gp;
  for (const uint i : jet.genparticles_indices()) {
    gp.push_back(&(genparticles->at(i)));
  }
  return gp;
}


/**
 * Get the collection of PFParticle*s for a given Jet
 */
std::vector<PFParticle*> QGAnalysisHists::get_jet_pfparticles(const Jet & jet, std::vector<PFParticle>* pfparticles) {
  std::vector<PFParticle*> pf;
  for (const uint i : jet.daughterIndices()) {
    pf.push_back(&(pfparticles->at(i)));
  }
  return pf;
}


void QGAnalysisHists::shift_neutral_hadron_pfparticles(std::vector<PFParticle*> pfparticles, float direction, float rel_shift) {
  for (auto & itr : pfparticles) {
    if (itr->particleID() == PFParticle::eH0) {
      itr->set_pt(itr->pt() * (1 + (direction * rel_shift)));
    }
  }
}


std::vector<PFParticle*> QGAnalysisHists::create_copy(std::vector<PFParticle*> pfparticles) {
  std::vector<PFParticle*> pfCopy;
  for (const auto & itr: pfparticles) {
    PFParticle * newPf = new PFParticle();
    newPf->set_particleID(itr->particleID());
    newPf->set_puppiWeight(itr->puppiWeight());
    newPf->set_v4(itr->v4());
    newPf->set_charge(itr->charge());
    pfCopy.push_back(newPf);
  }
  return pfCopy;
}

/**
 *
 */
// int QGAnalysisHists::get_jet_flavour(const Jet & jet, std::vector<GenParticle>* genparticles) {
//   cout << "jet:" << endl;
//   cout << jet.eta() << " : " << jet.phi() << endl;
//   if (jet.genParton() != nullptr) {
//     cout << "gp: " << jet.genParton()->eta() << " : " << jet.genParton()->phi() << endl;
//   }
//   for (const auto& gp : *genparticles) {

//     if (deltaR(jet.v4(), gp.v4()) < 0.1) cout << "*** ";
//     cout << gp.eta() << " : " << gp.phi() << " : " << gp.pdgId() << " : " << gp.status() << " : " << deltaR(jet.v4(), gp.v4()) << endl;
//     // if (gp.status)
//   }
//   return 1;
// }

QGAnalysisHists::~QGAnalysisHists(){}


QGJetTrigHists::QGJetTrigHists(Context & ctx, const string & dirname, const std::vector<std::string> & trigNames_):
Hists(ctx, dirname),
trigNames(trigNames_)
{
  int nPtBins = 800;
  float ptMin(0.), ptMax(2000.);
  nPtBins = (int)ptMax;
  int nEtaBins = 50;
  float etaMin(-5), etaMax(5);
  for (const auto & trig: trigNames) {
    hTrigs.push_back(book<TH2F>("pt_vs_eta_" + trig, trig+";Jet p_{T} [GeV];Jet |#eta|", nPtBins, ptMin, ptMax, nEtaBins, etaMin, etaMax));
    trigSels.push_back(TriggerSelection(trig));
  }
  hAll = book<TH2F>("pt_vs_eta_all", "All;Jet p_{T} [GeV];Jet |#eta|", nPtBins, ptMin, ptMax, nEtaBins, etaMin, etaMax);
}

void QGJetTrigHists::fill(const uhh2::Event & event) {
  if (event.jets->size()==0) return;

  auto leadingJet = event.jets->at(0);
  for (uint i=0; i < trigNames.size(); i++) {
    // Add check first to ensure trigger in event, otherwise throws
    auto ti = event.get_trigger_index(trigNames[i]);
    if (event.lookup_trigger_index(ti) && trigSels[i].passes(event)) {
      // cout << "fill: Fired: " << trigNames[i] << endl;
      hTrigs[i]->Fill(leadingJet.pt(), leadingJet.eta(), event.weight);
    }
  }
  hAll->Fill(leadingJet.pt(), leadingJet.eta(), event.weight);

}

QGJetTrigHists::~QGJetTrigHists(){}
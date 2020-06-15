#include "UHH2/QGAnalysis/include/QGAnalysisHists.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisHists::QGAnalysisHists(Context & ctx, const string & dirname,
                                 int useNJets,
                                 bool doGroomed,
                                 const string & selection,
                                 const string & reco_sel_handle_name, const string & gen_sel_handle_name,
                                 const string & reco_jetlambda_handle_name, const string & gen_jetlambda_handle_name
                                 ):
  Hists(ctx, dirname),
  dirname_(dirname),
  useNJets_(useNJets),
  doGroomed_(doGroomed),
  rsp_lowPt_cut_(30.),
  rsp_midPt_cut_(100.),
  rsp_highPt_cut_(250.),
  recoDauPtCut_(1.)
  {

  is_mc_ = ctx.get("dataset_type") == "MC";
  doPuppi_ = (ctx.get("PURemoval") == "PUPPI");

  if (useNJets_ < 0) useNJets_ = 99999; // Do them all
  else if (useNJets_ == 0) throw runtime_error("useNJets should be > 0, or < 0 to use all jets in the event");

  if (selection != "dijet" && selection != "zplusjets") {
    throw runtime_error("selection must be dijet or zplusjets");
  }

  gen_weight_handle = ctx.get_handle<double>("gen_weight");

  // book all histograms here
  int nWeightBins = 25;
  double weightBins [nWeightBins+1] = {1E-4, 1E-3, 1E-2, 1E-1, 1, 10, 100, 1000, 1E4, 2E4, 5E4, 7E4, 1E5, 2E5, 5E5, 7E5, 1E6, 2E6, 3E6, 4E6, 5E6, 6E6, 7E6, 8E6, 9E6, 1E7};

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
  h_jet_eta = book<TH1F>("jet_eta", ";y^{j};", nEtaBins, etaMin, etaMax);

  h_jet_flavour = book<TH1F>("jet_flavour", "jet flavour;PDGID;", 23, -0.5, 22.5);

  int nMultBins = 150;
  h_jet_multiplicity = book<TH1F>("jet_multiplicity", ";# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_jet_multiplicity_charged = book<TH1F>("jet_multiplicity_charged", ";# of charged constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity = book<TH1F>("jet_puppiMultiplicity", ";# of PUPPI-weighted constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_charged = book<TH1F>("jet_puppiMultiplicity_charged", ";# of PUPPI-weighted charged constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);

  int nBins = 100;
  h_jet_LHA = book<TH1F>("jet_LHA", ";LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_jet_pTD = book<TH1F>("jet_pTD", ";p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_jet_width = book<TH1F>("jet_width", ";Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_jet_thrust = book<TH1F>("jet_thrust", ";Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  h_jet_LHA_charged = book<TH1F>("jet_LHA_charged", ";LHA charged (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_jet_pTD_charged = book<TH1F>("jet_pTD_charged", ";p_{T}^{D} charged (#lambda_{0}^{2});", nBins, 0, 1);
  h_jet_width_charged = book<TH1F>("jet_width_charged", ";Width charged (#lambda_{1}^{1});", nBins, 0, 1);
  h_jet_thrust_charged = book<TH1F>("jet_thrust_charged", ";Thrust charged (#lambda_{2}^{1});", nBins, 0, 1);

  // gen-reco response hists
  float rsp_max = 5.;
  int nBinsNormRsp = 500;
  // hist names must always end in _response or _rel_response
  h_jet_multiplicity_response = book<TH2F>("jet_multiplicity_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_charged_response = book<TH2F>("jet_multiplicity_charged_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO), charged only", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_rel_response = book<TH2F>("jet_multiplicity_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_multiplicity_charged_rel_response = book<TH2F>("jet_multiplicity_charged_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);

  h_jet_multiplicity_lowPt_response = book<TH2F>("jet_multiplicity_lowPt_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_charged_lowPt_response = book<TH2F>("jet_multiplicity_charged_lowPt_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO), charged only", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_lowPt_rel_response = book<TH2F>("jet_multiplicity_lowPt_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_multiplicity_charged_lowPt_rel_response = book<TH2F>("jet_multiplicity_charged_lowPt_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);

  h_jet_multiplicity_midPt_response = book<TH2F>("jet_multiplicity_midPt_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_charged_midPt_response = book<TH2F>("jet_multiplicity_charged_midPt_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO), charged only", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_midPt_rel_response = book<TH2F>("jet_multiplicity_midPt_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_multiplicity_charged_midPt_rel_response = book<TH2F>("jet_multiplicity_charged_midPt_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);

  h_jet_multiplicity_highPt_response = book<TH2F>("jet_multiplicity_highPt_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_charged_highPt_response = book<TH2F>("jet_multiplicity_charged_highPt_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO), charged only", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_multiplicity_highPt_rel_response = book<TH2F>("jet_multiplicity_highPt_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN);# of constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_multiplicity_charged_highPt_rel_response = book<TH2F>("jet_multiplicity_charged_highPt_rel_response", ";# of constituents (#lambda_{0}^{0}) (GEN, charged only);# of constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);

  // PUPPI multiplicity
  h_jet_puppiMultiplicity_response = book<TH2F>("jet_puppiMultiplicity_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_charged_response = book<TH2F>("jet_puppiMultiplicity_charged_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO, charged only)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_rel_response = book<TH2F>("jet_puppiMultiplicity_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_puppiMultiplicity_charged_rel_response = book<TH2F>("jet_puppiMultiplicity_charged_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);

  h_jet_puppiMultiplicity_lowPt_response = book<TH2F>("jet_puppiMultiplicity_lowPt_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_charged_lowPt_response = book<TH2F>("jet_puppiMultiplicity_charged_lowPt_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO, charged only)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_lowPt_rel_response = book<TH2F>("jet_puppiMultiplicity_lowPt_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_puppiMultiplicity_charged_lowPt_rel_response = book<TH2F>("jet_puppiMultiplicity_charged_lowPt_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);

  h_jet_puppiMultiplicity_midPt_response = book<TH2F>("jet_puppiMultiplicity_midPt_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_charged_midPt_response = book<TH2F>("jet_puppiMultiplicity_charged_midPt_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO, charged only)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_midPt_rel_response = book<TH2F>("jet_puppiMultiplicity_midPt_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_puppiMultiplicity_charged_midPt_rel_response = book<TH2F>("jet_puppiMultiplicity_charged_midPt_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);

  h_jet_puppiMultiplicity_highPt_response = book<TH2F>("jet_puppiMultiplicity_highPt_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_charged_highPt_response = book<TH2F>("jet_puppiMultiplicity_charged_highPt_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO, charged only)", nMultBins, 0, nMultBins, nMultBins, 0, nMultBins);
  h_jet_puppiMultiplicity_highPt_rel_response = book<TH2F>("jet_puppiMultiplicity_highPt_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);
  h_jet_puppiMultiplicity_charged_highPt_rel_response = book<TH2F>("jet_puppiMultiplicity_charged_highPt_rel_response", ";# of PUPPI-weighted constituents (#lambda_{0}^{0}) (GEN, charged only);# of PUPPI-weighted constituents (#lambda_{0}^{0}) (RECO / GEN, charged only)", nMultBins, 0, nMultBins, nBinsNormRsp, 0, rsp_max);

  // LHA
  h_jet_LHA_response = book<TH2F>("jet_LHA_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_charged_response = book<TH2F>("jet_LHA_charged_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_rel_response = book<TH2F>("jet_LHA_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_LHA_charged_rel_response = book<TH2F>("jet_LHA_charged_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_LHA_lowPt_response = book<TH2F>("jet_LHA_lowPt_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_charged_lowPt_response = book<TH2F>("jet_LHA_charged_lowPt_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_lowPt_rel_response = book<TH2F>("jet_LHA_lowPt_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_LHA_charged_lowPt_rel_response = book<TH2F>("jet_LHA_charged_lowPt_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_LHA_midPt_response = book<TH2F>("jet_LHA_midPt_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_charged_midPt_response = book<TH2F>("jet_LHA_charged_midPt_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_midPt_rel_response = book<TH2F>("jet_LHA_midPt_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_LHA_charged_midPt_rel_response = book<TH2F>("jet_LHA_charged_midPt_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_LHA_highPt_response = book<TH2F>("jet_LHA_highPt_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_charged_highPt_response = book<TH2F>("jet_LHA_charged_highPt_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_LHA_highPt_rel_response = book<TH2F>("jet_LHA_highPt_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN);LHA (#lambda_{0.5}^{1}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_LHA_charged_highPt_rel_response = book<TH2F>("jet_LHA_charged_highPt_rel_response", ";LHA (#lambda_{0.5}^{1}) (GEN, charged only);LHA (#lambda_{0.5}^{1}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);

  // pTD
  h_jet_pTD_response = book<TH2F>("jet_pTD_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_charged_response = book<TH2F>("jet_pTD_charged_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN, charged);p_{T}^{D} (#lambda_{0}^{2}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_rel_response = book<TH2F>("jet_pTD_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_pTD_charged_rel_response = book<TH2F>("jet_pTD_charged_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_pTD_lowPt_response = book<TH2F>("jet_pTD_lowPt_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_charged_lowPt_response = book<TH2F>("jet_pTD_charged_lowPt_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN, charged);p_{T}^{D} (#lambda_{0}^{2}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_lowPt_rel_response = book<TH2F>("jet_pTD_lowPt_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_pTD_charged_lowPt_rel_response = book<TH2F>("jet_pTD_charged_lowPt_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_pTD_midPt_response = book<TH2F>("jet_pTD_midPt_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_charged_midPt_response = book<TH2F>("jet_pTD_charged_midPt_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN, charged);p_{T}^{D} (#lambda_{0}^{2}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_midPt_rel_response = book<TH2F>("jet_pTD_midPt_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_pTD_charged_midPt_rel_response = book<TH2F>("jet_pTD_charged_midPt_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_pTD_highPt_response = book<TH2F>("jet_pTD_highPt_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_charged_highPt_response = book<TH2F>("jet_pTD_charged_highPt_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN, charged);p_{T}^{D} (#lambda_{0}^{2}) (RECO, charged only)", nBins, 0, 1, nBins, 0, 1);
  h_jet_pTD_highPt_rel_response = book<TH2F>("jet_pTD_highPt_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_pTD_charged_highPt_rel_response = book<TH2F>("jet_pTD_charged_highPt_rel_response", ";p_{T}^{D} (#lambda_{0}^{2}) (GEN);p_{T}^{D} (#lambda_{0}^{2}) (RECO / GEN, charged only)", nBins, 0, 1, nBinsNormRsp, 0, rsp_max);

  // width
  // finer binning for these 2 variables as they can basically just peak @ 0.
  int nBinsFine = nBins*2;
  h_jet_width_response = book<TH2F>("jet_width_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_width_charged_response = book<TH2F>("jet_width_charged_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO, charged only)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_width_rel_response = book<TH2F>("jet_width_rel_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO / GEN)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_width_charged_rel_response = book<TH2F>("jet_width_charged_rel_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO / GEN, charged only)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_width_lowPt_response = book<TH2F>("jet_width_lowPt_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_width_charged_lowPt_response = book<TH2F>("jet_width_charged_lowPt_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO, charged only)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_width_lowPt_rel_response = book<TH2F>("jet_width_lowPt_rel_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO / GEN)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_width_charged_lowPt_rel_response = book<TH2F>("jet_width_charged_lowPt_rel_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO / GEN, charged only)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_width_midPt_response = book<TH2F>("jet_width_midPt_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_width_charged_midPt_response = book<TH2F>("jet_width_charged_midPt_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO, charged only)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_width_midPt_rel_response = book<TH2F>("jet_width_midPt_rel_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO / GEN)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_width_charged_midPt_rel_response = book<TH2F>("jet_width_charged_midPt_rel_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO / GEN, charged only)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_width_highPt_response = book<TH2F>("jet_width_highPt_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_width_charged_highPt_response = book<TH2F>("jet_width_charged_highPt_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO, charged only)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_width_highPt_rel_response = book<TH2F>("jet_width_highPt_rel_response", ";Width (#lambda_{1}^{1}) (GEN);Width (#lambda_{1}^{1}) (RECO / GEN)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_width_charged_highPt_rel_response = book<TH2F>("jet_width_charged_highPt_rel_response", ";Width (#lambda_{1}^{1}) (GEN, charged);Width (#lambda_{1}^{1}) (RECO / GEN, charged only)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);

  // thrust
  h_jet_thrust_response = book<TH2F>("jet_thrust_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_thrust_charged_response = book<TH2F>("jet_thrust_charged_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO, charged only)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_thrust_rel_response = book<TH2F>("jet_thrust_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO / GEN)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_thrust_charged_rel_response = book<TH2F>("jet_thrust_charged_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO / GEN, charged only)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_thrust_lowPt_response = book<TH2F>("jet_thrust_lowPt_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_thrust_charged_lowPt_response = book<TH2F>("jet_thrust_charged_lowPt_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO, charged only)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_thrust_lowPt_rel_response = book<TH2F>("jet_thrust_lowPt_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO / GEN)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_thrust_charged_lowPt_rel_response = book<TH2F>("jet_thrust_charged_lowPt_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO / GEN, charged only)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_thrust_midPt_response = book<TH2F>("jet_thrust_midPt_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_thrust_charged_midPt_response = book<TH2F>("jet_thrust_charged_midPt_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO, charged only)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_thrust_midPt_rel_response = book<TH2F>("jet_thrust_midPt_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO / GEN)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_thrust_charged_midPt_rel_response = book<TH2F>("jet_thrust_charged_midPt_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO / GEN, charged only)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);

  h_jet_thrust_highPt_response = book<TH2F>("jet_thrust_highPt_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_thrust_charged_highPt_response = book<TH2F>("jet_thrust_charged_highPt_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO, charged only)", nBinsFine, 0, 1, nBinsFine, 0, 1);
  h_jet_thrust_highPt_rel_response = book<TH2F>("jet_thrust_highPt_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN);Thrust (#lambda_{2}^{1}) (RECO / GEN)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);
  h_jet_thrust_charged_highPt_rel_response = book<TH2F>("jet_thrust_charged_highPt_rel_response", ";Thrust (#lambda_{2}^{1}) (GEN, charged);Thrust (#lambda_{2}^{1}) (RECO / GEN, charged only)", nBinsFine, 0, 1, nBinsNormRsp, 0, rsp_max);

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

  h_jet_multiplicity_charged_vs_pt = book<TH2F>("jet_multiplicity_charged_vs_pt", ";# of charged constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_jet_LHA_charged_vs_pt = book<TH2F>("jet_LHA_charged_vs_pt", ";LHA charged (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_pTD_charged_vs_pt = book<TH2F>("jet_pTD_charged_vs_pt", ";p_{T}^{D} charged (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_width_charged_vs_pt = book<TH2F>("jet_width_charged_vs_pt", ";Width charged (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_thrust_charged_vs_pt = book<TH2F>("jet_thrust_charged_vs_pt", ";Thrust charged (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  int rspNbins = 500;
  float rspMax = 5;
  h_jet_response_vs_genjet_pt = book<TH2F>("jet_response_vs_genjet_pt", ";response;GenJet p_{T} [GeV]", rspNbins, 0, rspMax, nPtBins, ptMin, ptMax);

  h_jet_flavour_vs_pt = book<TH2F>("jet_flavour_vs_pt", "jet flavour;PDGID;Jet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet1_flavour_vs_pt = book<TH2F>("jet1_flavour_vs_pt", "jet1 flavour;PDGID;Jet1 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet2_flavour_vs_pt = book<TH2F>("jet2_flavour_vs_pt", "jet2 flavour;PDGID;Jet2 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);

  h_genjet_flavour_vs_pt = book<TH2F>("genjet_flavour_vs_pt", "genjet flavour;PDGID;GenJet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_genjet1_flavour_vs_pt = book<TH2F>("genjet1_flavour_vs_pt", "genjet1 flavour;PDGID;GenJet1 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_genjet2_flavour_vs_pt = book<TH2F>("genjet2_flavour_vs_pt", "genjet2 flavour;PDGID;GenJet2 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);

  h_jet_flavour_vs_eta = book<TH2F>("jet_flavour_vs_eta", "jet flavour;PDGID;Jet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);

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
  h_jet_puppiMultiplicity_charged_vs_pt = book<TH2F>("jet_puppiMultiplicity_charged_vs_pt", ";# of charged constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
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
  h_genjet_eta = book<TH1F>("genjet_eta", ";y^{j};", nEtaBins, etaMin, etaMax);

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

  if (is_mc_) {
    // genJets_handle = ctx.get_handle< std::vector<GenJet> > ("GoodGenJets");
    genJetsLambda_handle = ctx.get_handle< std::vector<GenJetLambdaBundle> > (gen_jetlambda_handle_name);
    pass_gen_handle = ctx.get_handle<bool> (gen_sel_handle_name);
  }

  jetsLambda_handle = ctx.get_handle< std::vector<JetLambdaBundle> > (reco_jetlambda_handle_name);

  pass_reco_handle = ctx.get_handle<bool> (reco_sel_handle_name);
}


void QGAnalysisHists::fill(const Event & event){
  double weight = event.weight;

  // extract the separate gen & reco weight components, needed for TUnfold
  double gen_weight = event.get(gen_weight_handle);
  double reco_weight = weight / gen_weight;

  h_weights->Fill(weight);

  // const std::vector<GenJet> * genjets = nullptr;
  if (is_mc_) {
    if (event.genInfo->binningValues().size() > 0)
      h_pthat_vs_weight->Fill(weight, event.genInfo->binningValues()[0]);
    // genjets = &event.get(genJets_handle);
  }
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

  std::vector<JetLambdaBundle> jetLambdas = event.get(jetsLambda_handle);

  // Fill reco jet hists
  // ---------------------------------------------------------------------------
  // At this point, all jet filtering etc should have already been performed
  if (passReco){
    if (useNJets_ > (int) jetLambdas.size()) {
      cout << "useNJets_: " << useNJets_ << endl;
      cout << "jetLambdas: " << jetLambdas.size() << endl;
      throw runtime_error("useNJets_ > jetLambdas.size()");
    }

    for (int i = 0; i < useNJets_; i++) {
      const Jet & thisjet = jetLambdas.at(i).jet;
      LambdaCalculator<PFParticle> recoJetCalc = jetLambdas.at(i).getLambdaCalculator(false, doGroomed_);
      LambdaCalculator<PFParticle> recoJetCalcCharged = jetLambdas.at(i).getLambdaCalculator(true, doGroomed_);

      // To account for Lambda Calculators with 1 constituent from charged-only or grooming,
      // (which isn't really a jet) we test per jet, and treat it otherwise as a fail
      // TODO: should this be done outside?
      bool thisPassReco = (passReco && recoJetCalc.constits().size() > 1);
      bool thisPassRecoCharged = (passReco && recoJetCalcCharged.constits().size() > 1);
      if (!thisPassReco && !thisPassRecoCharged) continue;

      float jet_pt = thisjet.pt();
      h_weights_vs_pt->Fill(weight, jet_pt);
      if (is_mc_) {
        if (event.genInfo->binningValues().size() > 0)
          h_pthat_vs_jet_pt->Fill(thisjet.pt(), event.genInfo->binningValues()[0]);
      }
      h_jet_pt_unweighted->Fill(jet_pt);
      h_jet_pt->Fill(jet_pt, weight);
      h_jet_eta->Fill(thisjet.Rapidity(), weight);

      float puppiMult(0.), mult(0.), lha(0.), ptd(0.), width(0.), thrust(0.);

      if (thisPassReco) {
        puppiMult = recoJetCalc.getLambda(0, 0);
        mult = recoJetCalc.constits().size();
        lha = recoJetCalc.getLambda(1, 0.5);
        ptd = recoJetCalc.getLambda(2, 0);
        width = recoJetCalc.getLambda(1, 1);
        thrust = recoJetCalc.getLambda(1, 2);

        if (ptd > 1){
          cout << "ptd > 1: " << ptd << endl;
          cout << "ptSum: " << recoJetCalc.getPtSum() << endl;
          cout << " jet " << thisjet.pt() << " : " << thisjet.Rapidity() << " : " << thisjet.phi() << endl;
          cout << "calc constits: " << recoJetCalc.constits().size() << endl;
          cout << "constits: " << endl;
          float myptd = 0;
          float ptsum = 0;
          for (const auto & cind : thisjet.pfcand_indexs()) {
            const PFParticle & pf = event.pfparticles->at(cind);
            if (pf.pt() < 1) continue;
            cout << pf.pt() << " : " << pf.Rapidity() << " : " << pf.phi() << " : " << pf.puppiWeight() << endl;
            ptsum += (pf.pt() * pf.puppiWeight());
            myptd += pow(pf.pt() * pf.puppiWeight(), 2);
          }
          myptd /= pow(ptsum, 2);
          cout << "myptd: " << myptd << endl;
          cout << "myptsum: " << ptsum << endl;
        }

        h_jet_puppiMultiplicity->Fill(puppiMult, weight);
        h_jet_multiplicity->Fill(mult, weight);
        h_jet_LHA->Fill(lha, weight);
        h_jet_pTD->Fill(ptd, weight);
        h_jet_width->Fill(width, weight);
        h_jet_thrust->Fill(thrust, weight);

        h_jet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);
        h_jet_multiplicity_vs_pt->Fill(mult, jet_pt, weight);
        h_jet_LHA_vs_pt->Fill(lha, jet_pt, weight);
        h_jet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
        h_jet_width_vs_pt->Fill(width, jet_pt, weight);
        h_jet_thrust_vs_pt->Fill(thrust, jet_pt, weight);
      }

      float puppiMult_charged(0.), mult_charged(0.), lha_charged(0.), ptd_charged(0.), width_charged(0.), thrust_charged(0.);

      if (thisPassRecoCharged) {
        puppiMult_charged = recoJetCalcCharged.getLambda(0, 0);
        mult_charged = recoJetCalcCharged.constits().size();
        lha_charged = recoJetCalcCharged.getLambda(1, 0.5);
        ptd_charged = recoJetCalcCharged.getLambda(2, 0);
        width_charged = recoJetCalcCharged.getLambda(1, 1);
        thrust_charged = recoJetCalcCharged.getLambda(1, 2);

        h_jet_puppiMultiplicity_charged->Fill(puppiMult_charged, weight);
        h_jet_multiplicity_charged->Fill(mult_charged, weight);
        h_jet_LHA_charged->Fill(lha_charged, weight);
        h_jet_pTD_charged->Fill(ptd_charged, weight);
        h_jet_width_charged->Fill(width_charged, weight);
        h_jet_thrust_charged->Fill(thrust_charged, weight);

        h_jet_puppiMultiplicity_charged_vs_pt->Fill(puppiMult_charged, jet_pt, weight);
        h_jet_multiplicity_charged_vs_pt->Fill(mult_charged, jet_pt, weight);
        h_jet_LHA_charged_vs_pt->Fill(lha_charged, jet_pt, weight);
        h_jet_pTD_charged_vs_pt->Fill(ptd_charged, jet_pt, weight);
        h_jet_width_charged_vs_pt->Fill(width_charged, jet_pt, weight);
        h_jet_thrust_charged_vs_pt->Fill(thrust_charged, jet_pt, weight);
      }

      if (is_mc_) { // TODO: should there be a && passGen?
        // Store variables for matched GenJet
        bool matchedGenJet = (thisjet.genjet_index() > -1 && thisjet.genjet_index() < 5); // put upper limit to avoid weird matches
        float genjet_pt = -1.;
        float genjet_flav = -1.;
        float response = -1.;

        if (matchedGenJet) { // we don't care about passing the GEN selection, just looking for a match
          int thisInd = thisjet.genjet_index();
          if (thisInd >= (int) genjetLambdas->size()) {
            cout << "WARNING: wanted genjet_index " << thisInd << " but only have " << genjetLambdas->size() << " in genjetLambdas" << endl;
          }
          const GenJet & genjet = genjetLambdas->at(thisInd).jet;
          LambdaCalculator<GenParticle> matchedGenJetCalc = genjetLambdas->at(thisInd).getLambdaCalculator(false, doGroomed_);
          LambdaCalculator<GenParticle> matchedGenJetCalcCharged = genjetLambdas->at(thisInd).getLambdaCalculator(true, doGroomed_);

          bool thisPassGen = passGen && (matchedGenJetCalc.constits().size() > 1);
          bool thisPassGenCharged = passGen && (matchedGenJetCalcCharged.constits().size() > 1);

          genjet_pt = genjet.pt();
          genjet_flav = abs(genjet.partonFlavour());
          response = jet_pt/genjet_pt;
          h_jet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);

          if (thisPassGen && thisPassReco) {
            // Response hists over all pt, and high/low pt split ones
            float gen_mult = matchedGenJetCalc.getLambda(0, 0);
            fill_lambda_rsp_hists(mult, gen_mult, weight,
              h_jet_multiplicity_response, h_jet_multiplicity_rel_response,
              jet_pt,
              h_jet_multiplicity_lowPt_response, h_jet_multiplicity_midPt_response, h_jet_multiplicity_highPt_response,
              h_jet_multiplicity_lowPt_rel_response, h_jet_multiplicity_midPt_rel_response, h_jet_multiplicity_highPt_rel_response);

            fill_lambda_rsp_hists(puppiMult, gen_mult, weight,
              h_jet_puppiMultiplicity_response, h_jet_puppiMultiplicity_rel_response,
              jet_pt,
              h_jet_puppiMultiplicity_lowPt_response, h_jet_puppiMultiplicity_midPt_response, h_jet_puppiMultiplicity_highPt_response,
              h_jet_puppiMultiplicity_lowPt_rel_response, h_jet_puppiMultiplicity_midPt_rel_response, h_jet_puppiMultiplicity_highPt_rel_response);

            float gen_lha = matchedGenJetCalc.getLambda(1, 0.5);
            fill_lambda_rsp_hists(lha, gen_lha, weight,
              h_jet_LHA_response, h_jet_LHA_rel_response,
              jet_pt,
              h_jet_LHA_lowPt_response, h_jet_LHA_midPt_response, h_jet_LHA_highPt_response,
              h_jet_LHA_lowPt_rel_response, h_jet_LHA_midPt_rel_response, h_jet_LHA_highPt_rel_response);

            float gen_ptd = matchedGenJetCalc.getLambda(2, 0);
            fill_lambda_rsp_hists(ptd, gen_ptd, weight,
              h_jet_pTD_response, h_jet_pTD_rel_response,
              jet_pt,
              h_jet_pTD_lowPt_response, h_jet_pTD_midPt_response, h_jet_pTD_highPt_response,
              h_jet_pTD_lowPt_rel_response, h_jet_pTD_midPt_rel_response, h_jet_pTD_highPt_rel_response);

            float gen_width = matchedGenJetCalc.getLambda(1, 1);
            fill_lambda_rsp_hists(width, gen_width, weight,
              h_jet_width_response, h_jet_width_rel_response,
              jet_pt,
              h_jet_width_lowPt_response, h_jet_width_midPt_response, h_jet_width_highPt_response,
              h_jet_width_lowPt_rel_response, h_jet_width_midPt_rel_response, h_jet_width_highPt_rel_response);

            float gen_thrust = matchedGenJetCalc.getLambda(1, 2);
            fill_lambda_rsp_hists(thrust, gen_thrust, weight,
              h_jet_thrust_response, h_jet_thrust_rel_response,
              jet_pt,
              h_jet_thrust_lowPt_response, h_jet_thrust_midPt_response, h_jet_thrust_highPt_response,
              h_jet_thrust_lowPt_rel_response, h_jet_thrust_midPt_rel_response, h_jet_thrust_highPt_rel_response);
          }

          if (thisPassGenCharged && thisPassRecoCharged) {
            // Do charged-only daughters version
            float gen_mult_charged = matchedGenJetCalcCharged.getLambda(0, 0);
            fill_lambda_rsp_hists(mult_charged, gen_mult_charged, weight,
              h_jet_multiplicity_charged_response, h_jet_multiplicity_charged_rel_response,
              jet_pt,
              h_jet_multiplicity_charged_lowPt_response, h_jet_multiplicity_charged_midPt_response, h_jet_multiplicity_charged_highPt_response,
              h_jet_multiplicity_charged_lowPt_rel_response, h_jet_multiplicity_charged_midPt_rel_response, h_jet_multiplicity_charged_highPt_rel_response);

            fill_lambda_rsp_hists(puppiMult_charged, gen_mult_charged, weight,
              h_jet_puppiMultiplicity_charged_response, h_jet_puppiMultiplicity_charged_rel_response,
              jet_pt,
              h_jet_puppiMultiplicity_charged_lowPt_response, h_jet_puppiMultiplicity_charged_midPt_response, h_jet_puppiMultiplicity_charged_highPt_response,
              h_jet_puppiMultiplicity_charged_lowPt_rel_response, h_jet_puppiMultiplicity_charged_midPt_rel_response, h_jet_puppiMultiplicity_charged_highPt_rel_response);

            float gen_lha_charged = matchedGenJetCalcCharged.getLambda(1, 0.5);
            fill_lambda_rsp_hists(lha_charged, gen_lha_charged, weight,
              h_jet_LHA_charged_response, h_jet_LHA_charged_rel_response,
              jet_pt,
              h_jet_LHA_charged_lowPt_response, h_jet_LHA_charged_midPt_response, h_jet_LHA_charged_highPt_response,
              h_jet_LHA_charged_lowPt_rel_response, h_jet_LHA_charged_midPt_rel_response, h_jet_LHA_charged_highPt_rel_response);

            float gen_ptd_charged = matchedGenJetCalcCharged.getLambda(2, 0);
            fill_lambda_rsp_hists(ptd_charged, gen_ptd_charged, weight,
              h_jet_pTD_charged_response, h_jet_pTD_charged_rel_response,
              jet_pt,
              h_jet_pTD_charged_lowPt_response, h_jet_pTD_charged_midPt_response, h_jet_pTD_charged_highPt_response,
              h_jet_pTD_charged_lowPt_rel_response, h_jet_pTD_charged_midPt_rel_response, h_jet_pTD_charged_highPt_rel_response);

            float gen_width_charged = matchedGenJetCalcCharged.getLambda(1, 1);
            fill_lambda_rsp_hists(width_charged, gen_width_charged, weight,
              h_jet_width_charged_response, h_jet_width_charged_rel_response,
              jet_pt,
              h_jet_width_charged_lowPt_response, h_jet_width_charged_midPt_response, h_jet_width_charged_highPt_response,
              h_jet_width_charged_lowPt_rel_response, h_jet_width_charged_midPt_rel_response, h_jet_width_charged_highPt_rel_response);

            float gen_thrust_charged = matchedGenJetCalcCharged.getLambda(1, 2);
            fill_lambda_rsp_hists(thrust_charged, gen_thrust_charged, weight,
              h_jet_thrust_charged_response, h_jet_thrust_charged_rel_response,
              jet_pt,
              h_jet_thrust_charged_lowPt_response, h_jet_thrust_charged_midPt_response, h_jet_thrust_charged_highPt_response,
              h_jet_thrust_charged_lowPt_rel_response, h_jet_thrust_charged_midPt_rel_response, h_jet_thrust_charged_highPt_rel_response);
          }
        } // end if matchedgenjet

        // int jet_flav = get_jet_flavour(thisjet, event.genparticles);
        int jet_flav = abs(thisjet.partonFlavour());

        // Split by actual jet flavour - these only make sense for MC
        if (thisPassReco) {
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
            if (matchedGenJet) h_gjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
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
            if (matchedGenJet) h_qjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
            h_qjet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);
          }

          h_jet_flavour->Fill(jet_flav, weight);
          h_jet_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
          if (i == 0) {
            h_jet1_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
          }
          else if (i == 1) {
            h_jet2_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
          }
          if (matchedGenJet) {
            h_genjet_flavour_vs_pt->Fill(genjet_flav, genjet_pt, weight);
            if (i == 0) {
              h_genjet1_flavour_vs_pt->Fill(genjet_flav, genjet_pt, weight);
            }
            else if (i == 1) {
              h_genjet2_flavour_vs_pt->Fill(genjet_flav, genjet_pt, weight);
            }
          }

          h_jet_flavour_vs_eta->Fill(jet_flav, thisjet.Rapidity(), weight);
        }
      } // end is_mc_

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
    } // end loop over reco jets
  } //end if passReco

  // Fill GenJet hists
  // ---------------------------------------------------------------------------
  if (is_mc_ && passGen) { // note that there may be an implicit passReco as well from where this is called in the main module
    if (useNJets_ > (int) genjetLambdas->size()) {
      cout << "useNJets_: " << useNJets_ << endl;
      cout << "genjetLambdas: " << genjetLambdas->size() << endl;
      throw runtime_error("useNJets_ > genjetLambdas.size()");
    }

    for (int i = 0; i < useNJets_; i++) {
      const GenJet & thisjet = genjetLambdas->at(i).jet;
      LambdaCalculator<GenParticle> genJetCalc = genjetLambdas->at(i).getLambdaCalculator(false, doGroomed_);
      // FIXME check this corresponds to same jet as normal lambdas?
      LambdaCalculator<GenParticle> genJetCalcCharged = genjetLambdas->at(i).getLambdaCalculator(true, doGroomed_);

      float genjet_pt = thisjet.pt();

      h_genjet_pt->Fill(genjet_pt, gen_weight);
      h_genjet_eta->Fill(thisjet.Rapidity(), gen_weight);

      float gen_mult = genJetCalc.getLambda(0, 0);
      float gen_lha = genJetCalc.getLambda(1, 0.5);
      float gen_ptd = genJetCalc.getLambda(2, 0);
      float gen_width = genJetCalc.getLambda(1, 1);
      float gen_thrust = genJetCalc.getLambda(1, 2);

      h_genjet_multiplicity->Fill(gen_mult, gen_weight);
      h_genjet_LHA->Fill(gen_lha, gen_weight);
      h_genjet_pTD->Fill(gen_ptd, gen_weight);
      h_genjet_width->Fill(gen_width, gen_weight);
      h_genjet_thrust->Fill(gen_thrust, gen_weight);

      h_genjet_multiplicity_vs_pt->Fill(gen_mult, genjet_pt, gen_weight);
      h_genjet_LHA_vs_pt->Fill(gen_lha, genjet_pt, gen_weight);
      h_genjet_pTD_vs_pt->Fill(gen_ptd, genjet_pt, gen_weight);
      h_genjet_width_vs_pt->Fill(gen_width, genjet_pt, gen_weight);
      h_genjet_thrust_vs_pt->Fill(gen_thrust, genjet_pt, gen_weight);

      // if (thisjet.partonFlavour() == PDGID::GLUON) { // gluon jets
      //   h_ggenjet_multiplicity->Fill(gen_mult, gen_weight);
      //   h_ggenjet_LHA->Fill(gen_lha, gen_weight);
      //   h_ggenjet_pTD->Fill(gen_ptd, gen_weight);
      //   h_ggenjet_width->Fill(gen_width, gen_weight);
      //   h_ggenjet_thrust->Fill(gen_thrust, gen_weight);
      // } else if ((abs(thisjet.partonFlavour()) <= PDGID::STRANGE_QUARK) && (abs(thisjet.partonFlavour()) > PDGID::UNKNOWN)){ // uds jets
      //   h_qgenjet_multiplicity->Fill(gen_mult, gen_weight);
      //   h_qgenjet_LHA->Fill(gen_lha, gen_weight);
      //   h_qgenjet_pTD->Fill(gen_ptd, gen_weight);
      //   h_qgenjet_width->Fill(gen_width, gen_weight);
      //   h_qgenjet_thrust->Fill(gen_thrust, gen_weight);
      // }

    }
  }
}


void QGAnalysisHists::fill_lambda_rsp_hists(float reco_val, float gen_val, float weight,
                                            TH2F * response, TH2F * rel_response,
                                            float jet_pt,
                                            TH2F * response_lowPt, TH2F * response_midPt, TH2F * response_highPt,
                                            TH2F * rel_response_lowPt, TH2F * rel_response_midPt, TH2F * rel_response_highPt) {
  float rsp = reco_val / gen_val; // TODO handle infs
  response->Fill(gen_val, reco_val, weight);
  rel_response->Fill(gen_val, rsp, weight);
  if (response_lowPt != nullptr && response_midPt != nullptr && response_highPt != nullptr) {
    if (jet_pt > rsp_highPt_cut_) {
      response_highPt->Fill(gen_val, reco_val, weight);
    } else if (jet_pt > rsp_midPt_cut_) {
      response_midPt->Fill(gen_val, reco_val, weight);
    } else if (jet_pt > rsp_lowPt_cut_) {
      response_lowPt->Fill(gen_val, reco_val, weight);
    }
  }
  if (rel_response_lowPt != nullptr && rel_response_midPt != nullptr && rel_response_highPt != nullptr) {
    if (jet_pt > rsp_highPt_cut_) {
      rel_response_highPt->Fill(gen_val, rsp, weight);
    } else if (jet_pt > rsp_midPt_cut_) {
      rel_response_midPt->Fill(gen_val, rsp, weight);
    } else if (jet_pt > rsp_lowPt_cut_){
      rel_response_lowPt->Fill(gen_val, rsp, weight);
    }
  }
}


/**
 * Get the collection of GenParticle*s for a given GenJet
 */
// std::vector<GenParticle*> QGAnalysisHists::get_genjet_genparticles(const GenJet & jet, std::vector<GenParticle>* genparticles) {
//   std::vector<GenParticle*> gp;
//   for (const uint i : jet.genparticles_indices()) {
//     gp.push_back(&(genparticles->at(i)));
//   }
//   return gp;
// }


/**
 * Get the collection of PFParticle*s for a given Jet
 */
std::vector<PFParticle*> QGAnalysisHists::get_jet_pfparticles(const Jet & jet, std::vector<PFParticle>* pfparticles) {
  std::vector<PFParticle*> pf;
  for (const uint i : jet.pfcand_indexs()) {
    pf.push_back(&(pfparticles->at(i)));
  }
  return pf;
}


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
    hTrigs.push_back(book<TH2F>("pt_vs_eta_" + trig, trig+";Jet p_{T} [GeV];Jet |y|", nPtBins, ptMin, ptMax, nEtaBins, etaMin, etaMax));
    trigSels.push_back(TriggerSelection(trig));
  }
  hAll = book<TH2F>("pt_vs_eta_all", "All;Jet p_{T} [GeV];Jet |y|", nPtBins, ptMin, ptMax, nEtaBins, etaMin, etaMax);
}

void QGJetTrigHists::fill(const uhh2::Event & event) {
  if (event.jets->size()==0) return;

  auto leadingJet = event.jets->at(0);
  for (uint i=0; i < trigNames.size(); i++) {
    // Add check first to ensure trigger in event, otherwise throws
    auto ti = event.get_trigger_index(trigNames[i]);
    if (event.lookup_trigger_index(ti) && trigSels[i].passes(event)) {
      // cout << "fill: Fired: " << trigNames[i] << endl;
      hTrigs[i]->Fill(leadingJet.pt(), leadingJet.Rapidity(), event.weight);
    }
  }
  hAll->Fill(leadingJet.pt(), leadingJet.Rapidity(), event.weight);

}

QGJetTrigHists::~QGJetTrigHists(){}
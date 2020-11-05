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
                                 bool useStatus23Flavour,
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
  rsp_highPt2_cut_(500.),
  useStatus23Flavour_(useStatus23Flavour),
  N_PARTONS_MAX(4)
  {

  string jetCone = ctx.get("JetCone", "AK4");
  string pu_removal = ctx.get("PURemoval", "CHS");
  if (pu_removal != "CHS" && pu_removal != "PUPPI") {
      throw runtime_error("Only PURemoval == CHS, PUPPI supported for now");
  }
  jetRadius = get_jet_radius(jetCone);

  is_mc_ = ctx.get("dataset_type") == "MC";
  doPuppi_ = (ctx.get("PURemoval") == "PUPPI");

  if (useNJets_ < 0) useNJets_ = 99999; // Do them all
  else if (useNJets_ == 0) throw runtime_error("useNJets should be > 0, or < 0 to use all jets in the event");

  if (selection != "dijet" && selection != "zplusjets") {
    throw runtime_error("selection must be dijet or zplusjets");
  }

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
  float etaMin(-2.5), etaMax(2.5);
  h_jet_y = book<TH1F>("jet_y", ";y^{j};", nEtaBins, etaMin, etaMax);
  h_jet_eta = book<TH1F>("jet_eta", ";#eta^{j};", nEtaBins, etaMin, etaMax);

  h_jet_flavour = book<TH1F>("jet_flavour", "jet flavour;PDGID;", 23, -0.5, 22.5);

  int nMultBins = 150;

  int nBins = 100;

  if (is_mc_) {
    // gen-reco response hists
    float rsp_max = 5.;
    int nBinsNormRsp = 500;
    // hist names must always end in _response or _rel_response

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
  }

  // 2D versions vs PT
  // -----------------
  // All flavs
  h_jet_puppiMultiplicity_vs_pt = book<TH2F>("jet_puppiMultiplicity_vs_pt", ";# of constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_jet_LHA_vs_pt = book<TH2F>("jet_LHA_vs_pt", ";LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_pTD_vs_pt = book<TH2F>("jet_pTD_vs_pt", ";p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_width_vs_pt = book<TH2F>("jet_width_vs_pt", ";Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_thrust_vs_pt = book<TH2F>("jet_thrust_vs_pt", ";Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_jet_puppiMultiplicity_charged_vs_pt = book<TH2F>("jet_puppiMultiplicity_charged_vs_pt", ";# of charged constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_jet_LHA_charged_vs_pt = book<TH2F>("jet_LHA_charged_vs_pt", ";LHA charged (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_pTD_charged_vs_pt = book<TH2F>("jet_pTD_charged_vs_pt", ";p_{T}^{D} charged (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_width_charged_vs_pt = book<TH2F>("jet_width_charged_vs_pt", ";Width charged (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_thrust_charged_vs_pt = book<TH2F>("jet_thrust_charged_vs_pt", ";Thrust charged (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  if (is_mc_) {
    int rspNbins = 500;
    float rspMax = 5;
    h_jet_response_vs_genjet_pt = book<TH2F>("jet_response_vs_genjet_pt", ";response;GenJet p_{T} [GeV]", rspNbins, 0, rspMax, nPtBins, ptMin, ptMax);
    // with final unfolding binning
    h_jet_pt_vs_genjet_pt = book<TH2F>("jet_pt_vs_genjet_pt", ";Jet p_{T} [GeV];GenJet p_{T} [GeV]", Binning::nbins_pt_gen, &Binning::pt_bin_edges_gen[0], Binning::nbins_pt_gen, &Binning::pt_bin_edges_gen[0]);

    // deactivate for now as jets dont have proper flavour calculation
    // h_jet_flavour_vs_pt = book<TH2F>("jet_flavour_vs_pt", "jet flavour;PDGID;Jet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
    // h_jet1_flavour_vs_pt = book<TH2F>("jet1_flavour_vs_pt", "jet1 flavour;PDGID;Jet1 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
    // h_jet2_flavour_vs_pt = book<TH2F>("jet2_flavour_vs_pt", "jet2 flavour;PDGID;Jet2 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
    // h_jet_flavour_vs_eta = book<TH2F>("jet_flavour_vs_eta", "jet flavour;PDGID;Jet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);

    h_genjet_flavour_vs_pt = book<TH2F>("genjet_flavour_vs_pt", "genjet flavour;PDGID;GenJet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
    h_genjet1_flavour_vs_pt = book<TH2F>("genjet1_flavour_vs_pt", "genjet1 flavour;PDGID;GenJet1 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
    h_genjet2_flavour_vs_pt = book<TH2F>("genjet2_flavour_vs_pt", "genjet2 flavour;PDGID;GenJet2 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);

    h_genjet_flavour_vs_pt_nPartons.resize(N_PARTONS_MAX+1);
    h_genjet1_flavour_vs_pt_nPartons.resize(N_PARTONS_MAX+1);
    h_genjet2_flavour_vs_pt_nPartons.resize(N_PARTONS_MAX+1);

    for (uint n=0; n <= N_PARTONS_MAX; n++) {
      h_genjet_flavour_vs_pt_nPartons.at(n) = book<TH2F>(TString::Format("genjet_flavour_vs_pt_npartons_%d", n), TString::Format("genjet flavour for %d parton;PDGID;GenJet p_{T} [GeV]", n), 23, -0.5, 22.5, nPtBins, ptMin, ptMax);;
      h_genjet1_flavour_vs_pt_nPartons.at(n) = book<TH2F>(TString::Format("genjet1_flavour_vs_pt_npartons_%d", n), TString::Format("genjet1 flavour for %d parton;PDGID;GenJet1 p_{T} [GeV]", n), 23, -0.5, 22.5, nPtBins, ptMin, ptMax);;
      h_genjet2_flavour_vs_pt_nPartons.at(n) = book<TH2F>(TString::Format("genjet2_flavour_vs_pt_npartons_%d", n), TString::Format("genjet2 flavour for %d parton;PDGID;GenJet2 p_{T} [GeV]", n), 23, -0.5, 22.5, nPtBins, ptMin, ptMax);;
    }

    h_genjet_flavour_vs_eta_lowPt = book<TH2F>("genjet_flavour_vs_eta_lowPt", "genjet flavour;PDGID;GenJet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);
    h_genjet_flavour_vs_eta_midPt = book<TH2F>("genjet_flavour_vs_eta_midPt", "genjet flavour;PDGID;GenJet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);
    h_genjet_flavour_vs_eta_highPt = book<TH2F>("genjet_flavour_vs_eta_highPt", "genjet flavour;PDGID;GenJet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);
    h_genjet_flavour_vs_eta_highPt2 = book<TH2F>("genjet_flavour_vs_eta_highPt2", "genjet flavour;PDGID;GenJet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);

    // g jet only
    h_gjet_puppiMultiplicity_vs_pt = book<TH2F>("gjet_puppiMultiplicity_vs_pt", "g-flavour;# of constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
    h_gjet_LHA_vs_pt = book<TH2F>("gjet_LHA_vs_pt", "g-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_gjet_pTD_vs_pt = book<TH2F>("gjet_pTD_vs_pt", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_gjet_width_vs_pt = book<TH2F>("gjet_width_vs_pt", "g-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_gjet_thrust_vs_pt = book<TH2F>("gjet_thrust_vs_pt", "g-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

    h_gjet_puppiMultiplicity_charged_vs_pt = book<TH2F>("gjet_puppiMultiplicity_charged_vs_pt", "g-flavour;# of charged constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
    h_gjet_LHA_charged_vs_pt = book<TH2F>("gjet_LHA_charged_vs_pt", "g-flavour;LHA charged-only (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_gjet_pTD_charged_vs_pt = book<TH2F>("gjet_pTD_charged_vs_pt", "g-flavour;p_{T}^{D} charged-only (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_gjet_width_charged_vs_pt = book<TH2F>("gjet_width_charged_vs_pt", "g-flavour;Width charged-only (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_gjet_thrust_charged_vs_pt = book<TH2F>("gjet_thrust_charged_vs_pt", "g-flavour jet;Thrust charged-only (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

    h_gjet_response_vs_genjet_pt = book<TH2F>("gjet_response_vs_genjet_pt", ";g-flavour response;GenJet p_{T} [GeV]", rspNbins, 0, rspMax, nPtBins, ptMin, ptMax);

    // q jet only
    h_qjet_puppiMultiplicity_vs_pt = book<TH2F>("qjet_puppiMultiplicity_vs_pt", "uds-flavour;# of constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
    h_qjet_LHA_vs_pt = book<TH2F>("qjet_LHA_vs_pt", "uds-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_qjet_pTD_vs_pt = book<TH2F>("qjet_pTD_vs_pt", "uds-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_qjet_width_vs_pt = book<TH2F>("qjet_width_vs_pt", "uds-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_qjet_thrust_vs_pt = book<TH2F>("qjet_thrust_vs_pt", "uds-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

    h_qjet_puppiMultiplicity_charged_vs_pt = book<TH2F>("qjet_puppiMultiplicity_charged_vs_pt", "uds-flavour;# of charged constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
    h_qjet_LHA_charged_vs_pt = book<TH2F>("qjet_LHA_charged_vs_pt", "uds-flavour;LHA charged-only (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_qjet_pTD_charged_vs_pt = book<TH2F>("qjet_pTD_charged_vs_pt", "uds-flavour;p_{T}^{D} charged-only (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_qjet_width_charged_vs_pt = book<TH2F>("qjet_width_charged_vs_pt", "uds-flavour;Width charged-only (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_qjet_thrust_charged_vs_pt = book<TH2F>("qjet_thrust_charged_vs_pt", "uds-flavour jet;Thrust charged-only (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

    h_qjet_response_vs_genjet_pt = book<TH2F>("qjet_response_vs_genjet_pt", ";uds-flavour response;GenJet p_{T} [GeV]", rspNbins, 0, rspMax, nPtBins, ptMin, ptMax);

    // b/c jet only
    h_bcjet_puppiMultiplicity_vs_pt = book<TH2F>("bcjet_puppiMultiplicity_vs_pt", "b/c-flavour;# of constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
    h_bcjet_LHA_vs_pt = book<TH2F>("bcjet_LHA_vs_pt", "b/c-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_bcjet_pTD_vs_pt = book<TH2F>("bcjet_pTD_vs_pt", "b/c-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_bcjet_width_vs_pt = book<TH2F>("bcjet_width_vs_pt", "b/c-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_bcjet_thrust_vs_pt = book<TH2F>("bcjet_thrust_vs_pt", "b/c-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

    h_bcjet_puppiMultiplicity_charged_vs_pt = book<TH2F>("bcjet_puppiMultiplicity_charged_vs_pt", "b/c-flavour;# of charged constituents, PUPPI weighted (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
    h_bcjet_LHA_charged_vs_pt = book<TH2F>("bcjet_LHA_charged_vs_pt", "b/c-flavour;LHA charged-only (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_bcjet_pTD_charged_vs_pt = book<TH2F>("bcjet_pTD_charged_vs_pt", "b/c-flavour;p_{T}^{D} charged-only (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_bcjet_width_charged_vs_pt = book<TH2F>("bcjet_width_charged_vs_pt", "b/c-flavour;Width charged-only (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
    h_bcjet_thrust_charged_vs_pt = book<TH2F>("bcjet_thrust_charged_vs_pt", "b/c-flavour jet;Thrust charged-only (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

    h_bcjet_response_vs_genjet_pt = book<TH2F>("bcjet_response_vs_genjet_pt", ";b/c-flavour response;GenJet p_{T} [GeV]", rspNbins, 0, rspMax, nPtBins, ptMin, ptMax);
  }

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
      const LambdaCalculator<PFParticle> & recoJetCalc = jetLambdas.at(i).getLambdaCalculator(false, doGroomed_);
      const LambdaCalculator<PFParticle> & recoJetCalcCharged = jetLambdas.at(i).getLambdaCalculator(true, doGroomed_);

      float jet_pt = thisjet.pt();
      h_weights_vs_pt->Fill(weight, jet_pt);
      if (is_mc_) {
        if (event.genInfo->binningValues().size() > 0)
          h_pthat_vs_jet_pt->Fill(thisjet.pt(), event.genInfo->binningValues()[0]);
      }
      h_jet_pt_unweighted->Fill(jet_pt);
      h_jet_pt->Fill(jet_pt, weight);
      h_jet_y->Fill(thisjet.Rapidity(), weight);
      h_jet_eta->Fill(thisjet.eta(), weight);

      float puppiMult(0.), mult(0.), lha(0.), ptd(0.), width(0.), thrust(0.);
      float puppiMult_charged(0.), mult_charged(0.), lha_charged(0.), ptd_charged(0.), width_charged(0.), thrust_charged(0.);

      puppiMult = recoJetCalc.getLambda(Cuts::mult_args);
      mult = recoJetCalc.getLambda(Cuts::mult_args);
      lha = recoJetCalc.getLambda(Cuts::lha_args);
      ptd = recoJetCalc.getLambda(Cuts::pTD_args);
      width = recoJetCalc.getLambda(Cuts::width_args);
      thrust = recoJetCalc.getLambda(Cuts::thrust_args);

      h_jet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);
      h_jet_LHA_vs_pt->Fill(lha, jet_pt, weight);
      h_jet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
      h_jet_width_vs_pt->Fill(width, jet_pt, weight);
      h_jet_thrust_vs_pt->Fill(thrust, jet_pt, weight);

      // if (ptd > 0.5 && ptd < 0.501){
      //   cout << "+++++ ptd 0.5: " << ptd << endl;
      //   cout << "isgroomed?: " << doGroomed_ << endl;
      //   // cout << "ptSum: " << recoJetCalc.getPtSum() << endl;
      //   cout << " jet " << thisjet.pt() << " : " << thisjet.Rapidity() << " : " << thisjet.phi() << endl;
      //   cout << "calc constits: " << recoJetCalc.constits().size() << endl;
      //   for (const auto & constit : recoJetCalc.constits()) {
      //     cout << constit.pt() << " : " << constit.Rapidity() << " : " << constit.phi() << " : " << constit.puppiWeight() << endl;
      //   }

      //   cout << "jet constits: " << endl;
      //   for (const auto & cind : thisjet.pfcand_indexs()) {
      //     const PFParticle & pf = event.pfparticles->at(cind);
      //     cout << pf.pt() << " : " << pf.Rapidity() << " : " << pf.phi() << " : " << pf.puppiWeight() << endl;
      //   }
      // }
      // if (ptd > 0.33 && ptd < 0.335){
      //   cout << "***** ptd 0.33: " << ptd << endl;
      //   cout << "isgroomed?: " << doGroomed_ << endl;
      //   // cout << "ptSum: " << recoJetCalc.getPtSum() << endl;
      //   cout << " jet " << thisjet.pt() << " : " << thisjet.Rapidity() << " : " << thisjet.phi() << endl;
      //   cout << "calc constits: " << recoJetCalc.constits().size() << endl;
      //   for (const auto & constit : recoJetCalc.constits()) {
      //     cout << constit.pt() << " : " << constit.Rapidity() << " : " << constit.phi() << " : " << constit.puppiWeight() << endl;
      //   }

      //   cout << "jet constits: " << endl;
      //   for (const auto & cind : thisjet.pfcand_indexs()) {
      //     const PFParticle & pf = event.pfparticles->at(cind);
      //     cout << pf.pt() << " : " << pf.Rapidity() << " : " << pf.phi() << " : " << pf.puppiWeight() << endl;
      //   }
      // }

      puppiMult_charged = recoJetCalcCharged.getLambda(Cuts::mult_args);
      mult_charged = recoJetCalcCharged.getLambda(Cuts::mult_args);
      lha_charged = recoJetCalcCharged.getLambda(Cuts::lha_args);
      ptd_charged = recoJetCalcCharged.getLambda(Cuts::pTD_args);
      width_charged = recoJetCalcCharged.getLambda(Cuts::width_args);
      thrust_charged = recoJetCalcCharged.getLambda(Cuts::thrust_args);

      h_jet_puppiMultiplicity_charged_vs_pt->Fill(puppiMult_charged, jet_pt, weight);
      h_jet_LHA_charged_vs_pt->Fill(lha_charged, jet_pt, weight);
      h_jet_pTD_charged_vs_pt->Fill(ptd_charged, jet_pt, weight);
      h_jet_width_charged_vs_pt->Fill(width_charged, jet_pt, weight);
      h_jet_thrust_charged_vs_pt->Fill(thrust_charged, jet_pt, weight);

      // if (ptd_charged > 0.5 && ptd_charged < 0.501){
      //   cout << "+++++ ptd_charged 0.5: " << ptd_charged << endl;
      //   cout << "isgroomed?: " << doGroomed_ << endl;
      //   // cout << "ptSum: " << recoJetCalc.getPtSum() << endl;
      //   cout << " jet " << thisjet.pt() << " : " << thisjet.Rapidity() << " : " << thisjet.phi() << endl;
      //   cout << "calc constits: " << recoJetCalc.constits().size() << endl;
      //   for (const auto & constit : recoJetCalc.constits()) {
      //     cout << constit.pt() << " : " << constit.Rapidity() << " : " << constit.phi() << " : " << constit.puppiWeight() << " : " << constit.charge() << endl;
      //   }

      //   cout << "jet constits: " << endl;
      //   for (const auto & cind : thisjet.pfcand_indexs()) {
      //     const PFParticle & pf = event.pfparticles->at(cind);
      //     cout << pf.pt() << " : " << pf.Rapidity() << " : " << pf.phi() << " : " << pf.puppiWeight() << endl;
      //   }
      // }
      // if (ptd_charged > 0.33 && ptd_charged < 0.335){
      //   cout << "***** ptd_charged 0.33: " << ptd_charged << endl;
      //   cout << "isgroomed?: " << doGroomed_ << endl;
      //   // cout << "ptSum: " << recoJetCalc.getPtSum() << endl;
      //   cout << " jet " << thisjet.pt() << " : " << thisjet.Rapidity() << " : " << thisjet.phi() << endl;
      //   cout << "calc constits: " << recoJetCalc.constits().size() << endl;
      //   for (const auto & constit : recoJetCalc.constits()) {
      //     cout << constit.pt() << " : " << constit.Rapidity() << " : " << constit.phi() << " : " << constit.puppiWeight() << " : " << constit.charge() << endl;
      //   }

      //   cout << "jet constits: " << endl;
      //   for (const auto & cind : thisjet.pfcand_indexs()) {
      //     const PFParticle & pf = event.pfparticles->at(cind);
      //     cout << pf.pt() << " : " << pf.Rapidity() << " : " << pf.phi() << " : " << pf.puppiWeight() << " : " << pf.charge() <<  endl;
      //   }
      // }

      if (is_mc_ && passGen) {
        // Store variables for matched GenJet
        bool matchedGenJet = (thisjet.genjet_index() > -1 && thisjet.genjet_index() < 3); // put upper limit to avoid weird matches
        float genjet_pt = -1.;
        float genjet_y = 0.;
        float genjet_flav = PDGID::UNKNOWN;
        float response = -1.;

        if (matchedGenJet) {
          int thisInd = thisjet.genjet_index();
          if (thisInd >= (int) genjetLambdas->size()) {
            cout << "WARNING: wanted genjet_index " << thisInd << " but only have " << genjetLambdas->size() << " in genjetLambdas" << endl;
          }
          const GenJet & genjet = genjetLambdas->at(thisInd).jet;
          const LambdaCalculator<GenParticle> & matchedGenJetCalc = genjetLambdas->at(thisInd).getLambdaCalculator(false, doGroomed_);
          const LambdaCalculator<GenParticle> & matchedGenJetCalcCharged = genjetLambdas->at(thisInd).getLambdaCalculator(true, doGroomed_);

          genjet_pt = genjet.pt();
          genjet_y = genjet.Rapidity();
          genjet_flav = (useStatus23Flavour_) ? get_jet_flavour(genjet, event.genparticles, jetRadius/2., true) : genjet.partonFlavour();
          genjet_flav = abs(genjet_flav);
          if (genjet_flav > 100) genjet_flav = PDGID::UNKNOWN;

          response = jet_pt/genjet_pt;
          h_jet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
          h_jet_pt_vs_genjet_pt->Fill(jet_pt, genjet_pt, weight);

          // Response hists over all pt, and high/low pt split ones
          float gen_mult = matchedGenJetCalc.getLambda(Cuts::mult_args);
          fill_lambda_rsp_hists(puppiMult, gen_mult, weight,
            h_jet_puppiMultiplicity_response, h_jet_puppiMultiplicity_rel_response,
            jet_pt,
            h_jet_puppiMultiplicity_lowPt_response, h_jet_puppiMultiplicity_midPt_response, h_jet_puppiMultiplicity_highPt_response,
            h_jet_puppiMultiplicity_lowPt_rel_response, h_jet_puppiMultiplicity_midPt_rel_response, h_jet_puppiMultiplicity_highPt_rel_response);

          float gen_lha = matchedGenJetCalc.getLambda(Cuts::lha_args);
          fill_lambda_rsp_hists(lha, gen_lha, weight,
            h_jet_LHA_response, h_jet_LHA_rel_response,
            jet_pt,
            h_jet_LHA_lowPt_response, h_jet_LHA_midPt_response, h_jet_LHA_highPt_response,
            h_jet_LHA_lowPt_rel_response, h_jet_LHA_midPt_rel_response, h_jet_LHA_highPt_rel_response);

          float gen_ptd = matchedGenJetCalc.getLambda(Cuts::pTD_args);
          fill_lambda_rsp_hists(ptd, gen_ptd, weight,
            h_jet_pTD_response, h_jet_pTD_rel_response,
            jet_pt,
            h_jet_pTD_lowPt_response, h_jet_pTD_midPt_response, h_jet_pTD_highPt_response,
            h_jet_pTD_lowPt_rel_response, h_jet_pTD_midPt_rel_response, h_jet_pTD_highPt_rel_response);

          float gen_width = matchedGenJetCalc.getLambda(Cuts::width_args);
          fill_lambda_rsp_hists(width, gen_width, weight,
            h_jet_width_response, h_jet_width_rel_response,
            jet_pt,
            h_jet_width_lowPt_response, h_jet_width_midPt_response, h_jet_width_highPt_response,
            h_jet_width_lowPt_rel_response, h_jet_width_midPt_rel_response, h_jet_width_highPt_rel_response);

          float gen_thrust = matchedGenJetCalc.getLambda(Cuts::thrust_args);
          fill_lambda_rsp_hists(thrust, gen_thrust, weight,
            h_jet_thrust_response, h_jet_thrust_rel_response,
            jet_pt,
            h_jet_thrust_lowPt_response, h_jet_thrust_midPt_response, h_jet_thrust_highPt_response,
            h_jet_thrust_lowPt_rel_response, h_jet_thrust_midPt_rel_response, h_jet_thrust_highPt_rel_response);

          // Do charged-only daughters version
          float gen_mult_charged = matchedGenJetCalcCharged.getLambda(Cuts::mult_args);
          fill_lambda_rsp_hists(puppiMult_charged, gen_mult_charged, weight,
            h_jet_puppiMultiplicity_charged_response, h_jet_puppiMultiplicity_charged_rel_response,
            jet_pt,
            h_jet_puppiMultiplicity_charged_lowPt_response, h_jet_puppiMultiplicity_charged_midPt_response, h_jet_puppiMultiplicity_charged_highPt_response,
            h_jet_puppiMultiplicity_charged_lowPt_rel_response, h_jet_puppiMultiplicity_charged_midPt_rel_response, h_jet_puppiMultiplicity_charged_highPt_rel_response);

          float gen_lha_charged = matchedGenJetCalcCharged.getLambda(Cuts::lha_args);
          fill_lambda_rsp_hists(lha_charged, gen_lha_charged, weight,
            h_jet_LHA_charged_response, h_jet_LHA_charged_rel_response,
            jet_pt,
            h_jet_LHA_charged_lowPt_response, h_jet_LHA_charged_midPt_response, h_jet_LHA_charged_highPt_response,
            h_jet_LHA_charged_lowPt_rel_response, h_jet_LHA_charged_midPt_rel_response, h_jet_LHA_charged_highPt_rel_response);

          float gen_ptd_charged = matchedGenJetCalcCharged.getLambda(Cuts::pTD_args);
          fill_lambda_rsp_hists(ptd_charged, gen_ptd_charged, weight,
            h_jet_pTD_charged_response, h_jet_pTD_charged_rel_response,
            jet_pt,
            h_jet_pTD_charged_lowPt_response, h_jet_pTD_charged_midPt_response, h_jet_pTD_charged_highPt_response,
            h_jet_pTD_charged_lowPt_rel_response, h_jet_pTD_charged_midPt_rel_response, h_jet_pTD_charged_highPt_rel_response);

          float gen_width_charged = matchedGenJetCalcCharged.getLambda(Cuts::width_args);
          fill_lambda_rsp_hists(width_charged, gen_width_charged, weight,
            h_jet_width_charged_response, h_jet_width_charged_rel_response,
            jet_pt,
            h_jet_width_charged_lowPt_response, h_jet_width_charged_midPt_response, h_jet_width_charged_highPt_response,
            h_jet_width_charged_lowPt_rel_response, h_jet_width_charged_midPt_rel_response, h_jet_width_charged_highPt_rel_response);

          float gen_thrust_charged = matchedGenJetCalcCharged.getLambda(Cuts::thrust_args);
          fill_lambda_rsp_hists(thrust_charged, gen_thrust_charged, weight,
            h_jet_thrust_charged_response, h_jet_thrust_charged_rel_response,
            jet_pt,
            h_jet_thrust_charged_lowPt_response, h_jet_thrust_charged_midPt_response, h_jet_thrust_charged_highPt_response,
            h_jet_thrust_charged_lowPt_rel_response, h_jet_thrust_charged_midPt_rel_response, h_jet_thrust_charged_highPt_rel_response);

          h_genjet_flavour_vs_pt->Fill(genjet_flav, genjet_pt, weight);

          if (genjet_pt > rsp_highPt2_cut_)
            h_genjet_flavour_vs_eta_highPt2->Fill(genjet_flav, genjet_y, weight);
          else if (genjet_pt > rsp_highPt_cut_)
            h_genjet_flavour_vs_eta_highPt->Fill(genjet_flav, genjet_y, weight);
          else if (genjet_pt > rsp_midPt_cut_)
            h_genjet_flavour_vs_eta_midPt->Fill(genjet_flav, genjet_y, weight);
          else if (genjet_pt > rsp_lowPt_cut_)
            h_genjet_flavour_vs_eta_lowPt->Fill(genjet_flav, genjet_y, weight);

          if (i == 0) {
            h_genjet1_flavour_vs_pt->Fill(genjet_flav, genjet_pt, weight);
          } else if (i == 1) {
            h_genjet2_flavour_vs_pt->Fill(genjet_flav, genjet_pt, weight);
          }

          uint nPartons = get_num_outgoing_partons(*event.genparticles);
          if (nPartons > N_PARTONS_MAX) throw std::runtime_error("Too many partons " + nPartons);
          h_genjet_flavour_vs_pt_nPartons.at(nPartons)->Fill(genjet_flav, genjet_pt, weight);
          if (i == 0) {
            h_genjet1_flavour_vs_pt_nPartons.at(nPartons)->Fill(genjet_flav, genjet_pt, weight);
          } else if (i == 1) {
            h_genjet2_flavour_vs_pt_nPartons.at(nPartons)->Fill(genjet_flav, genjet_pt, weight);
          }

          if (genjet_flav == PDGID::GLUON) {
            h_gjet_puppiMultiplicity_vs_pt->Fill(mult, genjet_pt, weight);
            h_gjet_LHA_vs_pt->Fill(lha, genjet_pt, weight);
            h_gjet_pTD_vs_pt->Fill(ptd, genjet_pt, weight);
            h_gjet_width_vs_pt->Fill(width, genjet_pt, weight);
            h_gjet_thrust_vs_pt->Fill(thrust, genjet_pt, weight);

            h_gjet_puppiMultiplicity_charged_vs_pt->Fill(mult_charged, genjet_pt, weight);
            h_gjet_LHA_charged_vs_pt->Fill(lha_charged, genjet_pt, weight);
            h_gjet_pTD_charged_vs_pt->Fill(ptd_charged, genjet_pt, weight);
            h_gjet_width_charged_vs_pt->Fill(width_charged, genjet_pt, weight);
            h_gjet_thrust_charged_vs_pt->Fill(thrust_charged, genjet_pt, weight);

            h_gjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);

          } else if (genjet_flav >= PDGID::DOWN_QUARK && genjet_flav <= PDGID::STRANGE_QUARK) {
            h_qjet_puppiMultiplicity_vs_pt->Fill(mult, genjet_pt, weight);
            h_qjet_LHA_vs_pt->Fill(lha, genjet_pt, weight);
            h_qjet_pTD_vs_pt->Fill(ptd, genjet_pt, weight);
            h_qjet_width_vs_pt->Fill(width, genjet_pt, weight);
            h_qjet_thrust_vs_pt->Fill(thrust, genjet_pt, weight);

            h_qjet_puppiMultiplicity_charged_vs_pt->Fill(mult_charged, genjet_pt, weight);
            h_qjet_LHA_charged_vs_pt->Fill(lha_charged, genjet_pt, weight);
            h_qjet_pTD_charged_vs_pt->Fill(ptd_charged, genjet_pt, weight);
            h_qjet_width_charged_vs_pt->Fill(width_charged, genjet_pt, weight);
            h_qjet_thrust_charged_vs_pt->Fill(thrust_charged, genjet_pt, weight);

            h_qjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);

          } else if (genjet_flav == PDGID::CHARM_QUARK || genjet_flav == PDGID::BOTTOM_QUARK) {
            h_bcjet_puppiMultiplicity_vs_pt->Fill(mult, genjet_pt, weight);
            h_bcjet_LHA_vs_pt->Fill(lha, genjet_pt, weight);
            h_bcjet_pTD_vs_pt->Fill(ptd, genjet_pt, weight);
            h_bcjet_width_vs_pt->Fill(width, genjet_pt, weight);
            h_bcjet_thrust_vs_pt->Fill(thrust, genjet_pt, weight);

            h_bcjet_puppiMultiplicity_charged_vs_pt->Fill(mult_charged, genjet_pt, weight);
            h_bcjet_LHA_charged_vs_pt->Fill(lha_charged, genjet_pt, weight);
            h_bcjet_pTD_charged_vs_pt->Fill(ptd_charged, genjet_pt, weight);
            h_bcjet_width_charged_vs_pt->Fill(width_charged, genjet_pt, weight);
            h_bcjet_thrust_charged_vs_pt->Fill(thrust_charged, genjet_pt, weight);

            h_bcjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
          }
        } // end if matchedgenjet

        // int jet_flav = (useStatus23Flavour_) ? get_jet_flavour(thisjet, event.genparticles, jetRadius/2., true) : thisjet.partonFlavour();
        // jet_flav = abs(jet_flav);
        // if (jet_flav > 100) jet_flav = PDGID::UNKNOWN;

        // // Split by actual jet flavour - these only make sense for MC
        // if (jet_flav == PDGID::GLUON) { // gluon jets
        //   h_gjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
        //   h_gjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
        //   h_gjet_width_vs_pt->Fill(width, jet_pt, weight);
        //   h_gjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);
        //   if (matchedGenJet) h_gjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
        //   h_gjet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);

        // } else if ((jet_flav <= PDGID::STRANGE_QUARK) && (jet_flav >= PDGID::DOWN_QUARK)){ // uds jets
        //   h_qjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
        //   h_qjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
        //   h_qjet_width_vs_pt->Fill(width, jet_pt, weight);
        //   h_qjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);
        //   if (matchedGenJet) h_qjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
        //   h_qjet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);

        // } else {
        //   h_bcjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
        //   h_bcjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
        //   h_bcjet_width_vs_pt->Fill(width, jet_pt, weight);
        //   h_bcjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);
        //   if (matchedGenJet) h_bcjet_response_vs_genjet_pt->Fill(response, genjet_pt, weight);
        //   h_bcjet_puppiMultiplicity_vs_pt->Fill(puppiMult, jet_pt, weight);
        // }

        // h_jet_flavour->Fill(jet_flav, weight);
        // h_jet_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
        // if (i == 0) {
        //   h_jet1_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
        // } else if (i == 1) {
        //   h_jet2_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
        // }

        // h_jet_flavour_vs_eta->Fill(jet_flav, thisjet.Rapidity(), weight);

      } // end is_mc_ && passGen
    } // end loop over reco jets
  } //end if passReco
}


void QGAnalysisHists::fill_lambda_rsp_hists(float reco_val, float gen_val, float weight,
                                            TH2F * response, TH2F * rel_response,
                                            float jet_pt,
                                            TH2F * response_lowPt, TH2F * response_midPt, TH2F * response_highPt,
                                            TH2F * rel_response_lowPt, TH2F * rel_response_midPt, TH2F * rel_response_highPt) {
  float rsp = (gen_val != 0) ? reco_val / gen_val : 999999; // TODO handle infs
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
 * Get object flavour based on closest genparticle
 * @param  obj          [description]
 * @param  genparticles [description]
 * @param  dr_max       [description]
 * @param  pythiaMode   [description]
 * @return              [description]
 */
int QGAnalysisHists::get_jet_flavour(const Particle & obj, std::vector<GenParticle>* genparticles, float dr_max, bool pythiaMode) {

  // cout << "obj:" << obj.pt() << " : " << obj.eta() << " : " << obj.phi() << endl;
  float smallest_dr = dr_max;
  int pdgid = 0;
  for (const auto& gp : *genparticles) {
    if (pythiaMode && abs(gp.status()) != 23) continue;
    int absid = abs(gp.pdgId());
    if (absid > PDGID::GLUON) continue;
    if (absid > PDGID::TOP_QUARK && absid != PDGID::GLUON) continue;
    float dr = deltaR(obj.v4(), gp.v4());
    if (dr < smallest_dr) {
      pdgid = gp.pdgId();
      smallest_dr = dr;
      // cout << gp.pt() << " : " << gp.eta() << " : " << gp.phi() << " : " << gp.pdgId() << " : " << gp.status() << " : " << deltaR(obj.v4(), gp.v4()) << endl;
    }
  }
  return pdgid;
}


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

/**
 * Counter number of outgoing partons in Matrix Element.
 * Only works for Pythia-based status codes
 */
uint QGAnalysisHists::get_num_outgoing_partons(const std::vector<GenParticle> & genparticles) {
  uint counter = 0;
  for (const auto & gp : genparticles) {
    int pgd = abs(gp.pdgId());
    if ((gp.status() == 23) && isParton(gp.pdgId())) counter++;
  }
  return counter;
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
#include "UHH2/QGAnalysis/include/QGAnalysisGenHists.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisGenHists::QGAnalysisGenHists(Context & ctx,
                                       const string & dirname,
                                       int useNJets,
                                       bool doGroomed,
                                       bool useStatus23Flavour,
                                       const string & selection,
                                       const string & gen_sel_handle_name,
                                       const string & gen_jetlambda_handle_name
                                 ):
  Hists(ctx, dirname),
  dirname_(dirname),
  useNJets_(useNJets),
  doGroomed_(doGroomed),
  useStatus23Flavour_(useStatus23Flavour),
  N_PARTONS_MAX(4)
  {

  dataset_ = matchDatasetName(ctx.get("dataset_version"));
  string jetCone = ctx.get("JetCone", "AK4");
  jetRadius = get_jet_radius(jetCone);

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

  // 2D lambda vs PT
  // -----------------
  // All flavs
  h_jet_puppiMultiplicity_vs_pt = book<TH2F>("jet_puppiMultiplicity_vs_pt", ";# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_jet_LHA_vs_pt = book<TH2F>("jet_LHA_vs_pt", ";LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_pTD_vs_pt = book<TH2F>("jet_pTD_vs_pt", ";p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_width_vs_pt = book<TH2F>("jet_width_vs_pt", ";Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_thrust_vs_pt = book<TH2F>("jet_thrust_vs_pt", ";Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_jet_puppiMultiplicity_charged_vs_pt = book<TH2F>("jet_puppiMultiplicity_charged_vs_pt", ";# of charged constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_jet_LHA_charged_vs_pt = book<TH2F>("jet_LHA_charged_vs_pt", ";LHA charged (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_pTD_charged_vs_pt = book<TH2F>("jet_pTD_charged_vs_pt", ";p_{T}^{D} charged (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_width_charged_vs_pt = book<TH2F>("jet_width_charged_vs_pt", ";Width charged (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_thrust_charged_vs_pt = book<TH2F>("jet_thrust_charged_vs_pt", ";Thrust charged (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_jet_flavour_vs_pt = book<TH2F>("jet_flavour_vs_pt", "Genjet flavour;PDGID;GenJet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet1_flavour_vs_pt = book<TH2F>("jet1_flavour_vs_pt", "Genjet1 flavour;PDGID;GenJet1 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet2_flavour_vs_pt = book<TH2F>("jet2_flavour_vs_pt", "Genjet2 flavour;PDGID;GenJet2 p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);

  h_jet_flavour_vs_pt_nPartons.resize(N_PARTONS_MAX+1);
  h_jet1_flavour_vs_pt_nPartons.resize(N_PARTONS_MAX+1);
  h_jet2_flavour_vs_pt_nPartons.resize(N_PARTONS_MAX+1);

  for (uint n=0; n <= N_PARTONS_MAX; n++) {
    h_jet_flavour_vs_pt_nPartons.at(n) = book<TH2F>(TString::Format("jet_flavour_vs_pt_npartons_%d", n), TString::Format("Genjet flavour for %d parton;PDGID;GenJet p_{T} [GeV]", n), 23, -0.5, 22.5, nPtBins, ptMin, ptMax);;
    h_jet1_flavour_vs_pt_nPartons.at(n) = book<TH2F>(TString::Format("jet1_flavour_vs_pt_npartons_%d", n), TString::Format("Genjet1 flavour for %d parton;PDGID;GenJet1 p_{T} [GeV]", n), 23, -0.5, 22.5, nPtBins, ptMin, ptMax);;
    h_jet2_flavour_vs_pt_nPartons.at(n) = book<TH2F>(TString::Format("jet2_flavour_vs_pt_npartons_%d", n), TString::Format("Genjet2 flavour for %d parton;PDGID;GenJet2 p_{T} [GeV]", n), 23, -0.5, 22.5, nPtBins, ptMin, ptMax);;
  }

  h_jet_flavour_vs_eta_lowPt = book<TH2F>("jet_flavour_vs_eta_lowPt", "Genjet flavour;PDGID;Jet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);
  h_jet_flavour_vs_eta_midPt = book<TH2F>("jet_flavour_vs_eta_midPt", "Genjet flavour;PDGID;Jet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);
  h_jet_flavour_vs_eta_highPt = book<TH2F>("jet_flavour_vs_eta_highPt", "Genjet flavour;PDGID;Jet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);
  h_jet_flavour_vs_eta_highPt2 = book<TH2F>("jet_flavour_vs_eta_highPt2", "Genjet flavour;PDGID;Jet y", 23, -0.5, 22.5, nEtaBins, etaMin, etaMax);

  // 2D lambda vs pT, split by flavour groups
  // g jet only
  h_gjet_puppiMultiplicity_vs_pt = book<TH2F>("gjet_puppiMultiplicity_vs_pt", "g-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_gjet_LHA_vs_pt = book<TH2F>("gjet_LHA_vs_pt", "g-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_pTD_vs_pt = book<TH2F>("gjet_pTD_vs_pt", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_width_vs_pt = book<TH2F>("gjet_width_vs_pt", "g-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_thrust_vs_pt = book<TH2F>("gjet_thrust_vs_pt", "g-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_gjet_puppiMultiplicity_charged_vs_pt = book<TH2F>("gjet_puppiMultiplicity_charged_vs_pt", "g-flavour;# of charged constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_gjet_LHA_charged_vs_pt = book<TH2F>("gjet_LHA_charged_vs_pt", "g-flavour;LHA charged-only (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_pTD_charged_vs_pt = book<TH2F>("gjet_pTD_charged_vs_pt", "g-flavour;p_{T}^{D} charged-only (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_width_charged_vs_pt = book<TH2F>("gjet_width_charged_vs_pt", "g-flavour;Width charged-only (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_thrust_charged_vs_pt = book<TH2F>("gjet_thrust_charged_vs_pt", "g-flavour jet;Thrust charged-only (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  // uds jet only
  h_qjet_puppiMultiplicity_vs_pt = book<TH2F>("qjet_puppiMultiplicity_vs_pt", "uds-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_qjet_LHA_vs_pt = book<TH2F>("qjet_LHA_vs_pt", "uds-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_pTD_vs_pt = book<TH2F>("qjet_pTD_vs_pt", "uds-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_width_vs_pt = book<TH2F>("qjet_width_vs_pt", "uds-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_thrust_vs_pt = book<TH2F>("qjet_thrust_vs_pt", "uds-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_qjet_puppiMultiplicity_charged_vs_pt = book<TH2F>("qjet_puppiMultiplicity_charged_vs_pt", "uds-flavour;# of charged constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_qjet_LHA_charged_vs_pt = book<TH2F>("qjet_LHA_charged_vs_pt", "uds-flavour;LHA charged-only (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_pTD_charged_vs_pt = book<TH2F>("qjet_pTD_charged_vs_pt", "uds-flavour;p_{T}^{D} charged-only (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_width_charged_vs_pt = book<TH2F>("qjet_width_charged_vs_pt", "uds-flavour;Width charged-only (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_thrust_charged_vs_pt = book<TH2F>("qjet_thrust_charged_vs_pt", "uds-flavour jet;Thrust charged-only (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  // bc jet only
  h_bcjet_puppiMultiplicity_vs_pt = book<TH2F>("bcjet_puppiMultiplicity_vs_pt", "b/c-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_bcjet_LHA_vs_pt = book<TH2F>("bcjet_LHA_vs_pt", "b/c-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_bcjet_pTD_vs_pt = book<TH2F>("bcjet_pTD_vs_pt", "b/c-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_bcjet_width_vs_pt = book<TH2F>("bcjet_width_vs_pt", "b/c-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_bcjet_thrust_vs_pt = book<TH2F>("bcjet_thrust_vs_pt", "b/c-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_bcjet_puppiMultiplicity_charged_vs_pt = book<TH2F>("bcjet_puppiMultiplicity_charged_vs_pt", "b/c-flavour;# of charged constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_bcjet_LHA_charged_vs_pt = book<TH2F>("bcjet_LHA_charged_vs_pt", "b/c-flavour;LHA charged-only (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_bcjet_pTD_charged_vs_pt = book<TH2F>("bcjet_pTD_charged_vs_pt", "b/c-flavour;p_{T}^{D} charged-only (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_bcjet_width_charged_vs_pt = book<TH2F>("bcjet_width_charged_vs_pt", "b/c-flavour;Width charged-only (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_bcjet_thrust_charged_vs_pt = book<TH2F>("bcjet_thrust_charged_vs_pt", "b/c-flavour jet;Thrust charged-only (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  // genJets_handle = ctx.get_handle< std::vector<GenJet> > ("GoodGenJets");
  genJetsLambda_handle = ctx.get_handle< std::vector<GenJetLambdaBundle> > (gen_jetlambda_handle_name);
  pass_gen_handle = ctx.get_handle<bool> (gen_sel_handle_name);
}


void QGAnalysisGenHists::fill(const Event & event){
  if (event.isRealData) { return; }

  double weight = event.weight;

  // extract the separate gen & reco weight components, needed for TUnfold
  double gen_weight = event.get(gen_weight_handle);

  h_weights->Fill(weight);

  if (event.genInfo->binningValues().size() > 0)
    h_pthat_vs_weight->Fill(weight, event.genInfo->binningValues()[0]);

  const std::vector<GenJetLambdaBundle> & genjetLambdas = event.get(genJetsLambda_handle);;
  bool passGen = event.get(pass_gen_handle);
  if (!passGen) { return; }

  // Fill jet hists
  // ---------------------------------------------------------------------------
  // At this point, all jet filtering etc should have already been performed
  if (useNJets_ > (int) genjetLambdas.size()) {
    cout << "useNJets_: " << useNJets_ << endl;
    cout << "genjetLambdas: " << genjetLambdas.size() << endl;
    throw runtime_error("useNJets_ > genjetLambdas.size()");
  }

  for (int i = 0; i < useNJets_; i++) {
    const GenJet & thisjet = genjetLambdas.at(i).jet;
    const LambdaCalculator<GenParticle> & genJetCalc = genjetLambdas.at(i).getLambdaCalculator(false, doGroomed_);
    const LambdaCalculator<GenParticle> & genJetCalcCharged = genjetLambdas.at(i).getLambdaCalculator(true, doGroomed_);

    float jet_pt = thisjet.pt();
    float jet_y = thisjet.Rapidity();
    h_weights_vs_pt->Fill(weight, jet_pt);
    h_jet_pt_unweighted->Fill(jet_pt);
    h_jet_pt->Fill(jet_pt, weight);
    h_jet_y->Fill(jet_y, weight);
    h_jet_eta->Fill(thisjet.eta(), weight);

    float puppiMult(0.), mult(0.), lha(0.), ptd(0.), width(0.), thrust(0.);
    float puppiMult_charged(0.), mult_charged(0.), lha_charged(0.), ptd_charged(0.), width_charged(0.), thrust_charged(0.);

    mult = genJetCalc.getLambda(Cuts::mult_args);
    lha = genJetCalc.getLambda(Cuts::lha_args);
    ptd = genJetCalc.getLambda(Cuts::pTD_args);
    width = genJetCalc.getLambda(Cuts::width_args);
    thrust = genJetCalc.getLambda(Cuts::thrust_args);

    h_jet_puppiMultiplicity_vs_pt->Fill(mult, jet_pt, weight);
    h_jet_LHA_vs_pt->Fill(lha, jet_pt, weight);
    h_jet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
    h_jet_width_vs_pt->Fill(width, jet_pt, weight);
    h_jet_thrust_vs_pt->Fill(thrust, jet_pt, weight);

    mult_charged = genJetCalcCharged.getLambda(Cuts::mult_args);
    lha_charged = genJetCalcCharged.getLambda(Cuts::lha_args);
    ptd_charged = genJetCalcCharged.getLambda(Cuts::pTD_args);
    width_charged = genJetCalcCharged.getLambda(Cuts::width_args);
    thrust_charged = genJetCalcCharged.getLambda(Cuts::thrust_args);

    h_jet_puppiMultiplicity_charged_vs_pt->Fill(mult_charged, jet_pt, weight);
    h_jet_LHA_charged_vs_pt->Fill(lha_charged, jet_pt, weight);
    h_jet_pTD_charged_vs_pt->Fill(ptd_charged, jet_pt, weight);
    h_jet_width_charged_vs_pt->Fill(width_charged, jet_pt, weight);
    h_jet_thrust_charged_vs_pt->Fill(thrust_charged, jet_pt, weight);

    int jet_flav = (useStatus23Flavour_) ? get_jet_flavour(thisjet, event.genparticles, jetRadius/2., true) : thisjet.partonFlavour();
    jet_flav = abs(jet_flav);
    if (jet_flav > 100) jet_flav = PDGID::UNKNOWN;

    h_jet_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);

    if (jet_pt > 500)
      h_jet_flavour_vs_eta_highPt->Fill(jet_flav, jet_y, weight);
    else if (jet_pt > 250)
      h_jet_flavour_vs_eta_highPt->Fill(jet_flav, jet_y, weight);
    else if (jet_pt > 100)
      h_jet_flavour_vs_eta_midPt->Fill(jet_flav, jet_y, weight);
    else if (jet_pt > 30)
      h_jet_flavour_vs_eta_lowPt->Fill(jet_flav, jet_y, weight);

    if (i == 0) {
      h_jet1_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
    } else if (i == 1) {
      h_jet2_flavour_vs_pt->Fill(jet_flav, jet_pt, weight);
    }

    uint nPartons = get_num_outgoing_partons(*event.genparticles);
    if (nPartons > N_PARTONS_MAX) throw std::runtime_error("Too many partons " + nPartons);
    h_jet_flavour_vs_pt_nPartons.at(nPartons)->Fill(jet_flav, jet_pt, weight);
    if (i == 0) {
      h_jet1_flavour_vs_pt_nPartons.at(nPartons)->Fill(jet_flav, jet_pt, weight);
    } else if (i == 1) {
      h_jet2_flavour_vs_pt_nPartons.at(nPartons)->Fill(jet_flav, jet_pt, weight);
    }

    if (jet_flav == PDGID::GLUON) {
      h_gjet_puppiMultiplicity_vs_pt->Fill(mult, jet_pt, weight);
      h_gjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
      h_gjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
      h_gjet_width_vs_pt->Fill(width, jet_pt, weight);
      h_gjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);

      h_gjet_puppiMultiplicity_charged_vs_pt->Fill(mult_charged, jet_pt, weight);
      h_gjet_LHA_charged_vs_pt->Fill(lha_charged, jet_pt, weight);
      h_gjet_pTD_charged_vs_pt->Fill(ptd_charged, jet_pt, weight);
      h_gjet_width_charged_vs_pt->Fill(width_charged, jet_pt, weight);
      h_gjet_thrust_charged_vs_pt->Fill(thrust_charged, jet_pt, weight);

    } else if (jet_flav >= PDGID::DOWN_QUARK && jet_flav <= PDGID::STRANGE_QUARK) {
      h_qjet_puppiMultiplicity_vs_pt->Fill(mult, jet_pt, weight);
      h_qjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
      h_qjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
      h_qjet_width_vs_pt->Fill(width, jet_pt, weight);
      h_qjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);

      h_qjet_puppiMultiplicity_charged_vs_pt->Fill(mult_charged, jet_pt, weight);
      h_qjet_LHA_charged_vs_pt->Fill(lha_charged, jet_pt, weight);
      h_qjet_pTD_charged_vs_pt->Fill(ptd_charged, jet_pt, weight);
      h_qjet_width_charged_vs_pt->Fill(width_charged, jet_pt, weight);
      h_qjet_thrust_charged_vs_pt->Fill(thrust_charged, jet_pt, weight);

    } else if (jet_flav == PDGID::CHARM_QUARK || jet_flav == PDGID::BOTTOM_QUARK) {
      h_bcjet_puppiMultiplicity_vs_pt->Fill(mult, jet_pt, weight);
      h_bcjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
      h_bcjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
      h_bcjet_width_vs_pt->Fill(width, jet_pt, weight);
      h_bcjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);

      h_bcjet_puppiMultiplicity_charged_vs_pt->Fill(mult_charged, jet_pt, weight);
      h_bcjet_LHA_charged_vs_pt->Fill(lha_charged, jet_pt, weight);
      h_bcjet_pTD_charged_vs_pt->Fill(ptd_charged, jet_pt, weight);
      h_bcjet_width_charged_vs_pt->Fill(width_charged, jet_pt, weight);
      h_bcjet_thrust_charged_vs_pt->Fill(thrust_charged, jet_pt, weight);
    }
  } // end loop over jets
}


/**
 * Get the collection of GenParticle*s for a given GenJet
 */
// std::vector<GenParticle*> QGAnalysisGenHists::get_genjet_genparticles(const GenJet & jet, std::vector<GenParticle>* genparticles) {
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
int QGAnalysisGenHists::get_jet_flavour(const Particle & obj, std::vector<GenParticle>* genparticles, float dr_max, bool pythiaMode) {

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
 * Counter number of outgoing partons in Matrix Element.
 * Only works for Pythia-based status codes
 */
uint QGAnalysisGenHists::get_num_outgoing_partons(const std::vector<GenParticle> & genparticles) {
  uint counter = 0;
  for (const auto & gp : genparticles) {
    int pgd = abs(gp.pdgId());
    if ((gp.status() == 23) && isParton(gp.pdgId())) counter++;
  }
  return counter;
}

QGAnalysisGenHists::~QGAnalysisGenHists(){}

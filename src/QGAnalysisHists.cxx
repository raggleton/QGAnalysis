#include "UHH2/QGAnalysis/include/QGAnalysisHists.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisHists::QGAnalysisHists(Context & ctx, const string & dirname, int useNJets, const string & selection): Hists(ctx, dirname), useNJets_(useNJets) {
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
    reweightHist->SetDirectory(0);
  }

  // book all histograms here
  int nWeightBins = 11;
  double weightBins [nWeightBins+1] = {1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 10};

  h_weights = book<TH1F>("weights", ";weight;", nWeightBins, weightBins);
  h_weights_vs_pt = book<TH2F>("weights_vs_pt", ";weight;", nWeightBins, weightBins, 200, 0., 2000.);
  h_PU_pT_hat_max_vs_weight = book<TH2F>("PU_pT_hat_max_vs_pt", ";weight;PU_pT_hat_max", nWeightBins, weightBins, 200, 0, 200);
  h_qscale_vs_weight = book<TH2F>("qscale_vs_weight", ";weight;qScale", nWeightBins, weightBins, 200, 0, 2000);
  h_pthat_vs_weight = book<TH2F>("pthat_vs_weight", ";weight;ptHat", nWeightBins, weightBins, 200, 0, 2000);
  h_pthat_vs_jet_pt = book<TH2F>("pthat_vs_jet_pt", ";p_{T}^{leading jet};ptHat", 200, 0, 2000, 200, 0, 2000);

  // RECO jet hists
  // --------------
  // For all jets
  int nPtBins = 100;
  float ptMin(0.), ptMax(2000.);
  h_jet_pt = book<TH1F>("jet_pt", ";p_{T}^{j} [GeV];", 2*nPtBins, ptMin, ptMax); // use as a catch-all to see we haven't missed highPT jets
  h_jet_pt_unweighted = book<TH1F>("jet_pt_unweighted", ";p_{T}^{j} [GeV];", 2*nPtBins, ptMin, ptMax); // use as a catch-all to see we haven't missed highPT jets

  int nEtaBins = 50;
  float etaMin(-2.5), etaMax(2.5);
  h_jet_eta = book<TH1F>("jet_eta", ";#eta^{j};", nEtaBins, etaMin, etaMax);

  h_jet_flavour = book<TH1F>("jet_flavour", "jet flavour;PDGID;", 23, -0.5, 22.5);
  h_jet_genParton_flavour = book<TH1F>("jet_genParton_flavour", "jet flavour (genParton);PDGID;", 23, -0.5, 22.5);

  int nMultBins = 100;
  h_jet_multiplicity = book<TH1F>("jet_multiplicity", ";# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);

  int nBins = 100;
  h_jet_LHA = book<TH1F>("jet_LHA", ";LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_jet_pTD = book<TH1F>("jet_pTD", ";p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_jet_width = book<TH1F>("jet_width", ";Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_jet_thrust = book<TH1F>("jet_thrust", ";Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // Flavour-tagged
  h_qjet_multiplicity = book<TH1F>("qjet_multiplicity", "q-flavour;# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_qjet_LHA = book<TH1F>("qjet_LHA", "q-flavour;LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_qjet_pTD = book<TH1F>("qjet_pTD", "q-flavour;p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_qjet_width = book<TH1F>("qjet_width", "q-flavour;Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_qjet_thrust = book<TH1F>("qjet_thrust", "q-flavour jet;Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  h_gjet_multiplicity = book<TH1F>("gjet_multiplicity", "g-flavour;# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_gjet_LHA = book<TH1F>("gjet_LHA", "g-flavour;LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_gjet_pTD = book<TH1F>("gjet_pTD", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_gjet_width = book<TH1F>("gjet_width", "g-flavour;Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_gjet_thrust = book<TH1F>("gjet_thrust", "g-flavour jet;Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // 2D versions vs PT
  h_jet_multiplicity_vs_pt = book<TH2F>("jet_multiplicity_vs_pt", ";# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_jet_LHA_vs_pt = book<TH2F>("jet_LHA_vs_pt", ";LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_pTD_vs_pt = book<TH2F>("jet_pTD_vs_pt", ";p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_width_vs_pt = book<TH2F>("jet_width_vs_pt", ";Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_jet_thrust_vs_pt = book<TH2F>("jet_thrust_vs_pt", ";Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_jet_flavour_vs_pt = book<TH2F>("jet_flavour_vs_pt", "jet flavour;PDGID;Jet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);
  h_jet_genParton_flavour_vs_pt = book<TH2F>("jet_genParton_flavour_vs_pt", "jet flavour;PDGID;Jet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);

  h_qjet_multiplicity_vs_pt = book<TH2F>("qjet_multiplicity_vs_pt", "q-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_qjet_LHA_vs_pt = book<TH2F>("qjet_LHA_vs_pt", "q-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_pTD_vs_pt = book<TH2F>("qjet_pTD_vs_pt", "q-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_width_vs_pt = book<TH2F>("qjet_width_vs_pt", "q-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qjet_thrust_vs_pt = book<TH2F>("qjet_thrust_vs_pt", "q-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_gjet_multiplicity_vs_pt = book<TH2F>("gjet_multiplicity_vs_pt", "g-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_gjet_LHA_vs_pt = book<TH2F>("gjet_LHA_vs_pt", "g-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_pTD_vs_pt = book<TH2F>("gjet_pTD_vs_pt", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_width_vs_pt = book<TH2F>("gjet_width_vs_pt", "g-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_gjet_thrust_vs_pt = book<TH2F>("gjet_thrust_vs_pt", "g-flavour jet;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

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
  h_PU_pT_hat_max_vs_weight->Fill(weight, event.genInfo->PU_pT_hat_max());
  h_qscale_vs_weight->Fill(weight, event.genInfo->qScale());
  if (event.genInfo->binningValues().size() > 0)
    h_pthat_vs_weight->Fill(weight, event.genInfo->binningValues()[0]);

  // Fill reco jet hists
  for (int i = 0; i < useNJets_; i++) {
    const Jet & thisjet = jets->at(i);
    int mult = thisjet.numberOfDaughters();
    float lha = thisjet.LHA() / LHA_rescale;
    float ptd = thisjet.pTD();
    float width = thisjet.width() / width_rescale;
    float thrust = thisjet.thrust() / thrust_rescale;

    float jet_pt = thisjet.pt();
    h_weights_vs_pt->Fill(weight, jet_pt);
    if (i == 0) {
      h_weights_vs_pt->Fill(weight, thisjet.pt());
      if (event.genInfo->binningValues().size() > 0)
        h_pthat_vs_jet_pt->Fill(thisjet.pt(), event.genInfo->binningValues()[0]);
    }
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

    // int jet_flav = get_jet_flavour(thisjet, event.genparticles);
    int jet_flav = abs(thisjet.genPartonFlavor());

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
    }

    h_jet_flavour->Fill(abs(thisjet.flavor()), weight);
    h_jet_genParton_flavour->Fill(abs(thisjet.genPartonFlavor()), weight);
    h_jet_flavour_vs_pt->Fill(abs(thisjet.flavor()), jet_pt, weight);
    h_jet_genParton_flavour_vs_pt->Fill(abs(thisjet.genPartonFlavor()), jet_pt, weight);
  }

  // Fill GenJet hists
  std::vector<GenJetWithParts>* genjets = event.genjets;
  std::vector<GenParticle>* genparticles = event.genparticles;

  // for (int i = 0; i < useNJets_; i++) {
  //   const GenJetWithParts & thisjet = genjets->at(i);
  int counter = 0;
  for (const auto & thisjet : *genjets) {
    if (!(thisjet.pt() > 100 && fabs(thisjet.eta()) < 1.5))
      continue;

    counter++;
    if (counter > useNJets_)
      break;

    std::vector<GenParticle*> daughters = get_genjet_genparticles(thisjet, genparticles);

    h_genjet_pt->Fill(thisjet.pt(), weight);
    h_genjet_eta->Fill(thisjet.eta(), weight);

    // do special vars according to 1704.03878
    uint mult = 0;
    float ptd = 0, lha = 0, width = 0, thrust = 0;
    float pt_sum = 0;

    for (auto dtr : daughters) {
      pt_sum += dtr->pt();
      mult += 1;
    }

    for (auto dtr : daughters) {
      float z = dtr->pt() / pt_sum;
      float theta = deltaR(dtr->v4(), thisjet.v4()) / jetRadius;
      ptd += pow(z, 2);
      lha += z * pow(theta, 0.5);
      width += z*theta;
      thrust += z * pow(theta, 2);
    }

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

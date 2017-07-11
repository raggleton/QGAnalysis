#include "UHH2/QGAnalysis/include/QGAnalysisTheoryHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisTheoryHists::QGAnalysisTheoryHists(Context & ctx, const string & dirname, int useNJets):
  Hists(ctx, dirname),
  useNJets_(useNJets)
  {

  if (useNJets_ < 0) useNJets_ = 99999; // Do them all

  // book all histograms here
  int nPtBins = 100;
  float ptMin(0.), ptMax(2000.);

  int nEtaBins = 50;
  float etaMin(-2.5), etaMax(2.5);

  int nMultBins = 100;

  int nBins = 100;

  int nWeightBins = 11;
  double weightBins [nWeightBins+1] = {1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 10};
  h_weights = book<TH1F>("weights", ";weight;", nWeightBins, weightBins);
  h_weights_vs_pt = book<TH2F>("weights_vs_pt", ";weight;", nWeightBins, weightBins, 200, 0., 2000.);

  // GENJET hists
  // ------------
  // For all jets
  h_genjet_pt = book<TH1F>("genjet_pt", ";p_{T}^{j} [GeV];", 2*nPtBins, ptMin, 2*ptMax);
  h_genjet_pt_unweighted = book<TH1F>("genjet_pt_unweighted", ";p_{T}^{j} [GeV];", 2*nPtBins, ptMin, 2*ptMax);
  h_genjet_pt_all = book<TH1F>("genjet_pt_all", ";p_{T}^{j} [GeV];", 2*nPtBins, ptMin, 2*ptMax);
  h_genjet_ht = book<TH1F>("genjet_ht", ";H_{T} [GeV];", 2*nPtBins, ptMin, 2*ptMax);
  h_genjet_eta = book<TH1F>("genjet_eta", ";#eta^{j};", nEtaBins, etaMin, etaMax);
  h_genjet_flavour = book<TH1F>("genjet_flavour", "genjet flavour;PDGID;", 23, -0.5, 22.5);
  h_genjet_flavour_vs_pt = book<TH2F>("genjet_flavour_vs_pt", "genjet flavour;PDGID;Jet p_{T} [GeV]", 23, -0.5, 22.5, nPtBins, ptMin, ptMax);

  h_genjet_multiplicity = book<TH1F>("genjet_multiplicity", ";# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_genjet_LHA = book<TH1F>("genjet_LHA", ";LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_genjet_pTD = book<TH1F>("genjet_pTD", ";p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_genjet_width = book<TH1F>("genjet_width", ";Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_genjet_thrust = book<TH1F>("genjet_thrust", ";Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // Flavour-tagged
  // h_qgenjet_multiplicity = book<TH1F>("qgenjet_multiplicity", "q-flavour;# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  // h_qgenjet_LHA = book<TH1F>("qgenjet_LHA", "q-flavour;LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  // h_qgenjet_pTD = book<TH1F>("qgenjet_pTD", "q-flavour;p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  // h_qgenjet_width = book<TH1F>("qgenjet_width", "q-flavour;Width (#lambda_{1}^{1});", nBins, 0, 1);
  // h_qgenjet_thrust = book<TH1F>("qgenjet_thrust", "q-flavour jet;Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // h_ggenjet_multiplicity = book<TH1F>("ggenjet_multiplicity", "g-flavour;# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  // h_ggenjet_LHA = book<TH1F>("ggenjet_LHA", "g-flavour;LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  // h_ggenjet_pTD = book<TH1F>("ggenjet_pTD", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  // h_ggenjet_width = book<TH1F>("ggenjet_width", "g-flavour;Width (#lambda_{1}^{1});", nBins, 0, 1);
  // h_ggenjet_thrust = book<TH1F>("ggenjet_thrust", "g-flavour jet;Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  // 2D versions vs PT
  h_genjet_multiplicity_vs_pt = book<TH2F>("genjet_multiplicity_vs_pt", ";# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_genjet_LHA_vs_pt = book<TH2F>("genjet_LHA_vs_pt", ";LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_genjet_pTD_vs_pt = book<TH2F>("genjet_pTD_vs_pt", ";p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_genjet_width_vs_pt = book<TH2F>("genjet_width_vs_pt", ";Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_genjet_thrust_vs_pt = book<TH2F>("genjet_thrust_vs_pt", ";Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_qgenjet_multiplicity_vs_pt = book<TH2F>("qgenjet_multiplicity_vs_pt", "q-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_qgenjet_LHA_vs_pt = book<TH2F>("qgenjet_LHA_vs_pt", "q-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qgenjet_pTD_vs_pt = book<TH2F>("qgenjet_pTD_vs_pt", "q-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qgenjet_width_vs_pt = book<TH2F>("qgenjet_width_vs_pt", "q-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_qgenjet_thrust_vs_pt = book<TH2F>("qgenjet_thrust_vs_pt", "q-flavour;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  h_ggenjet_multiplicity_vs_pt = book<TH2F>("ggenjet_multiplicity_vs_pt", "g-flavour;# of constituents (#lambda_{0}^{0});Jet p_{T} [GeV]", nMultBins, 0, nMultBins, nPtBins, ptMin, ptMax);
  h_ggenjet_LHA_vs_pt = book<TH2F>("ggenjet_LHA_vs_pt", "g-flavour;LHA (#lambda_{0.5}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_ggenjet_pTD_vs_pt = book<TH2F>("ggenjet_pTD_vs_pt", "g-flavour;p_{T}^{D} (#lambda_{0}^{2});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_ggenjet_width_vs_pt = book<TH2F>("ggenjet_width_vs_pt", "g-flavour;Width (#lambda_{1}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);
  h_ggenjet_thrust_vs_pt = book<TH2F>("ggenjet_thrust_vs_pt", "g-flavour;Thrust (#lambda_{2}^{1});Jet p_{T} [GeV]", nBins, 0, 1, nPtBins, ptMin, ptMax);

  // Some normalisations that weren't done in the initial calculation
  std::string jet_cone = ctx.get("JetCone", "AK4");
  if (jet_cone.find("AK4") != string::npos)
    jetRadius = 0.4;
  else if (jet_cone.find("AK8") != string::npos)
    jetRadius = 0.8;
  else if (jet_cone.find("ca15") != string::npos)
    jetRadius = 1.5;
  else
    throw runtime_error("Cannot determine jetRadius in QGAnalysisTheoryHists");

  genJets_handle = ctx.declare_event_output< std::vector<GenJetWithParts> > ("GoodGenJets");
}


void QGAnalysisTheoryHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'\

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Fill GenJet hists
  // std::vector<GenJetWithParts>* genjets = event.genjets;
  const auto & genjets = event.get(genJets_handle);

  std::vector<GenParticle>* genparticles = event.genparticles;

  float ht = 0.;
  for (const auto & thisjet : genjets) {
    auto pt = thisjet.pt();
    h_genjet_pt_all->Fill(pt, weight);
    ht += pt;
  }
  h_genjet_ht->Fill(ht, weight);

  int counter = 0;
  for (const auto & thisjet : genjets) {
    counter++;
    if (counter > useNJets_)
      break;

    std::vector<GenParticle*> daughters = get_genjet_genparticles(thisjet, genparticles);
    h_weights->Fill(weight);
    h_weights_vs_pt->Fill(weight, thisjet.pt());

    h_genjet_pt->Fill(thisjet.pt(), weight);
    h_genjet_pt_unweighted->Fill(thisjet.pt());
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

    int flav = get_jet_flavour(thisjet, genparticles, jetRadius);
    h_genjet_flavour->Fill(abs(flav), weight);
    h_genjet_flavour_vs_pt->Fill(abs(flav), jet_pt, weight);

    if (abs(flav) == 21) { // gluon jets
      h_ggenjet_multiplicity_vs_pt->Fill(mult, jet_pt, weight);
      h_ggenjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
      h_ggenjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
      h_ggenjet_width_vs_pt->Fill(width, jet_pt, weight);
      h_ggenjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);
    } else if ((abs(flav) <= 3) && (abs(flav) > 0)){ // uds jets
      h_qgenjet_multiplicity_vs_pt->Fill(mult, jet_pt, weight);
      h_qgenjet_LHA_vs_pt->Fill(lha, jet_pt, weight);
      h_qgenjet_pTD_vs_pt->Fill(ptd, jet_pt, weight);
      h_qgenjet_width_vs_pt->Fill(width, jet_pt, weight);
      h_qgenjet_thrust_vs_pt->Fill(thrust, jet_pt, weight);
    }
  }

}


/**
 * Get the collection of GenParticle*s for a given GenJet
 */
std::vector<GenParticle*> QGAnalysisTheoryHists::get_genjet_genparticles(const GenJetWithParts & jet, std::vector<GenParticle>* genparticles) {
  std::vector<GenParticle*> gp;
  for (const uint i : jet.genparticles_indices()) {
    gp.push_back(&(genparticles->at(i)));
  }
  if (gp.size() < 6) {
    std::cout << "low mult jet" << std::endl;
    for (const auto *itr : gp) {std::cout << itr->pdgId() << " : " << itr->pt() << " : " << deltaR(jet.v4(), itr->v4()) << std::endl;}
    // for(const auto & itr : *genparticles) {std::cout << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << deltaR(itr.v4(), jet.v4()) << std::endl;}
  }
  return gp;
}

/**
 *
 */
int QGAnalysisTheoryHists::get_jet_flavour(const GenJetWithParts & jet, std::vector<GenParticle>* genparticles, float dr_max) {
  // cout << "jet:" << endl;
  // cout << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << endl;
  float smallest_dr = dr_max;
  int pdgid = 0;
  for (const auto& gp : *genparticles) {
    if (abs(gp.status()) != 23) continue;
    float dr = deltaR(jet.v4(), gp.v4());
    if (dr < smallest_dr) {
      pdgid = gp.pdgId();
      smallest_dr = dr;
    }
    // cout << gp.pt() << " : " << gp.eta() << " : " << gp.phi() << " : " << gp.pdgId() << " : " << gp.status() << " : " << deltaR(jet.v4(), gp.v4()) << endl;
  }
  return pdgid;
}

QGAnalysisTheoryHists::~QGAnalysisTheoryHists(){}

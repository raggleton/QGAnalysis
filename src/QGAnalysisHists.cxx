#include "UHH2/QGAnalysis/include/QGAnalysisHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisHists::QGAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  int nMultBins = 50;

  // RECO jet hists
  h_jet_multiplicity = book<TH1F>("jet_multiplicity", ";# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  int nBins = 50;
  h_jet_LHA = book<TH1F>("jet_LHA", ";LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_jet_pTD = book<TH1F>("jet_pTD", ";p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_jet_width = book<TH1F>("jet_width", ";Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_jet_thrust = book<TH1F>("jet_thrust", ";Thrust (#lambda_{2}^{1});", nBins, 0, 1);

  h_jet_flavour = book<TH1F>("jet_flavour", "jet flavour;PDGID;", 23, -0.5, 22.5);

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

  // GENJET hists
  h_genjet_multiplicity = book<TH1F>("genjet_multiplicity", ";# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  h_genjet_LHA = book<TH1F>("genjet_LHA", ";LHA (#lambda_{0.5}^{1});", nBins, 0, 1);
  h_genjet_pTD = book<TH1F>("genjet_pTD", ";p_{T}^{D} (#lambda_{0}^{2});", nBins, 0, 1);
  h_genjet_width = book<TH1F>("genjet_width", ";Width (#lambda_{1}^{1});", nBins, 0, 1);
  h_genjet_thrust = book<TH1F>("genjet_thrust", ";Thrust (#lambda_{2}^{1});", nBins, 0, 1);

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
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'\

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Fill reco jet hists
  std::vector<Jet>* jets = event.jets;
  for (const Jet & thisjet : *jets) {
    h_jet_multiplicity->Fill(thisjet.numberOfDaughters(), weight);
    h_jet_LHA->Fill(thisjet.LHA() / LHA_rescale, weight);
    h_jet_pTD->Fill(thisjet.pTD(), weight);
    h_jet_width->Fill(thisjet.width() / width_rescale, weight);
    h_jet_thrust->Fill(thisjet.thrust() / thrust_rescale, weight);

    // Split by actual jet flavour - these only make sense for MC
    if (thisjet.flavor() == 21) { // gluon jets
      h_gjet_multiplicity->Fill(thisjet.numberOfDaughters(), weight);
      h_gjet_LHA->Fill(thisjet.LHA() / LHA_rescale, weight);
      h_gjet_pTD->Fill(thisjet.pTD(), weight);
      h_gjet_width->Fill(thisjet.width() / width_rescale, weight);
      h_gjet_thrust->Fill(thisjet.thrust() / thrust_rescale, weight);
    } else if ((abs(thisjet.flavor()) <= 3) && (abs(thisjet.flavor()) > 0)){ // uds jets
      h_qjet_multiplicity->Fill(thisjet.numberOfDaughters(), weight);
      h_qjet_LHA->Fill(thisjet.LHA() / LHA_rescale, weight);
      h_qjet_pTD->Fill(thisjet.pTD(), weight);
      h_qjet_width->Fill(thisjet.width() / width_rescale, weight);
      h_qjet_thrust->Fill(thisjet.thrust() / thrust_rescale, weight);
    }
    h_jet_flavour->Fill(abs(thisjet.flavor()), weight);
  }

  // Fill GenJet hists
  std::vector<GenJetWithParts>* genjets = event.genjets;
  std::vector< GenParticle >* genparticles = event.genparticles;
  for (const GenJetWithParts & thisjet : * genjets) {

    std::vector<GenParticle*> daughters = get_genjet_genparticles(thisjet, genparticles);

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

float QGAnalysisHists::deltaR(const LorentzVector & a, const LorentzVector & b) {
  // Why am I writing this myself?
  float e1 = a.eta();
  float e2 = b.eta();
  float p1 = a.phi();
  float p2 = b.phi();
  float deta = std::abs(e1-e2);
  float dphi = p1-p2;
  if(dphi>M_PI) dphi -= (2*M_PI);
  return std::hypot(deta, dphi);
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

QGAnalysisHists::~QGAnalysisHists(){}

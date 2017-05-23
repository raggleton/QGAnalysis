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
  h_jet_multiplicity = book<TH1F>("jet_multiplicity", ";# of constituents (#lambda_{0}^{0});", nMultBins, 0, nMultBins);
  int nBins = 25;
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

  std::vector<Jet>* jets = event.jets;
  for (const Jet & thisjet : *jets) {
    h_jet_multiplicity->Fill(thisjet.numberOfDaughters(), weight);
    h_jet_LHA->Fill(thisjet.LHA() / LHA_rescale, weight);
    h_jet_pTD->Fill(thisjet.pTD(), weight);
    h_jet_width->Fill(thisjet.width() / width_rescale, weight);
    h_jet_thrust->Fill(thisjet.thrust() / thrust_rescale, weight);
    if (thisjet.flavor() == 22) {
      h_gjet_multiplicity->Fill(thisjet.numberOfDaughters(), weight);
      h_gjet_LHA->Fill(thisjet.LHA() / LHA_rescale, weight);
      h_gjet_pTD->Fill(thisjet.pTD(), weight);
      h_gjet_width->Fill(thisjet.width() / width_rescale, weight);
      h_gjet_thrust->Fill(thisjet.thrust() / thrust_rescale, weight);
    } else if ((abs(thisjet.flavor()) <= 3) && (abs(thisjet.flavor()) > 0)){
      h_qjet_multiplicity->Fill(thisjet.numberOfDaughters(), weight);
      h_qjet_LHA->Fill(thisjet.LHA() / LHA_rescale, weight);
      h_qjet_pTD->Fill(thisjet.pTD(), weight);
      h_qjet_width->Fill(thisjet.width() / width_rescale, weight);
      h_qjet_thrust->Fill(thisjet.thrust() / thrust_rescale, weight);
    }
    h_jet_flavour->Fill(abs(thisjet.flavor()), weight);
  }

}

QGAnalysisHists::~QGAnalysisHists(){}

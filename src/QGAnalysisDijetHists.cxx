#include "UHH2/QGAnalysis/include/QGAnalysisDijetHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisDijetHists::QGAnalysisDijetHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  n_jets = book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);

  pt_jet1 = book<TH1F>("pt_jet1", "p_{T}^{jet 1} [GeV/c]", 40, 0, 200);
  eta_jet1 = book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  phi_jet1 = book<TH1F>("phi_jet1", "#phi^{jet 1}", 60, -TMath::Pi(), TMath::Pi());

  pt_jet2 = book<TH1F>("pt_jet2", "p_{T}^{jet 2} [GeV/c]", 40, 0, 200);
  eta_jet2 = book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  phi_jet2 = book<TH1F>("phi_jet2", "#phi^{jet 2}", 60, -TMath::Pi(), TMath::Pi());

  m_jj = book<TH1F>("m_jj", "m_{jj} [GeV]", 40, 0, 400);

  deta_dphi_jj = book<TH2F>("deta_dphi_jj", ";#Delta #eta_{jj};#Delta #phi_{jj}", 60, 0, 6, 50, 0, TMath::Pi());

  // primary vertices
  book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);
}


void QGAnalysisDijetHists::fill(const Event & event){
  // fill the histograms. Please note better to
  // use histogram pointers as members instead of hist("name")

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  n_jets->Fill(Njets, weight);

  assert(Njets >= 2);

  Jet jet1 = jets->at(0);
  pt_jet1->Fill(jet1.pt(), weight);
  eta_jet1->Fill(jet1.eta(), weight);
  phi_jet1->Fill(jet1.phi(), weight);

  Jet jet2 = jets->at(1);
  pt_jet2->Fill(jet2.pt(), weight);
  eta_jet2->Fill(jet2.eta(), weight);
  phi_jet2->Fill(jet2.phi(), weight);

  double mass = (jet1.v4() + jet2.v4()).M();
  m_jj->Fill(mass, weight);

  auto diff = jet1.v4() - jet2.v4();
  deta_dphi_jj->Fill(fabs(diff.eta()), fabs(diff.phi()), weight);

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
}

QGAnalysisDijetHists::~QGAnalysisDijetHists(){}

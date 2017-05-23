#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisZPlusJetsHists::QGAnalysisZPlusJetsHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  n_jets = book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);

  // muons
  n_mu = book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);

  int nbins_pt = 100;
  int nbins_eta = 40;
  int nbins_reliso = 50;
  pt_mu1 = book<TH1F>("pt_mu1", "p_{T}^{#mu1} [GeV/c]", nbins_pt, 0, 200);
  eta_mu1 = book<TH1F>("eta_mu1", "#eta^{#mu1}", nbins_eta, -2.1, 2.1);
  reliso_mu1 = book<TH1F>("reliso_mu1", "#mu1 rel. Iso", nbins_reliso, 0, 0.5);

  pt_mu2 = book<TH1F>("pt_mu2", "p_{T}^{#mu2} [GeV/c]", nbins_pt, 0, 200);
  eta_mu2 = book<TH1F>("eta_mu2", "#eta^{#mu2}", nbins_eta, -2.1, 2.1);
  reliso_mu2 = book<TH1F>("reliso_mu2", "#mu2 rel. Iso", nbins_reliso, 0, 0.5);

  m_mumu = book<TH1F>("m_mumu", "m_{#mu#mu} [GeV]", 100, 90-50, 90+50);

  deta_dphi_mumu = book<TH2F>("deta_dphi_mumu", ";#Delta #eta_{#mu#mu};#Delta #phi_{#mu#mu}", 60, 0, 6, 60, 0, TMath::Pi());

  // primary vertices
  book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);

}


void QGAnalysisZPlusJetsHists::fill(const Event & event){
  // fill the histograms. Please note better to
  // use histogram pointers as members instead of hist("name")

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  std::vector<Jet> * jets = event.jets;
  n_jets->Fill(jets->size(), weight);

  std::vector<Muon> * muons = event.muons;
  int Nmuons = muons->size();
  n_mu->Fill(Nmuons, weight);

  assert(Nmuons >= 2);

  Muon mu1 = muons->at(0);
  pt_mu1->Fill(mu1.pt(), weight);
  eta_mu1->Fill(mu1.eta(), weight);
  reliso_mu1->Fill(mu1.relIso(), weight);

  Muon mu2 = muons->at(1);
  pt_mu2->Fill(mu2.pt(), weight);
  eta_mu2->Fill(mu2.eta(), weight);
  reliso_mu2->Fill(mu2.relIso(), weight);

  double mass = (mu1.v4() + mu2.v4()).M();
  m_mumu->Fill(mass, weight);

  auto diff = mu1.v4() - mu2.v4();
  deta_dphi_mumu->Fill(fabs(diff.eta()), fabs(diff.phi()), weight);

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
}

QGAnalysisZPlusJetsHists::~QGAnalysisZPlusJetsHists(){}

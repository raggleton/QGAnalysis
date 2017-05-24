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
  int nbins_pt = 50;
  float pt_max = 1000;
  int nbins_eta = 40;
  float eta_max = 3;
  n_jets = book<TH1F>("N_jets", "N_{jets}", 10, 0, 10);
  pt_jet1 = book<TH1F>("pt_jet1", "p_{T}^{jet 1} [GeV/c]", nbins_pt, 0, pt_max);
  eta_jet1 = book<TH1F>("eta_jet1", "#eta^{jet 1}", nbins_eta, -eta_max, eta_max);
  pt_jet1_z_ratio = book<TH1F>("pt_jet1_z_ratio", "p_{T}^{jet 1} / p_{T}^{Z}", 50, 0, 5);

  pt_max = 500;
  pt_jet2 = book<TH1F>("pt_jet2", "p_{T}^{jet 2} [GeV/c]", nbins_pt, 0, pt_max);
  eta_jet2 = book<TH1F>("eta_jet2", "#eta^{jet 2}", nbins_eta, -eta_max, eta_max);
  pt_jet2_z_ratio = book<TH1F>("pt_jet2_z_ratio", "p_{T}^{jet 2} / p_{T}^{Z}", 50, 0, 1);

  // muons
  n_mu = book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  int nbins_reliso = 60;
  float reliso_max = 0.3;
  pt_mu1 = book<TH1F>("pt_mu1", "p_{T}^{#mu1} [GeV/c]", nbins_pt, 0, pt_max);
  eta_mu1 = book<TH1F>("eta_mu1", "#eta^{#mu1}", nbins_eta, -eta_max, eta_max);
  reliso_mu1 = book<TH1F>("reliso_mu1", "#mu1 rel. Iso", nbins_reliso, 0, reliso_max);

  pt_mu2 = book<TH1F>("pt_mu2", "p_{T}^{#mu2} [GeV/c]", nbins_pt, 0, pt_max);
  eta_mu2 = book<TH1F>("eta_mu2", "#eta^{#mu2}", nbins_eta, -eta_max, eta_max);
  reliso_mu2 = book<TH1F>("reliso_mu2", "#mu2 rel. Iso", nbins_reliso, 0, reliso_max);

  m_mumu = book<TH1F>("m_mumu", "m_{#mu#mu} [GeV]", 80, 90-40, 90+40);
  pt_mumu = book<TH1F>("pt_mumu", "p_{T, #mu#mu} [GeV]", nbins_pt, 0, 750);

  deta_dphi_mumu = book<TH2F>("deta_dphi_mumu", ";#Delta #eta_{#mu#mu};#Delta #phi_{#mu#mu}", 60, 0, 6, 60, 0, TMath::Pi());
  deta_dphi_mumu_jet1 = book<TH2F>("deta_dphi_mumu_jet1", ";#Delta #eta_{#mu#mu, jet 1};#Delta #phi_{#mu#mu, jet1}", 60, 0, 6, 60, 0, TMath::Pi());

  // primary vertices
  book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);

}


void QGAnalysisZPlusJetsHists::fill(const Event & event){
  // fill the histograms. Please note better to
  // use histogram pointers as members instead of hist("name")

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Muons
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

  LorentzVector z_cand = mu1.v4() + mu2.v4();
  m_mumu->Fill(z_cand.M(), weight);
  pt_mumu->Fill(z_cand.pt(), weight);

  auto diff = mu1.v4() - mu2.v4();
  deta_dphi_mumu->Fill(fabs(diff.eta()), fabs(diff.phi()), weight);

  // Jets
  std::vector<Jet> * jets = event.jets;
  int Njets = jets->size();
  n_jets->Fill(Njets, weight);
  assert(Njets >= 1);

  Jet jet1 = jets->at(0);
  pt_jet1->Fill(jet1.pt(), weight);
  eta_jet1->Fill(jet1.eta(), weight);

  pt_jet1_z_ratio->Fill(jet1.pt() / z_cand.pt(), weight);
  diff = z_cand - jet1.v4();
  deta_dphi_mumu_jet1->Fill(fabs(diff.eta()), fabs(diff.phi()), weight);

  if (Njets >= 2) {
    Jet jet2 = jets->at(1);
    pt_jet2->Fill(jet2.pt(), weight);
    eta_jet2->Fill(jet2.eta(), weight);
    pt_jet2_z_ratio->Fill(jet2.pt() / z_cand.pt(), weight);
  }

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
}

QGAnalysisZPlusJetsHists::~QGAnalysisZPlusJetsHists(){}

#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsHists.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisZPlusJetsHists::QGAnalysisZPlusJetsHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  doHerwigReweighting = ctx.get("herwig_reweight_file", "") != "";
  if (doHerwigReweighting) {
    TFile f_weight(ctx.get("herwig_reweight_file", "").c_str());
    reweightHist = (TH1F*) f_weight.Get("zpj_reco");
    if (reweightHist == nullptr) {
      doHerwigReweighting = false;
      cout << "WARNING: could not find zpj_reco reweight hist - not reweighting ZPlusJetsHists!" << endl;
    } else {
      reweightHist->SetDirectory(0);
    }
  }

  // book all histograms here
  // jets
  int nbins_pt = 200;
  float pt_max = 2000;

  int nbins_eta = 200;
  float eta_max = 5;

  int nbins_phi = 128;
  float phi_max = 3.2;

  n_jets_vs_pt_jet1 = book<TH2F>("n_jets_vs_pt_jet1", ";N_{jets};p_{T}^{jet 1}", 10, 0, 10, nbins_pt, 0, pt_max);
  pt_jet1 = book<TH1F>("pt_jet1", ";p_{T}^{jet 1}", nbins_pt, 0, pt_max);
  eta_jet1_vs_pt_jet1 = book<TH2F>("eta_jet1_vs_pt_jet1", ";#eta^{jet 1};p_{T}^{jet 1}", nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  pt_jet1_z_ratio_vs_pt_jet1 = book<TH2F>("pt_jet1_z_ratio_vs_pt_jet1", ";p_{T}^{jet 1} / p_{T}^{Z};p_{T}^{jet 1}", 50, 0, 5, nbins_pt, 0, pt_max);

  pt_jet2_vs_pt_jet1 = book<TH2F>("pt_jet2_vs_pt_jet1", ";p_{T}^{jet 2};p_{T}^{jet 1}", nbins_pt, 0, pt_max, nbins_pt, 0, pt_max);
  eta_jet2_vs_pt_jet1 = book<TH2F>("eta_jet2_vs_pt_jet1", ";#eta^{jet 2};p_{T}^{jet 1}", nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  pt_jet2_z_ratio_vs_pt_jet1 = book<TH2F>("pt_jet2_z_ratio_vs_pt_jet1", ";p_{T}^{jet 2} / p_{T}^{Z};p_{T}^{jet 1}", 60, 0, 3, nbins_pt, 0, pt_max);
  pt_jet1_z_pt_jet2_z_ratio = book<TH2F>("pt_jet1_z_pt_jet2_z_ratio", ";p_{T}^{jet 1} / p_{T}^{Z};p_{T}^{jet 2} / p_{T}^{Z}", 50, 0, 5, 60, 0, 3);

  // muons
  n_mu_vs_pt_jet1 = book<TH2F>("n_mu_vs_pt_jet1", ";N^{#mu};p_{T}^{jet 1}", 10, 0, 10, nbins_pt, 0, pt_max);

  int nbins_reliso = 60;
  float reliso_max = 0.3;

  float mu_pt_max = 1000;

  pt_mu1_vs_pt_jet1 = book<TH2F>("pt_mu1_vs_pt_jet1", ";p_{T}^{#mu1};p_{T}^{jet 1}", nbins_pt, 0, mu_pt_max, nbins_pt, 0, pt_max);
  eta_mu1_vs_pt_jet1 = book<TH2F>("eta_mu1_vs_pt_jet1", ";#eta^{#mu1};p_{T}^{jet 1}", nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  reliso_mu1_vs_pt_jet1 = book<TH2F>("reliso_mu1_vs_pt_jet1", ";#mu1 rel. Iso;p_{T}^{jet 1}", nbins_reliso, 0, reliso_max, nbins_pt, 0, pt_max);

  pt_mu2_vs_pt_jet1 = book<TH2F>("pt_mu2_vs_pt_jet1", ";p_{T}^{#mu2};p_{T}^{jet 1}", nbins_pt, 0, mu_pt_max, nbins_pt, 0, pt_max);
  eta_mu2_vs_pt_jet1 = book<TH2F>("eta_mu2_vs_pt_jet1", ";#eta^{#mu2};p_{T}^{jet 1}", nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  reliso_mu2_vs_pt_jet1 = book<TH2F>("reliso_mu2_vs_pt_jet1", ";#mu2 rel. Iso;p_{T}^{jet 1}", nbins_reliso, 0, reliso_max, nbins_pt, 0, pt_max);

  m_mumu_vs_pt_jet1 = book<TH2F>("m_mumu_vs_pt_jet1", ";m_{#mu#mu} [GeV];p_{T}^{jet 1}", 80, 90-40, 90+40, nbins_pt, 0, pt_max);
  pt_mumu_vs_pt_jet1 = book<TH2F>("pt_mumu_vs_pt_jet1", ";p_{T, #mu#mu} [GeV];p_{T}^{jet 1}", nbins_pt, 0, 750, nbins_pt, 0, pt_max);

  dphi_j_z_vs_pt_jet1 = book<TH2F>("dphi_jet1_z_vs_pt_jet1", ";#Delta #phi_{#mu#mu, jet1};p_{T}^{jet 1}", 60, 0, 6, nbins_pt, 0, pt_max);

  deta_mumu_vs_pt_jet1 = book<TH2F>("deta_mumu_vs_pt_jet1", ";#Delta #eta_{#mu#mu};p_{T}^{jet 1}", nbins_eta, 0, 2*eta_max, nbins_pt, 0, pt_max);
  dphi_mumu_vs_pt_jet1 = book<TH2F>("dphi_mumu_vs_pt_jet1", ";#Delta #phi_{#mu#mu};p_{T}^{jet 1}", nbins_phi, 0, phi_max, nbins_pt, 0, pt_max);
  deta_mumu_jet1_vs_pt_jet1 = book<TH2F>("deta_mumu_jet1_vs_pt_jet1", ";#Delta #eta_{#mu#mu, jet 1};p_{T}^{jet 1}", nbins_eta, 0, 2*eta_max, nbins_pt, 0, pt_max);
  dphi_mumu_jet1_vs_pt_jet1 = book<TH2F>("dphi_mumu_jet1_vs_pt_jet1", ";#Delta #phi_{#mu#mu, jet1};p_{T}^{jet 1}", nbins_phi, 0, phi_max, nbins_pt, 0, pt_max);

  // primary vertices
  book<TH1F>("N_pv", ";N^{PV}", 50, 0, 50);

}


void QGAnalysisZPlusJetsHists::fill(const Event & event){
  std::vector<Jet> * jets = event.jets;
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

  // fill the histograms. Please note better to
  // use histogram pointers as members instead of hist("name")

  // Don't forget to always use the weight when filling.
  double weight = event.weight * herwig_weight;

  float jet1_pt = 0;
  if (Njets >= 1) {
    jet1_pt = jets->at(0).pt();
  }

  // Muons
  std::vector<Muon> * muons = event.muons;
  int Nmuons = muons->size();
  n_mu_vs_pt_jet1->Fill(Nmuons, jet1_pt, weight);

  if (Nmuons < 2) return;

  Muon mu1 = muons->at(0);
  pt_mu1_vs_pt_jet1->Fill(mu1.pt(), jet1_pt, weight);
  eta_mu1_vs_pt_jet1->Fill(mu1.eta(), jet1_pt, weight);
  reliso_mu1_vs_pt_jet1->Fill(mu1.relIso(), jet1_pt, weight);

  Muon mu2 = muons->at(1);
  pt_mu2_vs_pt_jet1->Fill(mu2.pt(), jet1_pt, weight);
  eta_mu2_vs_pt_jet1->Fill(mu2.eta(), jet1_pt, weight);
  reliso_mu2_vs_pt_jet1->Fill(mu2.relIso(), jet1_pt, weight);

  LorentzVector z_cand = mu1.v4() + mu2.v4();
  m_mumu_vs_pt_jet1->Fill(z_cand.M(), jet1_pt, weight);
  pt_mumu_vs_pt_jet1->Fill(z_cand.pt(), jet1_pt, weight);

  auto diff = mu1.v4() - mu2.v4();
  deta_mumu_vs_pt_jet1->Fill(fabs(diff.eta()), jet1_pt,  weight);
  dphi_mumu_vs_pt_jet1->Fill(fabs(diff.phi()), jet1_pt, weight);

  // Jets
  n_jets_vs_pt_jet1->Fill(Njets, jet1_pt, weight);

  if (Njets < 1) return;

  Jet jet1 = jets->at(0);
  pt_jet1->Fill(jet1_pt, weight);
  eta_jet1_vs_pt_jet1->Fill(jet1.eta(), jet1_pt, weight);

  pt_jet1_z_ratio_vs_pt_jet1->Fill(jet1.pt() / z_cand.pt(), jet1_pt, weight);
  diff = z_cand - jet1.v4();
  deta_mumu_jet1_vs_pt_jet1->Fill(fabs(diff.eta()), jet1_pt, weight);
  dphi_mumu_jet1_vs_pt_jet1->Fill(fabs(diff.phi()), jet1_pt, weight);

  dphi_j_z_vs_pt_jet1->Fill(deltaPhi(z_cand, jet1), jet1_pt, weight);

  if (Njets >= 2) {
    Jet jet2 = jets->at(1);
    pt_jet2_vs_pt_jet1->Fill(jet2.pt(), jet1_pt, weight);
    eta_jet2_vs_pt_jet1->Fill(jet2.eta(), jet1_pt, weight);
    pt_jet2_z_ratio_vs_pt_jet1->Fill(jet2.pt() / z_cand.pt(), jet1_pt, weight);
    pt_jet1_z_pt_jet2_z_ratio->Fill(jet1.pt() / z_cand.pt(), jet2.pt() / z_cand.pt(), weight);
  }

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
}

QGAnalysisZPlusJetsHists::~QGAnalysisZPlusJetsHists(){}

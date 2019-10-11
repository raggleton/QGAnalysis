#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsGenHists.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/QGAnalysis/include/QGAnalysisPrinters.h"

#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisZPlusJetsGenHists::QGAnalysisZPlusJetsGenHists(Context & ctx, const string & dirname):
Hists(ctx, dirname)
{
  if (ctx.get("dataset_type") == "MC") {
    genJets_handle = ctx.get_handle< std::vector<GenJetWithParts> > ("GoodGenJets");
    genMuons_handle = ctx.get_handle<std::vector<GenParticle>>("GoodGenMuons");

    int nbins_pt = 500;
    float pt_max = 2000;

    int nbins_eta = 200;
    float eta_max = 5;

    int nbins_phi = 128;
    float phi_max = 3.2;

    float mu_pt_max = 1000;

    // TString zName = "Z";
    TString zName = "#mu#mu";
    TString binByVarLabel = TString::Format("p_{T}^{%s} [GeV]", zName.Data());
    TString binByVar = "pt_z";

    deta_mumu = book<TH1F>("deta_mumu", TString::Format(";|#Delta#eta_{%s}|;", zName.Data()), nbins_eta, 0, 2*eta_max);
    deta_mumu_jet1 = book<TH1F>("deta_mumu_jet1", TString::Format(";|#Delta#eta_{%s, genjet 1}|;", zName.Data()), nbins_eta, 0, 2*eta_max);
    dphi_mumu = book<TH1F>("dphi_mumu", TString::Format(";|#Delta#phi_{%s}|;", zName.Data()), nbins_phi, 0, phi_max);
    dphi_mumu_jet1 = book<TH1F>("dphi_mumu_jet1", TString::Format(";|#Delta#phi_{%s, genjet 1}|;", zName.Data()), nbins_phi, 0, phi_max);
    eta_jet1 = book<TH1F>("eta_jet1", ";#eta^{genjet 1};", nbins_eta, -eta_max, eta_max);
    eta_jet2 = book<TH1F>("eta_jet2", ";#eta^{genjet 2};", nbins_eta, -eta_max, eta_max);
    eta_mu1 = book<TH1F>("eta_mu1", ";#eta^{#mu1};", nbins_eta, -eta_max, eta_max);
    eta_mu2 = book<TH1F>("eta_mu2", ";#eta^{#mu2};", nbins_eta, -eta_max, eta_max);
    eta_mumu = book<TH1F>("eta_mumu", TString::Format(";#eta^{%s};", zName.Data()), nbins_eta, -eta_max, eta_max);
    gen_ht = book<TH1F>("gen_ht", ";H_{T}^{Gen} [GeV];", 500, 0, 5000);
    pt_hat = book<TH1F>("ptHat", ";#hat{p}_{T} [GeV];", 2500, 0, 5000);
    m_mumu = book<TH1F>("m_mumu", TString::Format(";m_{%s} [GeV];", zName.Data()), 80, 90-40, 90+40);
    n_jets = book<TH1F>("n_jets", ";N_{genjets};", 10, 0, 10);
    n_mu = book<TH1F>("n_mu", ";N_{#mu};", 10, 0, 10);
    pt_jet1 = book<TH1F>("pt_jet1", ";p_{T}^{genjet 1} [GeV];", nbins_pt, 0, pt_max);
    pt_jet1_pt_jet2_ratio = book<TH1F>("pt_jet1_pt_jet2_ratio", ";p_{T}^{genjet 1} / p_{T}^{genjet 2};", 200, 0, 10);
    pt_jet1_z_ratio = book<TH1F>("pt_jet1_z_ratio", TString::Format(";p_{T}^{genjet 1}/ p_{T}^{%s};", zName.Data()), 50, 0, 5);
    pt_jet2 = book<TH1F>("pt_jet2", ";p_{T}^{genjet 2} [GeV];", nbins_pt, 0, pt_max);
    pt_jet2_z_ratio = book<TH1F>("pt_jet2_z_ratio", TString::Format(";p_{T}^{genjet 2} / p_{T}^{%s};", zName.Data()), 60, 0, 3);
    pt_jet_genHT_ratio = book<TH1F>("pt_jet_genHT_ratio", ";p_{T}^{genjet 1}/GenHT;", 250, 0, 2.5);
    pt_mu1 = book<TH1F>("pt_mu1", ";p_{T}^{#mu1};", nbins_pt, 0, mu_pt_max);
    pt_mu2 = book<TH1F>("pt_mu2", ";p_{T}^{#mu2};", nbins_pt, 0, mu_pt_max);
    pt_mumu = book<TH1F>("pt_mumu", TString::Format(";p_{T}^{%s} [GeV];", zName.Data()), nbins_pt, 0, mu_pt_max);
    q_scale = book<TH1F>("q_scale", ";q scale [GeV];", 250, 0, 500);
  }

}


void QGAnalysisZPlusJetsGenHists::fill(const Event & event){
  if (event.isRealData) return;
  double weight = event.weight;

  const std::vector<GenJetWithParts> * genjets = &event.get(genJets_handle);
  int Njets = genjets->size();
  n_jets->Fill(Njets, weight);

  const std::vector<GenParticle> * genmuons = &event.get(genMuons_handle);
  int Nmuons = genmuons->size();
  n_mu->Fill(Nmuons, weight);

  if (event.genInfo->binningValues().size() > 0) {
    double ptHat = event.genInfo->binningValues().at(0); // yes this is correct. no idea why
    pt_hat->Fill(ptHat, weight);
  }
  float qScale = event.genInfo->qScale();
  q_scale->Fill(qScale, weight);

  if (Nmuons < 2) return;
  GenParticle mu1 = genmuons->at(0);
  GenParticle mu2 = genmuons->at(1);
  LorentzVector z_cand = mu1.v4() + mu2.v4();
  float z_pt = z_cand.pt();
  pt_mumu->Fill(z_pt, weight);
  eta_mumu->Fill(z_cand.eta(), weight);

  pt_mu1->Fill(mu1.pt(), weight);
  eta_mu1->Fill(mu1.eta(), weight);

  pt_mu2->Fill(mu2.pt(), weight);
  eta_mu2->Fill(mu2.eta(), weight);

  m_mumu->Fill(z_cand.M(), weight);

  deta_mumu->Fill(fabs(mu1.eta() - mu2.eta()), weight);
  dphi_mumu->Fill(fabs(deltaPhi(mu1, mu2)), weight);

  // Jets
  if (Njets < 1) return;

  auto jet1 = genjets->at(0);
  float jet1_pt = jet1.pt();
  pt_jet1->Fill(jet1_pt, weight);
  eta_jet1->Fill(jet1.eta(), weight);

  float genHT = calcGenHT(*(event.genparticles));
  gen_ht->Fill(genHT, weight);
  pt_jet_genHT_ratio->Fill(jet1_pt / genHT, weight);

  pt_jet1_z_ratio->Fill(jet1_pt / z_pt, weight);
  deta_mumu_jet1->Fill(fabs(z_cand.eta() - jet1.eta()), weight);
  dphi_mumu_jet1->Fill(fabs(deltaPhi(jet1, z_cand)), weight);

  if (Njets >= 2) {
    auto jet2 = genjets->at(1);
    pt_jet2->Fill(jet2.pt(), weight);
    eta_jet2->Fill(jet2.eta(), weight);
    pt_jet2_z_ratio->Fill(jet2.pt() / z_pt, weight);
    pt_jet1_pt_jet2_ratio->Fill(jet1_pt / jet2.pt(), weight);
  }

}

QGAnalysisZPlusJetsGenHists::~QGAnalysisZPlusJetsGenHists(){}

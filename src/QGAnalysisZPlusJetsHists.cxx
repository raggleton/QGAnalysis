#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsHists.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisZPlusJetsHists::QGAnalysisZPlusJetsHists(Context & ctx, const string & dirname, const std::string & zLabel_): 
Hists(ctx, dirname),
hndlZ(ctx.get_handle<std::vector<Muon>>(zLabel_))
{
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

  TString zName = "Z";
  // TString zName = "#mu#mu";
  TString binByVarLabel = TString::Format("p_{T}^{%s} [GeV]", zName.Data());
  TString binByVar = "pt_z";

  n_jets_vs_pt = book<TH2F>(TString::Format("n_jets_vs_%s", binByVar.Data()), TString::Format(";N_{jets};%s", binByVarLabel.Data()), 10, 0, 10, nbins_pt, 0, pt_max);
  pt_jet1 = book<TH1F>("pt_jet1", ";p_{T}^{jet 1} [GeV];", nbins_pt, 0, pt_max);
  eta_jet1_vs_pt = book<TH2F>(TString::Format("eta_jet1_vs_%s", binByVar.Data()), TString::Format(";#eta^{jet 1};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  pt_jet1_z_ratio_vs_pt = book<TH2F>(TString::Format("pt_jet1_z_ratio_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 1}/ p_{T}^{%s};%s", zName.Data(), binByVarLabel.Data()), 50, 0, 5, nbins_pt, 0, pt_max);

  gen_ht = book<TH1F>("gen_ht", ";H_{T}^{Gen} [GeV]", 500, 0, 5000);

  pt_jet2_vs_pt = book<TH2F>(TString::Format("pt_jet2_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 2} [GeV];%s", binByVarLabel.Data()), nbins_pt, 0, pt_max, nbins_pt, 0, pt_max);
  eta_jet2_vs_pt = book<TH2F>(TString::Format("eta_jet2_vs_%s", binByVar.Data()), TString::Format(";#eta^{jet 2};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  pt_jet2_z_ratio_vs_pt = book<TH2F>(TString::Format("pt_jet2_z_ratio_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 2} / p_{T}^{%s};%s", zName.Data(), binByVarLabel.Data()), 60, 0, 3, nbins_pt, 0, pt_max);
  pt_jet1_pt_jet2_ratio_vs_pt = book<TH2F>(TString::Format("pt_jet1_pt_jet2_ratio_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 1} / p_{T}^{jet 2};%s", binByVarLabel.Data()), 100, 0, 5, nbins_pt, 0, pt_max);

  // muons
  n_mu_vs_pt = book<TH2F>(TString::Format("n_mu_vs_%s", binByVar.Data()), TString::Format(";N_{#mu};%s", binByVarLabel.Data()), 10, 0, 10, nbins_pt, 0, pt_max);

  int nbins_reliso = 100;
  float reliso_max = 1;

  float mu_pt_max = 1000;

  pt_mu1_vs_pt = book<TH2F>(TString::Format("pt_mu1_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{#mu1};%s", binByVarLabel.Data()), nbins_pt, 0, mu_pt_max, nbins_pt, 0, pt_max);
  eta_mu1_vs_pt = book<TH2F>(TString::Format("eta_mu1_vs_%s", binByVar.Data()), TString::Format(";#eta^{#mu1};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  reliso_mu1_vs_pt = book<TH2F>(TString::Format("reliso_mu1_vs_%s", binByVar.Data()), TString::Format(";#mu1 rel. Iso;%s", binByVarLabel.Data()), nbins_reliso, 0, reliso_max, nbins_pt, 0, pt_max);

  pt_mu2_vs_pt = book<TH2F>(TString::Format("pt_mu2_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{#mu2};%s", binByVarLabel.Data()), nbins_pt, 0, mu_pt_max, nbins_pt, 0, pt_max);
  eta_mu2_vs_pt = book<TH2F>(TString::Format("eta_mu2_vs_%s", binByVar.Data()), TString::Format(";#eta^{#mu2};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  reliso_mu2_vs_pt = book<TH2F>(TString::Format("reliso_mu2_vs_%s", binByVar.Data()), TString::Format(";#mu2 rel. Iso;%s", binByVarLabel.Data()), nbins_reliso, 0, reliso_max, nbins_pt, 0, pt_max);

  m_mumu_vs_pt = book<TH2F>(TString::Format("m_mumu_vs_%s", binByVar.Data()), TString::Format(";m_{%s} [GeV];%s", zName.Data(), binByVarLabel.Data()), 80, 90-40, 90+40, nbins_pt, 0, pt_max);
  pt_jet1_vs_pt = book<TH2F>(TString::Format("pt_jet1_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 1} [GeV];%s", binByVarLabel.Data()), 2*nbins_pt, 0, pt_max, 2*nbins_pt, 0, 2*pt_max);
  pt_mumu = book<TH1F>("pt_mumu", TString::Format(";p_{T}^{%s} [GeV];", zName.Data()), nbins_pt, 0, pt_max);

  dphi_j_z_vs_pt = book<TH2F>(TString::Format("dphi_jet1_z_vs_%s", binByVar.Data()), TString::Format(";|#Delta #phi_{%s, jet 1}|;%s", zName.Data(), binByVarLabel.Data()), 60, 0, 6, nbins_pt, 0, pt_max);

  deta_mumu_vs_pt = book<TH2F>(TString::Format("deta_mumu_vs_%s", binByVar.Data()), TString::Format(";|#Delta #eta_{%s}|;%s", zName.Data(), binByVarLabel.Data()), nbins_eta, 0, 2*eta_max, nbins_pt, 0, pt_max);
  dphi_mumu_vs_pt = book<TH2F>(TString::Format("dphi_mumu_vs_%s", binByVar.Data()), TString::Format(";|#Delta #phi_{%s}|;%s", zName.Data(), binByVarLabel.Data()), nbins_phi, 0, phi_max, nbins_pt, 0, pt_max);
  deta_mumu_jet1_vs_pt = book<TH2F>(TString::Format("deta_mumu_jet1_vs_%s", binByVar.Data()), TString::Format(";|#Delta #eta_{%s, jet 1}|;%s", zName.Data(), binByVarLabel.Data()), nbins_eta, 0, 2*eta_max, nbins_pt, 0, pt_max);
  dphi_mumu_jet1_vs_pt = book<TH2F>(TString::Format("dphi_mumu_jet1_vs_%s", binByVar.Data()), TString::Format(";|#Delta #phi_{%s, jet 1}|;%s", zName.Data(), binByVarLabel.Data()), nbins_phi, 0, phi_max, nbins_pt, 0, pt_max);

  // primary vertices
  n_pv = book<TH1F>("N_pv", ";N_{PV};", 50, 0, 50);

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
  auto & muons = event.get(hndlZ);
  Muon mu1 = muons.at(0);
  Muon mu2 = muons.at(1);
  LorentzVector z_cand = mu1.v4() + mu2.v4();
  float z_pt = z_cand.pt();
  pt_mumu->Fill(z_pt, weight);

  float binPt = z_pt;
  
  int Nmuons = event.muons->size();
  n_mu_vs_pt->Fill(Nmuons, binPt, weight);

  pt_mu1_vs_pt->Fill(mu1.pt(), binPt, weight);
  eta_mu1_vs_pt->Fill(mu1.eta(), binPt, weight);
  reliso_mu1_vs_pt->Fill(mu1.relIso(), binPt, weight);

  pt_mu2_vs_pt->Fill(mu2.pt(), binPt, weight);
  eta_mu2_vs_pt->Fill(mu2.eta(), binPt, weight);
  reliso_mu2_vs_pt->Fill(mu2.relIso(), binPt, weight);

  m_mumu_vs_pt->Fill(z_cand.M(), binPt, weight);
  pt_jet1_vs_pt->Fill(jet1_pt, binPt, weight);

  auto diff = mu1.v4() - mu2.v4();
  deta_mumu_vs_pt->Fill(fabs(diff.eta()), binPt,  weight);
  dphi_mumu_vs_pt->Fill(fabs(diff.phi()), binPt, weight);

  // Jets
  n_jets_vs_pt->Fill(Njets, binPt, weight);

  if (Njets < 1) return;

  Jet & jet1 = jets->at(0);
  pt_jet1->Fill(jet1_pt, weight);
  eta_jet1_vs_pt->Fill(jet1.eta(), binPt, weight);

  if (!event.isRealData) gen_ht->Fill(calcGenHT(*(event.genparticles)), weight);

  pt_jet1_z_ratio_vs_pt->Fill(jet1_pt / z_pt, binPt, weight);
  diff = z_cand - jet1.v4();
  deta_mumu_jet1_vs_pt->Fill(fabs(diff.eta()), binPt, weight);
  dphi_mumu_jet1_vs_pt->Fill(fabs(diff.phi()), binPt, weight);

  dphi_j_z_vs_pt->Fill(deltaPhi(z_cand, jet1), binPt, weight);

  if (Njets >= 2) {
    Jet jet2 = jets->at(1);
    pt_jet2_vs_pt->Fill(jet2.pt(), binPt, weight);
    eta_jet2_vs_pt->Fill(jet2.eta(), binPt, weight);
    pt_jet2_z_ratio_vs_pt->Fill(jet2.pt() / z_pt, binPt, weight);
    pt_jet1_pt_jet2_ratio_vs_pt->Fill(jet1_pt / jet2.pt(), binPt, weight);
  }

  int Npvs = event.pvs->size();
  n_pv->Fill(Npvs, weight);
}

QGAnalysisZPlusJetsHists::~QGAnalysisZPlusJetsHists(){}

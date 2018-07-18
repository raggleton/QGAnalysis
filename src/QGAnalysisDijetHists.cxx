#include "UHH2/QGAnalysis/include/QGAnalysisDijetHists.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisDijetHists::QGAnalysisDijetHists(Context & ctx, const string & dirname, const std::string & binning_): 
  Hists(ctx, dirname), 
  binning(binning_)
{
  doHerwigReweighting = ctx.get("herwig_reweight_file", "") != "";
  if (doHerwigReweighting) {
    TFile f_weight(ctx.get("herwig_reweight_file", "").c_str());
    reweightHist = (TH1F*) f_weight.Get("dijet_reco");
    if (reweightHist == nullptr) {
      doHerwigReweighting = false;
      cout << "WARNING: could not find dijet_reco reweight hist - not reweighting DijetHists!" << endl;
    } else {
      reweightHist->SetDirectory(0);
    }
  }

  if (binning != "ave" && binning != "tnp" && binning != "leading") {
    throw std::runtime_error("Binning should be 'ave', 'tnp', or 'leading'");
  }

  // book all histograms here
  // jets
  int nbins_pt = 400;
  float pt_max = 2000;

  int nbins_eta = 200;
  float eta_max = 5;

  int nbins_phi = 128;
  float phi_max = 3.2;

  TString binByVarLabel = "p_{T}^{jet 1} [GeV]";
  if (binning == "ave") {
    binByVarLabel = "#langle p_{T}^{jet 1, 2} #rangle [GeV]";
  } else if (binning == "tnp") {
    binByVarLabel = "p_{T}^{jet} [GeV]";
  }

  n_jets_vs_pt_jet = book<TH2F>("n_jets_vs_pt_jet", TString::Format(";N_{jets};%s", binByVarLabel.Data()), 10, 0, 10, nbins_pt, 0, pt_max);

  pt_jet = book<TH1F>("pt_jet1", binByVarLabel, nbins_pt, 0, pt_max);
  eta_jet1_vs_pt_jet = book<TH2F>("eta_jet1_vs_pt_jet", TString::Format(";#eta^{jet 1};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  phi_jet1_vs_pt_jet = book<TH2F>("phi_jet1_vs_pt_jet", TString::Format(";#phi^{jet 1};%s", binByVarLabel.Data()), nbins_phi, -phi_max, phi_max, nbins_pt, 0, pt_max);

  eta_jet1_vs_eta_jet2 = book<TH2F>("eta_jet1_vs_eta_jet2", ";#eta^{jet 1};#eta^{jet 2}", nbins_eta, -eta_max, eta_max, nbins_eta,-eta_max, eta_max);
  pt_jet2_vs_pt_jet = book<TH2F>("pt_jet2_vs_pt_jet", TString::Format(";p_{T}^{jet 2};%s", binByVarLabel.Data()), nbins_pt, 0, pt_max, nbins_pt, 0, pt_max);
  eta_jet2_vs_pt_jet = book<TH2F>("eta_jet2_vs_pt_jet", TString::Format(";#eta^{jet 2};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  phi_jet2_vs_pt_jet = book<TH2F>("phi_jet2_vs_pt_jet", TString::Format(";#phi^{jet 2};%s", binByVarLabel.Data()), nbins_phi, -phi_max, phi_max, nbins_pt, 0, pt_max);

  flav_jet1_jet2 = book<TH2F>("flav_jet1_jet2", ";jet 1 flav;jet 2 flav;", 23, -0.5, 22.5, 23, -0.5, 22.5);
  genparton_flav_jet1_jet2 = book<TH2F>("genparton_flav_jet1_jet2", ";jet 1 flav;jet 2 flav;", 23, -0.5, 22.5, 23, -0.5, 22.5);

  pt_jet1_jet2_ratio_vs_pt_jet = book<TH2F>("pt_jet1_jet2_ratio_vs_pt_jet", TString::Format(";p_{T}^{jet 2} / p_{T}^{jet 1};%s", binByVarLabel.Data()), 50, 0, 1, nbins_pt, 0, pt_max);
  jet1_jet2_asym_vs_pt_jet = book<TH2F>("jet1_jet2_asym_vs_pt_jet", TString::Format(";p_{T}^{jet 1} - p_{T}^{jet 2}/p_{T}^{jet 1} + p_{T}^{jet 2};%s", binByVarLabel.Data()), 50, 0, 1, nbins_pt, 0, pt_max);

  m_jj_vs_pt_jet = book<TH2F>("m_jj_vs_pt_jet", TString::Format(";m_{jj} [GeV];%s", binByVarLabel.Data()), 200, 0, 4000, nbins_pt, 0, pt_max);

  deta_jj_vs_pt_jet = book<TH2F>("deta_jj_vs_pt_jet", TString::Format(";|#Delta#eta_{jj}|;%s", binByVarLabel.Data()), nbins_eta, 0, 2*eta_max, nbins_pt, 0, pt_max);
  dphi_jj_vs_pt_jet = book<TH2F>("dphi_jj_vs_pt_jet", TString::Format(";|#Delta#phi_{jj}|;%s", binByVarLabel.Data()), nbins_phi, 0, phi_max, nbins_pt, 0, pt_max);
  sumeta_jj_vs_pt_jet = book<TH2F>("sumeta_jj_vs_pt_jet", TString::Format(";#sum#eta_{jj};%s", binByVarLabel.Data()), 2*nbins_eta, -2*eta_max, 2*eta_max, nbins_pt, 0, pt_max);

  // Possible 3rd jet in the event
  pt_jet3_vs_pt_jet = book<TH2F>("pt_jet3_vs_pt_jet", TString::Format(";p_{T}^{jet 3};%s", binByVarLabel.Data()), nbins_pt, 0, 500, nbins_pt, 0, pt_max);
  eta_jet3_vs_pt_jet = book<TH2F>("eta_jet3_vs_pt_jet", TString::Format(";#eta^{jet 3};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  pt_jet3_frac_vs_pt_jet = book<TH2F>("pt_jet3_frac_vs_pt_jet", TString::Format(";p_{T}^{jet 3} / #LT p_{T}^{jet 1}, p_{T}^{jet 2} #GT;%s", binByVarLabel.Data()), 50, 0, 1, nbins_pt, 0, pt_max);

  // MET
  int nbins_metSig(50);
  float metSig_max(10.);
  met_sig_vs_pt_jet = book<TH2F>("met_sig_vs_pt_jet", TString::Format(";MET signif.;%s", binByVarLabel.Data()), nbins_metSig, 0, metSig_max, nbins_pt, 0, pt_max);

  // primary vertices
  n_pv = book<TH1F>("N_pv", ";N^{PV};", 50, 0, 50);
}


void QGAnalysisDijetHists::fill(const Event & event){
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  if (Njets < 2) return;

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

  Jet jet1 = jets->at(0);
  double jet1_pt = jet1.pt();
  Jet jet2 = jets->at(1);
  double jet2_pt = jet2.pt();

  double binByVal = 0.;
  if (binning == "ave") {
    binByVal = (jet1_pt + jet2_pt)/2.;
  } else if (binning == "leading") {
    binByVal = jet1_pt;
  }

  n_jets_vs_pt_jet->Fill(Njets, binByVal, weight);
  pt_jet->Fill(binByVal, weight);
  eta_jet1_vs_pt_jet->Fill(jet1.eta(), binByVal, weight);
  phi_jet1_vs_pt_jet->Fill(jet1.phi(), binByVal, weight);
  
  eta_jet1_vs_eta_jet2->Fill(jet1.eta(), jet2.eta(), weight);
  pt_jet2_vs_pt_jet->Fill(jet2_pt, binByVal, weight);
  eta_jet2_vs_pt_jet->Fill(jet2.eta(), binByVal, weight);
  phi_jet2_vs_pt_jet->Fill(jet2.phi(), binByVal, weight);

  pt_jet1_jet2_ratio_vs_pt_jet->Fill(jet2_pt / jet1_pt, binByVal, weight);
  jet1_jet2_asym_vs_pt_jet->Fill((jet1_pt - jet2_pt) / (jet1_pt + jet2_pt), binByVal, weight);

  flav_jet1_jet2->Fill(abs(jet1.flavor()), abs(jet2.flavor()), weight);
  genparton_flav_jet1_jet2->Fill(abs(jet1.genPartonFlavor()), abs(jet2.genPartonFlavor()), weight);

  double mass = (jet1.v4() + jet2.v4()).M();
  m_jj_vs_pt_jet->Fill(mass, binByVal, weight);

  double dEta = fabs(jet1.eta() - jet2.eta());
  double dPhi = fabs(deltaPhi(jet1, jet2));
  deta_jj_vs_pt_jet->Fill(dEta, binByVal, weight);
  dphi_jj_vs_pt_jet->Fill(dPhi, binByVal, weight);
  sumeta_jj_vs_pt_jet->Fill(jet1.eta() + jet2.eta(), binByVal, weight);

  met_sig_vs_pt_jet->Fill(event.met->mEtSig(), binByVal, weight);

  if (Njets >= 3) {
    Jet jet3 = jets->at(2);
    pt_jet3_vs_pt_jet->Fill(jet3.pt(), binByVal, weight);
    eta_jet3_vs_pt_jet->Fill(jet3.eta(), binByVal, weight);
    pt_jet3_frac_vs_pt_jet->Fill(jet3.pt() / (0.5 * (jet1_pt + jet2_pt)), binByVal, weight);
  } else {
    pt_jet3_frac_vs_pt_jet->Fill(0.0, binByVal, weight);
  }

  int Npvs = event.pvs->size();
  n_pv->Fill(Npvs, weight);
}

QGAnalysisDijetHists::~QGAnalysisDijetHists(){}

#include "UHH2/QGAnalysis/include/QGAnalysisDijetHists.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisDijetHists::QGAnalysisDijetHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
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

  // book all histograms here
  // jets
  int nbins_pt = 200;
  float pt_max = 2000;

  int nbins_eta = 200;
  float eta_max = 5;

  int nbins_phi = 128;
  float phi_max = 3.2;

  n_jets_vs_pt_jet1 = book<TH2F>("n_jets_vs_pt_jet1", ";N_{jets};p_{T}^{jet 1}", 10, 0, 10, nbins_pt, 0, pt_max);

  pt_jet1 = book<TH1F>("pt_jet1", ";p_{T}^{jet 1};", nbins_pt, 0, pt_max);
  eta_jet1_vs_pt_jet1 = book<TH2F>("eta_jet1_vs_pt_jet1", ";#eta^{jet 1};p_{T}^{jet 1}", nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  phi_jet1_vs_pt_jet1 = book<TH2F>("phi_jet1_vs_pt_jet1", ";#phi^{jet 1};p_{T}^{jet 1}", nbins_phi, -phi_max, phi_max, nbins_pt, 0, pt_max);

  eta_jet1_vs_eta_jet2 = book<TH2F>("eta_jet1_vs_eta_jet2", ";#eta^{jet 1};#eta^{jet 2}", nbins_eta, -eta_max, eta_max, nbins_eta,-eta_max, eta_max);
  pt_jet2_vs_pt_jet1 = book<TH2F>("pt_jet2_vs_pt_jet1", ";p_{T}^{jet 2};p_{T}^{jet 1}", nbins_pt, 0, pt_max, nbins_pt, 0, pt_max);
  eta_jet2_vs_pt_jet1 = book<TH2F>("eta_jet2_vs_pt_jet1", ";#eta^{jet 2};p_{T}^{jet 1}", nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  phi_jet2_vs_pt_jet1 = book<TH2F>("phi_jet2_vs_pt_jet1", ";#phi^{jet 2};p_{T}^{jet 1}", nbins_phi, -phi_max, phi_max, nbins_pt, 0, pt_max);

  flav_jet1_jet2 = book<TH2F>("flav_jet1_jet2", ";jet 1 flav;jet 2 flav;", 23, -0.5, 22.5, 23, -0.5, 22.5);
  genparton_flav_jet1_jet2 = book<TH2F>("genparton_flav_jet1_jet2", ";jet 1 flav;jet 2 flav;", 23, -0.5, 22.5, 23, -0.5, 22.5);

  pt_jet1_jet2_ratio_vs_pt_jet1 = book<TH2F>("pt_jet1_jet2_ratio_vs_pt_jet1", ";p_{T}^{jet 2} / p_{T}^{jet 1};p_{T}^{jet 1};p_{T}^{jet 1}", 50, 0, 1, nbins_pt, 0, pt_max);

  m_jj_vs_pt_jet1 = book<TH2F>("m_jj_vs_pt_jet1", ";m_{jj} [GeV];p_{T}^{jet 1}", 200, 0, 4000, nbins_pt, 0, pt_max);

  deta_jj_vs_pt_jet1 = book<TH2F>("deta_jj_vs_pt_jet1", ";|#Delta#eta_{jj}|;p_{T}^{jet 1}", nbins_eta, 0, 2*eta_max, nbins_pt, 0, pt_max);
  dphi_jj_vs_pt_jet1 = book<TH2F>("dphi_jj_vs_pt_jet1", ";|#Delta#phi_{jj}|;p_{T}^{jet 1}", nbins_phi, 0, phi_max, nbins_pt, 0, pt_max);
  sumeta_jj_vs_pt_jet1 = book<TH2F>("sumeta_jj_vs_pt_jet1", ";#Sum#eta_{jj};p_{T}^{jet 1}", 2*nbins_eta, -2*eta_max, 2*eta_max, nbins_pt, 0, pt_max);

  // Possible 3rd jet in the event
  pt_jet3_vs_pt_jet1 = book<TH2F>("pt_jet3_vs_pt_jet1", ";p_{T}^{jet 3};p_{T}^{jet 1}", nbins_pt, 0, 500, nbins_pt, 0, pt_max);
  eta_jet3_vs_pt_jet1 = book<TH2F>("eta_jet3_vs_pt_jet1", ";#eta^{jet 3};p_{T}^{jet 1}", nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  pt_jet3_frac_vs_pt_jet1 = book<TH2F>("pt_jet3_frac_vs_pt_jet1", ";p_{T}^{jet 3} / #LT p_{T}^{jet 1}, p_{T}^{jet 2} #GT;p_{T}^{jet 1}", 50, 0, 1, nbins_pt, 0, pt_max);

  // MET
  int nbins_metSig(50);
  float metSig_max(10.);
  met_sig_vs_pt_jet1 = book<TH2F>("met_sig_vs_pt_jet1", ";MET/sumET;p_{T}^{jet 1}", nbins_metSig, 0, metSig_max, nbins_pt, 0, pt_max);

  // primary vertices
  book<TH1F>("N_pv", ";N^{PV};", 50, 0, 50);
}


void QGAnalysisDijetHists::fill(const Event & event){
  std::vector<Jet>* jets = event.jets;
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

  Jet jet1 = jets->at(0);
  double jet1_pt = jet1.pt();

  n_jets_vs_pt_jet1->Fill(Njets, jet1_pt, weight);

  if (Njets < 2) return;

  pt_jet1->Fill(jet1.pt(), weight);
  eta_jet1_vs_pt_jet1->Fill(jet1.eta(), jet1_pt, weight);
  phi_jet1_vs_pt_jet1->Fill(jet1.phi(), jet1_pt, weight);
  
  Jet jet2 = jets->at(1);
  eta_jet1_vs_eta_jet2->Fill(jet1.eta(), jet2.eta(), weight);
  pt_jet2_vs_pt_jet1->Fill(jet2.pt(), jet1_pt, weight);
  eta_jet2_vs_pt_jet1->Fill(jet2.eta(), jet1_pt, weight);
  phi_jet2_vs_pt_jet1->Fill(jet2.phi(), jet1_pt, weight);

  pt_jet1_jet2_ratio_vs_pt_jet1->Fill(jet2.pt() / jet1.pt(), jet1_pt, weight);

  flav_jet1_jet2->Fill(abs(jet1.flavor()), abs(jet2.flavor()), weight);
  genparton_flav_jet1_jet2->Fill(abs(jet1.genPartonFlavor()), abs(jet2.genPartonFlavor()), weight);

  double mass = (jet1.v4() + jet2.v4()).M();
  m_jj_vs_pt_jet1->Fill(mass, jet1_pt, weight);

  double dEta = fabs(jet1.eta() - jet2.eta());
  double dPhi = fabs(deltaPhi(jet1, jet2));
  deta_jj_vs_pt_jet1->Fill(dEta, jet1_pt, weight);
  dphi_jj_vs_pt_jet1->Fill(dPhi, jet1_pt, weight);
  sumeta_jj_vs_pt_jet1->Fill(jet1.eta() + jet2.eta(), jet1_pt, weight);

  met_sig_vs_pt_jet1->Fill(event.met->mEtSig(), jet1_pt, weight);

  if (Njets >= 3) {
    Jet jet3 = jets->at(2);
    pt_jet3_vs_pt_jet1->Fill(jet3.pt(), jet1_pt, weight);
    eta_jet3_vs_pt_jet1->Fill(jet3.eta(), jet1_pt, weight);
    pt_jet3_frac_vs_pt_jet1->Fill(jet3.pt() / (0.5 * (jet1.pt() + jet2.pt())), jet1_pt, weight);
  } else {
    pt_jet3_frac_vs_pt_jet1->Fill(0.0, jet1_pt, weight);
  }

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
}

QGAnalysisDijetHists::~QGAnalysisDijetHists(){}

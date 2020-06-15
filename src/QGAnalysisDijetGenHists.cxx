#include "UHH2/QGAnalysis/include/QGAnalysisDijetGenHists.h"
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

QGAnalysisDijetGenHists::QGAnalysisDijetGenHists(Context & ctx, const string & dirname, const string & genjets_name):
Hists(ctx, dirname)
{
  if (ctx.get("dataset_type") == "MC") {
    genJets_handle = ctx.get_handle< std::vector<GenJet> > (genjets_name);
    gen_weight_handle = ctx.get_handle<double>("gen_weight");

    int nbins_pt = 500;
    float pt_max = 2000;

    int nbins_eta = 200;
    float eta_max = 5;

    int nbins_phi = 128;
    float phi_max = 3.2;

    int nbins_r = 300;
    float r_max = 6;

    deta_jet = book<TH1F>("deta_jet", ";|#Deltay_{jet, jet}|;", nbins_eta, 0, 2*eta_max);
    dphi_jet = book<TH1F>("dphi_jet", ";|#Delta#phi_{jet, jet}|;", nbins_phi, 0, phi_max);
    dr_jet = book<TH1F>("dr_jet", ";|#DeltaR_{jet, jet}|;", nbins_r, 0, r_max);
    eta_jet1 = book<TH1F>("eta_jet1", ";y^{genjet 1};", nbins_eta, -eta_max, eta_max);
    eta_jet2 = book<TH1F>("eta_jet2", ";y^{genjet 2};", nbins_eta, -eta_max, eta_max);
    gen_ht = book<TH1F>("gen_ht", ";H_{T}^{Gen} [GeV];", 500, 0, 5000);
    pt_hat = book<TH1F>("ptHat", ";#hat{p}_{T} [GeV];", 2500, 0, 5000);
    n_jets = book<TH1F>("n_jets", ";N_{genjets};", 10, 0, 10);
    pt_jet1 = book<TH1F>("pt_jet1", ";p_{T}^{genjet 1} [GeV];", nbins_pt, 0, pt_max);
    pt_jet1_pt_jet2_ratio = book<TH1F>("pt_jet1_pt_jet2_ratio", ";p_{T}^{genjet 1} / p_{T}^{genjet 2};", 200, 0, 10);
    pt_jet2 = book<TH1F>("pt_jet2", ";p_{T}^{genjet 2} [GeV];", nbins_pt, 0, pt_max);
    pt_dijet_ave = book<TH1F>("pt_dijet_ave", ";#LTp_{T}^{genjet 1,2}#GT [GeV];", nbins_pt, 0, pt_max);
    pt_jet_genHT_ratio = book<TH1F>("pt_jet_genHT_ratio", ";p_{T}^{genjet 1}/GenHT;", 250, 0, 2.5);
    pt_hat_pt_jet_ratio = book<TH1F>("pt_hat_pt_jet_ratio", ";p_{T}^{genjet 1}/#hat{p}_{T};", 250, 0, 2.5);
    pt_asym = book<TH1F>("pt_asym", ";(p_{T}^{genjet 1} - p_{T}^{genjet 2}) / (p_{T}^{genjet 1} + p_{T}^{genjet 2});", 200, 0, 1);
    q_scale = book<TH1F>("q_scale", ";q scale [GeV];", 250, 0, 500);
  }

}


void QGAnalysisDijetGenHists::fill(const Event & event){
  if (event.isRealData) return;
  double weight = event.get(gen_weight_handle);

  const std::vector<GenJet> * genjets = &event.get(genJets_handle);
  int Njets = genjets->size();
  n_jets->Fill(Njets, weight);

  if (event.genInfo->binningValues().size() > 0) {
    double ptHat = event.genInfo->binningValues().at(0); // sometimes this is pthat, sometimes it means the hardest outgoing partoneg in H++
    pt_hat->Fill(ptHat, weight);
  }

  float qScale = event.genInfo->qScale();
  q_scale->Fill(qScale, weight);

  if (Njets < 2) return;

  auto jet1 = genjets->at(0);
  auto jet2 = genjets->at(1);
  
  float jet1_pt = jet1.pt();
  float jet2_pt = jet2.pt();

  deta_jet->Fill(fabs(jet1.Rapidity() - jet2.Rapidity()), weight);
  dphi_jet->Fill(fabs(deltaPhi(jet1, jet2)), weight);
  dr_jet->Fill(uhh2::deltaR(jet1, jet2), weight);

  pt_jet1->Fill(jet1_pt, weight);
  eta_jet1->Fill(jet1.Rapidity(), weight);
  if (event.genInfo->binningValues().size() > 0) {
    double ptHat = event.genInfo->binningValues().at(0);
    pt_hat_pt_jet_ratio->Fill(jet1_pt / ptHat, weight);
  }
  pt_jet2->Fill(jet2_pt, weight);
  eta_jet2->Fill(jet2.Rapidity(), weight);
  pt_dijet_ave->Fill(0.5*(jet1_pt + jet2_pt), weight);

  float genHT = calcGenHT(*(event.genparticles));
  gen_ht->Fill(genHT, weight);
  pt_jet_genHT_ratio->Fill(jet1_pt / genHT, weight);

  pt_jet1_pt_jet2_ratio->Fill(jet1_pt / jet2.pt(), weight);
  pt_asym->Fill((jet1_pt-jet2_pt) / (jet1_pt+jet2_pt), weight);

}

QGAnalysisDijetGenHists::~QGAnalysisDijetGenHists(){}

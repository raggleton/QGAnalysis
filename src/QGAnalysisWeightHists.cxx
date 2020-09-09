#include "UHH2/QGAnalysis/include/QGAnalysisWeightHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisPrinters.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisWeightHists::QGAnalysisWeightHists(Context & ctx, const string & dirname):
  Hists(ctx, dirname)
{
  is_mc_ = ctx.get("dataset_type") == "MC";
  if (is_mc_) genJets_handle = ctx.get_handle< std::vector<GenJet> > ("GoodGenJets");

  // book all histograms here
  // std::vector<double> pt_bin_edges = Binning::pt_bin_edges_reco_all;
  // int nbins_pt = pt_bin_edges.size()-1;
  vector<double> pt_bin_edges;
  float ptUpperLim = 2000;
  float interval = 2;
  for (float pt=15; pt<=ptUpperLim; pt+=interval) {
    pt_bin_edges.push_back(pt);
  }
  int nbins_pt = pt_bin_edges.size()-1;


  // vector<double> weight_bin_edges = {-6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9};
  vector<double> weight_bin_edges;
  float weightUpperLim = 9;
  interval = 15;
  for (float w=-6; w<=weightUpperLim; w+=interval) {
    weight_bin_edges.push_back(w);
  }
  int nbins_weight = weight_bin_edges.size()-1;

  vector<double> ratio_bin_edges;
  // first add fine binning
  float upperLim = 10.;
  interval = 0.1;
  for (float r=0; r<=upperLim; r+=interval) {
    ratio_bin_edges.push_back(r);
  }
  // now add coarser binning
  upperLim = 100.;
  interval = 5;
  for (float r=ratio_bin_edges.back()+interval; r<=upperLim; r+=interval) {
    ratio_bin_edges.push_back(r);
  }
  int nbins_ratio = ratio_bin_edges.size()-1;

  // for (auto w : weight_bin_edges) {
  //   cout << "w: " << w << endl;
  // }
  // for (auto r : ratio_bin_edges) {
  //   cout << "r: " << r << endl;
  // }


  weight_vs_pt_vs_pt_jet_genHT_ratio = book<TH3F>("weight_vs_pt_vs_pt_jet_genHT_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_jet_genHT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_pt_jet_genHT_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_pt_jet_genHT_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_jet_genHT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  
  // plot ve pt (reco) and pt(gen), since both are of importance  
  weight_vs_pt_vs_pt_genjet_genHT_ratio = book<TH3F>("weight_vs_pt_vs_pt_genjet_genHT_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_genjet_genHT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio = book<TH3F>("weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio", ";log_{10}(weight);p^{GenJet}_{T} [GeV];pt_genjet_genHT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_pt_genjet_genHT_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_pt_genjet_genHT_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_genjet_genHT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio_unweighted = book<TH3F>("weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio_unweighted", ";log_{10}(weight);p^{GenJet}_{T} [GeV];pt_genjet_genHT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  
  weight_vs_pt_vs_pt_jet_qScale_ratio = book<TH3F>("weight_vs_pt_vs_pt_jet_qScale_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_jet_qScale_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_pt_jet_qScale_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_pt_jet_qScale_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_jet_qScale_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);

  // plot ve pt (reco) and pt(gen), since both are of importance  
  weight_vs_pt_vs_pt_genjet_qScale_ratio = book<TH3F>("weight_vs_pt_vs_pt_genjet_qScale_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_genjet_qScale_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio = book<TH3F>("weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio", ";log_{10}(weight);p^{GenJet}_{T} [GeV];pt_genjet_qScale_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_pt_genjet_qScale_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_pt_genjet_qScale_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_genjet_qScale_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio_unweighted = book<TH3F>("weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio_unweighted", ";log_{10}(weight);p^{GenJet}_{T} [GeV];pt_genjet_qScale_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  
  weight_vs_pt_vs_pt_jet_ptHat_ratio = book<TH3F>("weight_vs_pt_vs_pt_jet_ptHat_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_jet_ptHat_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_pt_jet_ptHat_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_pt_jet_ptHat_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_jet_ptHat_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);

  // plot ve pt (reco) and pt(gen), since both are of importance  
  weight_vs_pt_vs_pt_genjet_ptHat_ratio = book<TH3F>("weight_vs_pt_vs_pt_genjet_ptHat_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_genjet_ptHat_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio = book<TH3F>("weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio", ";log_{10}(weight);p^{GenJet}_{T} [GeV];pt_genjet_ptHat_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_pt_genjet_ptHat_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_pt_genjet_ptHat_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_genjet_ptHat_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio_unweighted = book<TH3F>("weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio_unweighted", ";log_{10}(weight);p^{GenJet}_{T} [GeV];pt_genjet_ptHat_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  
  weight_vs_pt_vs_PU_ptHat_genHT_ratio = book<TH3F>("weight_vs_pt_vs_PU_ptHat_genHT_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];PU_ptHat_genHT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_PU_ptHat_genHT_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_PU_ptHat_genHT_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];PU_ptHat_genHT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);

  weight_vs_pt_vs_PU_ptHat_ptHat_ratio = book<TH3F>("weight_vs_pt_vs_PU_ptHat_ptHat_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];PU_ptHat_ptHat_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_PU_ptHat_ptHat_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_PU_ptHat_ptHat_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];PU_ptHat_ptHat_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);

  weight_vs_pt_vs_pt_jet_jetkT_ratio = book<TH3F>("weight_vs_pt_vs_pt_jet_jetkT_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_jet_jetkT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_pt_jet_jetkT_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_pt_jet_jetkT_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];pt_jet_jetkT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);

  weight_vs_pt_vs_PU_ptHat_jetkT_ratio = book<TH3F>("weight_vs_pt_vs_PU_ptHat_jetkT_ratio", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];PU_ptHat_jetkT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);
  weight_vs_pt_vs_PU_ptHat_jetkT_ratio_unweighted = book<TH3F>("weight_vs_pt_vs_PU_ptHat_jetkT_ratio_unweighted", ";log_{10}(weight);p^{RecoJet}_{T} [GeV];PU_ptHat_jetkT_ratio", nbins_weight, &weight_bin_edges[0], nbins_pt, &pt_bin_edges[0], nbins_ratio, &ratio_bin_edges[0]);

}


void QGAnalysisWeightHists::fill(const Event & event){
  if (event.isRealData) return;
  float genHT = calcGenHT(*(event.genparticles));
  float qScale = event.genInfo->qScale();

  const auto & genjets = event.get(genJets_handle);

  bool hasRecoJets = event.jets->size() > 0;
  float jet_pt = hasRecoJets ? event.jets->at(0).pt() : -1;

  bool hasGenJets = genjets.size() > 0;
  float genjet_pt = hasGenJets ? genjets.at(0).pt() : -1;

  float weight = event.weight;
  float log_weight = TMath::Log10(weight);
  float PU_pThat = event.genInfo->PU_pT_hat_max();

  if (genHT > 0) {
    if (hasRecoJets) {
      weight_vs_pt_vs_pt_jet_genHT_ratio->Fill(log_weight, jet_pt, jet_pt / genHT, weight);
      weight_vs_pt_vs_pt_jet_genHT_ratio_unweighted->Fill(log_weight, jet_pt, jet_pt / genHT);

      weight_vs_pt_vs_PU_ptHat_genHT_ratio->Fill(log_weight, jet_pt, PU_pThat / genHT, weight);
      weight_vs_pt_vs_PU_ptHat_genHT_ratio_unweighted->Fill(log_weight, jet_pt, PU_pThat / genHT);
    }
    if (hasGenJets) {
      weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio->Fill(log_weight, genjet_pt, genjet_pt / genHT, weight);
      weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio_unweighted->Fill(log_weight, genjet_pt, genjet_pt / genHT);
    }
    if (hasGenJets && hasRecoJets) {
      weight_vs_pt_vs_pt_genjet_genHT_ratio->Fill(log_weight, jet_pt, genjet_pt / genHT, weight);
      weight_vs_pt_vs_pt_genjet_genHT_ratio_unweighted->Fill(log_weight, jet_pt, genjet_pt / genHT);
    }
  }

  if (hasRecoJets) {
    weight_vs_pt_vs_pt_jet_qScale_ratio->Fill(log_weight, jet_pt, jet_pt / qScale, weight);
    weight_vs_pt_vs_pt_jet_qScale_ratio_unweighted->Fill(log_weight, jet_pt, jet_pt / qScale);
  }
  if (hasGenJets) {
    weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio->Fill(log_weight, genjet_pt, genjet_pt / qScale, weight);
    weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio_unweighted->Fill(log_weight, genjet_pt, genjet_pt / qScale);
  }
  if (hasGenJets && hasRecoJets) {
    weight_vs_pt_vs_pt_genjet_qScale_ratio->Fill(log_weight, jet_pt, genjet_pt / qScale, weight);
    weight_vs_pt_vs_pt_genjet_qScale_ratio_unweighted->Fill(log_weight, jet_pt, genjet_pt / qScale);
  }

  if (event.genInfo->binningValues().size() > 0) {
    double ptHat = event.genInfo->binningValues().at(0); // yes this is correct. no idea why
    if (hasRecoJets) {
      weight_vs_pt_vs_pt_jet_ptHat_ratio->Fill(log_weight, jet_pt, jet_pt / ptHat, weight);
      weight_vs_pt_vs_pt_jet_ptHat_ratio_unweighted->Fill(log_weight, jet_pt, jet_pt / ptHat);

      weight_vs_pt_vs_PU_ptHat_ptHat_ratio->Fill(log_weight, jet_pt, PU_pThat / ptHat, weight);
      weight_vs_pt_vs_PU_ptHat_ptHat_ratio_unweighted->Fill(log_weight, jet_pt, PU_pThat / ptHat);
    }
    if (hasGenJets) {
      weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio->Fill(log_weight, genjet_pt, genjet_pt / ptHat, weight);
      weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio_unweighted->Fill(log_weight, genjet_pt, genjet_pt / ptHat);
    }
    if (hasGenJets && hasRecoJets) {
      weight_vs_pt_vs_pt_genjet_ptHat_ratio->Fill(log_weight, jet_pt, genjet_pt / ptHat, weight);
      weight_vs_pt_vs_pt_genjet_ptHat_ratio_unweighted->Fill(log_weight, jet_pt, genjet_pt / ptHat);
    }
  }

  float jetkT = calcJetKt(*event.genparticles);
  if (jetkT > 0) {
    weight_vs_pt_vs_pt_jet_jetkT_ratio->Fill(log_weight, jet_pt, jet_pt / jetkT, weight);
    weight_vs_pt_vs_pt_jet_jetkT_ratio_unweighted->Fill(log_weight, jet_pt, jet_pt / jetkT);
    weight_vs_pt_vs_PU_ptHat_jetkT_ratio->Fill(log_weight, jet_pt, PU_pThat / jetkT, weight);
    weight_vs_pt_vs_PU_ptHat_jetkT_ratio_unweighted->Fill(log_weight, jet_pt, PU_pThat / jetkT);
  }

}

QGAnalysisWeightHists::~QGAnalysisWeightHists(){}

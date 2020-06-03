#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsHists.h"
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

QGAnalysisZPlusJetsHists::QGAnalysisZPlusJetsHists(Context & ctx, const string & dirname, const string & zLabel_, const string & genjets_name):
Hists(ctx, dirname),
hndlZ(ctx.get_handle<std::vector<Muon>>(zLabel_))
{
  string jet_cone = ctx.get("JetCone", "AK4");
  jetRadius = get_jet_radius(jet_cone);

  if (ctx.get("dataset_type") == "MC") genJets_handle = ctx.get_handle< std::vector<GenJetWithParts> > (genjets_name);

  // book all histograms here
  // jets
  int nbins_pt = 2000;
  float pt_max = 2000;

  int nbins_eta = 200;
  float eta_max = 5;

  int nbins_phi = 128;
  float phi_max = 3.2;

  // TString zName = "Z";
  TString zName = "#mu#mu";
  TString binByVarLabel = TString::Format("p_{T}^{%s} [GeV]", zName.Data());
  TString binByVar = "pt_z";

  n_jets_vs_pt = book<TH2F>(TString::Format("n_jets_vs_%s", binByVar.Data()), TString::Format(";N_{jets};%s", binByVarLabel.Data()), 10, 0, 10, nbins_pt, 0, pt_max);
  pt_jet1 = book<TH1F>("pt_jet1", ";p_{T}^{jet 1} [GeV];", nbins_pt, 0, pt_max);
  pt_jet1_unweighted = book<TH1F>("pt_jet1_unweighted", ";p_{T}^{jet 1} [GeV];", nbins_pt, 0, pt_max);
  pt_jet_response_binning = book<TH1F>("pt_jet_response_binning", TString::Format(";%s;", binByVarLabel.Data()), Binning::nbins_pt_zpj_reco_all, &Binning::pt_bin_edges_zpj_reco_all[0]);
  pt_genjet_response_binning = book<TH1F>("pt_genjet_response_binning", TString::Format(";%s;", binByVarLabel.Data()), Binning::nbins_pt_zpj_reco_all, &Binning::pt_bin_edges_zpj_reco_all[0]);

  eta_jet1_vs_pt = book<TH2F>(TString::Format("eta_jet1_vs_%s", binByVar.Data()), TString::Format(";y^{jet 1};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);

  flav_jet1_vs_pt_jet1 = book<TH2F>("flav_jet1_vs_pt_jet1", ";jet 1 flav;p_{T}^{jet1};", 23, -0.5, 22.5, nbins_pt, 0, pt_max);

  float ratio_lim = 8.;
  int nbins_ratio = ratio_lim*20;
  pt_jet1_z_ratio_vs_pt = book<TH2F>(TString::Format("pt_jet1_z_ratio_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 1} / p_{T}^{%s};%s", zName.Data(), binByVarLabel.Data()), nbins_ratio, 0, ratio_lim, nbins_pt, 0, pt_max);
  pt_jet1_z_ratio_vs_pt_jet1 = book<TH2F>("pt_jet1_z_ratio_vs_pt_jet1", TString::Format(";p_{T}^{jet 1} / p_{T}^{%s};p_{T}^{jet 1} [GeV]", zName.Data()), nbins_ratio, 0, ratio_lim, nbins_pt, 0, pt_max);

  float asym_lim = 1.; // maximum is ±1 by construction
  int nbins_asym = 2*asym_lim*20;
  jet1_z_asym_vs_pt = book<TH2F>(TString::Format("jet1_z_asym_vs_%s", binByVar.Data()), TString::Format(";(p_{T}^{jet 1} - p_{T}^{%s}) /(p_{T}^{jet 1} + p_{T}^{%s});%s", zName.Data(), zName.Data(), binByVarLabel.Data()), nbins_asym, -asym_lim, asym_lim, nbins_pt, 0, pt_max);
  jet1_z_asym_vs_pt_jet1 = book<TH2F>("jet1_z_asym_vs_pt_jet1", TString::Format(";(p_{T}^{jet 1} - p_{T}^{%s}) /(p_{T}^{jet 1} + p_{T}^{%s});p_{T}^{jet 1} [GeV]", zName.Data(), zName.Data()), nbins_asym, -asym_lim, asym_lim, nbins_pt, 0, pt_max);

  pt_jet_genHT_ratio = book<TH1F>("pt_jet_genHT_ratio", ";p_{T}^{jet 1}/GenHT", 250, 0, 2.5);
  pt_jet_response_fine = book<TH2F>("pt_jet_response_fine", ";p_{T}^{jet 1} (GEN);p_{T}^{jet 1} (RECO)", nbins_pt, 0, pt_max, nbins_pt, 0, pt_max);
  pt_jet_response = book<TH2F>("pt_jet_response", ";p_{T}^{jet 1} (GEN);p_{T}^{jet 1} (RECO)", Binning::nbins_pt_zpj_reco_all, &Binning::pt_bin_edges_zpj_reco_all[0], Binning::nbins_pt_zpj_reco_all, &Binning::pt_bin_edges_zpj_reco_all[0]);
  eta_jet_response = book<TH2F>("eta_jet_response", ";y^{jet} (GEN);y^{jet} (RECO)", nbins_eta, -eta_max, eta_max, nbins_eta, -eta_max, eta_max);

  gen_ht = book<TH1F>("gen_ht", ";H_{T}^{Gen} [GeV]", 500, 0, 5000);
  genjet_kt = book<TH1F>("genjet_kt", ";k_{T}^{Gen} [GeV]", 500, 0, 1000);
  // genjet_kt_vs_weight = book<TH2F>("genjet_kt_vs_weight", ";Genjet k_{T} [GeV];weight", 500, 0, 1000, );
  pt_jet2_vs_pt = book<TH2F>(TString::Format("pt_jet2_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 2} [GeV];%s", binByVarLabel.Data()), nbins_pt, 0, pt_max, nbins_pt, 0, pt_max);
  eta_jet2_vs_pt = book<TH2F>(TString::Format("eta_jet2_vs_%s", binByVar.Data()), TString::Format(";y^{jet 2};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  pt_jet2_z_ratio_vs_pt = book<TH2F>(TString::Format("pt_jet2_z_ratio_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 2} / p_{T}^{%s};%s", zName.Data(), binByVarLabel.Data()), 60, 0, 3, nbins_pt, 0, pt_max);
  pt_jet1_pt_jet2_ratio_vs_pt = book<TH2F>(TString::Format("pt_jet1_pt_jet2_ratio_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 1} / p_{T}^{jet 2};%s", binByVarLabel.Data()), 100, 0, 5, nbins_pt, 0, pt_max);

  DO_MATCHING_INDS = false;
  if (DO_MATCHING_INDS) {
    int nbins_jet_ind = 5;
    genjet1_ind_vs_pt_jet1 = book<TH2F>("genjet1_ind_vs_pt_jet1", ";GenJet index;p_{T}^{jet 1} [GeV]", nbins_jet_ind, 0-0.5, nbins_jet_ind-0.5, nbins_pt, 0, pt_max);
    for (int i=0; i < Binning::nbins_pt_zpj_reco_all; i++) {
      TH2F * tmp = book<TH2F>(TString::Format("genjet_ind_recojet_ind_pt_%g_%g", Binning::pt_bin_edges_zpj_reco_all.at(i), Binning::pt_bin_edges_zpj_reco_all.at(i+1)),
                              TString::Format("%g < p_{T}^{Gen} < %g GeV;GenJet index; RecoJet index", Binning::pt_bin_edges_zpj_reco_all.at(i), Binning::pt_bin_edges_zpj_reco_all.at(i+1)),
                              nbins_jet_ind+1, 0-1.5, nbins_jet_ind-0.5, nbins_jet_ind+1, 0-1.5, nbins_jet_ind-0.5);
      genjet_recojet_ind_binned.push_back(tmp);
    }
  }

  // met
  met_vs_pt = book<TH2F>(TString::Format("met_vs_%s", binByVar.Data()) , TString::Format(";p_{T}^{miss} [GeV];%s", binByVarLabel.Data()), 200, 0, 400, nbins_pt, 0, pt_max);
  int nbins_metSig(50);
  float metSig_max(10.);
  met_sig_vs_pt = book<TH2F>(TString::Format("met_sig_vs_%s", binByVar.Data()), TString::Format(";MET signif.;%s", binByVarLabel.Data()), nbins_metSig, 0, metSig_max, nbins_pt, 0, pt_max);

  // muons
  n_mu_vs_pt = book<TH2F>(TString::Format("n_mu_vs_%s", binByVar.Data()), TString::Format(";N_{#mu};%s", binByVarLabel.Data()), 10, 0, 10, nbins_pt, 0, pt_max);

  int nbins_reliso = 500;
  float reliso_max = 10;

  float mu_pt_max = 1000;

  pt_mu1_vs_pt = book<TH2F>(TString::Format("pt_mu1_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{#mu1};%s", binByVarLabel.Data()), nbins_pt, 0, mu_pt_max, nbins_pt, 0, pt_max);
  eta_mu1_vs_pt = book<TH2F>(TString::Format("eta_mu1_vs_%s", binByVar.Data()), TString::Format(";y^{#mu1};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  reliso_mu1_vs_pt = book<TH2F>(TString::Format("reliso_mu1_vs_%s", binByVar.Data()), TString::Format(";#mu1 rel. Iso;%s", binByVarLabel.Data()), nbins_reliso, 0, reliso_max, nbins_pt, 0, pt_max);

  pt_mu2_vs_pt = book<TH2F>(TString::Format("pt_mu2_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{#mu2};%s", binByVarLabel.Data()), nbins_pt, 0, mu_pt_max, nbins_pt, 0, pt_max);
  eta_mu2_vs_pt = book<TH2F>(TString::Format("eta_mu2_vs_%s", binByVar.Data()), TString::Format(";y^{#mu2};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt, 0, pt_max);
  reliso_mu2_vs_pt = book<TH2F>(TString::Format("reliso_mu2_vs_%s", binByVar.Data()), TString::Format(";#mu2 rel. Iso;%s", binByVarLabel.Data()), nbins_reliso, 0, reliso_max, nbins_pt, 0, pt_max);

  m_mumu_vs_pt = book<TH2F>(TString::Format("m_mumu_vs_%s", binByVar.Data()), TString::Format(";m_{%s} [GeV];%s", zName.Data(), binByVarLabel.Data()), 80, 90-40, 90+40, nbins_pt, 0, pt_max);
  pt_jet1_vs_pt = book<TH2F>(TString::Format("pt_jet1_vs_%s", binByVar.Data()), TString::Format(";p_{T}^{jet 1} [GeV];%s", binByVarLabel.Data()), 2*nbins_pt, 0, pt_max, 2*nbins_pt, 0, 2*pt_max);
  pt_mumu = book<TH1F>("pt_mumu", TString::Format(";p_{T}^{%s} [GeV];", zName.Data()), nbins_pt, 0, pt_max);


  deta_mumu_vs_pt = book<TH2F>(TString::Format("deta_mumu_vs_%s", binByVar.Data()), TString::Format(";|#Deltay(#mu, #mu)|;%s", binByVarLabel.Data()), nbins_eta, 0, 2*eta_max, nbins_pt, 0, pt_max);
  dphi_mumu_vs_pt = book<TH2F>(TString::Format("dphi_mumu_vs_%s", binByVar.Data()), TString::Format(";|#Delta#phi(#mu, #mu)|;%s",binByVarLabel.Data()), nbins_phi, 0, phi_max, nbins_pt, 0, pt_max);
  deta_mumu_jet1_vs_pt = book<TH2F>(TString::Format("deta_mumu_jet1_vs_%s", binByVar.Data()), TString::Format(";|#Deltay(%s, jet 1)|;%s", zName.Data(), binByVarLabel.Data()), nbins_eta, 0, 2*eta_max, nbins_pt, 0, pt_max);
  dphi_mumu_jet1_vs_pt = book<TH2F>(TString::Format("dphi_mumu_jet1_vs_%s", binByVar.Data()), TString::Format(";|#Delta#phi(%s, jet 1)|;%s", zName.Data(), binByVarLabel.Data()), nbins_phi, 0, phi_max, nbins_pt, 0, pt_max);

  // primary vertices
  n_pv = book<TH1F>("N_pv", ";N_{PV};", 50, 0, 50);

}


void QGAnalysisZPlusJetsHists::fill(const Event & event){
  std::vector<Jet> * jets = event.jets;
  int Njets = jets->size();

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

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
  eta_mu1_vs_pt->Fill(mu1.Rapidity(), binPt, weight);
  reliso_mu1_vs_pt->Fill(mu1.relIso(), binPt, weight);

  pt_mu2_vs_pt->Fill(mu2.pt(), binPt, weight);
  eta_mu2_vs_pt->Fill(mu2.Rapidity(), binPt, weight);
  reliso_mu2_vs_pt->Fill(mu2.relIso(), binPt, weight);

  m_mumu_vs_pt->Fill(z_cand.M(), binPt, weight);
  pt_jet1_vs_pt->Fill(jet1_pt, binPt, weight);

  deta_mumu_vs_pt->Fill(fabs(mu1.Rapidity() - mu2.Rapidity()), binPt,  weight);
  dphi_mumu_vs_pt->Fill(fabs(deltaPhi(mu1.v4(), mu2.v4())), binPt, weight);

  // Jets
  n_jets_vs_pt->Fill(Njets, binPt, weight);

  if (Njets < 1) return;

  Jet & jet1 = jets->at(0);
  pt_jet1->Fill(jet1_pt, weight);
  pt_jet1_unweighted->Fill(jet1_pt);
  pt_jet_response_binning->Fill(jet1_pt, weight);
  eta_jet1_vs_pt->Fill(jet1.Rapidity(), binPt, weight);
  flav_jet1_vs_pt_jet1->Fill(abs(jet1.flavor()), jet1_pt, weight);

  if (!event.isRealData) {
    float genHT = calcGenHT(*(event.genparticles));
    gen_ht->Fill(genHT, weight);
    pt_jet_genHT_ratio->Fill(jet1_pt / genHT, weight);

    float jetKt = calcJetKt(*(event.genparticles));
    genjet_kt->Fill(jetKt, weight);

    const std::vector<GenJetWithParts> * genjets = &event.get(genJets_handle);
    int gj_ind1 = jet1.genjet_index();

    // if (gj_ind1 < 0) throw std::runtime_error("gj_ind1 is < 0");
    if (gj_ind1 >= int(genjets->size())) throw std::runtime_error("gj_ind1 is larger than genjet collection");


    if (gj_ind1 >=0 ) {
      auto genjet1 = genjets->at(gj_ind1);

      pt_jet_response->Fill(genjet1.pt(), jet1_pt, weight);
      pt_jet_response_fine->Fill(genjet1.pt(), jet1_pt, weight);
      eta_jet_response->Fill(genjet1.Rapidity(), jet1.Rapidity(), weight);
      pt_genjet_response_binning->Fill(genjet1.pt(), weight);
    }

    if (DO_MATCHING_INDS) {
      genjet1_ind_vs_pt_jet1->Fill(gj_ind1, jet1_pt);

      // plot indices of matching genjets/recojets
      int genjet_ind = 0;
      std::vector<int> matches;
      for (const auto & gjItr : *genjets) {
        if (genjet_ind > 5) break;

        // find matching recojet for this genjet
        float dr_min = 9999;
        int reco_ind = -1;
        for (uint rj_ind=0; rj_ind < jets->size(); rj_ind++) {
          // skip if already matched reco jet
          if (std::find(matches.begin(), matches.end(), rj_ind) != matches.end()) continue;

          float this_dr = deltaR(jets->at(rj_ind), gjItr);
          if (this_dr < dr_min && this_dr < (jetRadius/2.)) {
            reco_ind = rj_ind;
          }
        }
        if (reco_ind >= 0) matches.push_back(reco_ind);

        //  find which bin is suitable for this genjet.
        auto pos = std::lower_bound(Binning::pt_bin_edges_zpj_reco_all.begin(), Binning::pt_bin_edges_zpj_reco_all.end(), gjItr.pt()) - Binning::pt_bin_edges_zpj_reco_all.begin();
        if (pos > (int) genjet_recojet_ind_binned.size()+1) {
          cout << "pos = " << pos << endl;
          cout << gjItr.pt() << endl;
          for (auto b : Binning::pt_bin_edges_zpj_reco_all) {
            cout << b << " ";
          }
          cout << endl;
          throw std::runtime_error("Invalid pos");
        }
        if (pos == 0) {
          cout << "pos = 0" << endl;
          cout << gjItr.pt() << endl;
          for (auto b : Binning::pt_bin_edges_zpj_reco_all) {
            cout << b << " ";
          }
          cout << endl;
          throw std::runtime_error("Invalid pos 0");
        }

        // need -1 as it will find the upper edge
        genjet_recojet_ind_binned.at(pos-1)->Fill(genjet_ind, reco_ind, weight);

        genjet_ind++;
      }

      // now go through the recojets and see if there are any missing matches
      int recojet_ind = 0;
      for (const auto & rjItr : *event.jets) {
        if (rjItr.genjet_index() < 0) {
          auto match = std::find(matches.begin(), matches.end(), recojet_ind);
          if (match != matches.end()) {
            cout << "Found match for recojet " << recojet_ind << " with genjet " << match - matches.begin() << endl;
            printJets(*event.jets, "BadMatching");
            printGenJets(*genjets, "BadMatching");
            // throw std::runtime_error("recojet.genjet_index indicates no match, but genjet match found");
          }

          //  find which bin is suitable
          //  not really correct as we're using reco pt, but there's no genjet anyway
          auto pos = std::lower_bound(Binning::pt_bin_edges_zpj_reco_all.begin(), Binning::pt_bin_edges_zpj_reco_all.end(), rjItr.pt()) - Binning::pt_bin_edges_zpj_reco_all.begin();
          if (pos > (int) genjet_recojet_ind_binned.size()+1) {
            cout << pos << endl;
            cout << rjItr.pt() << endl;
            for (auto b : Binning::pt_bin_edges_zpj_reco_all) {
              cout << b << " ";
            }
            cout << endl;
            throw std::runtime_error("Invalid pos");
          }
          // need -1 as it will find the upper edge
          genjet_recojet_ind_binned.at(pos-1)->Fill(-1, recojet_ind, weight);
        }
        recojet_ind++;
      }
    }
  }

  met_vs_pt->Fill(event.met->pt(), binPt, weight);
  met_sig_vs_pt->Fill(event.met->mEtSig(), binPt, weight);

  pt_jet1_z_ratio_vs_pt->Fill(jet1_pt / z_pt, binPt, weight);
  pt_jet1_z_ratio_vs_pt_jet1->Fill(jet1_pt / z_pt, jet1_pt, weight);

  jet1_z_asym_vs_pt->Fill((jet1_pt - z_pt) / (jet1_pt + z_pt), binPt, weight);
  jet1_z_asym_vs_pt_jet1->Fill((jet1_pt - z_pt) / (jet1_pt + z_pt), jet1_pt, weight);

  deta_mumu_jet1_vs_pt->Fill(fabs(z_cand.Rapidity() - jet1.Rapidity()), binPt, weight);
  dphi_mumu_jet1_vs_pt->Fill(fabs(deltaPhi(jet1, z_cand)), binPt, weight);


  if (Njets >= 2) {
    Jet jet2 = jets->at(1);
    pt_jet2_vs_pt->Fill(jet2.pt(), binPt, weight);
    eta_jet2_vs_pt->Fill(jet2.Rapidity(), binPt, weight);
    pt_jet2_z_ratio_vs_pt->Fill(jet2.pt() / z_pt, binPt, weight);
    pt_jet1_pt_jet2_ratio_vs_pt->Fill(jet1_pt / jet2.pt(), binPt, weight);
  }

  int Npvs = event.pvs->size();
  n_pv->Fill(Npvs, weight);
}

QGAnalysisZPlusJetsHists::~QGAnalysisZPlusJetsHists(){}

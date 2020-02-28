#include "UHH2/QGAnalysis/include/QGAnalysisDijetHists.h"
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

QGAnalysisDijetHists::QGAnalysisDijetHists(Context & ctx, const string & dirname, const std::string & binning_):
  Hists(ctx, dirname),
  binning(binning_)
{
  string jet_cone = ctx.get("JetCone", "AK4");
  jetRadius = get_jet_radius(jet_cone);

  if (binning != "ave" && binning != "tnp" && binning != "leading") {
    throw std::runtime_error("Binning should be 'ave', 'tnp', or 'leading'");
  }
  is_mc_ = ctx.get("dataset_type") == "MC";
  if (is_mc_) genJets_handle = ctx.get_handle< std::vector<GenJetWithParts> > ("GoodGenJets");

  // book all histograms here
  // jets
  int nbins_pt_equal = 2000;
  float pt_max = 2000;

  int nbins_eta = 200;
  float eta_max = 5;

  int nbins_phi = 128;
  float phi_max = 3.2;

  TString binByVarLabel = "p_{T}^{jet 1} [GeV]";
  if (binning == "ave") {
    binByVarLabel = "#LT p_{T}^{jet 1, 2} #GT [GeV]";
  } else if (binning == "tnp") {
    binByVarLabel = "p_{T}^{jet} [GeV]";
  }

  pt_bin_edges = Binning::pt_bin_edges_reco_all;
  nbins_pt = pt_bin_edges.size() - 1;

  n_jets_vs_pt_jet = book<TH2F>("n_jets_vs_pt_jet", TString::Format(";N_{jets};%s", binByVarLabel.Data()), 10, 0, 10, nbins_pt_equal, 0, pt_max);

  pt_jet = book<TH1F>("pt_jet", TString::Format(";%s;", binByVarLabel.Data()), nbins_pt_equal, 0, pt_max);
  pt_jet1 = book<TH1F>("pt_jet1", ";p_{T}^{jet 1} [GeV];", nbins_pt_equal, 0, pt_max);
  pt_jet2 = book<TH1F>("pt_jet2", ";p_{T}^{jet 2} [GeV];", nbins_pt_equal, 0, pt_max);
  pt_jet_unweighted = book<TH1F>("pt_jet_unweighted", TString::Format(";%s;", binByVarLabel.Data()), nbins_pt_equal, 0, pt_max);
  pt_jet1_unweighted = book<TH1F>("pt_jet1_unweighted", ";p_{T}^{jet 1} [GeV];", nbins_pt_equal, 0, pt_max);
  pt_jet_response_binning = book<TH1F>("pt_jet_response_binning", TString::Format(";%s;", binByVarLabel.Data()), nbins_pt, &pt_bin_edges[0]);
  pt_genjet_response_binning = book<TH1F>("pt_genjet_response_binning", TString::Format(";%s;", binByVarLabel.Data()), nbins_pt, &pt_bin_edges[0]);

  pt_jet_qScale_ratio = book<TH1F>("pt_jet_qScale_ratio", ";p_{T}^{jet 1}/qScale", 300, 0, 15);
  pt_jet_genHT_ratio = book<TH1F>("pt_jet_genHT_ratio", ";p_{T}^{jet 1}/GenHT", 250, 0, 2.5);
  pt_jet_vs_pdf_scalePDF = book<TH2F>("pt_jet_vs_pdf_scalePDF", ";p_{T}^{jet 1};pdf_scalePDF", nbins_pt_equal, 0, pt_max, nbins_pt_equal, 0, pt_max);
  pt_jet_vs_genHT = book<TH2F>("pt_jet_vs_genHT", ";p_{T}^{jet 1};GenHT", nbins_pt_equal, 0, pt_max, 500, 0, 5000);

  int nWeightBins = 25;
  double weightBins [nWeightBins+1] = {1E-4, 1E-3, 1E-2, 1E-1, 1, 10, 100, 1000, 1E4, 2E4, 5E4, 7E4, 1E5, 2E5, 5E5, 7E5, 1E6, 2E6, 3E6, 4E6, 5E6, 6E6, 7E6, 8E6, 9E6, 1E7};
  weight_vs_puHat_genHT_ratio = book<TH2F>("weight_vs_puHat_genHT_ratio", ";Weight;PU #hat{p}_{T} / Gen HT", nWeightBins, weightBins, 100, 0, 2);

  // dont' reverse axis direction - it doens't like it
  pt_jet_response_fine = book<TH2F>("pt_jet_response_fine", TString::Format(";%s (GEN);%s (RECO)", binByVarLabel.Data(), binByVarLabel.Data()), nbins_pt_equal, 0, pt_max, nbins_pt_equal, 0, pt_max);
  pt_jet_response = book<TH2F>("pt_jet_response", TString::Format(";%s (GEN);%s (RECO)", binByVarLabel.Data(), binByVarLabel.Data()), nbins_pt, &pt_bin_edges[0], nbins_pt, &pt_bin_edges[0]);
  eta_jet1_vs_pt_jet = book<TH2F>("eta_jet1_vs_pt_jet", TString::Format(";#eta^{jet 1};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt_equal, 0, pt_max);
  eta_jet_response = book<TH2F>("eta_jet_response", ";#eta^{jet} (GEN);#eta^{jet} (RECO)", nbins_eta, -eta_max, eta_max, nbins_eta, -eta_max, eta_max);
  phi_jet1_vs_pt_jet = book<TH2F>("phi_jet1_vs_pt_jet", TString::Format(";#phi^{jet 1};%s", binByVarLabel.Data()), nbins_phi, -phi_max, phi_max, nbins_pt_equal, 0, pt_max);

  eta_jet1_vs_eta_jet2 = book<TH2F>("eta_jet1_vs_eta_jet2", ";#eta^{jet 1};#eta^{jet 2}", nbins_eta, -eta_max, eta_max, nbins_eta,-eta_max, eta_max);
  pt_jet2_vs_pt_jet = book<TH2F>("pt_jet2_vs_pt_jet", TString::Format(";p_{T}^{jet 2};%s", binByVarLabel.Data()), nbins_pt_equal, 0, pt_max, nbins_pt_equal, 0, pt_max);
  eta_jet2_vs_pt_jet = book<TH2F>("eta_jet2_vs_pt_jet", TString::Format(";#eta^{jet 2};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt_equal, 0, pt_max);
  phi_jet2_vs_pt_jet = book<TH2F>("phi_jet2_vs_pt_jet", TString::Format(";#phi^{jet 2};%s", binByVarLabel.Data()), nbins_phi, -phi_max, phi_max, nbins_pt_equal, 0, pt_max);

  DO_MATCHING_INDS = false;
  if (DO_MATCHING_INDS) {
    int nbins_jet_ind = 5;
    // -1 on each axis = no match
    genjet1_ind_vs_pt_jet1 = book<TH2F>("genjet1_ind_vs_pt_jet1", ";GenJet index;p_{T}^{jet 1} [GeV]", nbins_jet_ind+1, 0-1.5, nbins_jet_ind-0.5, nbins_pt_equal, 0, pt_max);
    genjet2_ind_vs_pt_jet2 = book<TH2F>("genjet2_ind_vs_pt_jet2", ";GenJet index;p_{T}^{jet 2} [GeV]", nbins_jet_ind+1, 0-1.5, nbins_jet_ind-0.5, nbins_pt_equal, 0, pt_max);

    for (int i=0; i < nbins_pt; i++) {
      TH2F * tmp = book<TH2F>(TString::Format("genjet_ind_recojet_ind_pt_%g_%g", pt_bin_edges.at(i), pt_bin_edges.at(i+1)),
                              TString::Format("%g < p_{T}^{Gen} < %g GeV;GenJet index; RecoJet index", pt_bin_edges.at(i), pt_bin_edges.at(i+1)),
                              nbins_jet_ind+1, 0-1.5, nbins_jet_ind-0.5, nbins_jet_ind+1, 0-1.5, nbins_jet_ind-0.5);
      genjet_recojet_ind_binned.push_back(tmp);
    }
  }

  flav_jet1_jet2 = book<TH2F>("flav_jet1_jet2", ";jet 1 flav;jet 2 flav;", 23, -0.5, 22.5, 23, -0.5, 22.5);
  flav_jet1_vs_pt_jet = book<TH2F>("flav_jet1_vs_pt_jet", ";jet 1 flav;p_{T}^{jet1};", 23, -0.5, 22.5, nbins_pt_equal, 0, pt_max);
  flav_jet2_vs_pt_jet = book<TH2F>("flav_jet2_vs_pt_jet", ";jet 2 flav;p_{T}^{jet2};", 23, -0.5, 22.5, nbins_pt_equal, 0, pt_max);

  pt_jet1_jet2_ratio_vs_pt_jet = book<TH2F>("pt_jet1_jet2_ratio_vs_pt_jet", TString::Format(";p_{T}^{jet 2} / p_{T}^{jet 1};%s", binByVarLabel.Data()), 50, 0, 1, nbins_pt_equal, 0, pt_max);
  jet1_jet2_asym_vs_pt_jet = book<TH2F>("jet1_jet2_asym_vs_pt_jet", TString::Format(";|p_{T}^{jet 1} - p_{T}^{jet 2}|/p_{T}^{jet 1} + p_{T}^{jet 2};%s", binByVarLabel.Data()), 50, 0, 1, nbins_pt_equal, 0, pt_max);
  // for studying cuts
  jet1_jet2_asym_vs_pt_jet1 = book<TH2F>("jet1_jet2_asym_vs_pt_jet1", ";|p_{T}^{jet 1} - p_{T}^{jet 2}|/p_{T}^{jet 1} + p_{T}^{jet 2};p_{T}^{jet 1} [GeV]", 50, 0, 1, nbins_pt_equal, 0, pt_max);
  jet1_jet2_asym_vs_pt_jet2 = book<TH2F>("jet1_jet2_asym_vs_pt_jet2", ";|p_{T}^{jet 1} - p_{T}^{jet 2}|/p_{T}^{jet 1} + p_{T}^{jet 2};p_{T}^{jet 2} [GeV]", 50, 0, 1, nbins_pt_equal, 0, pt_max);

  m_jj_vs_pt_jet = book<TH2F>("m_jj_vs_pt_jet", TString::Format(";m_{jj} [GeV];%s", binByVarLabel.Data()), 200, 0, 4000, nbins_pt_equal, 0, pt_max);

  deta_jj_vs_pt_jet = book<TH2F>("deta_jj_vs_pt_jet", TString::Format(";|#Delta#eta(jet1, jet2)|;%s", binByVarLabel.Data()), nbins_eta, 0, 2*eta_max, nbins_pt_equal, 0, pt_max);
  dphi_jj_vs_pt_jet = book<TH2F>("dphi_jj_vs_pt_jet", TString::Format(";|#Delta#phi(jet1, jet2)|;%s", binByVarLabel.Data()), nbins_phi, 0, phi_max, nbins_pt_equal, 0, pt_max);
  sumeta_jj_vs_pt_jet = book<TH2F>("sumeta_jj_vs_pt_jet", TString::Format(";#sum#eta(jet1, jet2)};%s", binByVarLabel.Data()), 2*nbins_eta, -2*eta_max, 2*eta_max, nbins_pt_equal, 0, pt_max);

  // Possible 3rd jet in the event
  pt_jet3_vs_pt_jet = book<TH2F>("pt_jet3_vs_pt_jet", TString::Format(";p_{T}^{jet 3};%s", binByVarLabel.Data()), nbins_pt_equal, 0, 500, nbins_pt_equal, 0, pt_max);
  eta_jet3_vs_pt_jet = book<TH2F>("eta_jet3_vs_pt_jet", TString::Format(";#eta^{jet 3};%s", binByVarLabel.Data()), nbins_eta, -eta_max, eta_max, nbins_pt_equal, 0, pt_max);
  pt_jet3_frac_vs_pt_jet = book<TH2F>("pt_jet3_frac_vs_pt_jet", TString::Format(";p_{T}^{jet 3} / #LT p_{T}^{jet 1}, p_{T}^{jet 2} #GT;%s", binByVarLabel.Data()), 50, 0, 1, nbins_pt_equal, 0, pt_max);

  // MET
  met_vs_pt_jet = book<TH2F>("met_vs_pt_jet", TString::Format(";p_{T}^{miss} [GeV];%s", binByVarLabel.Data()), 200, 0, 400, nbins_pt_equal, 0, pt_max);
  int nbins_metSig(50);
  float metSig_max(10.);
  met_sig_vs_pt_jet = book<TH2F>("met_sig_vs_pt_jet", TString::Format(";MET signif.;%s", binByVarLabel.Data()), nbins_metSig, 0, metSig_max, nbins_pt_equal, 0, pt_max);

  // primary vertices
  n_pv = book<TH1F>("N_pv", ";N^{PV};", 50, 0, 50);
}


void QGAnalysisDijetHists::fill(const Event & event){
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  if (Njets < 2) return;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  Jet jet1 = jets->at(0);
  double jet1_pt = jet1.pt();
  pt_jet1->Fill(jet1_pt, weight);
  pt_jet1_unweighted->Fill(jet1_pt);

  Jet jet2 = jets->at(1);
  double jet2_pt = jet2.pt();
  pt_jet2->Fill(jet2_pt, weight);

  double binByVal = 0.;
  double largest_jet_pt = max(jet1_pt, jet2_pt);
  if (binning == "ave") {
    binByVal = (jet1_pt + jet2_pt)/2.;
  } else if (binning == "leading") {
    binByVal = largest_jet_pt;
  }


  // cout << "weight: " << weight << endl;
  // cout << "qScale: " << event.genInfo->qScale() << endl;
  // cout << "pdf_x1: " << event.genInfo->pdf_x1() << endl;
  // cout << "pdf_x2: " << event.genInfo->pdf_x2() << endl;
  // cout << "pdf_scalePDF: " << event.genInfo->pdf_scalePDF() << endl;
  // if (event.genInfo->binningValues().size()>0) cout << "ptHAT: " << event.genInfo->binningValues()[0] << endl;
  // if (binByVal > 100) {
  //   cout << "EVENT: " << event.event << " : " << event.luminosityBlock << endl;
  //   cout << "LARGE PT: " << binByVal << endl;
  //   printJets(*event.jets, "Silly Dijet jets");
  //   printGenJetsWithParts(event.get(genJets_handle), event.genparticles);
  //   printGenParticles(*event.genparticles);
  // }

  if (is_mc_) {
    float genHT = calcGenHT(*event.genparticles);
    pt_jet_qScale_ratio->Fill(largest_jet_pt / event.genInfo->pdf_scalePDF(), weight);
    pt_jet_genHT_ratio->Fill(largest_jet_pt / genHT, weight);
    pt_jet_vs_pdf_scalePDF->Fill(largest_jet_pt, event.genInfo->pdf_scalePDF(), weight);
    pt_jet_vs_genHT->Fill(largest_jet_pt, genHT, weight);
    weight_vs_puHat_genHT_ratio->Fill(weight, event.genInfo->PU_pT_hat_max() / genHT);

    const std::vector<GenJetWithParts> * genjets = &event.get(genJets_handle);
    int gj_ind1 = jet1.genjet_index();
    int gj_ind2 = jet2.genjet_index();

    if (gj_ind1 >= int(genjets->size())) throw std::runtime_error("gj_ind1 is larger than genjet collection");
    if (gj_ind2 >= int(genjets->size())) throw std::runtime_error("gj_ind2 is larger than genjet collection");

    double genVal = -1.;
    if (gj_ind1>=0) {
      auto genjet1 = genjets->at(gj_ind1);
      eta_jet_response->Fill(genjet1.eta(), jet1.eta(), weight);
    }
    if (gj_ind2 >=0) {
      auto genjet2 = genjets->at(gj_ind2);
      eta_jet_response->Fill(genjet2.eta(), jet2.eta(), weight);
    }
    if (gj_ind1>=0 && gj_ind2 >=0) {
      auto genjet1 = genjets->at(gj_ind1);
      auto genjet2 = genjets->at(gj_ind2);
      genVal = (genjet1.pt() + genjet2.pt())/2.;
      pt_genjet_response_binning->Fill(genVal, weight);
    }

    pt_jet_response_fine->Fill(genVal, binByVal, weight);
    pt_jet_response->Fill(genVal, binByVal, weight);

    // cout << " -- Doing matching" << endl;
    // printJets(*jets, "reco jets for matching");
    // printGenJets(*genjets, "gen jets for matching");

    if (DO_MATCHING_INDS) {
      // plot indices of matching genjets/recojets
      genjet1_ind_vs_pt_jet1->Fill(gj_ind1, jet1_pt, weight);
      genjet2_ind_vs_pt_jet2->Fill(gj_ind2, jet2_pt, weight);

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
        if (reco_ind >= 0) {
          matches.push_back(reco_ind);
          // cout << "GenJet-reco match: " << endl;
          // cout << "  Gen: " << genjet_ind << " : " << gjItr.pt() << " : " << gjItr.eta() << " : " << gjItr.phi() << endl;
          // cout << "  Reco: " << reco_ind << " : " << jets->at(reco_ind).pt() << " : " << jets->at(reco_ind).eta() << " : " << jets->at(reco_ind).phi() << endl;
        }

        //  find which pt bin is suitable for this genjet
        auto pos = std::lower_bound(pt_bin_edges.begin(), pt_bin_edges.end(), gjItr.pt()) - pt_bin_edges.begin();
        if (pos > (int) genjet_recojet_ind_binned.size()+1) {
          cout << "pos: " << pos << endl;
          cout << "genjet pt: " << gjItr.pt() << endl;
          cout << "bins: " << endl;
          for (auto b : pt_bin_edges) {
            cout << b << " ";
          }
          cout << endl;
          throw std::runtime_error("Invalid pos");
        }
        if (pos == 0) {
          for (auto & ptv : pt_bin_edges) { cout << ptv << ", " ;}
          cout << endl;
          cout << "pos = 0, genJet pt: " << gjItr.pt() << endl;
          throw std::runtime_error("pos=0, this will fail to find a gen-reco bin; decrease the first bin edge");
        }
        if (pos == (int) genjet_recojet_ind_binned.size()+1) {
          for (auto & ptv : pt_bin_edges) { cout << ptv << ", " ;}
          cout << endl;
          cout << "pos = " << genjet_recojet_ind_binned.size() << ", genJet pt: " << gjItr.pt() << endl;
          throw std::runtime_error("pos=genjet_recojet_ind_binned.size(), this will fail to find a gen-reco bin; increase the last bin edge");
        }
        // need -1 as it will find the upper edge
        genjet_recojet_ind_binned.at(pos-1)->Fill(genjet_ind, reco_ind, weight);

        genjet_ind++;
      }


      // now go through the recojets and see if there are any missing matches
      int recojet_ind = 0;
      for (const auto & rjItr : *event.jets) {
        auto match = std::find(matches.begin(), matches.end(), recojet_ind);

        if (rjItr.genjet_index() < 0) {
          if (match != matches.end()) {
            cout << "Found match for recojet " << recojet_ind << " with genjet " << match - matches.begin() << " but there wasn't a match earlier" << endl;
            printJets(*event.jets, "BadMatching");
            printGenJets(*genjets, "BadMatching");
            // throw std::runtime_error("recojet.genjet_index indicates no match, but genjet match found");
          } else {
            // here we have reco jet with no matching genjet
            // find which bin is suitable
            // not really correct as we're using reco pt, but there's no genjet anyway
            auto pos = std::lower_bound(pt_bin_edges.begin(), pt_bin_edges.end(), rjItr.pt()) - pt_bin_edges.begin();
            if (pos > (int) genjet_recojet_ind_binned.size()+1) {
              cout << pos << endl;
              cout << rjItr.pt() << endl;
              for (auto b : pt_bin_edges) {
                cout << b << " ";
              }
              cout << endl;
              throw std::runtime_error("Invalid pos");
            }
            if (pos == 0) {
              for (auto & ptv : pt_bin_edges) { cout << ptv << ", " ;}
              cout << endl;
              cout << "pos = 0, recojet pt: " << rjItr.pt() << endl;
              throw std::runtime_error("pos=0, this will fail to find a gen-reco bin; decrease the first bin edge");
            }
            // need -1 as it will find the upper edge
            genjet_recojet_ind_binned.at(pos-1)->Fill(-1, recojet_ind, weight);
          }
        }

        // here we have a recojet with a matching genjet, but it wasn't used earlier?
        // if (match == matches.end()) {
        //   // what are we doing here?
        //   // find which bin is suitable
        //   // not really correct as we're using reco pt, but there's no genjet anyway
        //   auto pos = std::lower_bound(pt_bin_edges.begin(), pt_bin_edges.end(), rjItr.pt()) - pt_bin_edges.begin();
        //   if (pos > (int) genjet_recojet_ind_binned.size()+1) {
        //     cout << pos << endl;
        //     cout << rjItr.pt() << endl;
        //     for (auto b : pt_bin_edges) {
        //       cout << b << " ";
        //     }
        //     cout << endl;
        //     throw std::runtime_error("Invalid pos");
        //   }
        //   if (pos == 0) {
        //     for (auto & ptv : pt_bin_edges) { cout << ptv << ", " ;}
        //     cout << endl;
        //     cout << "pos = 0, recojet pt: " << rjItr.pt() << endl;
        //     throw std::runtime_error("pos=0, this will fail to find a gen-reco bin; decrease the first bin edge");
        //   }
        //   // need -1 as it will find the upper edge
        //   genjet_recojet_ind_binned.at(pos-1)->Fill(-1, recojet_ind, weight);
        // }
        recojet_ind++;
      }
    }
  }

  n_jets_vs_pt_jet->Fill(Njets, binByVal, weight);
  pt_jet->Fill(binByVal, weight);
  pt_jet_unweighted->Fill(binByVal);
  pt_jet_response_binning->Fill(binByVal, weight);
  eta_jet1_vs_pt_jet->Fill(jet1.eta(), binByVal, weight);
  phi_jet1_vs_pt_jet->Fill(jet1.phi(), binByVal, weight);

  eta_jet1_vs_eta_jet2->Fill(jet1.eta(), jet2.eta(), weight);
  pt_jet2_vs_pt_jet->Fill(jet2_pt, binByVal, weight);
  eta_jet2_vs_pt_jet->Fill(jet2.eta(), binByVal, weight);
  phi_jet2_vs_pt_jet->Fill(jet2.phi(), binByVal, weight);

  // ensure it's always smaller / bigger
  float ratio = (jet1_pt > jet2_pt) ? jet2_pt / jet1_pt : jet1_pt / jet2_pt;
  pt_jet1_jet2_ratio_vs_pt_jet->Fill(ratio, binByVal, weight);
  float jetAsym = fabs(jet1_pt - jet2_pt) / (jet1_pt + jet2_pt);
  jet1_jet2_asym_vs_pt_jet->Fill(jetAsym, binByVal, weight);
  jet1_jet2_asym_vs_pt_jet1->Fill(jetAsym, jet1_pt, weight);
  jet1_jet2_asym_vs_pt_jet2->Fill(jetAsym, jet2_pt, weight);

  flav_jet1_jet2->Fill(abs(jet1.flavor()), abs(jet2.flavor()), weight);
  flav_jet1_vs_pt_jet->Fill(abs(jet1.flavor()), weight);
  flav_jet2_vs_pt_jet->Fill(abs(jet2.flavor()), weight);

  double mass = (jet1.v4() + jet2.v4()).M();
  m_jj_vs_pt_jet->Fill(mass, binByVal, weight);

  double dEta = fabs(jet1.eta() - jet2.eta());
  double dPhi = fabs(deltaPhi(jet1, jet2));
  deta_jj_vs_pt_jet->Fill(dEta, binByVal, weight);
  dphi_jj_vs_pt_jet->Fill(dPhi, binByVal, weight);
  sumeta_jj_vs_pt_jet->Fill(jet1.eta() + jet2.eta(), binByVal, weight);

  met_vs_pt_jet->Fill(event.met->pt(), binByVal, weight);
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

#include <iostream>
#include <vector>
#include <iomanip>      // std::setprecision

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TTimeStamp.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

using namespace std;

/**
 * Run with
 * root -q -b herwigpp_test_spectra.C
 *
 * DO NOT compile, i.e. use herwig_test_spectra.C+
 * Otherwise you will get a mountain of errors that I can't fix
 */

struct HistPack
{
  TH1D * h;
  TH1D * h_weighted;
  TH2D * h_kt_z_ratio_vs_kt;
  TH2D * h_kt_z_ratio_vs_kt_weighted;
  TH2D * h_kt_z_ratio_vs_weight;
  TH2D * h_weight_vs_kt;
  TH2D * h_jet3_jet1_ratio_vs_weight;
  double sum;
  double xsec;
};

HistPack fill_histograms(TChain * ch, const string & append, double xsec, long long limit) {
  HistPack hp;
  int nBins = 100;
  float minX = 0;
  float maxX = 1000;
  hp.h = new TH1D(("jetkt" + append).c_str(), ";Gen Parton k_{T} [GeV];N (raw)", nBins, minX, maxX);
  hp.h_weighted = new TH1D(("jetkt_weighted" + append).c_str(), ";Gen Parton k_{T} [GeV];N (weighted)", nBins, minX, maxX);
  
  int nBinsRatio = 100;
  float ratioMin = 0;
  float ratioMax = 5;
  // vector<double> weightBins = {
  //   1E-2, 2.5E-2, 5E-2, 7.5E-2, 1E-1, 2.5E-1, 5E-1, 7.5E-1, 1, 2.5, 5, 7.5, 10, 25, 50, 75, 100, 250, 500, 750, 1000, 2500, 7500, 10000
  // };

  vector<double> weightBins;
  double base = 10;
  for (double i=-2; i<=4.1; i+=0.1) {
    weightBins.push_back(pow(base, i));
    // cout << weightBins.back() << endl;
  }

  hp.h_kt_z_ratio_vs_kt = new TH2D(("kt_z_ratio_vs_kt" + append).c_str(), ";Gen Parton k_{T} [GeV];p_{T}^{Z} / Gen Parton k_{T}", nBins, minX, maxX, nBinsRatio, ratioMin, ratioMax);
  hp.h_kt_z_ratio_vs_kt_weighted = new TH2D(("kt_z_ratio_vs_kt_weighted" + append).c_str(), ";Gen Parton k_{T} [GeV];p_{T}^{Z} / Gen Parton k_{T}", nBins, minX, maxX, nBinsRatio, ratioMin, ratioMax);
  hp.h_kt_z_ratio_vs_weight = new TH2D(("kt_z_ratio_vs_weight" + append).c_str(), ";weight;p_{T}^{Z} / Gen Parton k_{T}", (int)(weightBins.size() - 1), &weightBins[0], nBinsRatio, ratioMin, ratioMax);
  hp.h_jet3_jet1_ratio_vs_weight = new TH2D(("jet3_jet1_ratio_vs_weight" + append).c_str(), ";weight;p_{T}^{genjet3} / p_{T}^{genjet1}", (int)(weightBins.size() - 1), &weightBins[0], nBinsRatio, ratioMin, 1);
  
  hp.h_weight_vs_kt = new TH2D(("weight_vs_kt" + append).c_str(), ";Gen Parton k_{T} [GeV];Weight", nBins, minX, maxX, (int)(weightBins.size() - 1), &weightBins[0]);

  ch->Print();
  TTreeReader reader(ch);
  // yes you do need this weird branch name, otherise will say
  // Error in <TTreeReaderValueBase::CreateProxy()>: The tree does not have a branch called genParticles.
  // Or setup an Alias

  // For GEN only files:
  // TTreeReaderValue<GenEventInfoProduct> genInfo(reader, "GenEventInfoProduct_generator__GEN.obj");
  // TTreeReaderValue<vector<reco::GenParticle>> genParticles(reader, "recoGenParticles_genParticles__GEN.obj");
  // TTreeReaderValue<vector<reco::GenJet>> genJets(reader, "recoGenJets_ak4GenJets__GEN.obj");

  // for GEN-SIM files"
  TTreeReaderValue<GenEventInfoProduct> genInfo(reader, "GenEventInfoProduct_generator__SIM.obj");
  TTreeReaderValue<vector<reco::GenParticle>> genParticles(reader, "recoGenParticles_genParticles__SIM.obj");
  TTreeReaderValue<vector<reco::GenJet>> genJets(reader, "recoGenJets_ak4GenJets__SIM.obj");

  double sum = 0.; // MUST be double otherwise loses precision
  int counter = 0;
  // int limit = 50000;
  limit = min(limit, ch->GetEntries());
  // int mark = limit/100;
  limit = 30;
  int mark = 1;
  cout << "===== Doing sample " << append << " (" << limit << " events)" << endl;
  TTimeStamp startTime;
  while (reader.Next() && counter < limit) {
    if (counter % mark == 0) {
      // predict finish time
      TTimeStamp now;
      double diff = now - startTime; // in seconds
      double rate = diff / counter;
      double timeLeft = rate * (limit - counter);
      string units = " seconds";
      if (timeLeft < 0) cout << "Warning overflow in timeLeft..." << endl;
      if (timeLeft > 60) {
        timeLeft /= 60;
        units = " minutes";
      }
      if (timeLeft > 60) {
        timeLeft /= 60;
        units = " hours";
      }
      if (timeLeft > 24) {
        timeLeft /= 24;
        units = " days";
      }
      cout << "..." << counter << " ETA: " << timeLeft << units << endl;
    }
    counter++;
    // if (counter > 48600) continue;

    float weight = genInfo->weights()[0];
    sum += weight;
    double largestPt = -99;
    double zPt = -99;
    int gpId = 0;
    for (const auto & gp: *genParticles) {
      if (gpId == 100) break;
      uint absId = abs(gp.pdgId());
      string marker = "";
      // if (gp.status() == 11 && gp.numberOfDaughters() == 1 && (absId<6 || absId == 21)) marker = " ********";
      cout << gpId << " : " << gp.pt() << " : " << gp.px() << " : " << gp.py() << " : " << gp.pz() << " : " << gp.status() << " : " << gp.pdgId() << " : " << gp.numberOfMothers() << " : " << gp.numberOfDaughters() << marker << endl;
      // cout << gpId << " : " << gp.pt() << " : " << gp.eta() << " : " << gp.phi() << " : " << gp.status() << " : " << gp.pdgId() << " : " << gp.numberOfMothers() << " : " << gp.numberOfDaughters() << endl;
      // if (gp.numberOfMothers() > 0) {
      //   for (int im=0; im< gp.numberOfMothers(); im++) {
      //     auto mum = gp.mother(im);
      //     cout << ".... mum " << mum->pt() << " : " << mum->status() << " : " << mum->pdgId() << endl;
      //   }
      // }
      // if (gp.numberOfDaughters() > 0) {
      //   for (int im=0; im< gp.numberOfDaughters(); im++) {
      //     auto dau = gp.daughter(im);
      //     cout << ".... dau " << dau->pt() << " : " << dau->status() << " : " << dau->pdgId() << endl;
      //   }
      // }

      gpId++;
      if (gp.status() != 11) continue;
      if (absId < 6 || absId == 21) {
        largestPt = max(largestPt, gp.pt());
      }
      if (absId == 23) {
        // if (zPt > 0 && zPt != gp.pt() && ((fabs(zPt-gp.pt()) / zPt) > 0.1 )) {
        //   cout << "Already set zPt - but different by > 10%: " << zPt << " vs " << gp.pt() << endl;
        // }
        if (zPt<0) zPt = gp.pt();
      }
    }
    
    // if (largestPt < 590 || largestPt > 600) continue;
    // if (weight < 250) continue;
    gpId = 0;
    // cout << "Event " << counter << endl;
    // for (const auto & gp: *genParticles) {
    //   if (gpId == 100) break;
    //   uint absId = abs(gp.pdgId());
    //   string marker = "";
    //   // if (gp.status() == 11 && gp.numberOfDaughters() == 1 && (absId<6 || absId == 21)) marker = " ********";
    //   if (gp.status() == 11 && (absId<6 || absId == 21)) marker = " ********";
    //   // cout << gpId << " : " << gp.pt() << " : " << gp.px() << " : " << gp.py() << " : " << gp.pz() << " : " << gp.status() << " : " << gp.pdgId() << " : " << gp.numberOfMothers() << " : " << gp.numberOfDaughters() << marker << " : " << largestPt << " : " << weight << endl;
    //   cout << gpId << " : " << gp.pt() << " : " << gp.eta() << " : " << gp.phi() << " : " << gp.status() << " : " << gp.pdgId() << " : " << gp.numberOfMothers() << " : " << gp.numberOfDaughters() << marker << " : " << largestPt << " : " << weight << endl;
    //   gpId++;
    // }

    // for (const auto & jet: *genJets) {
    //   cout << "Jet: " << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << endl;
    //   // cout << "Jet: " << jet.pt() << " : " << jet.px() << " : " << jet.py() << " : " << jet.pz() << endl;
    // }

    if (largestPt>0) {
      hp.h->Fill(largestPt);
      hp.h_weighted->Fill(largestPt, weight);
      hp.h_weight_vs_kt->Fill(largestPt, weight);
      if (zPt > 0) {
        // cout << largestPt / zPt << endl;
        hp.h_kt_z_ratio_vs_kt->Fill(largestPt, largestPt / zPt);
        hp.h_kt_z_ratio_vs_kt_weighted->Fill(largestPt, largestPt / zPt, weight);
        hp.h_kt_z_ratio_vs_weight->Fill(weight, largestPt / zPt);
      }
      if (genJets->size() >2) {
        hp.h_jet3_jet1_ratio_vs_weight->Fill(weight, genJets->at(2).pt() / genJets->at(0).pt());
      }
    }
  }
  
  TTimeStamp now;
  double totalTime = now - startTime;
  string units = " seconds";
  if (totalTime > 60) {
    totalTime /= 60;
    units = " minutes";
  }
  if (totalTime > 60) {
    totalTime /= 60;
    units = " hours";
  }
  if (totalTime > 24) {
    totalTime /= 24;
    units = " days";
  }
  cout << "Took " << totalTime << units << endl;
  cout << setprecision(9) << "sum: " << sum << endl;
  cout << setprecision(9) << "h->Integral(): " << hp.h->Integral() << endl;
  cout << setprecision(9) << "h_weighted->Integral(): " << hp.h_weighted->Integral() << endl;
  hp.sum = sum;
  hp.xsec = xsec;
  double lumi = 35918;
  hp.h_weighted->Scale(xsec * lumi / sum);
  return hp;
}


void herwigpp_test_spectra() {

  // TChain * chHighPtFlat = new TChain("Events");
  // chHighPtFlat->Add("/nfs/dust/cms/user/aggleton/QG/CMSSW_7_1_22/src/herwigpp_zjet_jetktmin170_pow2p5/job_*_JME-RunIISummer15GS-00016-GEN.root");
  // chHighPtFlat->Add("/nfs/dust/cms/user/aggleton/QG/CMSSW_7_1_22/src/JME-RunIISummer15GS-00016-GEN.root");

  TChain * chHighPt = new TChain("Events");
  chHighPt->Add("/pnfs/desy.de/cms/tier2/store/user/raggleto/PrivateHerwigppDYJetsToLL/DYJetsToLL_M-50_JetKtMin170_TuneCUETHS1_13TeV-herwigpp/herwigpp_zjet_jetktmin170_physical_1M_2020_04_02_v2/200402_214650/0000/JME-RunIISummer15GS-00016_2.root");
  // chHighPt->Add("/nfs/dust/cms/user/aggleton/QG/CMSSW_7_1_22/src/herwigpp_zjet_jetktmin170_physical_gen/job_*_JME-RunIISummer15GS-00016-GEN.root");

  // TChain * chHighPt100 = new TChain("Events");
  // chHighPt100->Add("/nfs/dust/cms/user/aggleton/QG/CMSSW_7_1_22/src/herwigpp_zjet_jetktmin100_physical_gen/job_*_JME-RunIISummer15GS-00016-GEN.root");

  // TChain * chIncl = new TChain("Events");
  // chIncl->Add("/nfs/dust/cms/user/aggleton/QG/CMSSW_7_1_22/src/herwigpp_zjet_physical_gen/job_*_JME-RunIISummer15GS-00016-GEN.root");

  double xsecHighPt = 3.0760;
  double xsecHighPt100 = 1.730e+01;
  double xsecIncl = 379.2;

  // HistPack hpHighPtFlat = fill_histograms(chHighPtFlat, "_highPtFlat", xsecHighPt, 50000);
  HistPack hpHighPt = fill_histograms(chHighPt, "_highPt", xsecHighPt, 5000000);
  // HistPack hpHighPt100 = fill_histograms(chHighPt100, "_highPt100", xsecHighPt100, 5000000);
  // HistPack hpIncl = fill_histograms(chIncl, "_incl", xsecIncl, 2500000);

  // // Scale & Plot
  // TCanvas * canv = new TCanvas("c", "", 800, 600);
  // canv->SetTicks(1, 1);
  // canv->SetLogy();
  // hpIncl.h->SetLineColor(kBlue);
  // hpIncl.h->Draw("HIST E");
  // hpHighPt.h->SetLineColor(kGreen+2);
  // hpHighPt.h->Draw("HIST E SAME");
  // hpHighPtFlat.h->SetLineColor(kRed);
  // hpHighPtFlat.h->Draw("HIST E SAME");
  // TLegend * leg = new TLegend(0.6, 0.6, 0.88, 0.88);
  // leg->AddEntry(hpIncl.h, "Inclusive", "L");
  // leg->AddEntry(hpHighPt.h, "Jet k_{T} > 170 GeV", "L");
  // leg->AddEntry(hpHighPtFlat.h, "Jet k_{T} > 170 GeV, Flat (pow=2.5)", "L");
  // leg->Draw();
  // canv->SaveAs("hpp_zjet_genjetkt.pdf");

  // canv->Clear();
  // hpIncl.h_weighted->SetLineColor(kBlue);
  // hpIncl.h_weighted->Draw("HIST E");
  // hpIncl.h_weighted->SetMinimum(1);
  // hpHighPt.h_weighted->SetLineColor(kGreen+2);
  // hpHighPt.h_weighted->Draw("HIST E SAME");
  // hpHighPtFlat.h_weighted->SetLineColor(kRed);
  // hpHighPtFlat.h_weighted->Draw("HIST E SAME");
  // leg->Draw();
  // canv->SaveAs("hpp_zjet_genjetkt_weighted.pdf");

  // Save hists to file
  // TFile * outf = new TFile("hpp_zjet_genjetkt.root", "RECREATE");
  // TFile * outf = new TFile("hpp_zjet_genjetkt_highPtFlat.root", "RECREATE");
  // TFile * outf = new TFile("hpp_zjet_genjetkt_highPt_100.root", "RECREATE");
  // TFile * outf = new TFile("hpp_zjet_genjetkt_incl.root", "RECREATE");
  
  // hpHighPtFlat.h->Write();
  // hpHighPtFlat.h_weighted->Write();
  // hpHighPtFlat.h_weight_vs_kt->Write();
  // hpHighPtFlat.h_kt_z_ratio_vs_kt->Write();
  // hpHighPtFlat.h_kt_z_ratio_vs_kt_weighted->Write();
  // hpHighPtFlat.h_kt_z_ratio_vs_weight->Write();
  // hpHighPtFlat.h_jet3_jet1_ratio_vs_weight->Write();
 
  // hpHighPt.h->Write();
  // hpHighPt.h_weighted->Write();
  
  // hpHighPt100.h->Write();
  // hpHighPt100.h_weighted->Write();
 
  // hpIncl.h->Write();
  // hpIncl.h_weighted->Write();
 
  // outf->Close();

}

int main(int argc, char * argv[]){
  herwigpp_test_spectra();
}
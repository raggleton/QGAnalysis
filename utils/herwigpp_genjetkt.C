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
#include "TChain.h"

#include "UHH2/core/include/GenInfo.h"
#include "UHH2/core/include/GenParticle.h"

using namespace std;

void herwigpp_genjetkt() {

  // analyse my high pt file
  TChain * ch = new TChain("AnalysisTree");
  ch->Add("/pnfs/desy.de/cms/tier2/store/user/raggleto/DYJetsToLL_M-50_JetKtMin175_TuneCUETHS1_13TeV-herwigpp/crab_DYJetsToLL_M-50_JetKtMin175_15_Mar_20_newRecoJetFlav_fixPuppi_v4_2/200315_110409/0000/Ntuple_*.root");

  int nBins = 500;
  float minX = 50;
  float maxX = 1050;

  TH1D * h = new TH1D("jetkt", ";Gen Jet k_{T} [GeV];N (raw)", nBins, minX, maxX);
  TH1D * h_weighted = new TH1D("jetkt_weights", ";Gen Jet k_{T} [GeV];N (weighted)", nBins, minX, maxX);

  TTreeReader reader(ch);
  TTreeReaderValue<GenInfo> genInfo(reader, "genInfo");
  TTreeReaderValue<vector<GenParticle>> genParticles(reader, "GenParticles");
  double sum = 0.; // MUST be double otherwise loses precision
  int counter = 0;
  int limit = 500000;
  limit = ch->GetEntries();
  int mark = limit/100;
  cout << "===== Doing my high pT sample (" << limit << ")" << endl;
  while (reader.Next() && counter < limit) {
    if (counter % mark == 0) { cout << "..." << counter << endl; }
    // yes I should do a check here, but the ROOT example doesn't work so bah
    float weight = genInfo->weights()[0];
    sum += weight;
    float largestPt = -99;
    for (const auto & gp: *genParticles) {
      if (gp.status() != 11) continue;
      uint absId = abs(gp.pdgId());
      if (absId < 6 || absId == 21) {
        largestPt = max(largestPt, gp.pt());
      }
    }

    if (largestPt>0) {
      h->Fill(largestPt);
      h_weighted->Fill(largestPt, weight);
    }
    counter++;
  }
  cout << setprecision(9) << "sum: " << sum << endl;
  cout << setprecision(9) << "h->Integral(): " << h->Integral() << endl;
  cout << setprecision(9) << "h_weighted->Integral(): " << h_weighted->Integral() << endl;

  // Same for inclusive files
  TChain * chIncl = new TChain("AnalysisTree");
  chIncl->Add("/pnfs/desy.de/cms/tier2/store/user/raggleto/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/crab_DYJetsToLL_M-50_mg-herwig_10_Feb_18_newRecoJetFlav_fixPuppi/180210_102303/0000/Ntuple_*.root");
  chIncl->Add("/pnfs/desy.de/cms/tier2/store/user/raggleto/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/crab_DYJetsToLL_M-50_mg-herwig_10_Feb_18_newRecoJetFlav_fixPuppi/180210_102303/0001/Ntuple_*.root");

  TH1D * hIncl = new TH1D("jetkt_incl", ";Gen Jet k_{T} [GeV];N (raw)", nBins, minX, maxX);
  TH1D * hIncl_weighted = new TH1D("jetkt_weights_incl", ";Gen Jet k_{T} [GeV];N (weighted)", nBins, minX, maxX);

  TTreeReader readerIncl(chIncl);
  TTreeReaderValue<GenInfo> genInfoIncl(readerIncl, "genInfo");
  TTreeReaderValue<vector<GenParticle>> genParticlesIncl(readerIncl, "GenParticles");
  double sumIncl = 0.; // MUST be double otherwise loses precision
  counter = 0;
  limit = 500000;
  limit = chIncl->GetEntries();
  mark = limit/100;
  cout << "===== Doing inclusive sample (" << limit << ")" << endl;
  while (readerIncl.Next() && counter < limit) {
    if (counter % mark == 0) { cout << "..." << counter << endl; }
    // yes I should do a check here, but the ROOT example doesn't work so bah
    float weight = genInfoIncl->weights()[0];
    sumIncl += weight;
    float largestPt = -99;
    for (const auto & gp: *genParticlesIncl) {
      if (gp.status() != 11) continue;
      uint absId = abs(gp.pdgId());
      if (absId < 6 || absId == 21) {
        largestPt = max(largestPt, gp.pt());
      }
    }

    if (largestPt>0) {
      hIncl->Fill(largestPt);
      hIncl_weighted->Fill(largestPt, weight);
    } else {
      cout << "!!! largestPt < 0" << endl;
      for (const auto & gp: *genParticlesIncl) {
        cout << gp.pt() << " : " << gp.pdgId() << " : " << gp.status() << endl;
      }
    }
    counter++;
  }
  cout << setprecision(9) << "sumIncl: " << sumIncl << endl;
  cout << setprecision(9) << "hIncl->Integral(): " << hIncl->Integral() << endl;
  cout << setprecision(9) << "hIncl_weighted->Integral(): " << hIncl_weighted->Integral() << endl;

  float lumi = 35918;
  float xsecHighPt = 3.0760;
  // h->Scale(xsecHighPt*35918 / sum);
  h_weighted->Scale(xsecHighPt*35918 / sum);

  float xsecIncl = 379.2;
  // hIncl->Scale(xsecIncl*35918 / sumIncl);
  hIncl_weighted->Scale(xsecIncl*35918 / sumIncl);

  // Scale & Plot
  TCanvas * canv = new TCanvas("c", "", 800, 600);
  canv->SetTicks(1, 1);
  canv->SetLogy();
  hIncl->SetLineColor(kBlue);
  hIncl->Draw();
  h->SetLineColor(kRed);
  h->Draw("SAME");
  canv->SaveAs("ntuple_hpp_zjet_nocuts_genjetkt.pdf");

  canv->Clear();
  hIncl_weighted->SetLineColor(kBlue);
  hIncl_weighted->Draw();
  h_weighted->SetLineColor(kRed);
  h_weighted->Draw("SAME");
  canv->SaveAs("ntuple_hpp_zjet_nocuts_genjetkt_weighted.pdf");

  // Save hists to file
  TFile * outf = new TFile("hpp_zjet_genjetkt.root", "RECREATE");
  h->Write();
  h_weighted->Write();
  hIncl->Write();
  hIncl_weighted->Write();
  outf->Close();

}

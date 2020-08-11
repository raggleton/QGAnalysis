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

QGAnalysisZPlusJetsGenHists::QGAnalysisZPlusJetsGenHists(Context & ctx, const string & dirname, const string & channel):
Hists(ctx, dirname),
channel_(channel)
{
  if (ctx.get("dataset_type") == "MC") {
    genJets_handle = ctx.get_handle< std::vector<GenJet> > ("GoodGenJets");
    if (channel_ != "muon" && channel_ != "electron")
      throw runtime_error("QGAnalysisZPlusJetsGenHists: channel should be muon or electron");

    genLeptons_handle = ctx.get_handle<std::vector<GenParticle>>(channel_ == "muon" ? "GoodGenMuons" : "GoodGenElectrons");
    gen_weight_handle = ctx.get_handle<double>("gen_weight");

    int nbins_pt = 500;
    float pt_max = 2000;

    int nbins_eta = 200;
    float eta_max = 5;

    int nbins_phi = 128;
    float phi_max = 3.2;

    float lepton_pt_max = 1000;

    // TString zName = "Z";
    TString lName = (channel_ == "muon") ? "#mu" : "e";
    TString zName = lName + lName;
    TString binByVar = "pt_z";

    deta_ll = book<TH1F>("deta_ll", TString::Format(";|#Deltay_{%s}|;", zName.Data()), nbins_eta, 0, 2*eta_max);
    deta_ll_jet1 = book<TH1F>("deta_ll_jet1", TString::Format(";|#Deltay_{%s, genjet 1}|;", zName.Data()), nbins_eta, 0, 2*eta_max);
    dphi_ll = book<TH1F>("dphi_ll", TString::Format(";|#Delta#phi_{%s}|;", zName.Data()), nbins_phi, 0, phi_max);
    dphi_ll_jet1 = book<TH1F>("dphi_ll_jet1", TString::Format(";|#Delta#phi_{%s, genjet 1}|;", zName.Data()), nbins_phi, 0, phi_max);
    eta_jet1 = book<TH1F>("eta_jet1", ";y^{genjet 1};", nbins_eta, -eta_max, eta_max);
    eta_jet2 = book<TH1F>("eta_jet2", ";y^{genjet 2};", nbins_eta, -eta_max, eta_max);
    eta_lepton1 = book<TH1F>("eta_lepton1", TString::Format(";y^{%s1};", lName.Data()), nbins_eta, -eta_max, eta_max);
    eta_lepton2 = book<TH1F>("eta_lepton2", TString::Format(";y^{%s2};", lName.Data()), nbins_eta, -eta_max, eta_max);
    eta_ll = book<TH1F>("eta_ll", TString::Format(";y^{%s};", zName.Data()), nbins_eta, -eta_max, eta_max);
    // eta_z = book<TH1F>("eta_z", ";y^{Z};", nbins_eta, -eta_max, eta_max);
    gen_ht = book<TH1F>("gen_ht", ";H_{T}^{Gen} [GeV];", 500, 0, 5000);
    pt_hat = book<TH1F>("ptHat", ";#hat{p}_{T} [GeV];", 2500, 0, 5000);
    m_ll = book<TH1F>("m_ll", TString::Format(";m_{%s} [GeV];", zName.Data()), 80, 90-40, 90+40);
    n_jets = book<TH1F>("n_jets", ";N_{genjets};", 10, 0, 10);
    n_leptons = book<TH1F>("n_leptons", TString::Format(";N_{%s};", lName.Data()), 10, 0, 10);
    pt_jet1 = book<TH1F>("pt_jet1", ";p_{T}^{genjet 1} [GeV];", nbins_pt, 0, pt_max);
    pt_jet1_pt_jet2_ratio = book<TH1F>("pt_jet1_pt_jet2_ratio", ";p_{T}^{genjet 1} / p_{T}^{genjet 2};", 200, 0, 10);
    pt_jet1_z_ratio = book<TH1F>("pt_jet1_z_ratio", TString::Format(";p_{T}^{genjet 1}/ p_{T}^{%s};", zName.Data()), 50, 0, 5);
    pt_jet2 = book<TH1F>("pt_jet2", ";p_{T}^{genjet 2} [GeV];", nbins_pt, 0, pt_max);
    pt_jet2_z_ratio = book<TH1F>("pt_jet2_z_ratio", TString::Format(";p_{T}^{genjet 2} / p_{T}^{%s};", zName.Data()), 60, 0, 3);
    pt_jet_genHT_ratio = book<TH1F>("pt_jet_genHT_ratio", ";p_{T}^{genjet 1}/GenHT;", 250, 0, 2.5);
    pt_hat_pt_jet_ratio = book<TH1F>("pt_hat_pt_jet_ratio", ";p_{T}^{genjet 1}/#hat{p}_{T};", 500, 0, 10);
    pt_lepton1 = book<TH1F>("pt_lepton1", TString::Format(";p_{T}^{%s1};", lName.Data()), nbins_pt, 0, lepton_pt_max);
    pt_lepton2 = book<TH1F>("pt_lepton2", TString::Format(";p_{T}^{%s2};", lName.Data()), nbins_pt, 0, lepton_pt_max);
    pt_ll = book<TH1F>("pt_ll", TString::Format(";p_{T}^{%s} [GeV];", zName.Data()), nbins_pt, 0, lepton_pt_max);
    pt_z = book<TH1F>("pt_z", ";p_{T}^{Z} [GeV];", nbins_pt, 0, lepton_pt_max);
    q_scale = book<TH1F>("q_scale", ";q scale [GeV];", 250, 0, 500);
    jet_kt = book<TH1F>("jet_kt", ";Parton k_{T} [GeV];", nbins_pt, 0, lepton_pt_max);
    jet_kt_pt_z_ratio = book<TH1F>("jet_kt_pt_z_ratio", ";Parton k_{T} / p_{T}^{Z};", 250, 0, 2.5);
  }

}


void QGAnalysisZPlusJetsGenHists::fill(const Event & event){
  if (event.isRealData) return;
  double weight = event.get(gen_weight_handle);

  const std::vector<GenJet> * genjets = &event.get(genJets_handle);
  int Njets = genjets->size();
  n_jets->Fill(Njets, weight);

  const std::vector<GenParticle> * genleptons = &event.get(genLeptons_handle);
  int Nleptons = genleptons->size();
  n_leptons->Fill(Nleptons, weight);

  // This is the Z made from ME leptons, not final state ones
  GenParticle genZ = findMatrixElementZ(*event.genparticles);
  // cout << genZ.pdgId() << " : " << genZ.status() << " : " << genZ.pt() << endl;
  pt_z->Fill(genZ.pt(), weight);
  // eta_z->Fill(genZ.Rapidity(), weight);

  float jetKt = calcJetKt(*event.genparticles);
  jet_kt->Fill(jetKt, weight);

  float ratio = 0;
  if (genZ.pt() > 0 && jetKt > 0) ratio = jetKt / genZ.pt();
  jet_kt_pt_z_ratio->Fill(ratio, weight);

  if (event.genInfo->binningValues().size() > 0) {
    double ptHat = event.genInfo->binningValues().at(0); // sometimes this is pthat, sometimes it means the hardest outgoing parton eg in H++? TBC
    pt_hat->Fill(ptHat, weight);
  }

  float qScale = event.genInfo->qScale();
  q_scale->Fill(qScale, weight);

  if (Nleptons < 2) return;
  GenParticle l1 = genleptons->at(0);
  GenParticle l2 = genleptons->at(1);
  LorentzVector z_cand = l1.v4() + l2.v4();
  float z_pt = z_cand.pt();
  pt_ll->Fill(z_pt, weight);
  eta_ll->Fill(z_cand.Rapidity(), weight);

  pt_lepton1->Fill(l1.pt(), weight);
  eta_lepton1->Fill(l1.Rapidity(), weight);

  pt_lepton2->Fill(l2.pt(), weight);
  eta_lepton2->Fill(l2.Rapidity(), weight);

  m_ll->Fill(z_cand.M(), weight);

  deta_ll->Fill(fabs(l1.Rapidity() - l2.Rapidity()), weight);
  dphi_ll->Fill(fabs(deltaPhi(l1, l2)), weight);

  // Jets
  if (Njets < 1) return;

  auto jet1 = genjets->at(0);
  float jet1_pt = jet1.pt();
  pt_jet1->Fill(jet1_pt, weight);
  eta_jet1->Fill(jet1.Rapidity(), weight);
  if (event.genInfo->binningValues().size() > 0) {
    double ptHat = event.genInfo->binningValues().at(0);
    pt_hat_pt_jet_ratio->Fill(jet1_pt / ptHat, weight);
  }

  float genHT = calcGenHT(*(event.genparticles));
  gen_ht->Fill(genHT, weight);
  pt_jet_genHT_ratio->Fill(jet1_pt / genHT, weight);

  pt_jet1_z_ratio->Fill(jet1_pt / z_pt, weight);
  deta_ll_jet1->Fill(fabs(z_cand.Rapidity() - jet1.Rapidity()), weight);
  dphi_ll_jet1->Fill(fabs(deltaPhi(jet1, z_cand)), weight);

  if (Njets >= 2) {
    auto jet2 = genjets->at(1);
    pt_jet2->Fill(jet2.pt(), weight);
    eta_jet2->Fill(jet2.Rapidity(), weight);
    pt_jet2_z_ratio->Fill(jet2.pt() / z_pt, weight);
    pt_jet1_pt_jet2_ratio->Fill(jet1_pt / jet2.pt(), weight);
  }

}

/**
 * Find Matrix Element "Z" by reconstructing dileptons
 * Note that due to v.occasional tau decays, and how the decay tree gets smushed,
 * we just look for those PDGIDs explicitly.
 */
// this should probably be in the main MC module and set a handle
GenParticle QGAnalysisZPlusJetsGenHists::findMatrixElementZ(std::vector<GenParticle> & gps) {

  bool foundFirstLepton = false;
  GenParticle firstLepton, secondLepton;
  auto targetPDGID = (channel_ == "muon") ? PDGID::MUON : PDGID::ELECTRON;
  for (const auto & itr : gps) {
    if (abs(itr.pdgId()) == targetPDGID) {
      if (!foundFirstLepton) {
        foundFirstLepton = true;
        firstLepton.set_pdgId(itr.pdgId());
        firstLepton.set_v4(itr.v4());
        firstLepton.set_status(itr.status());
      } else if (itr.pdgId() == -1*firstLepton.pdgId()) {
        secondLepton.set_pdgId(itr.pdgId());
        secondLepton.set_v4(itr.v4());
        secondLepton.set_status(itr.status());
        break;
      }
    }
  }
  if (firstLepton.pt() != 0 && secondLepton.pt() != 0) {
    // int counter = 0;
    // for (const auto & itr : gps) {
    //   counter++;
    //   cout << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << itr.v4().px() << " : " << itr.v4().py() << " : " << itr.v4().pz() << endl;
    //   if (counter == 100) break;
    // }
    // cout << " 1st lepton: " << firstLepton.pdgId() << " : " << firstLepton.status() << " : " << firstLepton.pt() << " : " << firstLepton.v4().px() << " : " << firstLepton.v4().py() << " : " << firstLepton.v4().pz() << endl;
    // cout << " 2nd lepton: " << secondLepton.pdgId() << " : " << secondLepton.status() << " : " << secondLepton.pt() << " : " << secondLepton.v4().px() << " : " << secondLepton.v4().py() << " : " << secondLepton.v4().pz() << endl;
    GenParticle genZ;
    genZ.set_v4(firstLepton.v4() + secondLepton.v4());
    genZ.set_pdgId(PDGID::Z); // in reality could be a photon
    return genZ;
  }

  cout << "**** No Z event: ****" << endl;
  int counter = 0;
  for (const auto & itr : gps) {
    counter++;
    cout << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << itr.v4().px() << " : " << itr.v4().py() << " : " << itr.v4().pz() << endl;
    if (counter == 200) break;
  }
  throw runtime_error("Couldn't find gen Z");
}

QGAnalysisZPlusJetsGenHists::~QGAnalysisZPlusJetsGenHists(){}

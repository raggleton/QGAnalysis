#include "UHH2/QGAnalysis/include/QGAnalysisFlavCompHists.h"
#include "UHH2/core/include/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

QGAnalysisFlavCompHists::QGAnalysisFlavCompHists(Context & ctx, const string & dirname, int useNJets, const string & selection, float jetPtMin, float jetPtMax, float jetEtaMin, float jetEtaMax): 
  Hists(ctx, dirname), 
  useNJets_(useNJets),
  jetPtMin_(jetPtMin),
  jetPtMax_(jetPtMax),
  jetEtaMin_(jetEtaMin),
  jetEtaMax_(jetEtaMax)
{
  if (useNJets_ < 0) useNJets_ = 99999; // Do them all

  if (useNJets_ == 0) throw runtime_error("useNJets should be > 0, or < 0 to use all jets in the event");


  doHerwigReweighting = ctx.get("herwig_reweight_file", "") != "";
  if (doHerwigReweighting) {
    if (selection != "dijet" && selection != "zplusjets") {
      throw runtime_error("selection must be dijet or zplusjets");
    }
    TFile f_weight(ctx.get("herwig_reweight_file", "").c_str());
    if (selection == "dijet")
      reweightHist = (TH1F*) f_weight.Get("dijet_reco");
    else if (selection == "zplusjets")
      reweightHist = (TH1F*) f_weight.Get("zpj_reco");

    if (reweightHist == nullptr) {
      doHerwigReweighting = false;
      cout << "WARNING: could not find reweight hist - not reweighting AnalysisHists!" << endl;
    } else {
      reweightHist->SetDirectory(0);
    }
  }

  h_jet_flavour_vs_hadron_flavour = book<TH2F>("jet_flavour_vs_hadron_flavour", ";parton flavour;hadron flavour", 23, -0.5, 22.5, 6, -0.5, 5.5);
  
}


void QGAnalysisFlavCompHists::fill(const Event & event){
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

  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  // Don't forget to always use the weight when filling.
  double weight = event.weight * herwig_weight;


  // Fill reco jet hists
  for (int i = 0; i < useNJets_; i++) {
    const Jet & thisjet = jets->at(i);
    float pt = thisjet.pt();
    float absEta = fabs(thisjet.eta());
    if ((pt < jetPtMin_) || (pt > jetPtMax_) || (absEta < jetEtaMin_) || (absEta > jetEtaMax_))
      continue;
    int jet_flav = abs(thisjet.partonFlavour());
    int jet_hadronflav = abs(thisjet.hadronFlavour());

    h_jet_flavour_vs_hadron_flavour->Fill(jet_flav, jet_hadronflav, weight);
  }
}

QGAnalysisFlavCompHists::~QGAnalysisFlavCompHists(){}

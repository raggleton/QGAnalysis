#include "UHH2/QGAnalysis/include/LambdaHistsFiller.h"

LambdaHistsFiller::LambdaHistsFiller(const PFLambdaArgs & pfLambdaArgs,
                                     TUnfoldBinning * recoBinning,
                                     const GenLambdaArgs & genLambdaArgs,
                                     TUnfoldBinning * genBinning):
pfLambdaArgs_(pfLambdaArgs),
recoBinning_(recoBinning),
passReco_(false),
genLambdaArgs_(genLambdaArgs),
genBinning_(genBinning),
passGen_(false),
recoBin_(0),
recoBinGenBinning_(0),
genBin_(0),
recoJetPt_(0),
genJetPt_(0)
{}

void LambdaHistsFiller::setPassReco(bool passReco) {
  passReco_ = passReco;
  if (!passReco_) {
    // if fail, then reset bins
    recoBin_ = 0;
    recoBinGenBinning_ = 0;
  }
}

void LambdaHistsFiller::setupReco(float recoJetPt,
                                  const LambdaCalculator<PFParticle> & recoJetCalc,
                                  bool passReco) {
  // reset bin numbers, bools
  recoBin_ = 0;
  recoBinGenBinning_ = 0;
  setPassReco(passReco);
  recoJetPt_ = recoJetPt;

  // calculate reco lambda variable, if those selections passed
  // if value is < 0 then indicates not enough constituents, so fail that selection
  // Get the correct TUnfoldBinning node (ie if underflow or not)
  // also calculate global bin numbers for TUnfold hists
  if (passReco_) {
    double recoLambda = recoJetCalc.getLambda(pfLambdaArgs_.kappa, pfLambdaArgs_.beta, pfLambdaArgs_.id);
    if (recoLambda < 0) {
      setPassReco(false);
    } else {
      bool isUnderflow = recoJetPt_ < Binning::pt_bin_edges_reco[0];
      if (recoBinning_ == nullptr) throw std::runtime_error("LambdaHistsFiller::recoBinning is nullptr");
      std::string nodeName = isUnderflow ? "detector_underflow" : "detector";
      const TUnfoldBinning * thisRecoBinning = recoBinning_->FindNode(nodeName.c_str());
      if (thisRecoBinning == nullptr) throw std::runtime_error("Cannot find TUnfoldBinning node with name " + nodeName);
      recoBin_ = thisRecoBinning->GetGlobalBinNumber(recoLambda, recoJetPt_);

      nodeName = isUnderflow ? "signal_underflow" : "signal";
      if (genBinning_ == nullptr) throw std::runtime_error("LambdaHistsFiller::genBinning is nullptr");
      const TUnfoldBinning * thisGenBinning = genBinning_->FindNode(nodeName.c_str());
      if (thisGenBinning == nullptr) throw std::runtime_error("Cannot find TUnfoldBinning node with name " + nodeName);
      recoBinGenBinning_ = thisGenBinning->GetGlobalBinNumber(recoLambda, recoJetPt_);
    }
  }
}

void LambdaHistsFiller::setPassGen(bool passGen) {
  passGen_ = passGen;
  if (!passGen_) {
    // if fail, then reset bin
    genBin_ = 0;
  }
}

void LambdaHistsFiller::setupGen(float genJetPt,
                                 const LambdaCalculator<GenParticle> & genJetCalc,
                                 bool passGen) {
  // reset bin numbers, bools
  genBin_ = 0;
  setPassGen(passGen);
  genJetPt_ = genJetPt;

  // calculate gen lambda variable, if those selections passed
  // if value is < 0 then indicates not enough constituents, so fail that selection
  // Get the correct TUnfoldBinning node (ie if underflow or not)
  // also calculate global bin numbers for TUnfold hists
  if (passGen_) {
    double genLambda = genJetCalc.getLambda(genLambdaArgs_.kappa, genLambdaArgs_.beta, genLambdaArgs_.id);
    if (genLambda < 0) {
      setPassGen(false);
    } else {
      bool isUnderflow = genJetPt_ < Binning::pt_bin_edges_gen[0];
      std::string nodeName = isUnderflow ? "signal_underflow" : "signal";
      const TUnfoldBinning * thisGenBinning = genBinning_->FindNode(nodeName.c_str());
      if (thisGenBinning == nullptr) throw std::runtime_error("Cannot find TUnfoldBinning node with name " + nodeName);
      genBin_ = thisGenBinning->GetGlobalBinNumber(genLambda, genJetPt_);
    }
  }
}

void LambdaHistsFiller::assignRecoHists(TH1D * allReco,
                                        TH1D * splitReco,
                                        TH1D * allRecoGenBinning,
                                        TH1D * splitRecoGenBinning,
                                        TH1D * allRecoFakes,
                                        TH1D * splitRecoFakes,
                                        TH1D * allRecoFakesGenBinning,
                                        std::vector<TH1D*> * jackknifeVariations,
                                        std::vector<TH1D*> * PDFVariations) {
  recoHist_ = allReco;
  recoSplitHist_ = splitReco;
  recoHistGenBinning_ = allRecoGenBinning;
  recoSplitHistGenBinning_ = splitRecoGenBinning;
  recoHistFakes_ = allRecoFakes;
  recoSplitHistFakes_ = splitRecoFakes;
  recoHistFakesGenBinning_ = allRecoFakesGenBinning;
  recoJackknifeVariations_ = jackknifeVariations;
  recoPDFVariations_ = PDFVariations;
}

void LambdaHistsFiller::assignGenHists(TH1D * allGen,
                                       TH1D * splitGen,
                                       std::vector<TH1D*> * jackknifeVariations,
                                       std::vector<TH1D*> * PDFVariations) {
  genHist_ = allGen;
  genSplitHist_ = splitGen;
  genJackknifeVariations_ = jackknifeVariations;
  genPDFVariations_ = PDFVariations;
}

void LambdaHistsFiller::assignResponseHists(TH2D * allResponse,
                                            TH2D * splitResponse,
                                            std::vector<TH2D*> * jackknifeVariations,
                                            std::vector<TH2D*> * PDFVariations) {
  responseHist_ = allResponse;
  responseSplitHist_ = splitResponse;
  responseJackknifeVariations_ = jackknifeVariations;
  responsePDFVariations_ = PDFVariations;
}

void LambdaHistsFiller::fillRecoTH1(TH1 * h, double weight) {
  if (h != nullptr && passReco()) {
    h->Fill(recoBin_, weight);
  }
}

void LambdaHistsFiller::fillRecoTH1GenBinning(TH1 * h, double weight) {
  if (h != nullptr && passReco()) {
    h->Fill(recoBinGenBinning_, weight);
  }
}

void LambdaHistsFiller::fillGenTH1(TH1 * h, double weight) {
  if (h != nullptr && passGen()) {
    h->Fill(genBin_, weight);
  }
}

void LambdaHistsFiller::fillResponseTH2(TH2 * h, double recoWeight, double genWeight) {
  // How to handle the fact that the matching reco jet may not be the one stored?
  // Assume the user has re-called setupReco with the matched reco jet
  if (h != nullptr && passGen()) {
    // fill if we pass gen requirements, irrespective of reco requirements
    // (will be bin 0 if fail)
    // fill twice: once normally with total event weight,
    // then again but in global underflow with (gen weight - total event weight).
    // This ensures that the projection on the gen axis matches the 1D gen distribution
    // (which would only have the gen weight)
    double totalWeight = recoWeight * genWeight;
    h->Fill(genBin_, recoBin_, totalWeight);
    double corrWeight = genWeight - totalWeight;
    int underflowBin = 0;
    // this is the underflow bin, since the histogram's bin edges are the global bin numbers,
    // not physical values
    h->Fill(genBin_, underflowBin, corrWeight);
  }
}

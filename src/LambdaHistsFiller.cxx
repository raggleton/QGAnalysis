#include "UHH2/QGAnalysis/include/LambdaHistsFiller.h"

LambdaHistsFiller::LambdaHistsFiller(const LambdaArgs & pfLambdaArgs,
                                     TUnfoldBinning * recoBinning,
                                     const LambdaArgs & genLambdaArgs,
                                     TUnfoldBinning * genBinning,
                                     const std::string & lambdaName):
pfLambdaArgs_(pfLambdaArgs),
recoBinning_(recoBinning),
passReco_(false),
genLambdaArgs_(genLambdaArgs),
genBinning_(genBinning),
passGen_(false),
lambdaName_(lambdaName),
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
    double recoLambda = recoJetCalc.getLambda(pfLambdaArgs_);
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
    double genLambda = genJetCalc.getLambda(genLambdaArgs_);
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


void LambdaHistsFiller::setupRecoHists(Context & ctx, const std::string & dirname, bool doSplitHists, bool doFakeHists, int nJackknifeVariations, int nPDFVariations) {
  TH1 * h_tu_reco_tmp = recoBinning_->CreateHistogram("hist_reco"); // template hist, since copying is quicker
  recoHist_ = copy_book_th1d(ctx, dirname, h_tu_reco_tmp, strWithLambdaName("hist_%s_reco_all"));

  // reco variable, but gen binning
  TH1 * h_tu_gen_tmp = genBinning_->CreateHistogram("hist_reco_gen_binning");
  recoHistGenBinning_ = copy_book_th1d(ctx, dirname, h_tu_gen_tmp, strWithLambdaName("hist_%s_reco_gen_binning"));;

  if (doSplitHists) {
    recoSplitHist_ = copy_book_th1d(ctx, dirname, h_tu_reco_tmp, strWithLambdaName("hist_%s_reco_split"));
    recoSplitHistGenBinning_ = copy_book_th1d(ctx, dirname, h_tu_gen_tmp, strWithLambdaName("hist_%s_reco_gen_binning_split"));
  }

  // for fakes
  if (doFakeHists) {
    recoHistFakes_ = copy_book_th1d(ctx, dirname, h_tu_reco_tmp, strWithLambdaName("hist_%s_reco_fake_all"));
    recoHistFakesGenBinning_ = copy_book_th1d(ctx, dirname, h_tu_gen_tmp, strWithLambdaName("hist_%s_reco_fake_gen_binning"));

    if (doSplitHists) recoSplitHistFakes_ = copy_book_th1d(ctx, dirname, h_tu_reco_tmp, strWithLambdaName("hist_%s_reco_fake_split"));
  }

  for (int i=0; i < nJackknifeVariations; i++) {
    recoJackknifeVariations_.push_back(copy_book_th1d(ctx, dirname, h_tu_reco_tmp, TString::Format("hist_%s_reco_all_jackknife_%d", lambdaName_.c_str(), i)));
  }

  for (int i=0; i < nPDFVariations; i++) {
    recoPDFVariations_.push_back(copy_book_th1d(ctx, dirname, h_tu_reco_tmp, TString::Format("hist_%s_reco_all_PDF_%d", lambdaName_.c_str(), i)));
  }

  delete h_tu_reco_tmp;
  delete h_tu_gen_tmp;
}


void LambdaHistsFiller::setupGenHists(Context & ctx, const std::string & dirname, bool doSplitHists, int nJackknifeVariations, int nPDFVariations) {
  TH1 * h_tu_gen_tmp = genBinning_->CreateHistogram("hist_gen_binning");
  genHist_ = copy_book_th1d(ctx, dirname, h_tu_gen_tmp, strWithLambdaName("hist_%s_truth_all"));

  if (doSplitHists) {
    genSplitHist_ = copy_book_th1d(ctx, dirname, h_tu_gen_tmp, strWithLambdaName("hist_%s_truth_split"));
  }

  for (int i=0; i < nJackknifeVariations; i++) {
    genJackknifeVariations_.push_back(copy_book_th1d(ctx, dirname, h_tu_gen_tmp, TString::Format("hist_%s_truth_all_jackknife_%d", lambdaName_.c_str(), i)));
  }

  for (int i=0; i < nPDFVariations; i++) {
    genPDFVariations_.push_back(copy_book_th1d(ctx, dirname, h_tu_gen_tmp, TString::Format("hist_%s_truth_all_PDF_%d", lambdaName_.c_str(), i)));
  }

  delete h_tu_gen_tmp;
}


void LambdaHistsFiller::setupResponseHists(Context & ctx, const std::string & dirname, bool doSplitHists, int nJackknifeVariations, int nPDFVariations) {
  TH2 * h_tu_response_tmp = TUnfoldBinning::CreateHistogramOfMigrations(genBinning_, recoBinning_, "tu_GenReco");
  responseHist_ = copy_book_th2d(ctx, dirname, h_tu_response_tmp, strWithLambdaName("tu_%s_GenReco_all"));

  if (doSplitHists) {
    responseSplitHist_ = copy_book_th2d(ctx, dirname, h_tu_response_tmp, strWithLambdaName("tu_%s_GenReco_split"));
  }

  for (int i=0; i < nJackknifeVariations; i++) {
    responseJackknifeVariations_.push_back(copy_book_th2d(ctx, dirname, h_tu_response_tmp, TString::Format("tu_%s_GenReco_all_jackknife_%d", lambdaName_.c_str(), i)));
  }

  for (int i=0; i < nPDFVariations; i++) {
    responsePDFVariations_.push_back(copy_book_th2d(ctx, dirname, h_tu_response_tmp, TString::Format("tu_%s_GenReco_all_PDF_%d", lambdaName_.c_str(), i)));
  }

  delete h_tu_response_tmp;
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

TH1D * LambdaHistsFiller::copy_book_th1d(Context & ctx, const std::string & dirname, TH1 * h, const std::string & newName) {
  // DO NOT USE h->GetXaxis()->GetXbins()->GetArray() to clone bins, it just doesn't work
  return book<TH1D>(ctx,
                    dirname,
                    newName.c_str(),
                    h->GetTitle(),
                    h->GetNbinsX(),
                    h->GetXaxis()->GetXmin(),
                    h->GetXaxis()->GetXmax());
}


TH1D * LambdaHistsFiller::copy_book_th1d(Context & ctx, const std::string & dirname, TH1 * h, const TString & newName) {
  // DO NOT USE h->GetXaxis()->GetXbins()->GetArray() to clone bins, it just doesn't work
  return book<TH1D>(ctx,
                    dirname,
                    newName,
                    h->GetTitle(),
                    h->GetNbinsX(),
                    h->GetXaxis()->GetXmin(),
                    h->GetXaxis()->GetXmax());
}

TH2D * LambdaHistsFiller::copy_book_th2d(Context & ctx, const std::string & dirname, TH2 * h, const std::string & newName) {
  return book<TH2D>(ctx,
                    dirname,
                    newName.c_str(),
                    h->GetTitle(),
                    h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                    h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

}


TH2D * LambdaHistsFiller::copy_book_th2d(Context & ctx, const std::string & dirname, TH2 * h, const TString & newName) {
  return book<TH2D>(ctx,
                    dirname,
                    newName,
                    h->GetTitle(),
                    h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                    h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

}


TString LambdaHistsFiller::strWithLambdaName(const std::string & str) {
  return TString::Format(str.c_str(), lambdaName_.c_str());
}
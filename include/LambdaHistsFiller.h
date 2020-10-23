#pragma once

#include "UHH2/QGAnalysis/include/QGAddModules.h"

#include "TUnfoldBinning.h"
#include "TH1.h"
#include "TH2.h"

/**
 * Helper class to handle filling all the hists of one lambda variable configuration,
 * both gen & reco
 *
 * TODO: store/setup the hists (all/split/jackknife/PDF etc)
 */
class LambdaHistsFiller {
public:
    LambdaHistsFiller(const PFLambdaArgs & pfLambdaArgs,
                      TUnfoldBinning * recoBinning,
                      const GenLambdaArgs & genLambdaArgs,
                      TUnfoldBinning * genBinning);
    // setup for particular reco, gen jets - recalculate variable and TUnfold bin numbers
    void setPassReco(bool passReco);
    void setupReco(float recoJetPt,
                   const LambdaCalculator<PFParticle> & recoJetCalc,
                   bool passReco);
    void setPassGen(bool passGen);
    void setupGen(float genJetPt,
                  const LambdaCalculator<GenParticle> & genJetCalc,
                  bool passGen);
    // setup hists
    void assignRecoHists(TH1D * allReco,
                         TH1D * splitReco,
                         TH1D * allRecoGenBinning,
                         TH1D * splitRecoGenBinning,
                         TH1D * allRecoFakes,
                         TH1D * splitRecoFakes,
                         TH1D * allRecoFakesGenBinning,
                         std::vector<TH1D*> * jackknifeVariations,
                         std::vector<TH1D*> * PDFVariations);
    void assignGenHists(TH1D * allGen,
                        TH1D * splitGen,
                        std::vector<TH1D*> * jackknifeVariations,
                        std::vector<TH1D*> * PDFVariations);
    void assignResponseHists(TH2D * allResponse,
                             TH2D * splitResponse,
                             std::vector<TH2D*> * jackknifeVariations,
                             std::vector<TH2D*> * PDFVariations);
    // Generic hist filling methods
    // These automatically account for check on # constituents
    void fillRecoTH1(TH1 * h, double weight);
    void fillRecoTH1GenBinning(TH1 * h, double weight);
    void fillGenTH1(TH1 * h, double weight);
    void fillResponseTH2(TH2 * h, double recoWeight, double genWeight);

    // getters
    float recoJetPt() { return recoJetPt_; }
    bool passReco() { return passReco_; }

    float genJetPt() { return genJetPt_; }
    bool passGen() { return passGen_; }

    TH1D * recoHist() { return recoHist_; }
    TH1D * recoSplitHist() { return recoSplitHist_; }
    TH1D * recoHistGenBinning() { return recoHistGenBinning_; }
    TH1D * recoSplitHistGenBinning() { return recoSplitHistGenBinning_; }
    TH1D * recoHistFakes() { return recoHistFakes_; }
    TH1D * recoSplitHistFakes() { return recoSplitHistFakes_; }
    TH1D * recoHistFakesGenBinning() { return recoHistFakesGenBinning_; }
    std::vector<TH1D*> * recoJackknifeVariations() { return recoJackknifeVariations_; }
    std::vector<TH1D*> * recoPDFVariations() { return recoPDFVariations_; }

    TH1D *genHist() { return genHist_; }
    TH1D *genSplitHist() { return genSplitHist_; }
    std::vector<TH1D*> *genJackknifeVariations() { return genJackknifeVariations_; }
    std::vector<TH1D*> *genPDFVariations() { return genPDFVariations_; }

    TH2D * responseHist() { return responseHist_; }
    TH2D * responseSplitHist() { return responseSplitHist_; }
    std::vector<TH2D*> *responseJackknifeVariations() { return responseJackknifeVariations_; }
    std::vector<TH2D*> *responsePDFVariations() { return responsePDFVariations_; }

private:
    PFLambdaArgs pfLambdaArgs_;
    TUnfoldBinning *recoBinning_;
    bool passReco_;
    GenLambdaArgs genLambdaArgs_;
    TUnfoldBinning *genBinning_;
    bool passGen_;
    int recoBin_, recoBinGenBinning_, genBin_;
    float recoJetPt_, genJetPt_;
    TH1D *recoHist_, *recoSplitHist_,
         *recoHistGenBinning_, *recoSplitHistGenBinning_,
         *recoHistFakes_, *recoSplitHistFakes_, *recoHistFakesGenBinning_;
    std::vector<TH1D*> *recoJackknifeVariations_, *recoPDFVariations_;
    TH1D *genHist_, *genSplitHist_;
    std::vector<TH1D*> *genJackknifeVariations_, *genPDFVariations_;
    TH2D *responseHist_, *responseSplitHist_;
    std::vector<TH2D*> *responseJackknifeVariations_, *responsePDFVariations_;
};
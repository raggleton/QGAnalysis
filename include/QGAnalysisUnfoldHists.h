#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"

#include "TUnfoldBinning.h"
#include "TRandom3.h"


namespace uhh2examples {
class LambdaHistsFiller;

/**  \brief Book and fill main analysis histograms
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class QGAnalysisUnfoldHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisUnfoldHists(uhh2::Context & ctx, const std::string & dirname,
                          int useNJets, bool doGroomed,
                          const std::string & selection,
                          const std::string & reco_sel_handle_name, const std::string & gen_sel_handle_name,
                          const std::string & reco_jetlambda_handle_name, const std::string & gen_jetlambda_handle_name);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisUnfoldHists();
protected:
    TH1D * copy_book_th1d(TH1 * h, const std::string & newName);
    TH2D * copy_book_th2d(TH2 * h, const std::string & newName);

    TH1D *h_fake_counter_raw, *h_fake_counter_weighted;
    TH1D *h_fake_counter_charged_raw, *h_fake_counter_charged_weighted;

    // For pt-only unfolding
    TUnfoldBinning *detector_tu_binning_pt, *detector_distribution_pt, *detector_distribution_underflow_pt;
    TUnfoldBinning *generator_tu_binning_pt, *generator_distribution_pt, *generator_distribution_underflow_pt;
    TH2D *h_tu_response_pt, *h_tu_response_pt_split;
    TH1D *h_tu_reco_pt, *h_tu_reco_pt_split, *h_tu_reco_pt_fake, *h_tu_reco_pt_fake_split;
    TH1D *h_tu_gen_pt, *h_tu_gen_pt_split;
    std::vector<TH1D*> h_tu_reco_pt_PDF_variations, h_tu_gen_pt_PDF_variations;
    std::vector<TH2D*> h_tu_response_pt_PDF_variations;

    // reco jet hists
    // For unfolding with TUnfold
    TUnfoldBinning *detector_tu_binning_LHA, *detector_distribution_LHA, *detector_distribution_underflow_LHA;
    TUnfoldBinning *generator_tu_binning_LHA, *generator_distribution_LHA, *generator_distribution_underflow_LHA;
    TH2D *h_tu_response_LHA, *h_tu_response_LHA_split; // response matrices. *_split histograms use half the MC to act as a check using independent sets of events
    TH1D *h_tu_reco_LHA, *h_tu_reco_LHA_fake, *h_tu_reco_LHA_split, *h_tu_reco_LHA_fake_split;  // detector quantities with reco binning
    TH1D *h_tu_reco_LHA_gen_binning, *h_tu_reco_LHA_fake_gen_binning, *h_tu_reco_LHA_gen_binning_split; // detector quantities with gen binning
    TH1D *h_tu_gen_LHA, *h_tu_gen_LHA_split; // truth histograms
    std::vector<TH1D*> h_tu_reco_LHA_PDF_variations, h_tu_gen_LHA_PDF_variations;
    std::vector<TH2D*> h_tu_response_LHA_PDF_variations;
    std::vector<TH1D*> h_tu_reco_LHA_jackknife_variations, h_tu_gen_LHA_jackknife_variations;
    std::vector<TH2D*> h_tu_response_LHA_jackknife_variations;

    TUnfoldBinning *detector_tu_binning_puppiMultiplicity, *detector_distribution_puppiMultiplicity, *detector_distribution_underflow_puppiMultiplicity;
    TUnfoldBinning *generator_tu_binning_puppiMultiplicity, *generator_distribution_puppiMultiplicity, *generator_distribution_underflow_puppiMultiplicity;
    TH2D *h_tu_response_puppiMultiplicity, *h_tu_response_puppiMultiplicity_split;
    TH1D *h_tu_reco_puppiMultiplicity, *h_tu_reco_puppiMultiplicity_fake, *h_tu_reco_puppiMultiplicity_split, *h_tu_reco_puppiMultiplicity_fake_split;
    TH1D *h_tu_reco_puppiMultiplicity_gen_binning, *h_tu_reco_puppiMultiplicity_fake_gen_binning, *h_tu_reco_puppiMultiplicity_gen_binning_split;
    TH1D *h_tu_gen_puppiMultiplicity, *h_tu_gen_puppiMultiplicity_split;
    std::vector<TH1D*> h_tu_reco_puppiMultiplicity_PDF_variations, h_tu_gen_puppiMultiplicity_PDF_variations;
    std::vector<TH2D*> h_tu_response_puppiMultiplicity_PDF_variations;
    std::vector<TH1D*> h_tu_reco_puppiMultiplicity_jackknife_variations, h_tu_gen_puppiMultiplicity_jackknife_variations;
    std::vector<TH2D*> h_tu_response_puppiMultiplicity_jackknife_variations;

    TUnfoldBinning *detector_tu_binning_pTD, *detector_distribution_pTD, *detector_distribution_underflow_pTD;
    TUnfoldBinning *generator_tu_binning_pTD, *generator_distribution_pTD, *generator_distribution_underflow_pTD;
    TH2D *h_tu_response_pTD, *h_tu_response_pTD_split;
    TH1D *h_tu_reco_pTD, *h_tu_reco_pTD_fake, *h_tu_reco_pTD_split, *h_tu_reco_pTD_fake_split;
    TH1D *h_tu_reco_pTD_gen_binning, *h_tu_reco_pTD_fake_gen_binning, *h_tu_reco_pTD_gen_binning_split;
    TH1D *h_tu_gen_pTD, *h_tu_gen_pTD_split;
    std::vector<TH1D* > h_tu_reco_pTD_PDF_variations, h_tu_gen_pTD_PDF_variations;
    std::vector<TH2D*> h_tu_response_pTD_PDF_variations;
    std::vector<TH1D* > h_tu_reco_pTD_jackknife_variations, h_tu_gen_pTD_jackknife_variations;
    std::vector<TH2D*> h_tu_response_pTD_jackknife_variations;

    TUnfoldBinning *detector_tu_binning_thrust, *detector_distribution_thrust, *detector_distribution_underflow_thrust;
    TUnfoldBinning *generator_tu_binning_thrust, *generator_distribution_thrust, *generator_distribution_underflow_thrust;
    TH2D *h_tu_response_thrust, *h_tu_response_thrust_split;
    TH1D *h_tu_reco_thrust, *h_tu_reco_thrust_fake, *h_tu_reco_thrust_split, *h_tu_reco_thrust_fake_split;
    TH1D *h_tu_reco_thrust_gen_binning, *h_tu_reco_thrust_fake_gen_binning, *h_tu_reco_thrust_gen_binning_split;
    TH1D *h_tu_gen_thrust, *h_tu_gen_thrust_split;
    std::vector<TH1D* > h_tu_reco_thrust_PDF_variations, h_tu_gen_thrust_PDF_variations;
    std::vector<TH2D*> h_tu_response_thrust_PDF_variations;
    std::vector<TH1D* > h_tu_reco_thrust_jackknife_variations, h_tu_gen_thrust_jackknife_variations;
    std::vector<TH2D*> h_tu_response_thrust_jackknife_variations;

    TUnfoldBinning *detector_tu_binning_width, *detector_distribution_width, *detector_distribution_underflow_width;
    TUnfoldBinning *generator_tu_binning_width, *generator_distribution_width, *generator_distribution_underflow_width;
    TH2D *h_tu_response_width, *h_tu_response_width_split;
    TH1D *h_tu_reco_width, *h_tu_reco_width_fake, *h_tu_reco_width_split, *h_tu_reco_width_fake_split;
    TH1D *h_tu_reco_width_gen_binning, *h_tu_reco_width_fake_gen_binning, *h_tu_reco_width_gen_binning_split;
    TH1D *h_tu_gen_width, *h_tu_gen_width_split;
    std::vector<TH1D* > h_tu_reco_width_PDF_variations, h_tu_gen_width_PDF_variations;
    std::vector<TH2D*> h_tu_response_width_PDF_variations;
    std::vector<TH1D* > h_tu_reco_width_jackknife_variations, h_tu_gen_width_jackknife_variations;
    std::vector<TH2D*> h_tu_response_width_jackknife_variations;

    // charged-only versions:
    TUnfoldBinning *detector_tu_binning_LHA_charged, *detector_distribution_LHA_charged, *detector_distribution_underflow_LHA_charged;
    TUnfoldBinning *generator_tu_binning_LHA_charged, *generator_distribution_LHA_charged, *generator_distribution_underflow_LHA_charged;
    TH2D *h_tu_response_LHA_charged, *h_tu_response_LHA_charged_split;
    TH1D *h_tu_reco_LHA_charged, *h_tu_reco_LHA_charged_fake, *h_tu_reco_LHA_charged_split, *h_tu_reco_LHA_charged_fake_split;
    TH1D *h_tu_reco_LHA_charged_gen_binning, *h_tu_reco_LHA_charged_fake_gen_binning, *h_tu_reco_LHA_charged_gen_binning_split;
    TH1D *h_tu_gen_LHA_charged, *h_tu_gen_LHA_charged_split;
    std::vector<TH1D* > h_tu_reco_LHA_charged_PDF_variations, h_tu_gen_LHA_charged_PDF_variations;
    std::vector<TH2D*> h_tu_response_LHA_charged_PDF_variations;
    std::vector<TH1D* > h_tu_reco_LHA_charged_jackknife_variations, h_tu_gen_LHA_charged_jackknife_variations;
    std::vector<TH2D*> h_tu_response_LHA_charged_jackknife_variations;

    TUnfoldBinning *detector_tu_binning_puppiMultiplicity_charged, *detector_distribution_puppiMultiplicity_charged, *detector_distribution_underflow_puppiMultiplicity_charged;
    TUnfoldBinning *generator_tu_binning_puppiMultiplicity_charged, *generator_distribution_puppiMultiplicity_charged, *generator_distribution_underflow_puppiMultiplicity_charged;
    TH2D *h_tu_response_puppiMultiplicity_charged, *h_tu_response_puppiMultiplicity_charged_split;
    TH1D *h_tu_reco_puppiMultiplicity_charged, *h_tu_reco_puppiMultiplicity_charged_fake, *h_tu_reco_puppiMultiplicity_charged_split, *h_tu_reco_puppiMultiplicity_charged_fake_split;
    TH1D *h_tu_reco_puppiMultiplicity_charged_gen_binning, *h_tu_reco_puppiMultiplicity_charged_fake_gen_binning, *h_tu_reco_puppiMultiplicity_charged_gen_binning_split;
    TH1D *h_tu_gen_puppiMultiplicity_charged, *h_tu_gen_puppiMultiplicity_charged_split;
    std::vector<TH1D*> h_tu_reco_puppiMultiplicity_charged_PDF_variations, h_tu_gen_puppiMultiplicity_charged_PDF_variations;
    std::vector<TH2D*> h_tu_response_puppiMultiplicity_charged_PDF_variations;
    std::vector<TH1D*> h_tu_reco_puppiMultiplicity_charged_jackknife_variations, h_tu_gen_puppiMultiplicity_charged_jackknife_variations;
    std::vector<TH2D*> h_tu_response_puppiMultiplicity_charged_jackknife_variations;

    TUnfoldBinning *detector_tu_binning_pTD_charged, *detector_distribution_pTD_charged, *detector_distribution_underflow_pTD_charged;
    TUnfoldBinning *generator_tu_binning_pTD_charged, *generator_distribution_pTD_charged, *generator_distribution_underflow_pTD_charged;
    TH2D *h_tu_response_pTD_charged, *h_tu_response_pTD_charged_split;
    TH1D *h_tu_reco_pTD_charged, *h_tu_reco_pTD_charged_fake, *h_tu_reco_pTD_charged_split, *h_tu_reco_pTD_charged_fake_split;
    TH1D *h_tu_reco_pTD_charged_gen_binning, *h_tu_reco_pTD_charged_fake_gen_binning, *h_tu_reco_pTD_charged_gen_binning_split;
    TH1D *h_tu_gen_pTD_charged, *h_tu_gen_pTD_charged_split;
    std::vector<TH1D* > h_tu_reco_pTD_charged_PDF_variations, h_tu_gen_pTD_charged_PDF_variations;
    std::vector<TH2D*> h_tu_response_pTD_charged_PDF_variations;
    std::vector<TH1D* > h_tu_reco_pTD_charged_jackknife_variations, h_tu_gen_pTD_charged_jackknife_variations;
    std::vector<TH2D*> h_tu_response_pTD_charged_jackknife_variations;

    TUnfoldBinning *detector_tu_binning_thrust_charged, *detector_distribution_thrust_charged, *detector_distribution_underflow_thrust_charged;
    TUnfoldBinning *generator_tu_binning_thrust_charged, *generator_distribution_thrust_charged, *generator_distribution_underflow_thrust_charged;
    TH2D *h_tu_response_thrust_charged, *h_tu_response_thrust_charged_split;
    TH1D *h_tu_reco_thrust_charged, *h_tu_reco_thrust_charged_fake, *h_tu_reco_thrust_charged_split, *h_tu_reco_thrust_charged_fake_split;
    TH1D *h_tu_reco_thrust_charged_gen_binning, *h_tu_reco_thrust_charged_fake_gen_binning, *h_tu_reco_thrust_charged_gen_binning_split;
    TH1D *h_tu_gen_thrust_charged, *h_tu_gen_thrust_charged_split;
    std::vector<TH1D* > h_tu_reco_thrust_charged_PDF_variations, h_tu_gen_thrust_charged_PDF_variations;
    std::vector<TH2D*> h_tu_response_thrust_charged_PDF_variations;
    std::vector<TH1D* > h_tu_reco_thrust_charged_jackknife_variations, h_tu_gen_thrust_charged_jackknife_variations;
    std::vector<TH2D*> h_tu_response_thrust_charged_jackknife_variations;

    TUnfoldBinning *detector_tu_binning_width_charged, *detector_distribution_width_charged, *detector_distribution_underflow_width_charged;
    TUnfoldBinning *generator_tu_binning_width_charged, *generator_distribution_width_charged, *generator_distribution_underflow_width_charged;
    TH2D *h_tu_response_width_charged, *h_tu_response_width_charged_split;
    TH1D *h_tu_reco_width_charged, *h_tu_reco_width_charged_fake, *h_tu_reco_width_charged_split, *h_tu_reco_width_charged_fake_split;
    TH1D *h_tu_reco_width_charged_gen_binning, *h_tu_reco_width_charged_fake_gen_binning, *h_tu_reco_width_charged_gen_binning_split;
    TH1D *h_tu_gen_width_charged, *h_tu_gen_width_charged_split;
    std::vector<TH1D* > h_tu_reco_width_charged_PDF_variations, h_tu_gen_width_charged_PDF_variations;
    std::vector<TH2D*> h_tu_response_width_charged_PDF_variations;
    std::vector<TH1D* > h_tu_reco_width_charged_jackknife_variations, h_tu_gen_width_charged_jackknife_variations;
    std::vector<TH2D*> h_tu_response_width_charged_jackknife_variations;

    int useNJets_;
    bool doGroomed_;

    uhh2::Event::Handle<std::vector<GenJetLambdaBundle> > genJetsLambda_handle;
    uhh2::Event::Handle<std::vector<JetLambdaBundle> > jetsLambda_handle;
    uhh2::Event::Handle<double> gen_weight_handle, pt_binning_reco_handle, pt_binning_gen_handle;
    uhh2::Event::Handle<bool> pass_reco_handle;
    uhh2::Event::Handle<bool> pass_gen_handle;
    bool is_mc_;

    TRandom3 rand_;
    bool doMCsplit_; // to do independent response & 1D samples to test closure using same MC dataset
    bool doJackknifeVariations_; // do jackknife variations of response matrix
    bool doPDFvariations_;
    bool useBinningValue_;
    const uint N_JACKKNIFE_VARIATIONS;
    const int N_PDF_VARIATIONS;
    int eventCounter_;

std::unique_ptr<LambdaHistsFiller> LHA_hist_filler,
                                    puppiMultiplicity_hist_filler,
                                    pTD_hist_filler,
                                    thrust_hist_filler,
                                    width_hist_filler,
                                    LHA_charged_hist_filler,
                                    puppiMultiplicity_charged_hist_filler,
                                    pTD_charged_hist_filler,
                                    thrust_charged_hist_filler,
                                    width_charged_hist_filler;

    std::vector<LambdaHistsFiller*> chargedPlusNeutralHistFillers, chargedOnlyHistFillers, allHistFillers;
};

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
    void fillFakesTH1(TH1 * h, double weight);
    void fillFakesTH1GenBinning(TH1 * h, double weight);
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

} // end of namespace
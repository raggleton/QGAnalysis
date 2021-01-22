#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/QGAnalysis/include/LambdaHistsFiller.h"

#include "TUnfoldBinning.h"
#include "TRandom3.h"


namespace uhh2examples {
// class LambdaHistsFiller;

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

    // For unfolding with TUnfold
    // charged+neutral
    TUnfoldBinning *detector_tu_binning_LHA, *detector_distribution_LHA, *detector_distribution_underflow_LHA;
    TUnfoldBinning *generator_tu_binning_LHA, *generator_distribution_LHA, *generator_distribution_underflow_LHA;

    TUnfoldBinning *detector_tu_binning_puppiMultiplicity, *detector_distribution_puppiMultiplicity, *detector_distribution_underflow_puppiMultiplicity;
    TUnfoldBinning *generator_tu_binning_puppiMultiplicity, *generator_distribution_puppiMultiplicity, *generator_distribution_underflow_puppiMultiplicity;

    TUnfoldBinning *detector_tu_binning_pTD, *detector_distribution_pTD, *detector_distribution_underflow_pTD;
    TUnfoldBinning *generator_tu_binning_pTD, *generator_distribution_pTD, *generator_distribution_underflow_pTD;

    TUnfoldBinning *detector_tu_binning_thrust, *detector_distribution_thrust, *detector_distribution_underflow_thrust;
    TUnfoldBinning *generator_tu_binning_thrust, *generator_distribution_thrust, *generator_distribution_underflow_thrust;

    TUnfoldBinning *detector_tu_binning_width, *detector_distribution_width, *detector_distribution_underflow_width;
    TUnfoldBinning *generator_tu_binning_width, *generator_distribution_width, *generator_distribution_underflow_width;

    // charged-only versions:
    TUnfoldBinning *detector_tu_binning_LHA_charged, *detector_distribution_LHA_charged, *detector_distribution_underflow_LHA_charged;
    TUnfoldBinning *generator_tu_binning_LHA_charged, *generator_distribution_LHA_charged, *generator_distribution_underflow_LHA_charged;

    TUnfoldBinning *detector_tu_binning_puppiMultiplicity_charged, *detector_distribution_puppiMultiplicity_charged, *detector_distribution_underflow_puppiMultiplicity_charged;
    TUnfoldBinning *generator_tu_binning_puppiMultiplicity_charged, *generator_distribution_puppiMultiplicity_charged, *generator_distribution_underflow_puppiMultiplicity_charged;

    TUnfoldBinning *detector_tu_binning_pTD_charged, *detector_distribution_pTD_charged, *detector_distribution_underflow_pTD_charged;
    TUnfoldBinning *generator_tu_binning_pTD_charged, *generator_distribution_pTD_charged, *generator_distribution_underflow_pTD_charged;

    TUnfoldBinning *detector_tu_binning_thrust_charged, *detector_distribution_thrust_charged, *detector_distribution_underflow_thrust_charged;
    TUnfoldBinning *generator_tu_binning_thrust_charged, *generator_distribution_thrust_charged, *generator_distribution_underflow_thrust_charged;

    TUnfoldBinning *detector_tu_binning_width_charged, *detector_distribution_width_charged, *detector_distribution_underflow_width_charged;
    TUnfoldBinning *generator_tu_binning_width_charged, *generator_distribution_width_charged, *generator_distribution_underflow_width_charged;

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
    uint N_JACKKNIFE_VARIATIONS;
    int N_PDF_VARIATIONS;
    int eventCounter_;

    std::unique_ptr<LambdaHistsFiller>  LHA_hist_filler,
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

} // end of namespace
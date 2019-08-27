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
    TH1F * copy_book_th1f(TH1 * h, const std::string & newName);
    TH2F * copy_book_th2f(TH2 * h, const std::string & newName);

    TH1F *h_fake_counter_raw, *h_fake_counter_weighted;

    // For pt-only unfolding
    TUnfoldBinning *detector_tu_binning_pt, *detector_distribution_pt, *detector_distribution_underflow_pt;
    TUnfoldBinning *generator_tu_binning_pt, *generator_distribution_pt, *generator_distribution_underflow_pt;
    TH2F *h_tu_response_pt, *h_tu_response_pt_split;
    TH1F *h_tu_reco_pt, *h_tu_reco_pt_split, *h_tu_reco_pt_fake, *h_tu_reco_pt_fake_split;
    TH1F *h_tu_gen_pt, *h_tu_gen_pt_split;

    // reco jet hists
    // For unfolding with TUnfold
    TUnfoldBinning *detector_tu_binning_LHA, *detector_distribution_LHA, *detector_distribution_underflow_LHA;
    TUnfoldBinning *generator_tu_binning_LHA, *generator_distribution_LHA, *generator_distribution_underflow_LHA;
    TH2F *h_tu_response_LHA, *h_tu_response_LHA_split; // response matrices. *_split histograms use half the MC to act as a check using independent sets of events
    TH1F *h_tu_reco_LHA, *h_tu_reco_LHA_fake, *h_tu_reco_LHA_split, *h_tu_reco_LHA_fake_split;  // detector quantities with reco binning
    TH1F *h_tu_reco_LHA_gen_binning, *h_tu_reco_LHA_fake_gen_binning, *h_tu_reco_LHA_gen_binning_split, *h_tu_reco_LHA_fake_gen_binning_split; // detector quantities with gen binning
    TH1F *h_tu_gen_LHA, *h_tu_gen_LHA_split; // truth histograms

    TUnfoldBinning *detector_tu_binning_puppiMultiplicity, *detector_distribution_puppiMultiplicity, *detector_distribution_underflow_puppiMultiplicity;
    TUnfoldBinning *generator_tu_binning_puppiMultiplicity, *generator_distribution_puppiMultiplicity, *generator_distribution_underflow_puppiMultiplicity;
    TH2F *h_tu_response_puppiMultiplicity, *h_tu_response_puppiMultiplicity_split;
    TH1F *h_tu_reco_puppiMultiplicity, *h_tu_reco_puppiMultiplicity_fake, *h_tu_reco_puppiMultiplicity_split, *h_tu_reco_puppiMultiplicity_fake_split;
    TH1F *h_tu_reco_puppiMultiplicity_gen_binning, *h_tu_reco_puppiMultiplicity_fake_gen_binning, *h_tu_reco_puppiMultiplicity_gen_binning_split, *h_tu_reco_puppiMultiplicity_fake_gen_binning_split;
    TH1F *h_tu_gen_puppiMultiplicity, *h_tu_gen_puppiMultiplicity_split;

    TUnfoldBinning *detector_tu_binning_pTD, *detector_distribution_pTD, *detector_distribution_underflow_pTD;
    TUnfoldBinning *generator_tu_binning_pTD, *generator_distribution_pTD, *generator_distribution_underflow_pTD;
    TH2F *h_tu_response_pTD, *h_tu_response_pTD_split;
    TH1F *h_tu_reco_pTD, *h_tu_reco_pTD_fake, *h_tu_reco_pTD_split, *h_tu_reco_pTD_fake_split;
    TH1F *h_tu_reco_pTD_gen_binning, *h_tu_reco_pTD_fake_gen_binning, *h_tu_reco_pTD_gen_binning_split, *h_tu_reco_pTD_fake_gen_binning_split;
    TH1F *h_tu_gen_pTD, *h_tu_gen_pTD_split;

    TUnfoldBinning *detector_tu_binning_thrust, *detector_distribution_thrust, *detector_distribution_underflow_thrust;
    TUnfoldBinning *generator_tu_binning_thrust, *generator_distribution_thrust, *generator_distribution_underflow_thrust;
    TH2F *h_tu_response_thrust, *h_tu_response_thrust_split;
    TH1F *h_tu_reco_thrust, *h_tu_reco_thrust_fake, *h_tu_reco_thrust_split, *h_tu_reco_thrust_fake_split;
    TH1F *h_tu_reco_thrust_gen_binning, *h_tu_reco_thrust_fake_gen_binning, *h_tu_reco_thrust_gen_binning_split, *h_tu_reco_thrust_fake_gen_binning_split;
    TH1F *h_tu_gen_thrust, *h_tu_gen_thrust_split;

    TUnfoldBinning *detector_tu_binning_width, *detector_distribution_width, *detector_distribution_underflow_width;
    TUnfoldBinning *generator_tu_binning_width, *generator_distribution_width, *generator_distribution_underflow_width;
    TH2F *h_tu_response_width, *h_tu_response_width_split;
    TH1F *h_tu_reco_width, *h_tu_reco_width_fake, *h_tu_reco_width_split, *h_tu_reco_width_fake_split;
    TH1F *h_tu_reco_width_gen_binning, *h_tu_reco_width_fake_gen_binning, *h_tu_reco_width_gen_binning_split, *h_tu_reco_width_fake_gen_binning_split;
    TH1F *h_tu_gen_width, *h_tu_gen_width_split;

    // charged-only versions:
    TUnfoldBinning *detector_tu_binning_LHA_charged, *detector_distribution_LHA_charged, *detector_distribution_underflow_LHA_charged;
    TUnfoldBinning *generator_tu_binning_LHA_charged, *generator_distribution_LHA_charged, *generator_distribution_underflow_LHA_charged;
    TH2F *h_tu_response_LHA_charged, *h_tu_response_LHA_charged_split;
    TH1F *h_tu_reco_LHA_charged, *h_tu_reco_LHA_charged_fake, *h_tu_reco_LHA_charged_split, *h_tu_reco_LHA_charged_fake_split;
    TH1F *h_tu_reco_LHA_charged_gen_binning, *h_tu_reco_LHA_charged_fake_gen_binning, *h_tu_reco_LHA_charged_gen_binning_split, *h_tu_reco_LHA_charged_fake_gen_binning_split;
    TH1F *h_tu_gen_LHA_charged, *h_tu_gen_LHA_charged_split;

    TUnfoldBinning *detector_tu_binning_puppiMultiplicity_charged, *detector_distribution_puppiMultiplicity_charged, *detector_distribution_underflow_puppiMultiplicity_charged;
    TUnfoldBinning *generator_tu_binning_puppiMultiplicity_charged, *generator_distribution_puppiMultiplicity_charged, *generator_distribution_underflow_puppiMultiplicity_charged;
    TH2F *h_tu_response_puppiMultiplicity_charged, *h_tu_response_puppiMultiplicity_charged_split;
    TH1F *h_tu_reco_puppiMultiplicity_charged, *h_tu_reco_puppiMultiplicity_charged_fake, *h_tu_reco_puppiMultiplicity_charged_split, *h_tu_reco_puppiMultiplicity_charged_fake_split;
    TH1F *h_tu_reco_puppiMultiplicity_charged_gen_binning, *h_tu_reco_puppiMultiplicity_charged_fake_gen_binning, *h_tu_reco_puppiMultiplicity_charged_gen_binning_split, *h_tu_reco_puppiMultiplicity_charged_fake_gen_binning_split;
    TH1F *h_tu_gen_puppiMultiplicity_charged, *h_tu_gen_puppiMultiplicity_charged_split;

    TUnfoldBinning *detector_tu_binning_pTD_charged, *detector_distribution_pTD_charged, *detector_distribution_underflow_pTD_charged;
    TUnfoldBinning *generator_tu_binning_pTD_charged, *generator_distribution_pTD_charged, *generator_distribution_underflow_pTD_charged;
    TH2F *h_tu_response_pTD_charged, *h_tu_response_pTD_charged_split;
    TH1F *h_tu_reco_pTD_charged, *h_tu_reco_pTD_charged_fake, *h_tu_reco_pTD_charged_split, *h_tu_reco_pTD_charged_fake_split;
    TH1F *h_tu_reco_pTD_charged_gen_binning, *h_tu_reco_pTD_charged_fake_gen_binning, *h_tu_reco_pTD_charged_gen_binning_split, *h_tu_reco_pTD_charged_fake_gen_binning_split;
    TH1F *h_tu_gen_pTD_charged, *h_tu_gen_pTD_charged_split;

    TUnfoldBinning *detector_tu_binning_thrust_charged, *detector_distribution_thrust_charged, *detector_distribution_underflow_thrust_charged;
    TUnfoldBinning *generator_tu_binning_thrust_charged, *generator_distribution_thrust_charged, *generator_distribution_underflow_thrust_charged;
    TH2F *h_tu_response_thrust_charged, *h_tu_response_thrust_charged_split;
    TH1F *h_tu_reco_thrust_charged, *h_tu_reco_thrust_charged_fake, *h_tu_reco_thrust_charged_split, *h_tu_reco_thrust_charged_fake_split;
    TH1F *h_tu_reco_thrust_charged_gen_binning, *h_tu_reco_thrust_charged_fake_gen_binning, *h_tu_reco_thrust_charged_gen_binning_split, *h_tu_reco_thrust_charged_fake_gen_binning_split;
    TH1F *h_tu_gen_thrust_charged, *h_tu_gen_thrust_charged_split;

    TUnfoldBinning *detector_tu_binning_width_charged, *detector_distribution_width_charged, *detector_distribution_underflow_width_charged;
    TUnfoldBinning *generator_tu_binning_width_charged, *generator_distribution_width_charged, *generator_distribution_underflow_width_charged;
    TH2F *h_tu_response_width_charged, *h_tu_response_width_charged_split;
    TH1F *h_tu_reco_width_charged, *h_tu_reco_width_charged_fake, *h_tu_reco_width_charged_split, *h_tu_reco_width_charged_fake_split;
    TH1F *h_tu_reco_width_charged_gen_binning, *h_tu_reco_width_charged_fake_gen_binning, *h_tu_reco_width_charged_gen_binning_split, *h_tu_reco_width_charged_fake_gen_binning_split;
    TH1F *h_tu_gen_width_charged, *h_tu_gen_width_charged_split;

    int useNJets_;
    bool doGroomed_;
    bool doHerwigReweighting;
    TH1F * reweightHist;

    uhh2::Event::Handle<std::vector<GenJetLambdaBundle> > genJetsLambda_handle;
    uhh2::Event::Handle<std::vector<JetLambdaBundle> > jetsLambda_handle;
    uhh2::Event::Handle<double> gen_weight_handle;
    uhh2::Event::Handle<bool> pass_reco_handle;
    uhh2::Event::Handle<bool> pass_gen_handle;
    bool is_mc_;

    TRandom3 rand_;
    bool doMCsplit_;
};

}

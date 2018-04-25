#include <iostream>
#include <memory>
#include <algorithm>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/QGAnalysis/include/QGAnalysisHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisDijetHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisTheoryHists.h"

using namespace std;
using namespace uhh2;

// Analysis module to calculate jet trigger efficiency

namespace uhh2examples {

const bool PRINTOUT = false;

// Easy way to refer to PDGIDs
enum PDGID {
    UNKNOWN = 0,
    DOWN_QUARK = 1,
    UP_QUARK = 2,
    STRANGE_QUARK = 3,
    CHARM_QUARK = 4,
    BOTTOM_QUARK = 5,
    TOP_QUARK = 6,
    ELECTRON = 11,
    MUON = 13,
    TAU = 15,
    GLUON = 21
};

/** \brief Basic analysis preselection
 *
 */
class QGAnalysisJetTrigEffModule: public AnalysisModule {
public:

    explicit QGAnalysisJetTrigEffModule(Context & ctx);
    virtual bool process(Event & event) override;

private:

    std::unique_ptr<CommonModules> common;

    std::unique_ptr<JetCorrector> jet_corrector_MC, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
    std::unique_ptr<GenericJetResolutionSmearer> jet_resolution_smearer;
    std::unique_ptr<JetCleaner> jet_cleaner;
    std::unique_ptr<JetElectronOverlapRemoval> jet_ele_cleaner;
    std::unique_ptr<JetMuonOverlapRemoval> jet_mu_cleaner;

    // Reco selections/hists
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, dijet_sel;
    std::unique_ptr<OrSelection> zplusjets_trigger_sel, dijet_trigger_sel;
    std::vector<TriggerSelection> dijet_trigger_sels;

    std::unique_ptr<QGJetTrigHists> trigHists;
    std::vector<QGJetTrigHists> jetTrigHists;
    std::vector<TriggerSelection> jet_trig_sels;

    bool is_mc;
    float jetRadius;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;

    std::unique_ptr<EventNumberSelection> event_sel;

    const std::vector<std::string> zpj_trigger_names = {
        "HLT_IsoMu24_v*",
        "HLT_IsoTkMu24_v*",
    };
    const std::vector<std::string> dj_trigger_names = {
        "HLT_PFJet40_v*",
        "HLT_PFJet60_v*",
        "HLT_PFJet80_v*",
        "HLT_PFJet140_v*",
        "HLT_PFJet200_v*",
        "HLT_PFJet260_v*",
        "HLT_PFJet320_v*",
        "HLT_PFJet400_v*",
        "HLT_PFJet450_v*",
        "HLT_PFJet500_v*"
        // "HLT_AK4PFJet30_v*",
        // "HLT_AK4PFJet50_v*",
        // "HLT_AK4PFJet80_v*",
        // "HLT_AK4PFJet100_v*",
        // "HLT_AK8PFJet40_v*",
        // "HLT_AK8PFJet60_v*",
        // "HLT_AK8PFJet80_v*",
        // "HLT_AK8PFJet140_v*",
        // "HLT_AK8PFJet200_v*",
        // "HLT_AK8PFJet260_v*",
        // "HLT_AK8PFJet320_v*",
        // "HLT_AK8PFJet400_v*",
        // "HLT_AK8PFJet450_v*",
        // "HLT_AK8PFJet500_v*",
    };
};


QGAnalysisJetTrigEffModule::QGAnalysisJetTrigEffModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.

    cout << "Running analysis module" << endl;

    is_mc = ctx.get("dataset_type") == "MC";
    // useGenPartonFlav = ctx.get("useGenPartonFlav") == "true";

    common.reset(new CommonModules());
    common->disable_mcpileupreweight();
    common->disable_jersmear();
    common->disable_jec(); // do it manually below
    common->change_pf_id(JetPFID::wp::WP_LOOSE);
    // common->change_pf_id(JetPFID::wp::WP_TIGHT_LEPVETO);
    common->set_muon_id(AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(26.0, 2.4), MuonIso(0.25)));
    common->set_electron_id(AndId<Electron>(ElectronID_Spring16_medium, PtEtaCut(20.0, 2.5)));
    common->switch_jetPtSorter(false);
    common->switch_jetlepcleaner(false);

    common->init(ctx);

    // Cleaning/JEC modules
    // Do manually and not in CommonModules to select correct cone size etc
    string jet_cone = ctx.get("JetCone", "AK4");
    string pu_removal = ctx.get("PURemoval", "CHS");
    if (pu_removal != "CHS" && pu_removal != "PUPPI") {
        throw runtime_error("Only PURemoval == CHS, PUPPI supported for now");
    }

    if (jet_cone.find("AK4") != string::npos)
        jetRadius = 0.4;
    else if (jet_cone.find("AK8") != string::npos)
        jetRadius = 0.8;
    else if (jet_cone.find("ca15") != string::npos)
        jetRadius = 1.5;
    else
        throw runtime_error("Cannot determine jetRadius in QGAnalysisTheoryHists");

    cout << "Running with jet cone: " << jet_cone << endl;
    cout << "Running with PUS: " << pu_removal << endl;

    if (is_mc) {
        std::vector<std::string> JEC_MC;
        std::string resolutionFilename;
        if (pu_removal == "CHS") {
            if (jet_cone == "AK4") {
                JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
                resolutionFilename = "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt";
            } else if (jet_cone == "AK8") {
                JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK8PFchs_MC;
                resolutionFilename = "Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt";  // actually doesn't matter, they're all the same
            }
        } else if (pu_removal == "PUPPI") {
            if (jet_cone == "AK4") {
                JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFPuppi_MC;
                resolutionFilename = "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt";
            } else if (jet_cone == "AK8") {
                JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK8PFPuppi_MC;
                resolutionFilename = "Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt";
            }
        }

        jet_corrector_MC.reset(new JetCorrector(ctx, JEC_MC));
        // jet_resolution_smearer.reset(new JetResolutionSmearer(ctx));
        jet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets", true, JERSmearing::SF_13TeV_2016_03Feb2017, resolutionFilename));
    } else {
        std::vector<std::string> JEC_BCD, JEC_EFearly, JEC_FlateG, JEC_H;
        if (pu_removal == "CHS") {
            if (jet_cone == "AK4") {
                JEC_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA;
                JEC_EFearly = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA;
                JEC_FlateG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA;
                JEC_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA;
            } else if (jet_cone == "AK8") {
                JEC_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFchs_DATA;
                JEC_EFearly = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFchs_DATA;
                JEC_FlateG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFchs_DATA;
                JEC_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFchs_DATA;
            }
        } else if (pu_removal == "PUPPI") {
            if (jet_cone == "AK4") {
                JEC_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFPuppi_DATA;
                JEC_EFearly = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFPuppi_DATA;
                JEC_FlateG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFPuppi_DATA;
                JEC_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFPuppi_DATA;
            } else if (jet_cone == "AK8") {
                JEC_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFPuppi_DATA;
                JEC_EFearly = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFPuppi_DATA;
                JEC_FlateG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFPuppi_DATA;
                JEC_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFPuppi_DATA;
            }
        }

        jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_BCD));
        jet_corrector_EFearly.reset(new JetCorrector(ctx, JEC_EFearly));
        jet_corrector_FlateG.reset(new JetCorrector(ctx, JEC_FlateG));
        jet_corrector_H.reset(new JetCorrector(ctx, JEC_H));
    }

    jet_cleaner.reset(new JetCleaner(ctx, PtEtaCut(30.0, 2.4)));
    jet_ele_cleaner.reset(new JetElectronOverlapRemoval(jetRadius));
    jet_mu_cleaner.reset(new JetMuonOverlapRemoval(jetRadius));

    // Event Selections
    njet_sel.reset(new NJetSelection(1));

    // Triggers for data
    zplusjets_trigger_sel.reset(new OrSelection());
    for (const auto & trig: zpj_trigger_names) {
        zplusjets_trigger_sel->add<TriggerSelection>(trig);
    }

    trigHists.reset(new QGJetTrigHists(ctx, "SingleMuRef", dj_trigger_names));

    for (uint i=0; i < (dj_trigger_names.size()-1); i++) {
        jet_trig_sels.push_back(TriggerSelection(dj_trigger_names[i]));
        jetTrigHists.push_back(QGJetTrigHists(ctx, dj_trigger_names[i]+"Ref", std::vector<std::string>{dj_trigger_names[i+1]}));
    }
}


bool QGAnalysisJetTrigEffModule::process(Event & event) {
    // if (!event_sel->passes(event)) return false;

    if (PRINTOUT) {cout << "-- Event: " << event.event << endl;}

    // This is the main procedure, called for each event.
    if (!common->process(event)) {return false;}


    // Do additional cleaning & JEC manually to allow custom JEC (common only does AK4)
    // - Apply JEC
    // - correct MET
    // - Smear jets if MC
   if (is_mc) {
        jet_corrector_MC->process(event);
        jet_corrector_MC->correct_met(event);
        jet_resolution_smearer->process(event);
    } else {
        if (event.run <= runnr_BCD) {
            jet_corrector_BCD->process(event);
            jet_corrector_BCD->correct_met(event);
        } else if (event.run < runnr_EFearly) { //< is correct, not <=
            jet_corrector_EFearly->process(event);
            jet_corrector_EFearly->correct_met(event);
        } else if (event.run <= runnr_FlateG) {
            jet_corrector_FlateG->process(event);
            jet_corrector_FlateG->correct_met(event);
        } else if (event.run > runnr_FlateG) {
            jet_corrector_H->process(event);
            jet_corrector_H->correct_met(event);
        } else {
            throw runtime_error("CommonModules.cxx: run number not covered by if-statements in process-routine.");
        }
    }

    // Do jet cleaning
    // jet_cleaner->process(event);
    jet_ele_cleaner->process(event);
    jet_mu_cleaner->process(event);

    // Resort by pT
    sort_by_pt(*event.jets);

    // RECO PART
    // if (!njet_sel->passes(event)) return false;

    if (zplusjets_trigger_sel->passes(event)) trigHists->fill(event);

    for (uint i=0; i<jet_trig_sels.size(); i++) {
        if (jet_trig_sels.at(i).passes(event))
            jetTrigHists.at(i).fill(event);
    }


    return true;
}


// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisJetTrigEffModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisJetTrigEffModule)

}

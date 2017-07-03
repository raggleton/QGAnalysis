#include <iostream>
#include <memory>
#include <algorithm>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/QGAnalysis/include/QGAnalysisHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisDijetHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/Utils.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

/** \brief Basic analysis preselection
 *
 */
class QGAnalysisModule: public AnalysisModule {
public:

    explicit QGAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

private:

    std::unique_ptr<CommonModules> common;

    std::unique_ptr<JetCorrector> jet_corrector_MC, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
    std::unique_ptr<JetLeptonCleaner> JLC_MC, JLC_BCD, JLC_EFearly, JLC_FlateG, JLC_H;
    std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;
    std::unique_ptr<JetCleaner> jet_cleaner;

    // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
    // to avoid memory leaks.
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, dijet_sel, first_jet_qflav_sel, first_jet_gflav_sel;
    std::unique_ptr<Hists> zplusjets_hists_presel, zplusjets_hists, zplusjets_qg_hists, zplusjets_hists_q, zplusjets_hists_g;
    std::unique_ptr<Hists> dijet_hists_presel, dijet_hists, dijet_qg_hists, dijet_hists_q, dijet_hists_g;

    std::unique_ptr<Hists> zplusjets_hists, zplusjets_qg_hists;
    std::unique_ptr<Hists> dijet_hists, dijet_qg_hists;

    bool is_mc;
    float jetRadius;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;
};


QGAnalysisModule::QGAnalysisModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.

    cout << "Running analysis module" << endl;

    is_mc = ctx.get("dataset_type") == "MC";

    common.reset(new CommonModules());
    common->disable_mcpileupreweight();
    common->disable_mclumiweight();
    common->disable_jersmear();
    common->change_pf_id(JetPFID::wp::WP_LOOSE);
    common->set_muon_id(MuonIDMedium_ICHEP());
    common->disable_jec();
    common->switch_jetPtSorter(false);
    common->switch_jetlepcleaner(false);

    common->init(ctx);

    // Cleaning/JEC modules
    // Do manually and not in CommonModules to select correct cone size etc
    string jet_cone = ctx.get("JetCone", "AK4");
    string pu_removal = ctx.get("PURemoval", "CHS");
    if (pu_removal != "CHS") {
        throw runtime_error("Only PURemoval == CHS supported for now");
    }

    if (jet_cone.find("AK4") != string::npos)
        jetRadius = 0.4;
    else if (jet_cone.find("AK8") != string::npos)
        jetRadius = 0.8;
    else if (jet_cone.find("ca15") != string::npos)
        jetRadius = 1.5;
    else
        throw runtime_error("Cannot determine jetRadius in QGAnalysisTheoryHists");

    if (is_mc) {
        std::vector<std::string> JEC_MC;
        if (pu_removal == "CHS") {
            if (jet_cone == "AK4") {
                JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
            } else if (jet_cone == "AK8") {
                JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK8PFchs_MC;
            } else {
                throw runtime_error("Unsupported JetCone value");
            }
        } else {
            throw runtime_error("Unsupported PURemoval value");
        }

        jet_corrector_MC.reset(new JetCorrector(ctx, JEC_MC));

        JLC_MC.reset(new JetLeptonCleaner(ctx, JEC_MC));
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
            } else {
                throw runtime_error("Unsupported JetCone value");
            }
        } else {
            throw runtime_error("Unsupported PURemoval value");
        }

        jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_BCD));
        jet_corrector_EFearly.reset(new JetCorrector(ctx, JEC_EFearly));
        jet_corrector_FlateG.reset(new JetCorrector(ctx, JEC_FlateG));
        jet_corrector_H.reset(new JetCorrector(ctx, JEC_H));

        JLC_BCD.reset(new JetLeptonCleaner(ctx, JEC_BCD));
        JLC_EFearly.reset(new JetLeptonCleaner(ctx, JEC_EFearly));
        JLC_FlateG.reset(new JetLeptonCleaner(ctx, JEC_FlateG));
        JLC_H.reset(new JetLeptonCleaner(ctx, JEC_H));
    }

    jet_cleaner.reset(new JetCleaner(ctx, PtEtaCut(30.0, 2.4)));

    // Event Selections
    njet_sel.reset(new NJetSelection(1));
    zplusjets_sel.reset(new ZplusJetsSelection());
    dijet_sel.reset(new DijetSelection());
    std::vector<int> q = {1, 2, 3};
    first_jet_qflav_sel.reset(new JetFlavourSelection(q, 0));
    std::vector<int> g = {21};
    first_jet_gflav_sel.reset(new JetFlavourSelection(g, 0));

    // Hists
    zplusjets_hists_presel.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel"));
    zplusjets_hists.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets"));
    zplusjets_hists_q.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_q"));
    zplusjets_hists_g.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_g"));
    zplusjets_qg_hists.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG", 1));

    dijet_hists_presel.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel"));
    dijet_hists.reset(new QGAnalysisDijetHists(ctx, "Dijet"));
    dijet_hists_q.reset(new QGAnalysisDijetHists(ctx, "Dijet_q"));
    dijet_hists_g.reset(new QGAnalysisDijetHists(ctx, "Dijet_g"));
    dijet_qg_hists.reset(new QGAnalysisHists(ctx, "Dijet_QG", 2));
}


bool QGAnalysisModule::process(Event & event) {
    // This is the main procedure, called for each event.
    if (!common->process(event)) return false;

    // Do additional cleaning & JEC manually to allow custom JEC (common only does AK4)
    // 1) jet-lepton cleaning
    // 2) Apply JEC
    if (is_mc) {
        JLC_MC->process(event);
        jet_corrector_MC->process(event);
    } else {
        if (event.run <= runnr_BCD) {
            JLC_BCD->process(event);
            jet_corrector_BCD->process(event);
        } else if (event.run < runnr_EFearly) { //< is correct, not <=
            JLC_EFearly->process(event);
            jet_corrector_EFearly->process(event);
        } else if (event.run <= runnr_FlateG) {
            JLC_FlateG->process(event);
            jet_corrector_FlateG->process(event);
        } else if (event.run > runnr_FlateG) {
            JLC_H->process(event);
            jet_corrector_H->process(event);
        } else {
            throw runtime_error("CommonModules.cxx: run number not covered by if-statements in process-routine.");
        }
    }

    // NB: Need to correct MET as well if you use it!
    // e.g. jet_corrector_MC->correct_met(event);

    // 3) Do jet cleaner
    jet_cleaner->process(event);

    // 4) Resort by pT
    sort_by_pt(*event.jets);


    // Preselection hists
    zplusjets_hists_presel->fill(event);
    dijet_hists_presel->fill(event);

    if (!njet_sel->passes(event)) return false;

    if (first_jet_gflav_sel->passes(event)) {
        zplusjets_hists_g->fill(event);
        dijet_hists_g->fill(event);
    } else if (first_jet_qflav_sel->passes(event)) {
        zplusjets_hists_q->fill(event);
        dijet_hists_q->fill(event);
    }

    bool zpj = zplusjets_sel->passes(event);
    if (zpj) {
        zplusjets_hists->fill(event);
        zplusjets_qg_hists->fill(event);
    }

    bool dj = dijet_sel->passes(event);
    if (dj) {
        dijet_hists->fill(event);
        dijet_qg_hists->fill(event);
    }

    if (zpj && dj) {
        cout << "Warning: event (runid, eventid) = ("  << event.run << ", " << event.event << ") passes both Z+jets and Dijet criteria" << endl;
    }

    return zpj || dj;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisModule)

}

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
#include "UHH2/QGAnalysis/include/QGAddModules.h"

using namespace std;
using namespace uhh2;

// Analysis module to calculate jet trigger efficiency

namespace uhh2examples {

const bool PRINTOUT = false;

/** \brief Basic analysis preselection
 *
 */
class QGAnalysisJetTrigEffModule: public AnalysisModule {
public:

    explicit QGAnalysisJetTrigEffModule(Context & ctx);
    virtual bool process(Event & event) override;

private:

    std::unique_ptr<GeneralEventSetup> common_setup;
    std::unique_ptr<RecoJetSetup> recojet_setup;

    // Reco selections/hists
    std::unique_ptr<Selection> njet_sel;
    std::unique_ptr<OrSelection> zplusjets_trigger_sel;
    std::unique_ptr<TriggerSelection> zb_trigger_sel;

    std::unique_ptr<QGJetTrigHists> singleMuTrigHists;
    std::vector<QGJetTrigHists> jetTrigHists;
    std::unique_ptr<QGJetTrigHists> zbTrigHists;
    std::vector<TriggerSelection> jet_trig_sels;

    bool is_mc;

    std::unique_ptr<EventNumberSelection> event_sel;

    const std::vector<std::string> zpj_trigger_names = {
        "HLT_IsoMu24_v*",
        "HLT_IsoTkMu24_v*",
    };
    const std::vector<std::string> zb_trigger_names = {
        "HLT_ZeroBias_v*",
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
        "HLT_PFJet500_v*",
        // "HLT_AK4PFJet30_v*",
        // "HLT_AK4PFJet50_v*",
        // "HLT_AK4PFJet80_v*",
        // "HLT_AK4PFJet100_v*",
        "HLT_AK8PFJet40_v*",
        "HLT_AK8PFJet60_v*",
        "HLT_AK8PFJet80_v*",
        "HLT_AK8PFJet140_v*",
        "HLT_AK8PFJet200_v*",
        "HLT_AK8PFJet260_v*",
        "HLT_AK8PFJet320_v*",
        "HLT_AK8PFJet400_v*",
        "HLT_AK8PFJet450_v*",
        "HLT_AK8PFJet500_v*",
    };
};


QGAnalysisJetTrigEffModule::QGAnalysisJetTrigEffModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.

    cout << "Running analysis module" << endl;

    is_mc = ctx.get("dataset_type") == "MC";

    // Cleaning/JEC modules
    // Do manually and not in CommonModules to select correct cone size etc
    string jet_cone = ctx.get("JetCone", "AK4");
    string pu_removal = ctx.get("PURemoval", "CHS");
    if (pu_removal != "CHS" && pu_removal != "PUPPI") {
        throw runtime_error("Only PURemoval == CHS, PUPPI supported for now");
    }
    float jet_radius = get_jet_radius(jet_cone);
    common_setup.reset(new GeneralEventSetup(ctx));
    float jet_pt_min = 10.;
    float jet_y_max = Cuts::jet_y_max; // allow both AK4 & AK8 to fall inside tracker, and keeps both consistent
    bool jetId = true;
    recojet_setup.reset(new RecoJetSetup(ctx, pu_removal, jet_cone, jet_radius, jet_pt_min, jet_y_max, jetId));
    cout << "Running with jet cone: " << jet_cone << endl;
    cout << "Running with PUS: " << pu_removal << endl;

    // Event Selections
    njet_sel.reset(new NJetSelection(1));

    // Triggers for data
    zplusjets_trigger_sel.reset(new OrSelection());
    for (const auto & trig: zpj_trigger_names) {
        zplusjets_trigger_sel->add<TriggerSelection>(trig);
    }
    singleMuTrigHists.reset(new QGJetTrigHists(ctx, "SingleMuRef", dj_trigger_names));

    zb_trigger_sel.reset(new TriggerSelection(zb_trigger_names[0]));
    zbTrigHists.reset(new QGJetTrigHists(ctx, "ZeroBiasRef", dj_trigger_names));

    for (uint i=0; i < (dj_trigger_names.size()-1); i++) {
        jet_trig_sels.push_back(TriggerSelection(dj_trigger_names[i]));
        jetTrigHists.push_back(QGJetTrigHists(ctx, dj_trigger_names[i]+"Ref", std::vector<std::string>{}));
    }
}


bool QGAnalysisJetTrigEffModule::process(Event & event) {
    // if (!event_sel->passes(event)) return false;

    // if (PRINTOUT) {cout << "-- Event: " << event.event << endl;}

    // This is the main procedure, called for each event.
    if (!common_setup->process(event)) {
        if (PRINTOUT) cout << " failed commonmodules" << endl;
        return false;
    }
    recojet_setup->process(event);

    if (PRINTOUT) {
        auto trigNames = event.get_current_triggernames();
        for (auto titr : trigNames) {
            // if (titr.find("HLT_") != std::string::npos) {
            if (titr.find("HLT_PFJet") != std::string::npos) {
                auto ti = event.get_trigger_index(titr);
                if (event.passes_trigger(ti))
                    cout << "    " << titr << " : " << event.passes_trigger(ti) << endl;
            }
        }
    }

    // RECO PART
    // if (!njet_sel->passes(event)) return false;

    if (zb_trigger_sel->passes(event)) zbTrigHists->fill(event);

    if (zplusjets_trigger_sel->passes(event)) singleMuTrigHists->fill(event);

    for (uint i=0; i<jet_trig_sels.size(); i++) {
        // Add check first to ensure trigger in event, otherwise throws
        auto ti = event.get_trigger_index(dj_trigger_names[i]);
        if (event.lookup_trigger_index(ti) && jet_trig_sels.at(i).passes(event)) {
            // if (i==4) cout << "Firing: " << dj_trigger_names[i] << endl;
            jetTrigHists.at(i).fill(event);
        }
    }

    return true;
}


// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisJetTrigEffModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisJetTrigEffModule)

}

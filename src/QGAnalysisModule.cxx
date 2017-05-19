#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/JetHists.h"

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

    std::unique_ptr<JetCleaner> jetcleaner;

    // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
    // to avoid memory leaks.
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, dijet_sel;

    std::unique_ptr<Hists> zplusjets_muon_hists, zplusjets_jet_hists, dijet_jet_hists;
};


QGAnalysisModule::QGAnalysisModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.

    cout << "Running analysis module" << endl;

    common.reset(new CommonModules());
    common->disable_mcpileupreweight();
    common->disable_mclumiweight();
    common->disable_jec();
    common->disable_jersmear();
    common->change_pf_id(JetPFID::wp::WP_LOOSE);
    // slightly tighter Jet selection on top
    common->set_jet_id(PtEtaCut(30.0, 2.4));
    common->set_muon_id(MuonIDMedium_ICHEP());
    common->switch_jetlepcleaner(true);
    common->init(ctx);

    // Selections
    njet_sel.reset(new NJetSelection(1));
    zplusjets_sel.reset(new ZplusJetsSelection());
    dijet_sel.reset(new DijetSelection());

    // Hists
    zplusjets_muon_hists.reset(new MuonHists(ctx, "ZPlusJets_Muon"));
    zplusjets_jet_hists.reset(new JetHists(ctx, "ZPlusJets_Jet"));

    dijet_jet_hists.reset(new JetHists(ctx, "Dijet_Jet"));
}


bool QGAnalysisModule::process(Event & event) {
    // This is the main procedure, called for each event.
    if (!common->process(event)) return false;

    if (!njet_sel->passes(event)) return false;

    bool zpj = zplusjets_sel->passes(event);
    if (zpj) {
        zplusjets_muon_hists->fill(event);
        zplusjets_jet_hists->fill(event);
    }

    bool dj = dijet_sel->passes(event);
    if (dj) {
        dijet_jet_hists->fill(event);
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

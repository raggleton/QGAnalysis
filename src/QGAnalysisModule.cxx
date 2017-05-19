#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"

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


};


QGAnalysisModule::QGAnalysisModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.

    cout << "Running preselection module" << endl;

    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }

    common.reset(new CommonModules());
    common->disable_mcpileupreweight();
    common->disable_mclumiweight();
    common->disable_jec();
    common->disable_jersmear();

    JetPFID loose_jet(JetPFID::wp::WP_LOOSE);
    common->set_jet_id(loose_jet);

    common->init(ctx);

    // slightly tighter Jet selection
    jetcleaner.reset(new JetCleaner(ctx, 30.0, 2.4));

    njet_sel.reset(new NJetSelection(1)); // see common/include/NSelections.h
    zplusjets_sel.reset(new ZplusJetsSelection());
    dijet_sel.reset(new DijetSelection());
}


bool QGAnalysisModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.

    // cout << "QGAnalysisModule: (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

    common->process(event);
    jetcleaner->process(event);

    bool njet_selection = njet_sel->passes(event);

    return njet_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisModule)

}

#include <iostream>
#include <memory>
#include <algorithm>

#include <boost/lexical_cast.hpp>

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
#include "UHH2/QGAnalysis/include/QGAnalysisPrinters.h"

using namespace std;
using namespace uhh2;

namespace DATASET {
    enum Name {
        SingleMu=0,
        JetHT,
        ZeroBias
    };
};

namespace uhh2examples {

const bool PRINTOUT = false;


/**
 * Analysis module for data datasets
 */
class QGAnalysisDataModule: public AnalysisModule {
public:

    explicit QGAnalysisDataModule(Context & ctx);
    virtual bool process(Event & event) override;
    DATASET::Name matchDatasetName(const std::string & name);

private:

    std::unique_ptr<GeneralEventSetup> common_setup;

    // Reco selections/hists
    std::unique_ptr<ZFinder> zFinder;
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, zplusjets_presel, dijet_sel;
    std::unique_ptr<Selection> zerobias_trigger_sel;
    std::unique_ptr<OrSelection> zplusjets_trigger_sel;
    std::unique_ptr<DataJetSelection> dijet_trigger_sel;
    std::vector<std::string> dj_trig_names;
    std::vector<float> dj_trig_prescales, dj_trig_thresholds;
    float jetht_zb_pt_boundary;

    std::unique_ptr<Hists> zplusjets_hists_presel, zplusjets_hists, zplusjets_qg_hists;

    std::unique_ptr<Hists> dijet_hists_presel, dijet_hists, dijet_qg_hists;

    std::unique_ptr<EventNumberSelection> event_sel;
    DATASET::Name dataset;

    std::string zLabel;

};


QGAnalysisDataModule::QGAnalysisDataModule(Context & ctx){
    cout << "Running analysis module" << endl;

    string datasetStr = ctx.get("Dataset", "");
    dataset = matchDatasetName(datasetStr);
    cout << "Running over " << datasetStr << endl;

    string jet_cone = ctx.get("JetCone", "AK4");
    string pu_removal = ctx.get("PURemoval", "CHS");
    if (pu_removal != "CHS" && pu_removal != "PUPPI") {
        throw runtime_error("Only PURemoval == CHS, PUPPI supported for now");
    }
    float jet_radius = get_jet_radius(jet_cone);

    cout << "Running with jet cone: " << jet_cone << endl;
    cout << "Running with PUS: " << pu_removal << endl;

    common_setup.reset(new GeneralEventSetup(ctx, pu_removal, jet_cone, jet_radius));

    // Event Selections
    njet_sel.reset(new NJetSelection(1));

    zLabel = "zMuonCand";
    zFinder.reset(new ZFinder(ctx, "muons", zLabel));

    // Z+JETS selection
    float mu1_pt = 20.;
    float mu2_pt = 20.;
    float mZ_window = 20.;
    float dphi_jet_z_min = 2.0;
    float second_jet_frac_max_zpj = 0.3;
    zplusjets_sel.reset(new ZplusJetsSelection(ctx, zLabel, mu1_pt, mu2_pt, mZ_window, dphi_jet_z_min, second_jet_frac_max_zpj));

    // Preselection for Z+J - only 2 muons to reco Z
    dphi_jet_z_min = 0.;
    second_jet_frac_max_zpj = 999.;
    zplusjets_presel.reset(new ZplusJetsSelection(ctx, zLabel, mu1_pt, mu2_pt, mZ_window, dphi_jet_z_min, second_jet_frac_max_zpj));

    // DIJET selection
    float dphi_min = 2.;
    float second_jet_frac_max_dj = 0.94;
    float third_jet_frac_max = 0.3;
    bool ss_eta = false;
    float deta = 12;
    float sumEta = 10.;
    dijet_sel.reset(new DijetSelection(dphi_min, second_jet_frac_max_dj, third_jet_frac_max, ss_eta, deta, sumEta));

    // Triggers for data
    zerobias_trigger_sel.reset(new TriggerSelection("HLT_ZeroBias_v*"));

    std::vector<std::string> zpj_trigger_names = {
        "HLT_IsoMu24_v*",
        "HLT_IsoTkMu24_v*",
    };
    zplusjets_trigger_sel.reset(new OrSelection());
    for (const auto & trig: zpj_trigger_names) {
        zplusjets_trigger_sel->add<TriggerSelection>(trig);
    }

    // TODO find some better structure to hold this data together?
    dj_trig_names = {
        "HLT_PFJet60_v*",
        "HLT_PFJet80_v*",
        "HLT_PFJet140_v*",
        "HLT_PFJet200_v*",
        "HLT_PFJet260_v*",
        "HLT_PFJet320_v*",
        "HLT_PFJet400_v*",
        "HLT_PFJet450_v*"
    };

    dj_trig_thresholds = {
        80,
        103,
        191,
        258,
        332,
        362,
        457,
        492
    };

    dj_trig_prescales = {
        49891.9453547,
        13120.4895678,
        1496.44452961,
        348.686346954,
        61.0210313345,
        20.446914767,
        3*2.38456,
        1.00010464076
    };


    std::vector<std::pair<float, float>> dj_trigger_bins;
    for (uint i=0; i<dj_trig_thresholds.size(); i++) {
        float lowerPt = dj_trig_thresholds.at(i);
        float upperPt = (i == dj_trig_thresholds.size()-1) ? 999999 : dj_trig_thresholds.at(i+1);
        dj_trigger_bins.push_back(std::make_pair(lowerPt, upperPt));
    }
    dijet_trigger_sel.reset(new DataJetSelection(dj_trig_names, dj_trigger_bins));

    jetht_zb_pt_boundary = dj_trig_thresholds.at(0);

    // Hists
    zplusjets_hists_presel.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel"));
    zplusjets_hists.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets"));
    zplusjets_qg_hists.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG", 1, "zplusjets"));

    dijet_hists_presel.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel"));
    dijet_hists.reset(new QGAnalysisDijetHists(ctx, "Dijet"));
    dijet_qg_hists.reset(new QGAnalysisHists(ctx, "Dijet_QG", 2, "dijet"));

    // event_sel.reset(new EventNumberSelection({111}));
}


bool QGAnalysisDataModule::process(Event & event) {
    // if (!event_sel->passes(event)) return false;

    if (PRINTOUT) {cout << "-- Event: " << event.event << endl; }

    if (PRINTOUT) printMuons(*event.muons, "Precleaning");
    if (PRINTOUT) printElectrons(*event.electrons, "Precleaning");
    if (PRINTOUT) printJets(*event.jets, "Precleaning");

    if (!njet_sel->passes(event)) return false;

    // Check trigger
    bool passTrigger(false);
    int dj_trig_ind(0);
    if (dataset == DATASET::SingleMu) {
        passTrigger = zplusjets_trigger_sel->passes(event);
    } else if (dataset == DATASET::JetHT) {
        if ((event.jets->at(0).pt() < jetht_zb_pt_boundary)) return false;
        passTrigger = dijet_trigger_sel->passes(event);
        // have to do dijet bit before any jet ID to figure out which jet fired trigger
        dj_trig_ind = dijet_trigger_sel->passIndex();
    } else if (dataset == DATASET::ZeroBias) {
        // to ensure we don't doublecount events
        // that are also in JetHT, we ensure that the first jet has pT < the first
        // bin in JetHT. (inside dijet_trigger_sel is an equivalent check for JetHT)
        // passTrigger = true;
        passTrigger = (zerobias_trigger_sel->passes(event)) && (event.jets->at(0).pt() < jetht_zb_pt_boundary);
    }
    if (!passTrigger) return false;

    if (!common_setup->process(event)) return false;

    if (PRINTOUT) printMuons(*event.muons);
    if (PRINTOUT) printElectrons(*event.electrons);

    if (!njet_sel->passes(event)) return false; //redo 1 jet cut again after cleaning etc

    if (PRINTOUT) printJets(*event.jets);

    // Selection & hists
    bool selected(false);


    if (dataset == DATASET::SingleMu) {
        if (!zFinder->process(event))
            return false;
        if (zplusjets_presel->passes(event)) {
            zplusjets_hists_presel->fill(event);
            selected = zplusjets_sel->passes(event);
            if (selected) {
                zplusjets_hists->fill(event);
                zplusjets_qg_hists->fill(event);
            }
        }
    } else if (dataset == DATASET::JetHT || dataset == DATASET::ZeroBias) {
        dijet_hists_presel->fill(event);
        if (event.jets->size() > 1) {
            selected = dijet_sel->passes(event);
            if (selected) {
                if (dataset == DATASET::JetHT) event.weight *= dj_trig_prescales.at(dj_trig_ind);
                dijet_hists->fill(event);
                dijet_qg_hists->fill(event);
            }
        }
    }

    return selected;
}


DATASET::Name QGAnalysisDataModule::matchDatasetName(const std::string & name) {
    if (name == "SingleMu") {
        return DATASET::SingleMu;
    } else if (name == "JetHT") {
        return DATASET::JetHT;
    } else if (name == "ZeroBias") {
        return DATASET::ZeroBias;
    } else {
        throw std::runtime_error("Cannot understand dataset with name " + name);
    }
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisDataModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisDataModule)

}

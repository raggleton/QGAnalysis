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

using namespace std;
using namespace uhh2;

namespace uhh2examples {

const bool PRINTOUT = false;


/** \brief Basic analysis preselection
 *
 */
class QGAnalysisPreselDataModule: public AnalysisModule {
public:

    explicit QGAnalysisPreselDataModule(Context & ctx);
    virtual bool process(Event & event) override;
    void printJets(const std::vector<Jet> & jets, const std::string & info="", Color::Code color=Color::FG_GREEN);
    void printMuons(const std::vector<Muon> & muons, const std::string & info="", Color::Code color=Color::FG_RED);
    void printElectrons(const std::vector<Electron> & electrons, const std::string & info="", Color::Code color=Color::FG_YELLOW);

private:

    std::unique_ptr<CommonModules> common;

    std::unique_ptr<DataJetMetCorrector> jet_met_corrector;
    std::unique_ptr<JetCleaner> jet_cleaner_loose, jet_cleaner_tight;
    std::unique_ptr<JetElectronOverlapRemoval> jet_ele_cleaner;
    std::unique_ptr<JetMuonOverlapRemoval> jet_mu_cleaner;

    // Reco selections/hists
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, zplusjets_presel, dijet_sel;
    std::unique_ptr<OrSelection> zplusjets_trigger_sel;
    std::unique_ptr<DataJetSelection> dijet_trigger_sel;
    std::vector<std::string> dj_trig_names;
    std::vector<float> dj_trig_prescales, dj_trig_thresholds;

    std::unique_ptr<Hists> zplusjets_hists_presel, zplusjets_hists, zplusjets_qg_hists;

    std::unique_ptr<Hists> dijet_hists_presel, dijet_hists, dijet_qg_hists;

    std::unique_ptr<Hists> dijet_hists_presel_highPt, dijet_hists_highPt, dijet_qg_hists_highPt;

    // for sweeping over PU
    std::vector<std::pair<int, int>> pu_bins = {
        std::make_pair(5, 15),
        std::make_pair(20, 25),
        std::make_pair(30, 40)
    };
    std::vector< std::unique_ptr<Selection> > sel_pu_binned;
    std::vector< std::unique_ptr<Hists> > zplusjets_qg_hists_pu_binned;
    std::vector< std::unique_ptr<Hists> > dijet_qg_hists_pu_binned;

    float jetRadius;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;

    const bool DO_PU_BINNED_HISTS = false;

    std::unique_ptr<EventNumberSelection> event_sel;

};


QGAnalysisPreselDataModule::QGAnalysisPreselDataModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.

    cout << "Running analysis module" << endl;

    common.reset(new CommonModules());
    common->disable_mcpileupreweight();
    common->disable_jersmear();
    common->disable_jec(); // do it manually below
    common->change_pf_id(JetPFID::wp::WP_LOOSE);
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

    jet_met_corrector.reset(new DataJetMetCorrector(ctx, pu_removal, jet_cone));
    
    jet_cleaner_loose.reset(new JetCleaner(ctx, PtEtaCut(30.0, 4.7)));
    jet_cleaner_tight.reset(new JetCleaner(ctx, PtEtaCut(30.0, 2.4)));
    jet_ele_cleaner.reset(new JetElectronOverlapRemoval(jetRadius));
    jet_mu_cleaner.reset(new JetMuonOverlapRemoval(jetRadius));

    // Event Selections
    njet_sel.reset(new NJetSelection(1));

    // Z+JETS selection
    float mu1_pt = 20.;
    float mu2_pt = 20.;
    float mZ_window = 20.;
    float dphi_jet_z_min = 2.0;
    float second_jet_frac_max_zpj = 0.3;
    zplusjets_sel.reset(new ZplusJetsSelection(mu1_pt, mu2_pt, mZ_window, dphi_jet_z_min, second_jet_frac_max_zpj));

    // Preselection for Z+J - only 2 muons to reco Z
    dphi_jet_z_min = 0.;
    second_jet_frac_max_zpj = 999.;
    zplusjets_presel.reset(new ZplusJetsSelection(mu1_pt, mu2_pt, mZ_window, dphi_jet_z_min, second_jet_frac_max_zpj));

    // DIJET selection
    float dphi_min = 2.;
    float second_jet_frac_max_dj = 0.94;
    float third_jet_frac_max = 0.3;
    bool ss_eta = false;
    float deta = 12;
    float sumEta = 10.;
    dijet_sel.reset(new DijetSelection(dphi_min, second_jet_frac_max_dj, third_jet_frac_max, ss_eta, deta, sumEta));

    // Triggers for data
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

    // Hists
    zplusjets_hists_presel.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel"));
    // preselection hists, if jet is quark, or gluon
    zplusjets_hists.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets"));
    zplusjets_qg_hists.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG", 1, "zplusjets"));

    dijet_hists_presel.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel"));
    dijet_hists.reset(new QGAnalysisDijetHists(ctx, "Dijet"));
    dijet_qg_hists.reset(new QGAnalysisHists(ctx, "Dijet_QG", 2, "dijet"));

    dijet_hists_presel_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_highPt"));
    dijet_hists_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_highPt"));
    dijet_qg_hists_highPt.reset(new QGAnalysisHists(ctx, "Dijet_QG_highPt", 2, "dijet"));

    if (DO_PU_BINNED_HISTS) {
        for (auto puBin : pu_bins) {
            std::unique_ptr<Selection> pu_sel(new NPVSelection(puBin.first, puBin.second));
            sel_pu_binned.push_back(std::move(pu_sel));
            std::unique_ptr<QGAnalysisHists> zpj(new QGAnalysisHists(ctx, TString::Format("ZPlusJets_QG_PU_%d_to_%d", puBin.first, puBin.second).Data(), 1, "zplusjets"));
            zplusjets_qg_hists_pu_binned.push_back(std::move(zpj));
            std::unique_ptr<QGAnalysisHists> dj(new QGAnalysisHists(ctx, TString::Format("Dijet_QG_PU_%d_to_%d", puBin.first, puBin.second).Data(), 2, "dijet"));
            dijet_qg_hists_pu_binned.push_back(std::move(dj));
        }
    }

    // event_sel.reset(new EventNumberSelection({111}));
}


bool QGAnalysisPreselDataModule::process(Event & event) {
    // if (!event_sel->passes(event)) return false;

    if (PRINTOUT) {cout << "-- Event: " << event.event << endl; }

    printMuons(*event.muons, "Precleaning");
    printElectrons(*event.electrons, "Precleaning");
    printJets(*event.jets, "Precleaning");

    // Gen-level HT cut if necessary

    // This is the main procedure, called for each event.
    if (!common->process(event)) {return false;}

    printMuons(*event.muons);
    printElectrons(*event.electrons);

    bool pass_zpj_trig = zplusjets_trigger_sel->passes(event);
    // have to do dijet bit before any jet ID to figure out which jet fired trigger
    bool pass_dj_trig = dijet_trigger_sel->passes(event);
    int dj_trig_ind = dijet_trigger_sel->passIndex();

    if (!(pass_zpj_trig || pass_dj_trig)) return false;

    // Do jet cleaning/correcting
    jet_met_corrector->process(event);

    // jet_cleaner_loose->process(event);
    jet_cleaner_tight->process(event);
    jet_ele_cleaner->process(event);
    jet_mu_cleaner->process(event);

    // Resort by pT
    sort_by_pt(*event.jets);

    // RECO PART
    if (!njet_sel->passes(event)) return false;

    printJets(*event.jets);

    // Preselection hists
    zplusjets_hists_presel->fill(event);
    dijet_hists_presel->fill(event);

    // Full selection & hists
    bool zpj(false), dj(false), dj_highPt(false);

    // if (zplusjets_presel->passes(event)) {
    //     zpj = zplusjets_sel->passes(event);
    //     if (zpj && pass_zpj_trig) {
    //         zplusjets_hists->fill(event);
    //         zplusjets_qg_hists->fill(event);
    //     }
    // }

    if (event.jets->size() > 1) {
        dj = dijet_sel->passes(event);
        if (dj && pass_dj_trig) {
            event.weight *= dj_trig_prescales.at(dj_trig_ind);
            dijet_hists->fill(event);
            dijet_qg_hists->fill(event);
        }
    }

    // do pu-binned hists
    if (DO_PU_BINNED_HISTS) {
        for (uint i=0; i<sel_pu_binned.size(); i++) {
            if (sel_pu_binned.at(i)->passes(event)) {
                if (zpj) zplusjets_qg_hists_pu_binned.at(i)->fill(event);
                if (dj) dijet_qg_hists_pu_binned.at(i)->fill(event);
            }
        }
    }

    // do high pt jet version - both jets must pass much higher pt threshold
    // don't need to do a Z+jets version as only care about leading jet.
    // float ptCut = 500;
    // if ((event.jets->size() > 1) && (event.jets->at(0).pt() > ptCut) && (event.jets->at(1).pt() > ptCut)) {
    //     dijet_hists_presel_highPt->fill(event);

    //     dj_highPt = dijet_sel->passes(event);
    //     if (dj_highPt) { // flav2 only sensible if passed pt cut
    //         dijet_hists_highPt->fill(event);
    //         dijet_qg_hists_highPt->fill(event);
    //     }
    // }

    if (zpj && dj) {
        cout << "Warning: event (runid, eventid) = ("  << event.run << ", " << event.event << ") passes both Z+jets and Dijet criteria" << endl;
    }

    return zpj || dj;
    // return zpj || dj || dj_highPt;
}

void QGAnalysisPreselDataModule::printJets(const std::vector<Jet> & jets, const std::string & label, Color::Code color) {
    if (!PRINTOUT) return;
    for (auto & itr: jets) {
        cout << color << "jet";
        if (label != "") {
            cout << " [" << label << "]";
        }
        cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << endl;
    }
}

void QGAnalysisPreselDataModule::printMuons(const std::vector<Muon> & muons, const std::string & label, Color::Code color) {
    if (!PRINTOUT) return;
    for (auto & itr: muons) {
        cout << color << "muon";
        if (label != "") {
            cout << " [" << label << "]";
        }
        cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << endl;
    }
}

void QGAnalysisPreselDataModule::printElectrons(const std::vector<Electron> & electrons, const std::string & label, Color::Code color) {
    if (!PRINTOUT) return;
    for (auto & itr: electrons) {
        cout << color << "electron";
        if (label != "") {
            cout << " [" << label << "]";
        }
        cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << endl;
    }
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisPreselDataModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisPreselDataModule)

}

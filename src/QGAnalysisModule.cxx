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
#include "UHH2/QGAnalysis/include/QGAnalysisTheoryHists.h"
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
    std::vector<GenJetWithParts> getGenJets(std::vector<GenJetWithParts> * jets, std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=1.5, float lepton_overlap_dr=0.2);
    std::vector<GenParticle> getGenMuons(std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=2.5);
    std::vector<Jet> getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax=0.8);
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

    std::unique_ptr<Selection> zplusjets_theory_sel, dijet_theory_sel;
    std::unique_ptr<Hists> zplusjets_hists_theory, dijet_hists_theory;

    // for sweeping over ptMin
    std::vector<float> theory_pt_bins = {50, 100, 200, 400, 800};
    std::vector< std::unique_ptr<Selection> > zplusjets_theory_sel_pt_binned;
    std::vector< std::unique_ptr<Selection> > dijet_theory_sel_pt_binned;
    std::vector< std::unique_ptr<Hists> > zplusjets_hists_theory_pt_binned;
    std::vector< std::unique_ptr<Hists> > dijet_hists_theory_pt_binned;

    Event::Handle<std::vector<GenJetWithParts>> genjets_handle;
    Event::Handle<std::vector<GenParticle>> genmuons_handle;

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
    common->disable_jersmear();
    common->disable_jec(); // do it manually below
    common->change_pf_id(JetPFID::wp::WP_LOOSE);
    common->set_muon_id(MuonIDMedium_ICHEP());
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
        jet_resolution_smearer.reset(new JetResolutionSmearer(ctx));
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

    genjets_handle = ctx.declare_event_output< std::vector<GenJetWithParts> > ("GoodGenJets");
    genmuons_handle = ctx.declare_event_output< std::vector<GenParticle> > ("GoodGenMuons");

    // Event Selections
    njet_sel.reset(new NJetSelection(1));
    zplusjets_sel.reset(new ZplusJetsSelection());
    dijet_sel.reset(new DijetSelection());
    std::vector<int> q = {1, 2, 3};
    first_jet_qflav_sel.reset(new JetFlavourSelection(q, 0));
    std::vector<int> g = {21};
    first_jet_gflav_sel.reset(new JetFlavourSelection(g, 0));

    zplusjets_theory_sel.reset(new ZplusJetsTheorySelection(ctx));
    dijet_theory_sel.reset(new DijetTheorySelection(ctx));

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

    zplusjets_hists_theory.reset(new QGAnalysisTheoryHists(ctx, "ZPlusJets_genjet", 1));
    dijet_hists_theory.reset(new QGAnalysisTheoryHists(ctx, "Dijet_genjet", 2));

    for (auto pt : theory_pt_bins) {
        std::unique_ptr<Selection> A(new ZplusJetsTheorySelection(ctx, pt));
        zplusjets_theory_sel_pt_binned.push_back(std::move(A));
        std::unique_ptr<Selection> B(new DijetTheorySelection(ctx, pt));
        dijet_theory_sel_pt_binned.push_back(std::move(B));

        std::unique_ptr<QGAnalysisTheoryHists> a(new QGAnalysisTheoryHists(ctx, TString::Format("ZPlusJets_genjet_ptMin_%d", int(pt)).Data(), 1));
        zplusjets_hists_theory_pt_binned.push_back(std::move(a));
        std::unique_ptr<QGAnalysisTheoryHists> b(new QGAnalysisTheoryHists(ctx, TString::Format("Dijet_genjet_ptMin_%d", int(pt)).Data(), 2));
        dijet_hists_theory_pt_binned.push_back(std::move(b));
    }
}


bool QGAnalysisModule::process(Event & event) {
    // This is the main procedure, called for each event.
    if (!common->process(event)) {return false;}

    // Do additional cleaning & JEC manually to allow custom JEC (common only does AK4)
    // 1) jet-lepton cleaning
    // 2) Apply JEC
    // 3) correct MET
    // 4) Smear jets if MC
    if (is_mc) {
        JLC_MC->process(event);
        jet_corrector_MC->process(event);
        jet_corrector_MC->correct_met(event);
        jet_resolution_smearer->process(event);
    } else {
        if (event.run <= runnr_BCD) {
            JLC_BCD->process(event);
            jet_corrector_BCD->process(event);
            jet_corrector_BCD->correct_met(event);
        } else if (event.run < runnr_EFearly) { //< is correct, not <=
            JLC_EFearly->process(event);
            jet_corrector_EFearly->process(event);
            jet_corrector_EFearly->correct_met(event);
        } else if (event.run <= runnr_FlateG) {
            JLC_FlateG->process(event);
            jet_corrector_FlateG->process(event);
            jet_corrector_FlateG->correct_met(event);
        } else if (event.run > runnr_FlateG) {
            JLC_H->process(event);
            jet_corrector_H->process(event);
            jet_corrector_H->correct_met(event);
        } else {
            throw runtime_error("CommonModules.cxx: run number not covered by if-statements in process-routine.");
        }
    }

    // 3) Do jet cleaner
    jet_cleaner->process(event);

    // 4) Resort by pT
    sort_by_pt(*event.jets);

   // THEORY PART
    std::vector<GenJetWithParts> goodGenJets = getGenJets(event.genjets, event.genparticles, 5., 1.5, jetRadius);
    if (goodGenJets.size() < 1) {
        return false;
    }
    event.set(genjets_handle, std::move(goodGenJets));

    std::vector<GenParticle> goodGenMuons = getGenMuons(event.genparticles, 5., 2.5);
    event.set(genmuons_handle, std::move(goodGenMuons));

    bool zpj_th = zplusjets_theory_sel->passes(event);
    bool dj_th = dijet_theory_sel->passes(event);

    if (zpj_th && dj_th) {
        cout << "Warning: event (runid, eventid) = ("  << event.run << ", " << event.event << ") passes both Z+jets and Dijet theory criteria" << endl;
    }

    if (zpj_th) {
        zplusjets_hists_theory->fill(event);
    }

    if (dj_th) {
        dijet_hists_theory->fill(event);
    }

    // ptMin binned hists
    for (uint i=0; i < theory_pt_bins.size(); i++) {
        if (zplusjets_theory_sel_pt_binned.at(i)->passes(event)) {
            zplusjets_hists_theory_pt_binned.at(i)->fill(event);
        }
        if (dijet_theory_sel_pt_binned.at(i)->passes(event)) {
            dijet_hists_theory_pt_binned.at(i)->fill(event);
        }
    }

    // RECO PART
    // Only consider jets that have a matching GenJet - this is to avoid dijet events solely from PU,
    // which cause havoc if there is an event weight dervied from a significantly lower GenJet pT.
    if (is_mc) {
        std::vector<Jet> goodJets = getMatchedJets(event.jets, &event.get(genjets_handle), 0.8);
        std::swap(goodJets, *event.jets);
    }

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



    // For checking genparticle/jet assignments:
    // std::cout << "JETS" << std::endl;
    // for (const auto & itr : *event.jets) {
    //     std::cout << itr.eta() << " : " << itr.phi() << " : " << itr.genPartonFlavor() << std::endl;
    // }
    // std::cout << "GenPARTICLES" << std::endl;
    // for (const auto & itr : *event.genparticles) {
    //     if (abs(itr.status()) > 1)
    //     std::cout << itr.pdgId() << " : " << itr.status() << " : " << itr.eta() << " : " << itr.phi() << std::endl;
    // }

    return zpj || dj || zpj_th || dj_th;
}


/**
 * Get GenJets, ignoring those that are basically a lepton, and have some minimum pt and maximum eta.
 */
std::vector<GenJetWithParts> QGAnalysisModule::getGenJets(std::vector<GenJetWithParts> * jets, std::vector<GenParticle> * genparticles, float pt_min, float eta_max, float lepton_overlap_dr) {
    std::vector<GenJetWithParts> genjets;
    for (const auto itr : *jets) {
        bool found = (std::find(genjets.begin(), genjets.end(), itr) != genjets.end());
        // avoid jets that are just leptons + a few spurious gluons
        bool leptonOverlap = false;
        if (genparticles != nullptr) {
            for (const auto & ptr : *genparticles) {
                leptonOverlap = leptonOverlap || (((abs(ptr.pdgId()) == 13) || (abs(ptr.pdgId()) == 11) || (abs(ptr.pdgId()) == 15)) && (deltaR(ptr.v4(), itr.v4()) < lepton_overlap_dr));
            }
        }
        if ((itr.pt() > pt_min) && (fabs(itr.eta()) < eta_max) && !found && !leptonOverlap) genjets.push_back(itr);
    }
    sort_by_pt(genjets);
    return genjets;
}


/**
 * Select gen muons from all genparticles, that have some minimum pt and maximum eta
 */
std::vector<GenParticle> QGAnalysisModule::getGenMuons(std::vector<GenParticle> * genparticles, float pt_min, float eta_max) {
    std::vector<GenParticle> muons;
    for(const auto itr : *genparticles) {
        bool found = (std::find(muons.begin(), muons.end(), itr) != muons.end());
        if ((abs(itr.pdgId()) == 13) && (itr.status() == 1) && (itr.pt() > pt_min) && (fabs(itr.eta()) < eta_max) && !found) {
            muons.push_back(itr);
        }
    }
    sort_by_pt(muons);
    return muons;
}

/**
 * Only select reco jets that have a matching GenJet within some DR
 */
std::vector<Jet> QGAnalysisModule::getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax) {
    std::vector<Jet> goodJets;
    for (const auto jtr: *jets) {
        for (const auto genjtr: * genjets) {
            if (deltaR(jtr, genjtr) < drMax) {
                goodJets.push_back(jtr);
                break;
            }
        }
    }
    return goodJets;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisModule)

}

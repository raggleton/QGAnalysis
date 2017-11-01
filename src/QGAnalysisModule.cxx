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

namespace Color {
    enum Code {
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_YELLOW   = 33,
        FG_BLUE     = 34,
        FG_MAGENTA  = 35,
        FG_CYAN     = 36,
        FG_DEFAULT  = 39,
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };
    std::ostream& operator<<(std::ostream& os, Code code) {
        return os << "\033[" << static_cast<int>(code) << "m";
    }
}

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
class QGAnalysisModule: public AnalysisModule {
public:

    explicit QGAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    std::vector<GenJetWithParts> getGenJets(std::vector<GenJetWithParts> * jets, std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=1.5, float lepton_overlap_dr=0.2);
    std::vector<GenParticle> getGenMuons(std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=2.5);
    std::vector<Jet> getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax=0.8, bool uniqueMatch=true);
    void printGenParticles(const std::vector<GenParticle> & gps, const std::string & info="", Color::Code color=Color::FG_DEFAULT);
    void printGenJets(const std::vector<GenJetWithParts> & gps, const std::string & info="", Color::Code color=Color::FG_BLUE);
    void printJets(const std::vector<Jet> & jets, const std::string & info="", Color::Code color=Color::FG_GREEN);

private:

    std::unique_ptr<CommonModules> common;

    std::unique_ptr<JetCorrector> jet_corrector_MC, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
    std::unique_ptr<JetLeptonCleaner> JLC_MC, JLC_BCD, JLC_EFearly, JLC_FlateG, JLC_H;
    // std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;
    std::unique_ptr<GenericJetResolutionSmearer> jet_resolution_smearer;
    std::unique_ptr<JetCleaner> jet_cleaner;

    // Reco selections/hists
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, dijet_sel;
    std::unique_ptr<Hists> zplusjets_hists_presel, zplusjets_hists, zplusjets_qg_hists, zplusjets_hists_q, zplusjets_hists_g;
    std::unique_ptr<Hists> zplusjets_hists_presel_q, zplusjets_hists_presel_g, zplusjets_hists_presel_unknown;
    std::unique_ptr<Hists> dijet_hists_presel, dijet_hists, dijet_qg_hists, dijet_hists_q, dijet_hists_g;
    std::unique_ptr<Hists> dijet_hists_presel_gg, dijet_hists_presel_qg, dijet_hists_presel_qq;
    std::unique_ptr<Hists> dijet_hists_presel_q_unknown, dijet_hists_presel_g_unknown, dijet_hists_presel_unknown_unknown, dijet_hists_presel_unknown_q, dijet_hists_presel_unknown_g;

    std::unique_ptr<Hists> zplusjets_qg_hists_u, zplusjets_qg_hists_d, zplusjets_qg_hists_s;
    std::unique_ptr<Hists> dijet_qg_hists_u, dijet_qg_hists_d, dijet_qg_hists_s, dijet_qg_hists_c;

    // for sweeping over PU
    std::vector<std::pair<int, int>> pu_bins = {
        std::make_pair(5, 15),
        std::make_pair(20, 25),
        std::make_pair(30, 40)
    };
    std::vector< std::unique_ptr<Selection> > sel_pu_binned;
    std::vector< std::unique_ptr<Hists> > zplusjets_qg_hists_pu_binned;
    std::vector< std::unique_ptr<Hists> > dijet_qg_hists_pu_binned;

    // Theory selections/hists
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

        JLC_BCD.reset(new JetLeptonCleaner(ctx, JEC_BCD));
        JLC_EFearly.reset(new JetLeptonCleaner(ctx, JEC_EFearly));
        JLC_FlateG.reset(new JetLeptonCleaner(ctx, JEC_FlateG));
        JLC_H.reset(new JetLeptonCleaner(ctx, JEC_H));
    }

    jet_cleaner.reset(new JetCleaner(ctx, PtEtaCut(30.0, 2.5)));

    genjets_handle = ctx.declare_event_output< std::vector<GenJetWithParts> > ("GoodGenJets");
    genmuons_handle = ctx.declare_event_output< std::vector<GenParticle> > ("GoodGenMuons");

    // Event Selections
    njet_sel.reset(new NJetSelection(1));

    zplusjets_sel.reset(new ZplusJetsSelection());
    dijet_sel.reset(new DijetSelection());

    zplusjets_theory_sel.reset(new ZplusJetsTheorySelection(ctx));
    dijet_theory_sel.reset(new DijetTheorySelection(ctx));

    // Hists
    zplusjets_hists_presel.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel"));
    // preselection hists, if jet is quark, or gluon
    zplusjets_hists_presel_q.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_q"));
    zplusjets_hists_presel_g.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_g"));
    zplusjets_hists_presel_unknown.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_unknown"));
    zplusjets_hists.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets"));
    zplusjets_qg_hists.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG", 1, "zplusjets"));

    dijet_hists_presel.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel"));
    // preselection hiss, if both gluon jets, one gluon, or both quark, or one or both unknown
    dijet_hists_presel_gg.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_gg"));
    dijet_hists_presel_qg.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_qg"));
    dijet_hists_presel_qq.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_qq"));
    dijet_hists_presel_unknown_q.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_q"));
    dijet_hists_presel_unknown_g.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_g"));
    dijet_hists_presel_unknown_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_unknown"));
    dijet_hists_presel_q_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_q_unknown"));
    dijet_hists_presel_g_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_g_unknown"));
    dijet_hists.reset(new QGAnalysisDijetHists(ctx, "Dijet"));
    dijet_qg_hists.reset(new QGAnalysisHists(ctx, "Dijet_QG", 2, "dijet"));

    zplusjets_hists_theory.reset(new QGAnalysisTheoryHists(ctx, "ZPlusJets_genjet", 1, "zplusjets"));
    dijet_hists_theory.reset(new QGAnalysisTheoryHists(ctx, "Dijet_genjet", 2, "dijet"));

    for (auto pt : theory_pt_bins) {
        std::unique_ptr<Selection> A(new ZplusJetsTheorySelection(ctx, pt));
        zplusjets_theory_sel_pt_binned.push_back(std::move(A));
        std::unique_ptr<Selection> B(new DijetTheorySelection(ctx, pt));
        dijet_theory_sel_pt_binned.push_back(std::move(B));

        std::unique_ptr<QGAnalysisTheoryHists> a(new QGAnalysisTheoryHists(ctx, TString::Format("ZPlusJets_genjet_ptMin_%d", int(pt)).Data(), 1, "zplusjets"));
        zplusjets_hists_theory_pt_binned.push_back(std::move(a));
        std::unique_ptr<QGAnalysisTheoryHists> b(new QGAnalysisTheoryHists(ctx, TString::Format("Dijet_genjet_ptMin_%d", int(pt)).Data(), 2, "dijet"));
        dijet_hists_theory_pt_binned.push_back(std::move(b));
    }

    for (auto puBin : pu_bins) {
        std::unique_ptr<Selection> pu_sel(new NPVSelection(puBin.first, puBin.second));
        sel_pu_binned.push_back(std::move(pu_sel));
        std::unique_ptr<QGAnalysisHists> zpj(new QGAnalysisHists(ctx, TString::Format("ZPlusJets_QG_PU_%d_to_%d", puBin.first, puBin.second).Data(), 1, "zplusjets"));
        zplusjets_qg_hists_pu_binned.push_back(std::move(zpj));
        std::unique_ptr<QGAnalysisHists> dj(new QGAnalysisHists(ctx, TString::Format("Dijet_QG_PU_%d_to_%d", puBin.first, puBin.second).Data(), 2, "dijet"));
        dijet_qg_hists_pu_binned.push_back(std::move(dj));
    }
}


bool QGAnalysisModule::process(Event & event) {
    if (PRINTOUT) {cout << "-- Event: " << event.event << endl;}

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
    // printGenParticles(*event.genparticles);

    std::vector<GenJetWithParts> goodGenJets = getGenJets(event.genjets, event.genparticles, 5., 5., jetRadius);
    if (goodGenJets.size() == 0) return false;
    if (event.jets->size() == 0) return false;

    // Check event weight is sensible based on pthat
    if (event.genInfo->binningValues().size() > 0) {
        double ptHat = event.genInfo->binningValues().at(0); // yes this is correct. no idea why
        if (goodGenJets[0].pt() / ptHat > 2) return false;
        // Again, need to do this as sometimes reco dodgy but gen ok can end up with dodgy weight
        if (event.jets->at(0).pt() / ptHat > 2) return false;
    }

    event.set(genjets_handle, std::move(goodGenJets));

    // Determine if weight is good - weight calculated based on leading jet
    // Bad if it's a PU jet!
    if (is_mc) {
        bool goodEvent = false;
        for (const auto genjtr: event.get(genjets_handle)) {
            if (deltaR(event.jets->at(0), genjtr) < jetRadius){
                goodEvent = true;
                break;
            }
        }
        if (!goodEvent) return false;
    }

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
        std::vector<Jet> goodJets = getMatchedJets(event.jets, &event.get(genjets_handle), jetRadius/2.);
        std::swap(goodJets, *event.jets);
    }

    // Preselection hists
    zplusjets_hists_presel->fill(event);
    dijet_hists_presel->fill(event);

    // flav-specific preselection hists, useful for optimising selection
    uint flav1 = event.jets->at(0).genPartonFlavor();
    uint flav2(99999999);
    if (event.jets->size() > 1) {
        flav2 = event.jets->at(1).genPartonFlavor();
    }
    if (flav1 == PDGID::GLUON) {
        zplusjets_hists_presel_g->fill(event);
        if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
            dijet_hists_presel_qg->fill(event);
        } else if (flav2 == PDGID::GLUON) {
            dijet_hists_presel_gg->fill(event);
        } else if (flav2 == PDGID::UNKNOWN) {
            dijet_hists_presel_g_unknown->fill(event);
        }
    } else if (flav1 > PDGID::UNKNOWN && flav1 < PDGID::CHARM_QUARK) {
        zplusjets_hists_presel_q->fill(event);
        if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
            dijet_hists_presel_qq->fill(event);
        } else if (flav2 == PDGID::GLUON) {
            dijet_hists_presel_qg->fill(event);
        } else if (flav2 == PDGID::UNKNOWN) {
            dijet_hists_presel_q_unknown->fill(event);
        }
    } else if (flav1 == PDGID::UNKNOWN) {
        zplusjets_hists_presel_unknown->fill(event);
        if (flav2 == PDGID::GLUON) {
            dijet_hists_presel_unknown_g->fill(event);
        } else if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
            dijet_hists_presel_unknown_q->fill(event);
        } else if (flav2 == PDGID::UNKNOWN) {
            dijet_hists_presel_unknown_unknown->fill(event);
        }
    }

    // Do reco selection and hist filling
    if (!njet_sel->passes(event)) return false;

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

    // do pu-binned hists
    for (uint i=0; i<sel_pu_binned.size(); i++) {
        if (sel_pu_binned.at(i)->passes(event)) {
            if (zpj) zplusjets_qg_hists_pu_binned.at(i)->fill(event);
            if (dj) dijet_qg_hists_pu_binned.at(i)->fill(event);
        }
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
                leptonOverlap = leptonOverlap || (((abs(ptr.pdgId()) == PDGID::MUON) || (abs(ptr.pdgId()) == PDGID::ELECTRON) || (abs(ptr.pdgId()) == PDGID::TAU)) && (deltaR(ptr.v4(), itr.v4()) < lepton_overlap_dr));
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
    // Do in reverse order to pick up most evolved muons first
    for (auto itr = genparticles->rbegin(); itr != genparticles->rend(); ++itr){
        // We check to see if we already have a very similar, bu tno exact, muon
        // since the MC "evolves" the particle and slightly changes pt/eta/phi
        bool alreadyFound = std::any_of(muons.begin(), muons.end(), [&itr] (const GenParticle & mtr) { return deltaR(*itr, mtr) < 0.05 && itr->charge() == mtr.charge(); });
        if ((abs(itr->pdgId()) == PDGID::MUON) && (itr->status() == 1) && (itr->pt() > pt_min) && (fabs(itr->eta()) < eta_max) && !alreadyFound) {
            muons.push_back(*itr);
        }
    }
    // turn this off, because we want the latest muons, not similar ones
    sort_by_pt(muons);
    return muons;
}

/**
 * Select reco jets that have a matching GenJet within some DR
 *
 * Stores index of matching GenJet.
 * Will take the closest matching GenJet as the match, provided it is within drMax.
 * uniqueMatch controls whether matching GenJets must be unique (i.e 2 reco jets can't match the same GenJet)
 *
 */
std::vector<Jet> QGAnalysisModule::getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax, bool uniqueMatch) {
    std::vector<Jet> goodJets;
    std::vector<uint> matchedIndices;
    for (auto jtr: *jets) {
        double minDR = 9999.;
        int matchInd = -1; // sensible default - not 0!

        for (uint gjInd=0; gjInd < genjets->size(); gjInd++) {
            // If we want unique matches and we've already matched then skip this genjet
            if (uniqueMatch && std::find(matchedIndices.begin(), matchedIndices.end(), gjInd) != matchedIndices.end())
                continue;

            const auto genjtr = genjets->at(gjInd);
            auto thisDR = deltaR(jtr, genjtr);
            if (thisDR < drMax && thisDR < minDR) {
                matchInd = gjInd;
                minDR = thisDR;
            }
        }

        jtr.set_genjet_index(matchInd);
        if (matchInd > -1) {
            goodJets.push_back(jtr);
            matchedIndices.push_back(matchInd);
            // cout << "Found a match at index " << matchInd << " with dr " << minDR << endl;
            // cout << "RECO pt/eta/phi: " << jtr.pt() << " : " << jtr.eta() << " : " << jtr.phi() << endl;
            // cout << "GEN pt/eta/phi: " << genjets->at(matchInd).pt() << " : " << genjets->at(matchInd).eta() << " : " << genjets->at(matchInd).phi() << endl;
        }
    }
    return goodJets;
}

void QGAnalysisModule::printGenParticles(const std::vector<GenParticle> & gps, const std::string & label, Color::Code color) {
    if (!PRINTOUT) return;
    for (auto & itr: gps) {
        if (itr.status() != 1) continue;
        cout << color << "GP";
        if (label != "") {
            cout << " [" << label << "]";
        }
        cout << ": " << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << endl;
    }
}

void QGAnalysisModule::printGenJets(const std::vector<GenJetWithParts> & gps, const std::string & label, Color::Code color) {
    if (!PRINTOUT) return;
    for (auto & itr: gps) {
        cout << color << "GenJet";
        if (label != "") {
            cout << " [" << label << "]";
        }
        cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << endl;
    }
}

void QGAnalysisModule::printJets(const std::vector<Jet> & jets, const std::string & label, Color::Code color) {
    if (!PRINTOUT) return;
    for (auto & itr: jets) {
        cout << color << "jet";
        if (label != "") {
            cout << " [" << label << "]";
        }
        cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << endl;
    }
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisModule)

}

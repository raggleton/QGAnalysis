#include <iostream>
#include <memory>
#include <algorithm>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/QGAnalysis/include/QGAnalysisFlavCompHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/Utils.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

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
class QGAnalysisFlavModule: public AnalysisModule {
public:

    explicit QGAnalysisFlavModule(Context & ctx);
    virtual bool process(Event & event) override;
    std::vector<GenJetWithParts> getGenJets(std::vector<GenJetWithParts> * jets, std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=1.5, float lepton_overlap_dr=0.2);

private:

    std::unique_ptr<CommonModules> common;

    std::unique_ptr<JetCorrector> jet_corrector_MC, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
    std::unique_ptr<GenericJetResolutionSmearer> jet_resolution_smearer;
    std::unique_ptr<JetCleaner> jet_cleaner;
    std::unique_ptr<JetElectronOverlapRemoval> jet_ele_cleaner;
    std::unique_ptr<JetMuonOverlapRemoval> jet_mu_cleaner;

    // Reco selections/hists
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, dijet_sel;

    std::unique_ptr<Hists> dijet_flav_correlation_presel, zpj_flav_correlation_presel;

    Event::Handle<std::vector<GenJetWithParts>> genjets_handle;

    bool is_mc;
    float jetRadius;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;

    std::vector<float> pt_mins = {30, 50, 100, 200, 500, 1000};
    std::vector<std::unique_ptr<Hists>> zplusjets_hists_binned;
    std::vector<std::unique_ptr<Hists>> dijet_hists_binned;
};


QGAnalysisFlavModule::QGAnalysisFlavModule(Context & ctx){
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
    common->set_electron_id(ElectronID_Spring16_medium);
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

    genjets_handle = ctx.declare_event_output< std::vector<GenJetWithParts> > ("GoodGenJets");

    // Event Selections
    njet_sel.reset(new NJetSelection(1));

    zplusjets_sel.reset(new ZplusJetsSelection());
    float deta = 1.2;
    float sumEta = 10.;
    dijet_sel.reset(new DijetSelection(2, 0.94, 0.3, true, deta, sumEta));

    // Hists
    dijet_flav_correlation_presel.reset(new QGAnalysisFlavCompHists(ctx, "Dijet_presel", 2, "dijet", 0));
    zpj_flav_correlation_presel.reset(new QGAnalysisFlavCompHists(ctx, "ZPlusJets_presel", 1, "zplusjets", 0));

    for (auto ptMin : pt_mins ) {
        std::unique_ptr<QGAnalysisFlavCompHists> dj(new QGAnalysisFlavCompHists(ctx, TString::Format("Dijet_%d", int(ptMin)).Data(), 2, "dijet", ptMin));
        dijet_hists_binned.push_back(std::move(dj));
        std::unique_ptr<QGAnalysisFlavCompHists> zpj(new QGAnalysisFlavCompHists(ctx, TString::Format("ZPlusJets_%d", int(ptMin)).Data(), 1, "zplusjets", ptMin));
        zplusjets_hists_binned.push_back(std::move(zpj));
    }
}


bool QGAnalysisFlavModule::process(Event & event) {

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
    jet_cleaner->process(event);
    jet_ele_cleaner->process(event);
    jet_mu_cleaner->process(event);

    // Resort by pT
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

    // RECO PART
    if (!njet_sel->passes(event)) return false;
    zpj_flav_correlation_presel->fill(event);
    
    if (event.jets->size() >=2) dijet_flav_correlation_presel->fill(event);

    bool zpj = zplusjets_sel->passes(event);
    if (zpj) {
        for (uint i=0; i < pt_mins.size(); i++) {
            zplusjets_hists_binned.at(i)->fill(event);
        }
    }

    bool dj = dijet_sel->passes(event);
    if (dj) {
        for (uint i=0; i < pt_mins.size(); i++) {
            dijet_hists_binned.at(i)->fill(event);
        }
    }

    return zpj || dj;
}


/**
 * Get GenJets, ignoring those that are basically a lepton, and have some minimum pt and maximum eta.
 */
std::vector<GenJetWithParts> QGAnalysisFlavModule::getGenJets(std::vector<GenJetWithParts> * jets, std::vector<GenParticle> * genparticles, float pt_min, float eta_max, float lepton_overlap_dr) {
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

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisFlavModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisFlavModule)

}

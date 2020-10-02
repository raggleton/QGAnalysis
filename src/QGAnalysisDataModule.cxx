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
#include "UHH2/common/include/JetHists.h"

#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/QGAnalysis/include/QGAnalysisHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisUnfoldHists.h"
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
    virtual void endInputData() override;

private:

    std::unique_ptr<GeneralEventSetup> common_setup;
    std::unique_ptr<RecoJetSetup> recojet_setup;
    Event::Handle<double> gen_weight_handle, pt_binning_reco_handle, pt_binning_gen_handle;;

    // Reco selections/hists
    std::unique_ptr<ZFinder> zFinder;
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, zplusjets_presel, dijet_sel;
    std::unique_ptr<Selection> zerobias_trigger_sel;
    std::unique_ptr<OrSelection> zplusjets_trigger_sel;
    std::unique_ptr<DataJetSelection> dijet_trigger_sel;
    std::vector<std::string> dj_trig_names, ak4dj_trig_names, ak8dj_trig_names;
    std::vector<double> dj_trig_prescales, ak4dj_trig_prescales, ak8dj_trig_prescales, dj_trig_thresholds;
    float jetht_zb_pt_boundary;

    std::unique_ptr<QGAnalysisJetLambda> jetLambdaCreatorPtSorted;
    std::unique_ptr<JetLambdaCopier> jetLambdaCopierForward, jetLambdaCopierCentral;

    Event::Handle<bool> pass_zpj_sel_handle, pass_dj_sel_handle;
    Event::Handle<std::vector<Jet>> dijet_forward_handle, dijet_central_handle;

    std::unique_ptr<Hists> zplusjets_hists_presel, zplusjets_hists, zplusjets_qg_hists, zplusjets_qg_hists_groomed;
    std::unique_ptr<Hists> zplusjets_qg_unfold_hists, zplusjets_qg_unfold_hists_groomed;

    std::unique_ptr<Hists> dijet_hists_presel, dijet_hists, dijet_hists_eta_ordered;
    std::unique_ptr<Hists> dijet_qg_hists, dijet_qg_hists_central_tighter, dijet_qg_hists_forward_tighter;
    std::unique_ptr<Hists> dijet_qg_hists_central_tighter_groomed, dijet_qg_hists_forward_tighter_groomed;
    std::unique_ptr<Hists> dijet_qg_unfold_hists_central_tighter, dijet_qg_unfold_hists_forward_tighter;
    std::unique_ptr<Hists> dijet_qg_unfold_hists_central_tighter_groomed, dijet_qg_unfold_hists_forward_tighter_groomed;
    std::vector<JetHists> dijet_jet_hists, dijet_jet_hists_unweighted;

    std::unique_ptr<EventNumberSelection> event_sel;
    DATASET::Name dataset;
    bool isZPlusJets = false;

    bool DO_UNFOLD_HISTS = true;
    bool DO_KINEMATIC_HISTS = true;
    bool DO_LAMBDA_HISTS = true;

    std::string zLabel;

    double TOTAL_LUMI = 35918.219492947;
    double ZB_LUMI = 0.029048362;
    double zb_prescale = TOTAL_LUMI / ZB_LUMI;

    std::unique_ptr<TriggerSelection> trig320;

    int nFailTrigEvents, nSelectedEvents, nTotalEvents;
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
    float jetRadius = get_jet_radius(jet_cone);

    cout << "Running with jet cone: " << jet_cone << endl;
    cout << "Running with PUS: " << pu_removal << endl;
    cout << "Is Z+jets: " << isZPlusJets << endl;
    isZPlusJets = (dataset == DATASET::SingleMu);

    DO_UNFOLD_HISTS = string2bool(ctx.get("DO_UNFOLD_HISTS", "true"));
    DO_KINEMATIC_HISTS = string2bool(ctx.get("DO_KINEMATIC_HISTS", "true"));
    DO_LAMBDA_HISTS = string2bool(ctx.get("DO_LAMBDA_HISTS", "true"));

    common_setup.reset(new GeneralEventSetup(ctx));
    recojet_setup.reset(new RecoJetSetup(ctx, pu_removal, jet_cone, jetRadius, Cuts::reco_jet_pt_min, 5.)); // set y large here, do y selection as part of dijet selection
    gen_weight_handle = ctx.get_handle<double>("gen_weight");
    pt_binning_reco_handle = ctx.get_handle<double>("pt_binning_reco_value"); // the value to use for reco pt bin e.g dijet average
    pt_binning_gen_handle = ctx.get_handle<double>("pt_binning_gen_value"); // the value to use for gen pt bin e.g dijet average

    // Save explicitly the forward/central jets
    std::string dijet_forward_handle_name("DijetForwardJet"), dijet_central_handle_name("DijetCentralJet");
    dijet_forward_handle = ctx.get_handle<std::vector<Jet>>(dijet_forward_handle_name); // vector so easily interchangeable with event.jets etc
    dijet_central_handle = ctx.get_handle<std::vector<Jet>>(dijet_central_handle_name);

    // Save handles to global pass/fail selection bools
    // Gen ones are dummy
    std::string pass_zpj_sel_handle_name("ZPlusJetsSelection"), pass_zpj_gen_sel_handle_name("ZPlusJetsGenSelection");
    std::string pass_dj_sel_handle_name("DijetSelection"), pass_dj_gen_sel_handle_name("DijetGenSelection");
    pass_zpj_sel_handle = ctx.get_handle<bool> (pass_zpj_sel_handle_name);
    pass_dj_sel_handle = ctx.get_handle<bool> (pass_dj_sel_handle_name);

    // Event Selections
    int NJETS_ZPJ = 1;
    int NJETS_DIJET = 2;
    int minNJets = isZPlusJets ? NJETS_ZPJ : NJETS_DIJET;
    njet_sel.reset(new NJetSelection(minNJets));

    // Z+JETS selection
    if (isZPlusJets) {
        zLabel = "zMuonCand";
        zFinder.reset(new ZFinder(ctx, "muons", zLabel));

        float second_jet_frac_max_zpj = 1000.3;
        float z_jet_asym_max = 100.4;
        zplusjets_sel.reset(new ZplusJetsSelection(ctx, zLabel, Cuts::reco_jet_pt_min, Cuts::jet_y_max, Cuts::reco_muon_pt_min, Cuts::reco_muon_pt_min, Cuts::mZ_window, Cuts::dphi_jet_z_min, second_jet_frac_max_zpj, Cuts::z_pt_min, Cuts::z_asym_max, "ZPlusJetsSel"));

        // Preselection for Z+J - only 2 muons to reco Z, no other event cuts
        float dphi_jet_z_min = 0.;
        second_jet_frac_max_zpj = 999.;
        z_jet_asym_max = 1.;
        zplusjets_presel.reset(new ZplusJetsSelection(ctx, zLabel, Cuts::reco_jet_pt_min, Cuts::jet_y_max, Cuts::reco_muon_pt_min, Cuts::reco_muon_pt_min, Cuts::mZ_window, dphi_jet_z_min, second_jet_frac_max_zpj, Cuts::z_pt_min, z_jet_asym_max, "ZPlusJetsSel_presel"));
    } else {
        // DIJET selection
        float second_jet_frac_max_dj = 10.94;
        bool ss_eta = false;
        float deta = 12;
        float sumEta = 10.;
        dijet_sel.reset(new DijetSelection(ctx, Cuts::reco_jet_pt_min, Cuts::jet_y_max, Cuts::dijet_dphi_min, second_jet_frac_max_dj, Cuts::jet_asym_max, ss_eta, deta, sumEta, "DijetSelCutFlow"));
    }

    // Lambda calculators
    bool doPuppi = (pu_removal == "PUPPI");
    int maxNJets = minNJets;
    float recoConstitPtMin = Cuts::constit_pt_min;
    float recoConstitEtaMax = Cuts::constit_eta_max;
    // FIXME: get stuff from ctx not extra args?
    std::string reco_jetlambda_handle_name = "JetLambdas";
    jetLambdaCreatorPtSorted.reset(new QGAnalysisJetLambda(ctx, jetRadius, maxNJets, doPuppi,
                                                           PtEtaCut(recoConstitPtMin, recoConstitEtaMax),
                                                           "jets", reco_jetlambda_handle_name));

    // create forward & central (Gen)JetLambda bundles by copying from the main ones
    std::string reco_jetlambda_forward_handle_name = "JetLambdasForward";
    jetLambdaCopierForward.reset(new JetLambdaCopier(ctx,
                                                     dijet_forward_handle_name,
                                                     reco_jetlambda_handle_name,
                                                     reco_jetlambda_forward_handle_name));

    std::string reco_jetlambda_central_handle_name = "JetLambdasCentral";
    jetLambdaCopierCentral.reset(new JetLambdaCopier(ctx,
                                                     dijet_central_handle_name,
                                                     reco_jetlambda_handle_name,
                                                     reco_jetlambda_central_handle_name));

    // dummy values
    std::string gen_jetlambda_handle_name = "GoodGenJetLambdas";
    std::string gen_jetlambda_forward_handle_name = "GoodGenJetLambdas";
    std::string gen_jetlambda_central_handle_name = "GoodGenJetLambdas";

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
    ak4dj_trig_names = {
        "HLT_PFJet40_v*",
        "HLT_PFJet60_v*",
        "HLT_PFJet80_v*",
        "HLT_PFJet140_v*",
        "HLT_PFJet200_v*",
        "HLT_PFJet260_v*",
        "HLT_PFJet320_v*",
        "HLT_PFJet400_v*",
        "HLT_PFJet450_v*"
    };
    ak8dj_trig_names = {
        "HLT_AK8PFJet40_v*",
        "HLT_AK8PFJet60_v*",
        "HLT_AK8PFJet80_v*",
        "HLT_AK8PFJet140_v*",
        "HLT_AK8PFJet200_v*",
        "HLT_AK8PFJet260_v*",
        "HLT_AK8PFJet320_v*",
        "HLT_AK8PFJet400_v*",
        "HLT_AK8PFJet450_v*"
    };

    dj_trig_thresholds = {
        65,
        88,
        114,
        180,
        254,
        318,
        400,
        500,
        625,
        // 50,
        // 65,
        // 88,
        // 109,
        // 180,
        // 250,
        // 318,
        // 388,
        // 479,
        // 535,
    };

    ak4dj_trig_prescales = {
        135892.2580,
        49972.7348,
        13168.7348,
        1499.7422,
        349.8364,
        61.1339,
        20.4793,
        6.9885,
        1.0000,
    };

    // ak8dj_trig_prescales = {
    //     674753.7,
    //     102172.0,
    //     33363.6,
    //     3316.4,
    //     390.9,
    //     64.7,
    //     22.0,
    //     7.3,
    //     1.0,
    // };

    ak8dj_trig_prescales = {
        TOTAL_LUMI / 0.04917,
        TOTAL_LUMI / 0.32476,
        TOTAL_LUMI / 0.99456,
        TOTAL_LUMI / 10.00565,
        TOTAL_LUMI / 84.89303,
        TOTAL_LUMI / 512.84106,
        TOTAL_LUMI / 1510.15511,
        TOTAL_LUMI / 4544.78556,
        TOTAL_LUMI / 33182.26210,
    };

    trig320.reset(new TriggerSelection("HLT_AK8PFJet320_v*"));

    std::vector<std::pair<float, float>> dj_trigger_bins;
    for (uint i=0; i<dj_trig_thresholds.size(); i++) {
        float lowerPt = dj_trig_thresholds.at(i);
        float upperPt = (i == dj_trig_thresholds.size()-1) ? 999999 : dj_trig_thresholds.at(i+1);
        dj_trigger_bins.push_back(std::make_pair(lowerPt, upperPt));
    }
    if (jet_cone == "AK8") {
        dijet_trigger_sel.reset(new DataJetSelection(ak8dj_trig_names, dj_trigger_bins));
        dj_trig_prescales = ak8dj_trig_prescales;
    } else {
        dijet_trigger_sel.reset(new DataJetSelection(ak4dj_trig_names, dj_trigger_bins));
        dj_trig_prescales = ak4dj_trig_prescales;
    }

    jetht_zb_pt_boundary = dj_trig_thresholds.at(0);

    // only setup needed hists to avoid extra objects in output file
    if (isZPlusJets) {
        // Hists
        if (DO_KINEMATIC_HISTS) {
            zplusjets_hists_presel.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel", zLabel));
            zplusjets_hists.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets", zLabel));
        }
        std::string zpj_sel = "zplusjets";
        if (DO_LAMBDA_HISTS) {
            zplusjets_qg_hists.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG",
                                                         NJETS_ZPJ, false, false, zpj_sel,
                                                         pass_zpj_sel_handle_name, pass_zpj_gen_sel_handle_name,
                                                         reco_jetlambda_handle_name, gen_jetlambda_handle_name));
            zplusjets_qg_hists_groomed.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG_groomed",
                                                                 NJETS_ZPJ, true, false, zpj_sel,
                                                                 pass_zpj_sel_handle_name, pass_zpj_gen_sel_handle_name,
                                                                 reco_jetlambda_handle_name, gen_jetlambda_handle_name));
        }
        if (DO_UNFOLD_HISTS) {
            zplusjets_qg_unfold_hists.reset(new QGAnalysisUnfoldHists(ctx, "ZPlusJets_QG_Unfold",
                                                                      NJETS_ZPJ, false, zpj_sel,
                                                                      pass_zpj_sel_handle_name, pass_zpj_gen_sel_handle_name,
                                                                      reco_jetlambda_handle_name, gen_jetlambda_handle_name));
            zplusjets_qg_unfold_hists_groomed.reset(new QGAnalysisUnfoldHists(ctx, "ZPlusJets_QG_Unfold_groomed",
                                                                              NJETS_ZPJ, true, zpj_sel,
                                                                              pass_zpj_sel_handle_name, pass_zpj_gen_sel_handle_name,
                                                                              reco_jetlambda_handle_name, gen_jetlambda_handle_name));
        }

    } else {
        std::string binning_method = "ave";
        std::string dj_sel = "dijet";
        if (DO_KINEMATIC_HISTS) {
            dijet_hists_presel.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel", binning_method));
            dijet_hists.reset(new QGAnalysisDijetHists(ctx, "Dijet_tighter", binning_method));
            dijet_hists_eta_ordered.reset(new QGAnalysisDijetHists(ctx, "Dijet_eta_ordered", binning_method));
            // add in simple jet kinematic hists, for each trigger, so we can show a manually stitched plot
            int counter = 0;
            for (const auto & trg_pair : dj_trigger_bins) {
                float pt_min = trg_pair.first;
                float pt_max = trg_pair.second;
                dijet_jet_hists.push_back(JetHists(ctx, TString::Format("Dijet_jet_hist_%d", counter).Data(), 2));
                dijet_jet_hists_unweighted.push_back(JetHists(ctx, TString::Format("Dijet_jet_hist_unweighted_%d", counter).Data(), 2, "", false));
                counter++;
            }
        }

        if (DO_LAMBDA_HISTS) {
            dijet_qg_hists.reset(new QGAnalysisHists(ctx, "Dijet_QG_tighter",
                                                     NJETS_DIJET, false, false, dj_sel,
                                                     pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                     reco_jetlambda_handle_name, gen_jetlambda_handle_name));
            dijet_qg_hists_central_tighter.reset(new QGAnalysisHists(ctx, "Dijet_QG_central_tighter",
                                                                     1, false, false, dj_sel,
                                                                     pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                     reco_jetlambda_central_handle_name, gen_jetlambda_central_handle_name));
            dijet_qg_hists_forward_tighter.reset(new QGAnalysisHists(ctx, "Dijet_QG_forward_tighter",
                                                                     1, false, false, dj_sel,
                                                                     pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                     reco_jetlambda_forward_handle_name, gen_jetlambda_forward_handle_name));

            dijet_qg_hists_central_tighter_groomed.reset(new QGAnalysisHists(ctx, "Dijet_QG_central_tighter_groomed",
                                                                             1, true, false, dj_sel,
                                                                             pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                             reco_jetlambda_central_handle_name, gen_jetlambda_central_handle_name));
            dijet_qg_hists_forward_tighter_groomed.reset(new QGAnalysisHists(ctx, "Dijet_QG_forward_tighter_groomed",
                                                                             1, true, false, dj_sel,
                                                                             pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                             reco_jetlambda_forward_handle_name, gen_jetlambda_forward_handle_name));
        }

        if (DO_UNFOLD_HISTS) {
            // unfolding hists
            // note that each of these does neutral+charged, and charged-only
            dijet_qg_unfold_hists_central_tighter.reset(new QGAnalysisUnfoldHists(ctx, "Dijet_QG_Unfold_central_tighter",
                                                                                  1, false, dj_sel,
                                                                                  pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                                  reco_jetlambda_central_handle_name, gen_jetlambda_central_handle_name));
            dijet_qg_unfold_hists_forward_tighter.reset(new QGAnalysisUnfoldHists(ctx, "Dijet_QG_Unfold_forward_tighter",
                                                                                  1, false, dj_sel,
                                                                                  pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                                  reco_jetlambda_forward_handle_name, gen_jetlambda_forward_handle_name));

            dijet_qg_unfold_hists_central_tighter_groomed.reset(new QGAnalysisUnfoldHists(ctx, "Dijet_QG_Unfold_central_tighter_groomed",
                                                                                          1, true, dj_sel,
                                                                                          pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                                          reco_jetlambda_central_handle_name, gen_jetlambda_central_handle_name));
            dijet_qg_unfold_hists_forward_tighter_groomed.reset(new QGAnalysisUnfoldHists(ctx, "Dijet_QG_Unfold_forward_tighter_groomed",
                                                                                          1, true, dj_sel,
                                                                                          pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                                          reco_jetlambda_forward_handle_name, gen_jetlambda_forward_handle_name));
        }
    } // end if else on DATASET

    nFailTrigEvents = 0;
    nSelectedEvents = 0;
    nTotalEvents = 0;

    // event_sel.reset(new EventNumberSelection({111}));
}


bool QGAnalysisDataModule::process(Event & event) {
    double orig_weight = 1.;
    event.set(gen_weight_handle, orig_weight); // need to set this at the start
    event.set(pt_binning_reco_handle, 0); // need to set this at the start
    event.set(pt_binning_gen_handle, 0); // need to set this at the start

    // Setup selection handles
    // will be overwritten below
    // -------------------------------------------------------------------------
    bool selected(false);
    event.set(pass_zpj_sel_handle, false);
    event.set(pass_dj_sel_handle, false);

    // Setup dijet jet handles to blanks incase they don't get set later
    // -------------------------------------------------------------------------
    event.set(dijet_forward_handle, std::vector<Jet>());
    event.set(dijet_central_handle, std::vector<Jet>());

    nTotalEvents++;

    // if (!event_sel->passes(event)) return false;
    if (event.run < 274954) {
        if (trig320->passes(event)) {
            throw std::runtime_error("SHOULDN'T BE PASSING");
        }
    }

    if (PRINTOUT) cout << Color::BG_RED << "-- Event: " << event.event << " : " << event.luminosityBlock << " : " << event.run << Color::BG_DEFAULT << endl;

    if (PRINTOUT) printMuons(*event.muons, "Precleaning");
    if (PRINTOUT) printElectrons(*event.electrons, "Precleaning");
    if (PRINTOUT) printJets(*event.jets, "Precleaning");

    if (!common_setup->process(event)) {
        if (PRINTOUT) cout << Color::FG_BLUE << "Failed common setup" << Color::FG_DEFAULT << endl;
        return false;
    }

    if (!recojet_setup->process(event)) {
        if (PRINTOUT) cout << Color::FG_BLUE << "Failed recojet setup" << Color::FG_DEFAULT << endl;
        return false;
    }

    if (!njet_sel->passes(event)) {
        if (PRINTOUT) cout << Color::FG_BLUE << "Failed njet_sel" << Color::FG_DEFAULT << endl;
        return false;
    }

    if (PRINTOUT) printJets(*event.jets, "Postcleaning");

    // Check trigger
    // Note that the ZeroBias/JetHT gets split depending on leading jet pT
    // So do it after JEC/JER etc (if that's how the threshold was measured)
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
        // to ensure we don't double-count events
        // that are also in JetHT, we ensure that the first jet has pT < the first
        // bin in JetHT. (inside dijet_trigger_sel is an equivalent check for JetHT)
        // passTrigger = true;
        passTrigger = (zerobias_trigger_sel->passes(event)) && (event.jets->at(0).pt() < jetht_zb_pt_boundary);
    }
    if (!passTrigger) {
        nFailTrigEvents++;
        if (PRINTOUT) {
            if (event.jets->at(0).pt() > 100) {
                for (auto tname : event.get_current_triggernames()) {
                    if (tname.rfind("Flag_", 0) == 0) continue;
                    if (tname.rfind("muGT", 0) == 0) continue;
                    auto ti = event.get_trigger_index(tname);
                    if (event.passes_trigger(ti))
                        cout << "   " << tname << " : " << event.passes_trigger(ti) << endl;
                }
            }

            cout << Color::BG_YELLOW << Color::FG_BLUE << " ! Failed trigger" << Color::BG_DEFAULT << Color::FG_DEFAULT << endl;
        }
        return false;
    }

    if (PRINTOUT) printMuons(*event.muons);
    if (PRINTOUT) printElectrons(*event.electrons);
    if (PRINTOUT) printJetsWithParts(*event.jets, event.pfparticles);

    if (isZPlusJets) {
        if (!zFinder->process(event))
            return false;

        event.set(pt_binning_reco_handle, event.jets->at(0).pt());
        if (zplusjets_presel->passes(event)) {
            if (DO_KINEMATIC_HISTS) {
                zplusjets_hists_presel->fill(event);
            }
            selected = zplusjets_sel->passes(event);
            event.set(pass_zpj_sel_handle, selected);
            if (selected) {
                // Calculate lambda vars for jets (pt ordered)
                // These will be used in various histogram classes
                // At this point, all objects should have had all necessary corrections, filtering, etc
                jetLambdaCreatorPtSorted->process(event);

                if (DO_KINEMATIC_HISTS) {
                    zplusjets_hists->fill(event);
                }
                if (DO_LAMBDA_HISTS) {
                    zplusjets_qg_hists->fill(event);
                    zplusjets_qg_hists_groomed->fill(event);
                }
                if (DO_UNFOLD_HISTS) {
                    zplusjets_qg_unfold_hists->fill(event);
                    zplusjets_qg_unfold_hists_groomed->fill(event);
                }
            } else if (PRINTOUT) {
                cout << Color::BG_WHITE << Color::FG_RED << " ! Failed Z+J selection" << Color::BG_DEFAULT << Color::FG_DEFAULT << endl;
            }
        } else if (PRINTOUT) {
            cout << Color::BG_WHITE << Color::FG_RED << " ! Failed Z+J selection" << Color::BG_DEFAULT << Color::FG_DEFAULT << endl;
        }

    } else {
        dijet_hists_presel->fill(event);
        selected = dijet_sel->passes(event);
        event.set(pass_dj_sel_handle, selected);
        if (selected) {
            // Sort by eta & assign to handles
            // assign forward/central jets to handles for later use
            std::vector<Jet> leadingJets(event.jets->begin(), event.jets->begin()+2);
            if (leadingJets.size() != 2) {
                throw std::runtime_error("Slicing jets gone wrong!");
            }
            double ave_pt = 0.5*(leadingJets[0].pt() + leadingJets[1].pt());
            event.set(pt_binning_reco_handle, ave_pt);
            sort_by_y(leadingJets);
            std::vector<Jet> forwardJet = {leadingJets[0]};
            std::vector<Jet> centralJet = {leadingJets[1]};
            event.set(dijet_forward_handle, forwardJet);
            event.set(dijet_central_handle, centralJet);

            // Calculate lambda vars for jets (pt ordered)
            // These will be used in various histogram classes
            // At this point, all objects should have had all necessary corrections, filtering, etc
            jetLambdaCreatorPtSorted->process(event);

            // Calculate lambda vars for recojets for dijets
            // These will be used in various histogram classes
            // At this point, all objects should have had all necessary corrections, filtering, etc
            jetLambdaCopierForward->process(event);
            jetLambdaCopierCentral->process(event);

            // apply prescale factor
            if (dataset == DATASET::JetHT) event.weight *= dj_trig_prescales.at(dj_trig_ind);
            else if (dataset == DATASET::ZeroBias) event.weight *= zb_prescale;

            if (DO_KINEMATIC_HISTS) {
                dijet_hists->fill(event);
                dijet_jet_hists[dj_trig_ind].fill(event);
                dijet_jet_hists_unweighted[dj_trig_ind].fill(event);
            }
            if (DO_LAMBDA_HISTS) {
                dijet_qg_hists->fill(event);
                dijet_qg_hists_central_tighter->fill(event);
                dijet_qg_hists_forward_tighter->fill(event);
                dijet_qg_hists_central_tighter_groomed->fill(event);
                dijet_qg_hists_forward_tighter_groomed->fill(event);
            }
            if (DO_UNFOLD_HISTS) {
                dijet_qg_unfold_hists_central_tighter->fill(event);
                dijet_qg_unfold_hists_forward_tighter->fill(event);
                dijet_qg_unfold_hists_central_tighter_groomed->fill(event);
                dijet_qg_unfold_hists_forward_tighter_groomed->fill(event);
            }

            // do eta-sorted dijet hists (where we need both jets)
            if (DO_KINEMATIC_HISTS) {
                // do dijet hists but sorted by eta (only one that matters about eta-ordering)
                // get them from event.jets and not the central/forward handles,
                // since event.jets has correct genjet_index
                std::vector<Jet> leadingJets(event.jets->begin(), event.jets->begin()+2);
                if (leadingJets.size() != 2) {
                    throw std::runtime_error("Slicing jets gone wrong!");
                }
                sort_by_y(leadingJets);  // by descending eta, so jets[0] = fwd, jets[1] = cen
                std::swap(*event.jets, leadingJets);
                dijet_hists_eta_ordered->fill(event);
            }
        } else if (PRINTOUT) {
            cout << Color::BG_WHITE << Color::FG_RED << " ! Failed dijet selection" << Color::BG_DEFAULT << Color::FG_DEFAULT << endl;
        }
    }
    nSelectedEvents++;
    return selected;
}


void QGAnalysisDataModule::endInputData(){
    cout << "Summary stats: " << endl;
    cout << " # fail trigger: " << nFailTrigEvents << endl;
    cout << " # selected: " << nSelectedEvents << endl;
    cout << " # input: " << nTotalEvents << endl;
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

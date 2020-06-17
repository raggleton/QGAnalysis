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
#include "UHH2/common/include/GenJetsHists.h"

#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/QGAnalysis/include/QGAnalysisHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisUnfoldHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsGenHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisDijetGenHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisDijetHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisTheoryHists.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/QGAnalysis/include/QGAnalysisPrinters.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

// print out info about collections, only use for debugging
const bool PRINTOUT = false;


/**
 * Analysis module for MC datasets
 */
class QGAnalysisMCModule: public AnalysisModule {
public:

    explicit QGAnalysisMCModule(Context & ctx);
    virtual bool process(Event & event) override;
    std::vector<GenJetWithParts> getGenJets(std::vector<GenJetWithParts> * jets, std::vector<GenParticle> * genparticles, float pt_min=5., float y_max=1.5, float lepton_overlap_dr=0.2);
    std::vector<GenParticle> getGenMuons(std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=2.5);
    std::vector<Jet> getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax=0.8, bool uniqueMatch=true);
    virtual void endInputData() override;
private:
    std::unique_ptr<GenJetClusterer> genjet_cluster; // to cluster genjets - needed for AK8 since MiniAOD starts at 150
    Event::Handle<std::vector<GenJetWithParts>> new_genjets_handle;
    std::unique_ptr<GeneralEventSetup> common_setup;
    std::unique_ptr<RecoJetSetup> recojet_setup;
    std::unique_ptr<GenJetSelector> genJet_selector;
    std::unique_ptr<MCReweighting> mc_reweight;
    std::unique_ptr<TrackingEfficiency> tracking_eff;
    std::unique_ptr<JetCleaner> jet_pf_id;

    std::unique_ptr<QGAnalysisJetLambda> jetLambdaCreatorPtSorted, jetLambdaCreatorForward, jetLambdaCreatorCentral;
    std::unique_ptr<QGAnalysisGenJetLambda> genjetLambdaCreatorPtSorted, genjetLambdaCreatorForward, genjetLambdaCreatorCentral;

    // For doing reco/gen jet matching e.g. after forward/central eta split
    std::unique_ptr<JetMatcher> jetMatcherPtOrdered, jetMatcherForward, jetMatcherCentral;

    Event::Handle<std::vector<GenJetWithParts>> genjets_handle;
    Event::Handle<std::vector<GenParticle>> genmuons_handle;
    Event::Handle<std::vector<Jet>> dijet_forward_handle, dijet_central_handle;
    Event::Handle<std::vector<GenJetWithParts>> dijet_gen_forward_handle, dijet_gen_central_handle;

    // Reco selections/hists
    std::unique_ptr<ZFinder> zFinder;
    std::unique_ptr<Selection> njet_min_sel, njet_two_sel, ngenjet_min_sel, ngenjet_two_sel, ngenjet_good_sel, ngenjet_good_two_sel;
    std::unique_ptr<Selection> zplusjets_sel, zplusjets_sel_passGen, zplusjets_presel, dijet_sel, dijet_sel_tighter, dijet_sel_tighter_passGen;

    std::unique_ptr<Hists> zplusjets_gen_hists;
    std::unique_ptr<Hists> zplusjets_hists_presel, zplusjets_hists;
    std::unique_ptr<Hists> zplusjets_hists_presel_q, zplusjets_hists_presel_g, zplusjets_hists_presel_unknown;
    std::unique_ptr<Hists> zplusjets_hists_q, zplusjets_hists_g, zplusjets_hists_unknown;
    std::unique_ptr<Hists> zplusjets_qg_hists, zplusjets_qg_hists_groomed;
    std::unique_ptr<Hists> zplusjets_qg_unfold_hists, zplusjets_qg_unfold_hists_groomed;

    std::unique_ptr<Hists> dijet_gen_hists;
    std::unique_ptr<Hists> dijet_hists_presel, dijet_hists, dijet_hists_eta_ordered, dijet_hists_tighter;
    std::unique_ptr<Hists> dijet_hists_presel_gg, dijet_hists_presel_qg, dijet_hists_presel_gq, dijet_hists_presel_qq;
    std::unique_ptr<Hists> dijet_hists_eta_ordered_gg, dijet_hists_eta_ordered_qg, dijet_hists_eta_ordered_gq, dijet_hists_eta_ordered_qq;
    std::unique_ptr<Hists> dijet_hists_presel_q_unknown, dijet_hists_presel_g_unknown, dijet_hists_presel_unknown_unknown, dijet_hists_presel_unknown_q, dijet_hists_presel_unknown_g;
    std::unique_ptr<Hists> dijet_qg_hists, dijet_qg_hists_tighter, dijet_qg_hists_central_tighter, dijet_qg_hists_forward_tighter, dijet_qg_hists_central_tighter_groomed, dijet_qg_hists_forward_tighter_groomed;
    std::unique_ptr<Hists> dijet_qg_unfold_hists_central_tighter, dijet_qg_unfold_hists_forward_tighter;
    std::unique_ptr<Hists> dijet_qg_unfold_hists_central_tighter_groomed, dijet_qg_unfold_hists_forward_tighter_groomed;

    std::unique_ptr<Hists> genjet_hists, genjet_hists_passZpJReco, genjet_hists_passDijetReco;

    // high pT-specific hists for tests
    std::unique_ptr<Hists> dijet_hists_presel_highPt, dijet_hists_highPt;
    std::unique_ptr<Hists> dijet_hists_presel_gg_highPt, dijet_hists_presel_qg_highPt, dijet_hists_presel_gq_highPt, dijet_hists_presel_qq_highPt;
    std::unique_ptr<Hists> dijet_hists_presel_q_unknown_highPt, dijet_hists_presel_g_unknown_highPt, dijet_hists_presel_unknown_unknown_highPt, dijet_hists_presel_unknown_q_highPt, dijet_hists_presel_unknown_g_highPt;
    std::unique_ptr<Hists> dijet_qg_hists_highPt;

    // for sweeping over PU
    std::vector<std::pair<int, int>> pu_bins = {
        std::make_pair(5, 15),
        std::make_pair(20, 25),
        std::make_pair(30, 40)
    };
    std::vector< std::unique_ptr<Selection> > sel_pu_binned;
    std::vector< std::unique_ptr<Hists> > zplusjets_qg_hists_pu_binned, zplusjets_qg_hists_groomed_pu_binned;
    std::vector< std::unique_ptr<Hists> > dijet_qg_hists_central_pu_binned, dijet_qg_hists_central_groomed_pu_binned, dijet_qg_hists_forward_pu_binned, dijet_qg_hists_forward_groomed_pu_binned;

    // Gen selection
    std::unique_ptr<Selection> zplusjets_gen_sel, zplusjets_gen_sel_passReco, dijet_gen_sel, dijet_gen_sel_passReco;

    Event::Handle<double> gen_weight_handle, pt_binning_reco_handle, pt_binning_gen_handle;
    Event::Handle<bool> pass_zpj_sel_handle, pass_zpj_gen_sel_handle, pass_dj_sel_handle, pass_dj_gen_sel_handle;

    float jetRadius;
    string jetCone;
    float htMax, partonKtMin, partonKtMax;

    bool DO_PU_BINNED_HISTS;
    bool DO_UNFOLD_HISTS;
    bool DO_FLAVOUR_HISTS;
    bool DO_KINEMATIC_HISTS;
    bool DO_LAMBDA_HISTS;

    std::unique_ptr<EventNumberSelection> event_sel, event_sel_printout;

    int NJETS_ZPJ, NJETS_DIJET;
    uint nOverlapEvents, nZPJEvents, nDijetEvents, nPassEvents;

    bool isZPlusJets;
};


QGAnalysisMCModule::QGAnalysisMCModule(Context & ctx){
    cout << "Running analysis module" << endl;

    htMax = boost::lexical_cast<float>(ctx.get("maxHT", "-1"));
    partonKtMin = boost::lexical_cast<float>(ctx.get("partonKtMin", "-1"));
    partonKtMax = boost::lexical_cast<float>(ctx.get("partonKtMax", "-1")); // -ve disables cut

    jetCone = ctx.get("JetCone", "AK4");
    string pu_removal = ctx.get("PURemoval", "CHS");
    if (pu_removal != "CHS" && pu_removal != "PUPPI") {
        throw runtime_error("Only PURemoval == CHS, PUPPI supported for now");
    }
    jetRadius = get_jet_radius(jetCone);

    const std::string & datasetVersion = ctx.get("dataset_version");
    isZPlusJets = (isSubstring(datasetVersion, "DYJetsToLL", true) || isSubstring(datasetVersion, "SingleMu", true) || isSubstring(datasetVersion, "ZPJ", true));
    ctx.set("isZPlusJets", bool2string(isZPlusJets));

    DO_PU_BINNED_HISTS = string2bool(ctx.get("DO_PU_BINNED_HISTS", "false"));
    DO_UNFOLD_HISTS = string2bool(ctx.get("DO_UNFOLD_HISTS", "true"));
    DO_FLAVOUR_HISTS = string2bool(ctx.get("DO_FLAVOUR_HISTS", "false"));
    DO_KINEMATIC_HISTS = string2bool(ctx.get("DO_KINEMATIC_HISTS", "true"));
    DO_LAMBDA_HISTS = string2bool(ctx.get("DO_LAMBDA_HISTS", "true"));

    cout << "Running with jet cone: " << jetCone << endl;
    cout << "Running with PUS: " << pu_removal << endl;
    cout << "Is Z+jets: " << isZPlusJets << endl;
    cout << "jetRadius==0.8: " << (bool)(jetRadius == 0.8) << endl;
    if (jetCone == "AK8") {
        // FIXME add cut nonu, final state
        string new_genjet_name = "newgenjets";
        genjet_cluster.reset(new GenJetClusterer(ctx, new_genjet_name, 0.8, "genparticles"));
        new_genjets_handle = ctx.get_handle< std::vector<GenJetWithParts> > (new_genjet_name);
    }

    // FIXME: get everything from ctx not extra args
    common_setup.reset(new GeneralEventSetup(ctx));
    bool update4vec = false;
    tracking_eff.reset(new TrackingEfficiency(ctx, update4vec));
    bool doJetId = false; // do it AFTER the tracking SF and not before - could have some promoted particles
    recojet_setup.reset(new RecoJetSetup(ctx, pu_removal, jetCone, jetRadius, Cuts::reco_jet_pt_min, Cuts::jet_y_max, doJetId));

    // another jet ID check after tracking SFs applied (basically constituent check)
    jet_pf_id.reset(new JetCleaner(ctx, JetPFID(Cuts::RECO_JET_ID)));

    // Setup Gen objects
    std::string genjet_handle_name = "GoodGenJets";
    genJet_selector.reset(new GenJetSelector(ctx, Cuts::gen_jet_pt_min, Cuts::jet_y_max, jetRadius, "genjets", genjet_handle_name, "genparticles"));
    std::string genmuon_handle_name = "GoodGenMuons";
    genjets_handle = ctx.get_handle< std::vector<GenJetWithParts> > (genjet_handle_name);
    genmuons_handle = ctx.get_handle< std::vector<GenParticle> > (genmuon_handle_name);

    mc_reweight.reset(new MCReweighting(ctx, genjet_handle_name, genmuon_handle_name));

    gen_weight_handle = ctx.get_handle<double>("gen_weight");
    pt_binning_reco_handle = ctx.get_handle<double>("pt_binning_reco_value"); // the value to use for reco pt bin e.g dijet average
    pt_binning_gen_handle = ctx.get_handle<double>("pt_binning_gen_value"); // the value to use for gen pt bin e.g dijet average

    // Save explicitly the forward/central (gen) jets
    std::string dijet_forward_handle_name("DijetForwardJet"), dijet_central_handle_name("DijetCentralJet");
    std::string dijet_gen_forward_handle_name("DijetForwardGenJet"), dijet_gen_central_handle_name("DijetCentralGenJet");
    dijet_forward_handle = ctx.get_handle<std::vector<Jet>>(dijet_forward_handle_name); // vector so easily interchangeable with event.jets etc
    dijet_central_handle = ctx.get_handle<std::vector<Jet>>(dijet_central_handle_name);
    dijet_gen_forward_handle = ctx.get_handle<std::vector<GenJetWithParts>>(dijet_gen_forward_handle_name);
    dijet_gen_central_handle = ctx.get_handle<std::vector<GenJetWithParts>>(dijet_gen_central_handle_name);

    // Save handles to global pass/fail selection bools
    std::string pass_zpj_sel_handle_name("ZPlusJetsSelection"), pass_zpj_gen_sel_handle_name("ZPlusJetsGenSelection");
    std::string pass_dj_sel_handle_name("DijetSelection"), pass_dj_gen_sel_handle_name("DijetGenSelection");
    pass_zpj_sel_handle = ctx.get_handle<bool> (pass_zpj_sel_handle_name);
    pass_zpj_gen_sel_handle = ctx.get_handle<bool> (pass_zpj_gen_sel_handle_name);
    pass_dj_sel_handle = ctx.get_handle<bool> (pass_dj_sel_handle_name);
    pass_dj_gen_sel_handle = ctx.get_handle<bool> (pass_dj_gen_sel_handle_name);

    float drMatch = jetRadius/2.;
    bool uniqueMatch = true;
    jetMatcherPtOrdered.reset(new JetMatcher(ctx, "jets", genjet_handle_name, drMatch, uniqueMatch));
    jetMatcherForward.reset(new JetMatcher(ctx, dijet_forward_handle_name, dijet_gen_forward_handle_name, drMatch, uniqueMatch));
    jetMatcherCentral.reset(new JetMatcher(ctx, dijet_central_handle_name, dijet_gen_central_handle_name, drMatch, uniqueMatch));

    // Event Selections
    NJETS_ZPJ = 1;
    NJETS_DIJET = 2;
    int minNJets = isZPlusJets ? NJETS_ZPJ : NJETS_DIJET;
    njet_min_sel.reset(new NJetSelection(minNJets));
    // njet_two_sel.reset(new NJetSelection(NJETS_DIJET));
    ngenjet_min_sel.reset(new NGenJetWithPartsSelection(minNJets));
    // ngenjet_two_sel.reset(new NGenJetSelection(NJETS_DIJET)); //, -1, boost::none, genjet_handle_name));
    ngenjet_good_sel.reset(new NGenJetWithPartsSelection(minNJets, -1, boost::none, genjets_handle));
    // ngenjet_good_two_sel.reset(new NGenJetSelection(2, -1, boost::none, genjets_handle));

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

    float genConstitPtMin = Cuts::constit_pt_min;
    float genConstitEtaMax = Cuts::constit_eta_max;
    std::string gen_jetlambda_handle_name = "GoodGenJetLambdas";
    genjetLambdaCreatorPtSorted.reset(new QGAnalysisGenJetLambda(ctx, jetRadius, 5, // allow more jets for possible reco/gen matching outside of top 2
                                                                 PtEtaCut(genConstitPtMin, genConstitEtaMax),
                                                                 genjet_handle_name, gen_jetlambda_handle_name));

    std::string reco_jetlambda_forward_handle_name = "JetLambdasForward";
    jetLambdaCreatorForward.reset(new QGAnalysisJetLambda(ctx, jetRadius, 1, doPuppi,
                                                          PtEtaCut(recoConstitPtMin, recoConstitEtaMax),
                                                          dijet_forward_handle_name, reco_jetlambda_forward_handle_name));

    std::string gen_jetlambda_forward_handle_name = "GoodGenJetLambdasForward";
    genjetLambdaCreatorForward.reset(new QGAnalysisGenJetLambda(ctx, jetRadius, 1,
                                                                PtEtaCut(genConstitPtMin, genConstitEtaMax),
                                                                dijet_gen_forward_handle_name, gen_jetlambda_forward_handle_name));

    std::string reco_jetlambda_central_handle_name = "JetLambdasCentral";
    jetLambdaCreatorCentral.reset(new QGAnalysisJetLambda(ctx, jetRadius, 1, doPuppi,
                                                          PtEtaCut(recoConstitPtMin, recoConstitEtaMax),
                                                          dijet_central_handle_name, reco_jetlambda_central_handle_name));

    std::string gen_jetlambda_central_handle_name = "GoodGenJetLambdasCentral";
    genjetLambdaCreatorCentral.reset(new QGAnalysisGenJetLambda(ctx, jetRadius, 1,
                                                                PtEtaCut(genConstitPtMin, genConstitEtaMax),
                                                                dijet_gen_central_handle_name, gen_jetlambda_central_handle_name));

    // Setup for systematics
    // FIXME put all this inside the ctor as it has ctx!
    std::string chargedHadronShift = ctx.get("chargedHadronShift", "nominal");
    float chargedHadronShiftAmount = 0.01;
    if (chargedHadronShift == "nominal") {
        // pass
    } else if (chargedHadronShift == "up") {
        jetLambdaCreatorPtSorted->set_charged_hadron_shift(1, chargedHadronShiftAmount);
        jetLambdaCreatorCentral->set_charged_hadron_shift(1, chargedHadronShiftAmount);
        jetLambdaCreatorForward->set_charged_hadron_shift(1, chargedHadronShiftAmount);
    } else if (chargedHadronShift == "down") {
        jetLambdaCreatorPtSorted->set_charged_hadron_shift(-1, chargedHadronShiftAmount);
        jetLambdaCreatorCentral->set_charged_hadron_shift(-1, chargedHadronShiftAmount);
        jetLambdaCreatorForward->set_charged_hadron_shift(-1, chargedHadronShiftAmount);
    } else {
        throw runtime_error("chargedHadronShift must be nominal, up, or down");
    }

    std::string neutralHadronShift = ctx.get("neutralHadronShift", "nominal");
    float neutralHadronShiftAmount = 0.03;
    if (neutralHadronShift == "nominal") {
        // pass
    } else if (neutralHadronShift == "up") {
        jetLambdaCreatorPtSorted->set_neutral_hadron_shift(1, neutralHadronShiftAmount);
        jetLambdaCreatorCentral->set_neutral_hadron_shift(1, neutralHadronShiftAmount);
        jetLambdaCreatorForward->set_neutral_hadron_shift(1, neutralHadronShiftAmount);
    } else if (neutralHadronShift == "down") {
        jetLambdaCreatorPtSorted->set_neutral_hadron_shift(-1, neutralHadronShiftAmount);
        jetLambdaCreatorCentral->set_neutral_hadron_shift(-1, neutralHadronShiftAmount);
        jetLambdaCreatorForward->set_neutral_hadron_shift(-1, neutralHadronShiftAmount);
    } else {
        throw runtime_error("neutralHadronShift must be nominal, up, or down");
    }

    std::string photonShift = ctx.get("photonShift", "nominal");
    float photonShiftAmount = 0.01;
    if (photonShift == "nominal") {
        // pass
    } else if (photonShift == "up") {
        jetLambdaCreatorPtSorted->set_photon_shift(1, photonShiftAmount);
        jetLambdaCreatorCentral->set_photon_shift(1, photonShiftAmount);
        jetLambdaCreatorForward->set_photon_shift(1, photonShiftAmount);
    } else if (photonShift == "down") {
        jetLambdaCreatorPtSorted->set_photon_shift(-1, photonShiftAmount);
        jetLambdaCreatorCentral->set_photon_shift(-1, photonShiftAmount);
        jetLambdaCreatorForward->set_photon_shift(-1, photonShiftAmount);
    } else {
        throw runtime_error("photonShift must be nominal, up, or down");
    }

    std::string zLabel = "zMuonCand";



    // Note that Gen selections have basically same cuts as reco,
    // to avoid over-reliance on MC to extrapolate to wider gen phase space
    // A lot of variables are very well-measured anyway

    if (isZPlusJets) {
        zFinder.reset(new ZFinder(ctx, "muons", zLabel));

        // Z+JETS selection
        float mu1_pt = Cuts::reco_muon_pt_min;
        float mu2_pt = Cuts::reco_muon_pt_min;
        float second_jet_frac_max_zpj = 1000.3;
        float z_jet_asym_max = 100.;
        zplusjets_sel.reset(new ZplusJetsSelection(ctx, zLabel, mu1_pt, mu2_pt, Cuts::mZ_window, Cuts::dphi_jet_z_min, second_jet_frac_max_zpj, Cuts::z_pt_min, z_jet_asym_max, "ZPlusJetsSelCutFlow"));
        // just to plot cutflow only when passGen==true
        zplusjets_sel_passGen.reset(new ZplusJetsSelection(ctx, zLabel, mu1_pt, mu2_pt, Cuts::mZ_window, Cuts::dphi_jet_z_min, second_jet_frac_max_zpj, Cuts::z_pt_min, z_jet_asym_max, "ZPlusJetsSelPassGenCutFlow"));

        zplusjets_gen_sel.reset(new ZplusJetsGenSelection(ctx, mu1_pt, mu2_pt, Cuts::mZ_window, Cuts::dphi_jet_z_min, second_jet_frac_max_zpj, Cuts::z_pt_min, z_jet_asym_max,
                                                          "ZPlusJetsGenSelCutFlow", genjet_handle_name, genmuon_handle_name));
        // just to plot gen cutflow only when passReco==true
        zplusjets_gen_sel_passReco.reset(new ZplusJetsGenSelection(ctx, mu1_pt, mu2_pt, Cuts::mZ_window, Cuts::dphi_jet_z_min, second_jet_frac_max_zpj, Cuts::z_pt_min, z_jet_asym_max,
                                                                   "ZPlusJetsGenSelPassRecoCutFlow", genjet_handle_name, genmuon_handle_name));

        // Preselection for Z+J - only 2 muons to reco Z
        // dphi_jet_z_min = 0.;
        // second_jet_frac_max_zpj = 999.;
        // z_jet_asym_max = 1.;
        // zplusjets_presel.reset(new ZplusJetsSelection(ctx, zLabel, mu1_pt, mu2_pt, Cuts::mZ_window, dphi_jet_z_min, second_jet_frac_max_zpj, Cuts::z_pt_min, z_jet_asym_max));

    } else {

        // DIJET selection
        float second_jet_frac_max_dj = 10.94;
        bool ss_eta = false;
        float deta = 12;
        float sumEta = 10.;
        dijet_sel.reset(new DijetSelection(ctx, Cuts::dijet_dphi_min, second_jet_frac_max_dj, 1000, ss_eta, deta, sumEta, "DijetSelCutFlow")); // this is without any asymmetry cut
        dijet_sel_tighter.reset(new DijetSelection(ctx, Cuts::dijet_dphi_min, second_jet_frac_max_dj, Cuts::jet_asym_max, ss_eta, deta, sumEta, "DijetSelTighterCutFlow"));
        // just to plot cutflow only when passGen == true
        dijet_sel_tighter_passGen.reset(new DijetSelection(ctx, Cuts::dijet_dphi_min, second_jet_frac_max_dj, Cuts::jet_asym_max, ss_eta, deta, sumEta, "DijetSelTighterPassGenCutFlow"));

        dijet_gen_sel.reset(new DijetGenSelection(ctx, Cuts::dijet_dphi_min, second_jet_frac_max_dj, Cuts::jet_asym_max, ss_eta, deta, sumEta,
                                                  "DijetGenSelCutFlow", genjet_handle_name));
        // just to plot gen cutflow only when passReco==true
        dijet_gen_sel_passReco.reset(new DijetGenSelection(ctx, Cuts::dijet_dphi_min, second_jet_frac_max_dj, Cuts::jet_asym_max, ss_eta, deta, sumEta,
                                                           "DijetGenSelPassRecoCutFlow", genjet_handle_name));
    }

    if (DO_PU_BINNED_HISTS) {
        for (auto puBin : pu_bins) {
            std::unique_ptr<Selection> pu_sel(new NPVSelection(puBin.first, puBin.second));
            sel_pu_binned.push_back(std::move(pu_sel));
        }
    }

    // Hists
    // -------------------------------------------------------------------------
    std::string zpj_sel = "zplusjets";
    if (isZPlusJets) {
        // Z+JETS hists
        if (DO_KINEMATIC_HISTS) {
            zplusjets_gen_hists.reset(new QGAnalysisZPlusJetsGenHists(ctx, "ZPlusJets_gen"));
            zplusjets_hists_presel.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel", zLabel));
            zplusjets_hists.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets", zLabel));
        }

        if (DO_FLAVOUR_HISTS) {
            // preselection hists, if jet is quark, or gluon
            // zplusjets_hists_presel_q.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_q", zLabel));
            // zplusjets_hists_presel_g.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_g", zLabel));
            // zplusjets_hists_presel_unknown.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_unknown", zLabel));

            zplusjets_hists_q.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_q", zLabel));
            zplusjets_hists_g.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_g", zLabel));
            zplusjets_hists_unknown.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_unknown", zLabel));
        }

        // Lambda variables, used for e.g. response, determine binning
        if (DO_LAMBDA_HISTS) {
            // note that each of these does neutral+charged, and charged-only
            zplusjets_qg_hists.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG",
                                                         NJETS_ZPJ, false, zpj_sel,
                                                         pass_zpj_sel_handle_name, pass_zpj_gen_sel_handle_name,
                                                         reco_jetlambda_handle_name, gen_jetlambda_handle_name));
            zplusjets_qg_hists_groomed.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG_groomed",
                                                                 NJETS_ZPJ, true, zpj_sel,
                                                                 pass_zpj_sel_handle_name, pass_zpj_gen_sel_handle_name,
                                                                 reco_jetlambda_handle_name, gen_jetlambda_handle_name));
        }

        // With special binning for unfolding
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

        if (DO_PU_BINNED_HISTS) {
            for (auto puBin : pu_bins) {
                std::unique_ptr<Selection> pu_sel(new NPVSelection(puBin.first, puBin.second));
                sel_pu_binned.push_back(std::move(pu_sel));

                // ungroomed & groomed z+j
                std::unique_ptr<QGAnalysisHists> zpj(new QGAnalysisHists(ctx, TString::Format("ZPlusJets_QG_PU_%d_to_%d", puBin.first, puBin.second).Data(),
                                                                         NJETS_ZPJ, false, zpj_sel,
                                                                         pass_zpj_sel_handle_name, pass_zpj_gen_sel_handle_name,
                                                                         reco_jetlambda_handle_name, gen_jetlambda_handle_name));
                zplusjets_qg_hists_pu_binned.push_back(std::move(zpj));
                std::unique_ptr<QGAnalysisHists> zpjg(new QGAnalysisHists(ctx, TString::Format("ZPlusJets_QG_PU_%d_to_%d_groomed", puBin.first, puBin.second).Data(),
                                                                          NJETS_ZPJ, true, zpj_sel,
                                                                          pass_zpj_sel_handle_name, pass_zpj_gen_sel_handle_name,
                                                                          reco_jetlambda_handle_name, gen_jetlambda_handle_name));
                zplusjets_qg_hists_groomed_pu_binned.push_back(std::move(zpjg));
            }
        }

    } else {

        // DIJET hists
        std::string binning_method = "ave";
        if (DO_KINEMATIC_HISTS) {
            dijet_gen_hists.reset(new QGAnalysisDijetGenHists(ctx, "Dijet_gen"));
            dijet_hists_presel.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel", binning_method));
            dijet_hists.reset(new QGAnalysisDijetHists(ctx, "Dijet", binning_method));
            dijet_hists_eta_ordered.reset(new QGAnalysisDijetHists(ctx, "Dijet_eta_ordered", binning_method));
            dijet_hists_tighter.reset(new QGAnalysisDijetHists(ctx, "Dijet_tighter", binning_method));
        }

        if (DO_FLAVOUR_HISTS) {
            // preselection hiss, if both gluon jets, one gluon, or both quark, or one or both unknown
            // dijet_hists_presel_gg.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_gg", binning_method));
            // dijet_hists_presel_qg.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_qg", binning_method));
            // dijet_hists_presel_gq.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_gq", binning_method));
            // dijet_hists_presel_qq.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_qq", binning_method));
            // dijet_hists_presel_unknown_q.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_q", binning_method));
            // dijet_hists_presel_unknown_g.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_g", binning_method));
            // dijet_hists_presel_unknown_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_unknown", binning_method));
            // dijet_hists_presel_q_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_q_unknown", binning_method));
            // dijet_hists_presel_g_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_g_unknown", binning_method));

            // eta-ordered flavour-specific hists
            dijet_hists_eta_ordered_gg.reset(new QGAnalysisDijetHists(ctx, "Dijet_eta_ordered_gg", binning_method));
            dijet_hists_eta_ordered_qg.reset(new QGAnalysisDijetHists(ctx, "Dijet_eta_ordered_qg", binning_method));
            dijet_hists_eta_ordered_gq.reset(new QGAnalysisDijetHists(ctx, "Dijet_eta_ordered_gq", binning_method));
            dijet_hists_eta_ordered_qq.reset(new QGAnalysisDijetHists(ctx, "Dijet_eta_ordered_qq", binning_method));
        }

        // note that each of these does neutral+charged, and charged-only
        std::string dj_sel = "dijet";
        if (DO_LAMBDA_HISTS) {
            dijet_qg_hists.reset(new QGAnalysisHists(ctx, "Dijet_QG",
                                                     NJETS_DIJET, false, dj_sel,
                                                     pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                     reco_jetlambda_handle_name, gen_jetlambda_handle_name));
            dijet_qg_hists_tighter.reset(new QGAnalysisHists(ctx, "Dijet_QG_tighter",
                                                             NJETS_DIJET, false, dj_sel,
                                                             pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                             reco_jetlambda_handle_name, gen_jetlambda_handle_name));
            dijet_qg_hists_central_tighter.reset(new QGAnalysisHists(ctx, "Dijet_QG_central_tighter",
                                                                     1, false, dj_sel,
                                                                     pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                     reco_jetlambda_central_handle_name, gen_jetlambda_central_handle_name));
            dijet_qg_hists_forward_tighter.reset(new QGAnalysisHists(ctx, "Dijet_QG_forward_tighter",
                                                                     1, false, dj_sel,
                                                                     pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                     reco_jetlambda_forward_handle_name, gen_jetlambda_forward_handle_name));

            dijet_qg_hists_central_tighter_groomed.reset(new QGAnalysisHists(ctx, "Dijet_QG_central_tighter_groomed",
                                                                             1, true, dj_sel,
                                                                             pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                             reco_jetlambda_central_handle_name, gen_jetlambda_central_handle_name));
            dijet_qg_hists_forward_tighter_groomed.reset(new QGAnalysisHists(ctx, "Dijet_QG_forward_tighter_groomed",
                                                                             1, true, dj_sel,
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

        if (DO_PU_BINNED_HISTS) {
            for (auto puBin : pu_bins) {
                // ungroomed dijet
                std::unique_ptr<QGAnalysisHists> djc(new QGAnalysisHists(ctx, TString::Format("Dijet_QG_central_tighter_PU_%d_to_%d", puBin.first, puBin.second).Data(),
                                                                         1, false, dj_sel,
                                                                         pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                         reco_jetlambda_central_handle_name, gen_jetlambda_central_handle_name));
                dijet_qg_hists_central_pu_binned.push_back(std::move(djc));
                std::unique_ptr<QGAnalysisHists> djf(new QGAnalysisHists(ctx, TString::Format("Dijet_QG_forward_tighter_PU_%d_to_%d", puBin.first, puBin.second).Data(),
                                                                         1, false, dj_sel,
                                                                         pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                         reco_jetlambda_forward_handle_name, gen_jetlambda_forward_handle_name));
                dijet_qg_hists_forward_pu_binned.push_back(std::move(djf));

                // groomed dijet
                std::unique_ptr<QGAnalysisHists> djcg(new QGAnalysisHists(ctx, TString::Format("Dijet_QG_central_tighter_PU_%d_to_%d_groomed", puBin.first, puBin.second).Data(),
                                                                          1, true, dj_sel,
                                                                          pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                          reco_jetlambda_central_handle_name, gen_jetlambda_central_handle_name));
                dijet_qg_hists_central_groomed_pu_binned.push_back(std::move(djcg));
                std::unique_ptr<QGAnalysisHists> djfg(new QGAnalysisHists(ctx, TString::Format("Dijet_QG_forward_tighter_PU_%d_to_%d_groomed", puBin.first, puBin.second).Data(),
                                                                          1, true, dj_sel,
                                                                          pass_dj_sel_handle_name, pass_dj_gen_sel_handle_name,
                                                                          reco_jetlambda_forward_handle_name, gen_jetlambda_forward_handle_name));
                dijet_qg_hists_forward_groomed_pu_binned.push_back(std::move(djfg));
            }
        }
    }

    // dijet_hists_presel_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_highPt", binning_method));
    // preselection hiss, if both gluon jets, one gluon, or both quark, or one or both unknown
    // dijet_hists_presel_gg_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_gg_highPt", binning_method));
    // dijet_hists_presel_qg_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_qg_highPt", binning_method));
    // dijet_hists_presel_gq_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_gq_highPt", binning_method));
    // dijet_hists_presel_qq_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_qq_highPt", binning_method));
    // dijet_hists_presel_unknown_q_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_q_highPt", binning_method));
    // dijet_hists_presel_unknown_g_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_g_highPt", binning_method));
    // dijet_hists_presel_unknown_unknown_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_unknown_highPt", binning_method));
    // dijet_hists_presel_q_unknown_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_q_unknown_highPt", binning_method));
    // dijet_hists_presel_g_unknown_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_g_unknown_highPt", binning_method));
    // dijet_hists_highPt.reset(new QGAnalysisDijetHists(ctx, "Dijet_highPt", binning_method));
    // dijet_qg_hists_highPt.reset(new QGAnalysisHists(ctx, "Dijet_QG_highPt", NJETS_DIJET, "dijet"));

    bool useRapidity = true;
    genjet_hists.reset(new GenJetsHists(ctx, "GenJetsPresel", 3, genjet_handle_name, useRapidity));
    if (isZPlusJets) {
        genjet_hists_passZpJReco.reset(new GenJetsHists(ctx, "GenJetsPassZPlusJetsReco", 3, genjet_handle_name, useRapidity));
    } else {
        genjet_hists_passDijetReco.reset(new GenJetsHists(ctx, "GenJetsPassDijetReco", 3, genjet_handle_name, useRapidity));
    }

    // event_sel.reset(new EventNumberSelection({65406240}));
    event_sel_printout.reset(new EventNumberSelection({65406240, 98686924}));

    nOverlapEvents = 0;
    nZPJEvents = 0;
    nDijetEvents = 0;
    nPassEvents = 0;
}


bool QGAnalysisMCModule::process(Event & event) {

    bool printout_event_sel = event_sel_printout->passes(event);

    // if (!event_sel->passes(event)) return false;

    event.set(gen_weight_handle, event.weight); // need to set this at the start
    event.set(pt_binning_reco_handle, 0); // need to set this at the start
    event.set(pt_binning_gen_handle, 0); // need to set this at the start

    // Save selection flags incase they don't get set later
    // -------------------------------------------------------------------------
    bool pass_zpj_reco(false), pass_dj_reco(false);
    bool pass_zpj_gen(false), pass_dj_gen(false);
    // incase they doesn't get set later
    event.set(pass_zpj_sel_handle, pass_zpj_reco);
    event.set(pass_dj_sel_handle, pass_dj_reco);
    event.set(pass_zpj_gen_sel_handle, pass_zpj_gen);
    event.set(pass_dj_gen_sel_handle, pass_dj_gen);
    // bool pass_dj_highPt(false);

    // Setup dijet jet handles to blanks incase they don't get set later
    // -------------------------------------------------------------------------
    event.set(dijet_forward_handle, std::vector<Jet>());
    event.set(dijet_central_handle, std::vector<Jet>());
    event.set(dijet_gen_forward_handle, std::vector<GenJetWithParts>());
    event.set(dijet_gen_central_handle, std::vector<GenJetWithParts>());

    if (PRINTOUT || printout_event_sel) { cout << "-- Event: " << event.event << endl; }
    // cout << "-- Event: " << event.event << endl;

    // Redo AK8 genjet clustering
    // Do before any ngenjet cuts
    // This is so we can get all gen jets, since the ntuple only had pt > 170
    // -------------------------------------------------------------------------
    if (jetCone == "AK8") {
        // for (const auto & itr : *(event.genjets)) {
        //     cout << "Original gj: " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << " : " << itr.genparticles_indices().size() << endl;
        //     for (const auto & citr : itr.genparticles_indices()) {
        //         auto constit = event.genparticles->at(citr);
        //         cout << "    constit: " << constit.pt() << " : " << constit.eta() << " : " << constit.phi() << " : " << constit.pdgId() << " : " << constit.status() << " : " << constit.daughter1() << " : " << constit.daughter2() << endl;
        //     }
        // }
        genjet_cluster->process(event);
        std::vector<GenJetWithParts> newGenJets = event.get(new_genjets_handle);
        std::swap(*event.genjets, newGenJets);
        // for (const auto & itr : *(event.genjets)) {
        //     cout << "New gj: " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << " : " << itr.genparticles_indices().size() << endl;
        //     for (const auto & citr : itr.genparticles_indices()) {
        //         auto constit = event.genparticles->at(citr);
        //         cout << "    constit: " << constit.pt() << " : " << constit.eta() << " : " << constit.phi() << " : " << constit.pdgId() << " : " << constit.status() << " : " << constit.daughter1() << " : " << constit.daughter2() << endl;
        //     }
        // }

    }

    if (!(njet_min_sel->passes(event) || ngenjet_min_sel->passes(event))) return false;

    if (PRINTOUT || printout_event_sel) printMuons(*event.muons, "Precleaning");
    if (PRINTOUT || printout_event_sel) printElectrons(*event.electrons, "Precleaning");
    if (PRINTOUT || printout_event_sel) printJets(*event.jets, "Precleaning");

    // Gen-level HT cut if necessary
    // -------------------------------------------------------------------------
    float genHT = 0;
    if (!event.isRealData) {
        genHT = calcGenHT(*(event.genparticles));
    }
    if ((htMax > 0) && (genHT > htMax)) { return false; }

    // Gen-level cut (Herwig samples)
    // -------------------------------------------------------------------------
    float qScale = event.genInfo->qScale();
    // if ((qScaleMax>0) || (qScaleMin>0)) {
    //     if ((qScaleMax > 0) && (qScale > qScaleMax)) return false;
    //     if ((qScaleMin > 0) && (qScale < qScaleMin)) return false;
    // }

    float partonKt = calcJetKt(*event.genparticles);
    if ((partonKtMax>0) || (partonKtMin>0)) {
        if ((partonKtMax > 0) && (partonKt > partonKtMax)) return false;
        if ((partonKtMin > 0) && (partonKt < partonKtMin)) return false;
    }

    // Common things
    // -------------------------------------------------------------------------
    // Note that we only care about this for reco-specific bits,
    // not gen-specific (only false if fails MET filters)
    bool passCommonRecoSetup = common_setup->process(event);
    recojet_setup->process(event);

    if (PRINTOUT || printout_event_sel) printMuons(*event.muons);
    if (PRINTOUT || printout_event_sel) printElectrons(*event.electrons);
    // if (PRINTOUT || printout_event_sel) printGenParticles(*event.genparticles);

    // Get Gen muons
    // -------------------------------------------------------------------------
    std::vector<GenParticle> goodGenMuons = getGenMuons(event.genparticles, Cuts::gen_muon_pt_min, Cuts::muon_eta_max);
    // printGenParticles(goodGenMuons);
    event.set(genmuons_handle, std::move(goodGenMuons));

    // Get good GenJets, store in event
    // -------------------------------------------------------------------------
    genJet_selector->process(event);
    // Need these as loosest possible requirement to run reco- or gen-specific bits
    bool hasGenJets = ngenjet_good_sel->passes(event);

    if (!(njet_min_sel->passes(event) || hasGenJets)) return false;

    // MC-specific parts like reweighting for SF, for muR/F scale, etc
    mc_reweight->process(event); // also responsible for setting gen weight, so do after scale variations

    // Cuts to throw away high-weight events from lower pT bins
    // (e.g. where leading jet actually PU jet)
    // -------------------------------------------------------------------------
    // 1. Cut on pt/genHt to avoid weird events
    if (genHT > 0 && (njet_min_sel->passes(event) && ((event.jets->at(0).pt() / genHT) > 2))) { return false; }
    if (genHT > 0 && (hasGenJets && ((event.get(genjets_handle)[0].pt() / genHT) > 2))) { return false; }

    // if (jetKt > 0 && (hasGenJets && ((event.get(genjets_handle)[0].pt() / genHT) > 2))) { return false; }

    if (njet_min_sel->passes(event) && ((event.jets->at(0).pt() / qScale) > 2)) { return false; }

    float PU_pThat = event.genInfo->PU_pT_hat_max();
    if (genHT > 0 && njet_min_sel->passes(event) && ((PU_pThat / genHT) > 2)) { return false; }

    // cout << "*** EVENT:" << endl;
    // cout << "genHT: " << genHT << endl;
    // cout << "qScale: " << qScale << endl;
    // cout << "PU_pThat: " << PU_pThat << endl;
    // cout << "pdf_scalePDF: " << event.genInfo->pdf_scalePDF() << endl;
    // cout << "weight: " << event.weight << endl;
    // if (njet_min_sel->passes(event)) cout << "jet1pt: " << event.jets->at(0).pt() << endl;

    // 2. Check event weight is sensible based on pthat - but isn't always available
    if (event.genInfo->binningValues().size() > 0) {
        // double ptHat = event.genInfo->binningValues().at(0); // yes this is correct. no idea why
        // cout << "ptHAt: " << ptHat << endl;
        // if (hasRecoJets && (event.jets->at(0).pt() / ptHat > 2)) return false;
        // if (hasGenJets && (event.get(genjets_handle)[0].pt() / ptHat > 2)) return false;
    }

    genjet_hists->fill(event);

    // RECO PART
    if (PRINTOUT || printout_event_sel) printJets(*event.jets, "Original jets");

    // Get matching GenJets for reco jets
    // -------------------------------------------------------------------------
    jetMatcherPtOrdered->process(event);

    // if (PRINTOUT || printout_event_sel) printJets(*event.jets, "Matched Jets");
    // if (PRINTOUT || printout_event_sel) printGenJets(event.get(genjets_handle), "GoodGenJets");
    if (PRINTOUT || printout_event_sel) printJetsWithParts(*event.jets, event.pfparticles, "Matched Jets");
    if (PRINTOUT || printout_event_sel) printGenJetsWithParts(event.get(genjets_handle), event.genparticles, "GoodGenJets");

    // Apply tracking SFs, but only after JECs, etc applied
    // - we want to use the original jet pT
    tracking_eff->process(event);

    // Apply Jet PF ID since the jet constituents changes
    jet_pf_id->process(event);

    if (PRINTOUT || printout_event_sel) printJetsWithParts(*event.jets, event.pfparticles, "Matched Jets after tracking SF");

    // Do AFTER all things that could affect number of jets e.g. track SF, IDs
    bool hasRecoJets = njet_min_sel->passes(event) && passCommonRecoSetup; // commonReco bit here as common for all reco parts

    // We need recojets and/or genjets (want both fakes and miss-recos)
    if (!(hasRecoJets || hasGenJets)) return false;

    // Calculate lambda vars for reco jets for Z+jets
    // At this point, all objects should have had all necessary corrections, filtering, etc
    // -------------------------------------------------------------------------
    jetLambdaCreatorPtSorted->process(event);
    genjetLambdaCreatorPtSorted->process(event);

    // Do Z+Jet hists & selection
    // -------------------------------------------------------------------------
    if (isZPlusJets) {
        if (DO_KINEMATIC_HISTS) zplusjets_gen_hists->fill(event);

        pass_zpj_gen = zplusjets_gen_sel->passes(event);
        event.set(pass_zpj_gen_sel_handle, pass_zpj_gen);

        bool found_reco_z = zFinder->process(event);

        if (pass_zpj_gen) {
            event.set(pt_binning_gen_handle, event.get(genjets_handle)[0].pt());
            if (found_reco_z) zplusjets_sel_passGen->passes(event); // just to plot cutflow, need the if since it uses handle internally
        }

        if (hasRecoJets) {
            event.set(pt_binning_reco_handle, event.jets->at(0).pt());
            // flav-specific preselection hists, useful for optimising selection
            uint flav1 = event.jets->at(0).flavor();
            if (found_reco_z) {
                if (DO_KINEMATIC_HISTS) zplusjets_hists_presel->fill(event);
                // if (zplusjets_presel->passes(event)) {
                    // if (DO_FLAVOUR_HISTS) {
                    //     if (flav1 == PDGID::GLUON) {
                    //         zplusjets_hists_presel_g->fill(event);
                    //     } else if (flav1 > PDGID::UNKNOWN && flav1 < PDGID::CHARM_QUARK) {
                    //         zplusjets_hists_presel_q->fill(event);
                    //     } else if (flav1 == PDGID::UNKNOWN) {
                    //         zplusjets_hists_presel_unknown->fill(event);
                    //     }
                    // }

                    pass_zpj_reco = zplusjets_sel->passes(event);
                    event.set(pass_zpj_sel_handle, pass_zpj_reco);
                    if (pass_zpj_reco) {
                        genjet_hists_passZpJReco->fill(event);
                        if (DO_KINEMATIC_HISTS) {
                            zplusjets_gen_sel_passReco->passes(event); // this plots gen cutflow as well
                            zplusjets_hists->fill(event);
                        }

                        if (DO_LAMBDA_HISTS) {
                            zplusjets_qg_hists->fill(event);
                            zplusjets_qg_hists_groomed->fill(event);
                        }

                        if (DO_FLAVOUR_HISTS) {
                            if (flav1 == PDGID::GLUON) {
                                zplusjets_hists_g->fill(event);
                            } else if (flav1 > PDGID::UNKNOWN && flav1 < PDGID::CHARM_QUARK) {
                                zplusjets_hists_q->fill(event);
                            } else if (flav1 == PDGID::UNKNOWN) {
                                zplusjets_hists_unknown->fill(event);
                            }
                        }
                    }
                // }
            }
        }

        // Do unfolding hists
        // ---------------------------------------------------------------------
        if (DO_UNFOLD_HISTS && (pass_zpj_reco || pass_zpj_gen)) {
            zplusjets_qg_unfold_hists->fill(event);
            zplusjets_qg_unfold_hists_groomed->fill(event);
        }

        // Do pu-binned hists
        // ---------------------------------------------------------------------
        if (DO_PU_BINNED_HISTS && pass_zpj_reco) {
            for (uint i=0; i<sel_pu_binned.size(); i++) {
                if (sel_pu_binned.at(i)->passes(event)) {
                    zplusjets_qg_hists_pu_binned.at(i)->fill(event);
                    zplusjets_qg_hists_groomed_pu_binned.at(i)->fill(event);
                }
            }
        }

    } else {
        // Do DiJet hists & selection
        // ---------------------------------------------------------------------
        // For dijet, we sort our leading 2 jets by eta, and use the largest and
        // smallest abs(eta) jets separately for unfolding etc

        // This is pretty horrible - need to run the Lambda bundle creators NO MATTER
        // wether we have enough jets or not. This is because we need to call event.set()
        // anyway, otherwise e.g. unfolding Hists module dies
        //
        // So we do eta sorting, and jet lambda, then do all the hist filling,
        // since they might use these handles.

        if (DO_KINEMATIC_HISTS) dijet_gen_hists->fill(event);

        if (hasRecoJets) {
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

            // Do Reco selection
            pass_dj_reco = dijet_sel_tighter->passes(event);
            event.set(pass_dj_sel_handle, pass_dj_reco);

        }

        if (hasGenJets) {
            // Save forward/central genjets to own handles
            std::vector<GenJetWithParts> leadingGenJets(event.get(genjets_handle).begin(), event.get(genjets_handle).begin()+2);
            if (leadingGenJets.size() != 2) {
                throw std::runtime_error("Slicing genjets gone wrong!");
            }
            double ave_pt = 0.5*(leadingGenJets[0].pt() + leadingGenJets[1].pt());
            event.set(pt_binning_gen_handle, ave_pt);
            sort_by_y(leadingGenJets);
            std::vector<GenJetWithParts> forwardGenJet = {leadingGenJets[0]};
            std::vector<GenJetWithParts> centralGenJet = {leadingGenJets[1]};
            event.set(dijet_gen_forward_handle, forwardGenJet);
            event.set(dijet_gen_central_handle, centralGenJet);

            // Do Gen selection
            pass_dj_gen = dijet_gen_sel->passes(event);
            event.set(pass_dj_gen_sel_handle, pass_dj_gen);
        }

        // Do genjet matching specifically for forward/central
        jetMatcherForward->process(event);
        jetMatcherCentral->process(event);

        // Calculate lambda vars for genjets for dijets
        // These will be used in various histogram classes
        // At this point, all objects should have had all necessary corrections, filtering, etc
        // Have to do outside of any if(), because we always need it to run event.set()
        // otherwise unfolding module dies
        jetLambdaCreatorForward->process(event);
        jetLambdaCreatorCentral->process(event);

        genjetLambdaCreatorForward->process(event);
        genjetLambdaCreatorCentral->process(event);

        if (hasRecoJets) {
            // flav-specific preselection hists, useful for optimising selection
            uint flav1 = event.jets->at(0).flavor();
            uint flav2 = event.jets->at(1).flavor();

            // Fill hists
            if (DO_KINEMATIC_HISTS) dijet_hists_presel->fill(event);
            // if (DO_FLAVOUR_HISTS) {
            //     if (flav1 == PDGID::GLUON) {
            //         if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
            //             dijet_hists_presel_gq->fill(event);
            //         } else if (flav2 == PDGID::GLUON) {
            //             dijet_hists_presel_gg->fill(event);
            //         } else if (flav2 == PDGID::UNKNOWN) {
            //             dijet_hists_presel_g_unknown->fill(event);
            //         }
            //     } else if (flav1 > PDGID::UNKNOWN && flav1 < PDGID::CHARM_QUARK) {
            //         if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
            //             dijet_hists_presel_qq->fill(event);
            //         } else if (flav2 == PDGID::GLUON) {
            //             dijet_hists_presel_qg->fill(event);
            //         } else if (flav2 == PDGID::UNKNOWN) {
            //             dijet_hists_presel_q_unknown->fill(event);
            //         }
            //     } else if (flav1 == PDGID::UNKNOWN) {
            //         if (flav2 == PDGID::GLUON) {
            //             dijet_hists_presel_unknown_g->fill(event);
            //         } else if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
            //             dijet_hists_presel_unknown_q->fill(event);
            //         } else if (flav2 == PDGID::UNKNOWN) {
            //             dijet_hists_presel_unknown_unknown->fill(event);
            //         }
            //     }
            // }

            if (pass_dj_reco) {
                dijet_gen_sel_passReco->passes(event); // this plots cutflow as well - only if we pass dijet reco selection
                genjet_hists_passDijetReco->fill(event);
            }

            if (pass_dj_gen) {
                // to plot reco cutflow when passGen==true
                dijet_sel_tighter_passGen->passes(event);
            }

            bool standard_sel = dijet_sel->passes(event);
            bool tight_sel = dijet_sel_tighter->passes(event); // aka passReco
            if (DO_KINEMATIC_HISTS) {
                if (standard_sel) {
                    dijet_hists->fill(event);
                    dijet_qg_hists->fill(event);
                }
                if (tight_sel) {
                    dijet_hists_tighter->fill(event);
                    dijet_qg_hists_tighter->fill(event);
                    dijet_qg_hists_central_tighter->fill(event);
                    dijet_qg_hists_forward_tighter->fill(event);
                    dijet_qg_hists_central_tighter_groomed->fill(event);
                    dijet_qg_hists_forward_tighter_groomed->fill(event);
                }
            }

            if (DO_LAMBDA_HISTS) {
                if (standard_sel) {
                    dijet_qg_hists->fill(event);
                }
                if (tight_sel) {
                    dijet_qg_hists_tighter->fill(event);
                    dijet_qg_hists_central_tighter->fill(event);
                    dijet_qg_hists_forward_tighter->fill(event);
                    dijet_qg_hists_central_tighter_groomed->fill(event);
                    dijet_qg_hists_forward_tighter_groomed->fill(event);
                }
            }

            // do eta-sorted dijet hists (where we need both jets)

            if (DO_KINEMATIC_HISTS && standard_sel) {
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

                // flav-specific eta-orderd hists, useful for optimising selection
                flav1 = event.jets->at(0).flavor();
                flav2 = event.jets->at(1).flavor();

                // cout << "eta order plots: " << endl;
                // if (event.jets->size() > 0) cout << "jet[0].genjet_index: " << event.jets->at(0).genjet_index() << endl;
                // if (event.jets->size() > 1) cout << "jet[1].genjet_index: " << event.jets->at(1).genjet_index() << endl;

                if (DO_FLAVOUR_HISTS) {
                    if (flav1 == PDGID::GLUON) {
                        if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
                            dijet_hists_eta_ordered_gq->fill(event);
                        } else if (flav2 == PDGID::GLUON) {
                            dijet_hists_eta_ordered_gg->fill(event);
                        }
                    } else if (flav1 > PDGID::UNKNOWN && flav1 < PDGID::CHARM_QUARK) {
                        if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
                            dijet_hists_eta_ordered_qq->fill(event);
                        } else if (flav2 == PDGID::GLUON) {
                            dijet_hists_eta_ordered_qg->fill(event);
                        }
                    }
                }
            }
        } // end hasRecoJets

        // Do unfolding hists
        // ---------------------------------------------------------------------
        if (DO_UNFOLD_HISTS && (pass_dj_reco || pass_dj_gen)) {
            dijet_qg_unfold_hists_central_tighter->fill(event);
            dijet_qg_unfold_hists_forward_tighter->fill(event);
            dijet_qg_unfold_hists_central_tighter_groomed->fill(event);
            dijet_qg_unfold_hists_forward_tighter_groomed->fill(event);
        }

        // Do pu-binned hists
        // ---------------------------------------------------------------------
        if (DO_PU_BINNED_HISTS && pass_dj_reco) {
            for (uint i=0; i<sel_pu_binned.size(); i++) {
                if (sel_pu_binned.at(i)->passes(event)) {
                    dijet_qg_hists_central_pu_binned.at(i)->fill(event);
                    dijet_qg_hists_forward_pu_binned.at(i)->fill(event);
                    dijet_qg_hists_central_groomed_pu_binned.at(i)->fill(event);
                    dijet_qg_hists_forward_groomed_pu_binned.at(i)->fill(event);
                }
            }
        }

    /*
        // Do high pt jet version
        // ---------------------------------------------------------------------
        // Both jets must pass much higher pt threshold
        // don't need to do a Z+jets version as only care about leading jet.
        float ptCut = 500;
        if (event.jets->at(0).pt() < ptCut) return false;
        flav2 = 99999999;
        if ((event.jets->size() > 1) && (event.jets->at(1).pt() > ptCut)) {
            flav2 = event.jets->at(1).flavor();
            dijet_hists_presel_highPt->fill(event);

            // flav-specific preselection hists, useful for optimising selection
            if (flav1 == PDGID::GLUON) {
                if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
                    dijet_hists_presel_gq_highPt->fill(event);
                } else if (flav2 == PDGID::GLUON) {
                    dijet_hists_presel_gg_highPt->fill(event);
                } else if (flav2 == PDGID::UNKNOWN) {
                    dijet_hists_presel_g_unknown_highPt->fill(event);
                }
            } else if (flav1 > PDGID::UNKNOWN && flav1 < PDGID::CHARM_QUARK) {
                if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
                    dijet_hists_presel_qq_highPt->fill(event);
                } else if (flav2 == PDGID::GLUON) {
                    dijet_hists_presel_qg_highPt->fill(event);
                } else if (flav2 == PDGID::UNKNOWN) {
                    dijet_hists_presel_q_unknown_highPt->fill(event);
                }
            } else if (flav1 == PDGID::UNKNOWN) {
                if (flav2 == PDGID::GLUON) {
                    dijet_hists_presel_unknown_g_highPt->fill(event);
                } else if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
                    dijet_hists_presel_unknown_q_highPt->fill(event);
                } else if (flav2 == PDGID::UNKNOWN) {
                    dijet_hists_presel_unknown_unknown_highPt->fill(event);
                }
            }

            dj_highPt = dijet_sel->passes(event);
            if (dj_highPt && flav2 < 100) { // flav2 only sensible if passed pt cut
                dijet_hists_highPt->fill(event);
                dijet_qg_hists_highPt->fill(event);
            }
        }
    */
    } // end if isZPlusJets

    if (pass_zpj_reco && pass_dj_reco) {
        nOverlapEvents++;
        cout << "Warning: event (runid, eventid) = ("  << event.run << ", "
             << event.event << ") passes both Z+jets and Dijet criteria ("
             << nOverlapEvents << " total)" << endl;
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
    if (pass_zpj_reco || pass_zpj_gen) nZPJEvents++;
    if (pass_dj_reco || pass_dj_gen) nDijetEvents++;
    if (pass_zpj_reco || pass_dj_reco || pass_zpj_gen || pass_dj_gen) nPassEvents++;
    return pass_zpj_reco || pass_dj_reco || pass_zpj_gen || pass_dj_gen;
}


void QGAnalysisMCModule::endInputData(){
    cout << "Summary stats: " << endl;
    cout << " # dijet reco||gen: " << nDijetEvents << endl;
    cout << " # z+j reco||gen: " << nZPJEvents << endl;
    cout << " # dijet||z+j reco||gen: " << nPassEvents << endl;
    cout << " # dijet && z+j reco: " << nOverlapEvents << endl;
}



/**
 * Select gen muons from all genparticles, that have some minimum pt and maximum eta
 */
std::vector<GenParticle> QGAnalysisMCModule::getGenMuons(std::vector<GenParticle> * genparticles, float pt_min, float eta_max) {
    std::vector<GenParticle> muons;
    // Do in reverse order to pick up most evolved muons first
    for (auto itr = genparticles->rbegin(); itr != genparticles->rend(); ++itr){
        // We check to see if we already have a very similar, but not exact, muon
        // since the MC "evolves" the particle and slightly changes pt/eta/phi
        // status check may not be reliable for e.g. herwig
        bool alreadyFound = std::any_of(muons.begin(), muons.end(), [&itr] (const GenParticle & mtr) { return deltaR(*itr, mtr) < 0.05 && itr->charge() == mtr.charge(); });
        if ((abs(itr->pdgId()) == PDGID::MUON) && (itr->status() == 1) && (itr->pt() > pt_min) && (fabs(itr->eta()) < eta_max) && !alreadyFound) {
        // if ((abs(itr->pdgId()) == PDGID::MUON) && (itr->pt() > pt_min) && (fabs(itr->eta()) < eta_max)) {
            muons.push_back(*itr);
        }
    }
    sort_by_pt(muons);
    return muons;
}

/**
 * Select reco jets that have a matching GenJet within some DR
 * Also stores index of matching GenJets in the passed Jet collection
 * Will take the closest matching GenJet as the match, provided it is within drMax.
 * uniqueMatch controls whether matching GenJets must be unique (i.e 2 reco jets can't match the same GenJet)
 *
 */
std::vector<Jet> QGAnalysisMCModule::getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax, bool uniqueMatch) {
    std::vector<Jet> goodJets;
    std::vector<uint> matchedIndices;
    for (auto & jtr: *jets) {
        double minDR = 9999.;
        int matchInd = -1; // sensible default - not 0!

        for (uint gjInd=0; gjInd < genjets->size(); gjInd++) {
            // If we want unique matches and we've already matched then skip this genjet
            if (uniqueMatch && std::find(matchedIndices.begin(), matchedIndices.end(), gjInd) != matchedIndices.end())
                continue;

            const auto genjtr = genjets->at(gjInd);
            auto thisDR = deltaRUsingY(jtr, genjtr);
            if (thisDR < drMax && thisDR < minDR) {
                matchInd = gjInd;
                minDR = thisDR;
            }
        }

        jtr.set_genjet_index(matchInd);
        if (matchInd > -1) {
            jtr.set_genjet_index(matchInd);
            goodJets.push_back(jtr);
            matchedIndices.push_back(matchInd);
            // cout << "Found a match at index " << matchInd << " with dr " << minDR << endl;
            // cout << "RECO pt/eta/phi: " << jtr.pt() << " : " << jtr.eta() << " : " << jtr.phi() << endl;
            // cout << "GEN pt/eta/phi: " << genjets->at(matchInd).pt() << " : " << genjets->at(matchInd).eta() << " : " << genjets->at(matchInd).phi() << endl;
        } else {
            if (PRINTOUT) cout << "Cannot find match for jet " << jtr.pt() << " : " << jtr.eta() << " : " << jtr.phi() << endl;
        }
    }
    return goodJets;
}


// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisMCModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisMCModule)

}

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
#include "UHH2/QGAnalysis/include/QGAnalysisUnfoldHists.h"
#include "UHH2/QGAnalysis/include/QGAnalysisZPlusJetsHists.h"
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
    std::vector<GenJetWithParts> getGenJets(std::vector<GenJetWithParts> * jets, std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=1.5, float lepton_overlap_dr=0.2);
    std::vector<GenParticle> getGenMuons(std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=2.5);
    std::vector<Jet> getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax=0.8, bool uniqueMatch=true);

private:

    std::unique_ptr<GeneralEventSetup> common_setup;
    std::unique_ptr<MCReweighting> mc_reweight;
    std::unique_ptr<MCScaleVariation> mc_scalevar;

    std::unique_ptr<QGAnalysisJetLambda> jetLambdaCreator, jetChargedLambdaCreator;
    std::unique_ptr<QGAnalysisGenJetLambda> genjetLambdaCreator, genjetChargedLambdaCreator;

    // Reco selections/hists
    std::unique_ptr<ZFinder> zFinder;
    std::unique_ptr<Selection> njet_sel, ngenjet_sel, ngenjet_good_sel, zplusjets_sel, zplusjets_presel, dijet_sel, dijet_sel_tighter;

    std::unique_ptr<Hists> zplusjets_hists_presel, zplusjets_hists;
    std::unique_ptr<Hists> zplusjets_hists_presel_q, zplusjets_hists_presel_g, zplusjets_hists_presel_unknown;
    std::unique_ptr<Hists> zplusjets_qg_hists;
    std::unique_ptr<Hists> zplusjets_qg_unfold_hists;

    std::unique_ptr<Hists> dijet_hists_presel, dijet_hists, dijet_hists_tighter;
    std::unique_ptr<Hists> dijet_hists_presel_gg, dijet_hists_presel_qg, dijet_hists_presel_gq, dijet_hists_presel_qq;
    std::unique_ptr<Hists> dijet_hists_presel_q_unknown, dijet_hists_presel_g_unknown, dijet_hists_presel_unknown_unknown, dijet_hists_presel_unknown_q, dijet_hists_presel_unknown_g;
    std::unique_ptr<Hists> dijet_qg_hists, dijet_qg_hists_tighter;
    std::unique_ptr<Hists> dijet_qg_unfold_hists_tighter;

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
    std::vector< std::unique_ptr<Hists> > zplusjets_qg_hists_pu_binned;
    std::vector< std::unique_ptr<Hists> > dijet_qg_hists_pu_binned;

    // Gen selection
    std::unique_ptr<Selection> zplusjets_gen_sel, dijet_gen_sel;

    // Theory selections/hists
    std::unique_ptr<Selection> zplusjets_theory_sel, dijet_theory_sel;
    std::unique_ptr<Hists> zplusjets_hists_theory, dijet_hists_theory;

    // for sweeping over ptMin
    // std::vector<float> theory_pt_bins = {50, 100, 200, 400, 800};
    // std::vector< std::unique_ptr<Selection> > zplusjets_theory_sel_pt_binned;
    // std::vector< std::unique_ptr<Selection> > dijet_theory_sel_pt_binned;
    // std::vector< std::unique_ptr<Hists> > zplusjets_hists_theory_pt_binned;
    // std::vector< std::unique_ptr<Hists> > dijet_hists_theory_pt_binned;

    Event::Handle<std::vector<GenJetWithParts>> genjets_handle;
    Event::Handle<std::vector<GenParticle>> genmuons_handle;
    Event::Handle<double> gen_weight_handle;
    Event::Handle<bool> pass_zpj_sel_handle, pass_zpj_gen_sel_handle, pass_dj_sel_handle, pass_dj_gen_sel_handle;

    float jetRadius;
    float htMax;

    const bool DO_PU_BINNED_HISTS = true;

    std::unique_ptr<EventNumberSelection> event_sel;

    std::string zLabel;

    int NJETS_ZPJ, NJETS_DIJET;
    uint nOverlapEvents;
};


QGAnalysisMCModule::QGAnalysisMCModule(Context & ctx){
    cout << "Running analysis module" << endl;

    htMax =  boost::lexical_cast<float>(ctx.get("maxHT", "-1"));

    string jet_cone = ctx.get("JetCone", "AK4");
    string pu_removal = ctx.get("PURemoval", "CHS");
    if (pu_removal != "CHS" && pu_removal != "PUPPI") {
        throw runtime_error("Only PURemoval == CHS, PUPPI supported for now");
    }
    jetRadius = get_jet_radius(jet_cone);

    cout << "Running with jet cone: " << jet_cone << endl;
    cout << "Running with PUS: " << pu_removal << endl;

    // FIXME: get everything from ctx not extra args
    common_setup.reset(new GeneralEventSetup(ctx, pu_removal, jet_cone, jetRadius));
    mc_reweight.reset(new MCReweighting(ctx));
    mc_scalevar.reset(new MCScaleVariation(ctx));

    genjets_handle = ctx.declare_event_output< std::vector<GenJetWithParts> > ("GoodGenJets");
    genmuons_handle = ctx.declare_event_output< std::vector<GenParticle> > ("GoodGenMuons");
    gen_weight_handle = ctx.declare_event_output<double>("gen_weight");

    pass_zpj_sel_handle = ctx.declare_event_output<bool> ("ZPlusJetsSelection");
    pass_zpj_gen_sel_handle = ctx.declare_event_output<bool> ("ZPlusJetsGenSelection");
    pass_dj_sel_handle = ctx.declare_event_output<bool> ("DijetSelection");
    pass_dj_gen_sel_handle = ctx.declare_event_output<bool> ("DijetGenSelection");

    // Event Selections
    NJETS_ZPJ = 1;
    NJETS_DIJET = 2;
    njet_sel.reset(new NJetSelection(NJETS_ZPJ));
    ngenjet_sel.reset(new NGenJetWithPartsSelection(NJETS_ZPJ));
    ngenjet_good_sel.reset(new NGenJetWithPartsSelection(NJETS_ZPJ, 10000, boost::none, genjets_handle));

    // Lambda calculators
    bool doPuppi = (pu_removal == "PUPPI");
    int maxNJets = max(NJETS_ZPJ, NJETS_DIJET);
    float recoConstitPtMin = 1.;

    // FIXME: get stuff from ctx not extra args
    jetLambdaCreator.reset(new QGAnalysisJetLambda(ctx, jetRadius, maxNJets, doPuppi, PtEtaCut(recoConstitPtMin, 5.), "jets", "JetLambdas"));
    jetChargedLambdaCreator.reset(new QGAnalysisJetLambda(ctx, jetRadius, maxNJets, doPuppi, AndId<PFParticle>(PtEtaCut(recoConstitPtMin, 5.), ChargedCut()), "jets", "JetChargedLambdas"));
    // do Lambda for more genjets, since we sometimes are interested in reco-genjet matches for genjets >= #3
    genjetLambdaCreator.reset(new QGAnalysisGenJetLambda(ctx, jetRadius, 5, PtEtaCut(0., 5.), "GoodGenJets", "GoodGenJetLambdas"));
    genjetChargedLambdaCreator.reset(new QGAnalysisGenJetLambda(ctx, jetRadius, 5, ChargedCut(), "GoodGenJets", "GoodGenJetChargedLambdas"));

    // Setup for systematics
    // FIXME put all this inside the ctor as it has ctx!
    std::string neutralHadronShift = ctx.get("neutralHadronShift", "nominal");
    float neutralHadronShiftAmount = 0.1;
    if (neutralHadronShift == "nominal") {
        // pass
    } else if (neutralHadronShift == "up") {
        jetLambdaCreator->set_neutral_hadron_shift(1, neutralHadronShiftAmount);
        jetChargedLambdaCreator->set_neutral_hadron_shift(1, neutralHadronShiftAmount);
    } else if (neutralHadronShift == "down") {
        jetLambdaCreator->set_neutral_hadron_shift(-1, neutralHadronShiftAmount);
        jetChargedLambdaCreator->set_neutral_hadron_shift(-1, neutralHadronShiftAmount);
    } else {
        throw runtime_error("neutralHadronShift must be nominal, up, or down");
    }

    std::string photonShift = ctx.get("photonShift", "nominal");
    float photonShiftAmount = 0.01;
    if (photonShift == "nominal") {
        // pass
    } else if (photonShift == "up") {
        jetLambdaCreator->set_photon_shift(1, photonShiftAmount);
        jetChargedLambdaCreator->set_photon_shift(1, photonShiftAmount);
    } else if (photonShift == "down") {
        jetLambdaCreator->set_photon_shift(-1, photonShiftAmount);
        jetChargedLambdaCreator->set_photon_shift(-1, photonShiftAmount);
    } else {
        throw runtime_error("photonShift must be nominal, up, or down");
    }

    zLabel = "zMuonCand";
    zFinder.reset(new ZFinder(ctx, "muons", zLabel, ctx.get("z_reweight_file", "")));

    // For gen selection, do as per reco selection but looser by this factor
    float mcSelFactor = 1.25;

    // Z+JETS selection
    float mu1_pt = 26.; // muon pt cut comes from the cleaner in GeneralEventSetup
    float mu2_pt = 26.;
    float mZ_window = 20.;
    float dphi_jet_z_min = 2.0;
    float second_jet_frac_max_zpj = 0.3;
    zplusjets_sel.reset(new ZplusJetsSelection(ctx, zLabel, mu1_pt, mu2_pt, mZ_window, dphi_jet_z_min, second_jet_frac_max_zpj));

    zplusjets_gen_sel.reset(new ZplusJetsGenSelection(ctx, mu1_pt/mcSelFactor, mu2_pt/mcSelFactor, mZ_window*mcSelFactor, dphi_jet_z_min/mcSelFactor, second_jet_frac_max_zpj*mcSelFactor, "GoodGenJets", "GoodGenMuons"));

    // Preselection for Z+J - only 2 muons to reco Z
    dphi_jet_z_min = 0.;
    second_jet_frac_max_zpj = 999.;
    zplusjets_presel.reset(new ZplusJetsSelection(ctx, zLabel, mu1_pt, mu2_pt, mZ_window, dphi_jet_z_min, second_jet_frac_max_zpj));

    // DIJET selection
    float dphi_min = 2.;
    float second_jet_frac_max_dj = 10.94;
    float jet_asym_max = 0.3;
    bool ss_eta = false;
    float deta = 12;
    float sumEta = 10.;
    dijet_sel.reset(new DijetSelection(dphi_min, second_jet_frac_max_dj, 1000, ss_eta, deta, sumEta));
    dijet_sel_tighter.reset(new DijetSelection(dphi_min, second_jet_frac_max_dj, jet_asym_max, ss_eta, deta, sumEta));

    dijet_gen_sel.reset(new DijetGenSelection(ctx, dphi_min*mcSelFactor, second_jet_frac_max_dj*mcSelFactor, jet_asym_max*mcSelFactor, ss_eta, deta*mcSelFactor, sumEta*mcSelFactor, "GoodGenJets"));

    // zplusjets_theory_sel.reset(new ZplusJetsTheorySelection(ctx));
    // dijet_theory_sel.reset(new DijetTheorySelection(ctx));

    // Hists
    zplusjets_hists_presel.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel", zLabel));
    // preselection hists, if jet is quark, or gluon
    zplusjets_hists_presel_q.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_q", zLabel));
    zplusjets_hists_presel_g.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_g", zLabel));
    zplusjets_hists_presel_unknown.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets_Presel_unknown", zLabel));
    zplusjets_hists.reset(new QGAnalysisZPlusJetsHists(ctx, "ZPlusJets", zLabel));
    zplusjets_qg_hists.reset(new QGAnalysisHists(ctx, "ZPlusJets_QG", NJETS_ZPJ, "zplusjets", "ZPlusJetsSelection", "ZPlusJetsGenSelection"));
    zplusjets_qg_unfold_hists.reset(new QGAnalysisUnfoldHists(ctx, "ZPlusJets_QG_Unfold", NJETS_ZPJ, "zplusjets", "ZPlusJetsSelection", "ZPlusJetsGenSelection"));

    std::string binning_method = "ave";
    dijet_hists_presel.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel", binning_method));
    // preselection hiss, if both gluon jets, one gluon, or both quark, or one or both unknown
    dijet_hists_presel_gg.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_gg", binning_method));
    dijet_hists_presel_qg.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_qg", binning_method));
    dijet_hists_presel_gq.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_gq", binning_method));
    dijet_hists_presel_qq.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_qq", binning_method));
    dijet_hists_presel_unknown_q.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_q", binning_method));
    dijet_hists_presel_unknown_g.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_g", binning_method));
    dijet_hists_presel_unknown_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_unknown_unknown", binning_method));
    dijet_hists_presel_q_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_q_unknown", binning_method));
    dijet_hists_presel_g_unknown.reset(new QGAnalysisDijetHists(ctx, "Dijet_Presel_g_unknown", binning_method));
    dijet_hists.reset(new QGAnalysisDijetHists(ctx, "Dijet", binning_method));
    dijet_hists_tighter.reset(new QGAnalysisDijetHists(ctx, "Dijet_tighter", binning_method));
    dijet_qg_hists.reset(new QGAnalysisHists(ctx, "Dijet_QG", NJETS_DIJET, "dijet", "DijetSelection", "DijetGenSelection"));
    dijet_qg_hists_tighter.reset(new QGAnalysisHists(ctx, "Dijet_QG_tighter", NJETS_DIJET, "dijet", "DijetSelection", "DijetGenSelection"));
    dijet_qg_unfold_hists_tighter.reset(new QGAnalysisUnfoldHists(ctx, "Dijet_QG_Unfold_tighter", NJETS_DIJET, "dijet", "DijetSelection", "DijetGenSelection"));

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


    // zplusjets_hists_theory.reset(new QGAnalysisTheoryHists(ctx, "ZPlusJets_genjet", NJETS_ZPJ, "zplusjets"));
    // dijet_hists_theory.reset(new QGAnalysisTheoryHists(ctx, "Dijet_genjet", NJETS_DIJET, "dijet"));

    // for (auto pt : theory_pt_bins) {
    //     std::unique_ptr<Selection> A(new ZplusJetsTheorySelection(ctx, pt));
    //     zplusjets_theory_sel_pt_binned.push_back(std::move(A));
    //     std::unique_ptr<Selection> B(new DijetTheorySelection(ctx, pt));
    //     dijet_theory_sel_pt_binned.push_back(std::move(B));

    //     std::unique_ptr<QGAnalysisTheoryHists> a(new QGAnalysisTheoryHists(ctx, TString::Format("ZPlusJets_genjet_ptMin_%d", int(pt)).Data(), NJETS_ZPJ, "zplusjets"));
    //     zplusjets_hists_theory_pt_binned.push_back(std::move(a));
    //     std::unique_ptr<QGAnalysisTheoryHists> b(new QGAnalysisTheoryHists(ctx, TString::Format("Dijet_genjet_ptMin_%d", int(pt)).Data(), NJETS_DIJET, "dijet"));
    //     dijet_hists_theory_pt_binned.push_back(std::move(b));
    // }

    if (DO_PU_BINNED_HISTS) {
        for (auto puBin : pu_bins) {
            std::unique_ptr<Selection> pu_sel(new NPVSelection(puBin.first, puBin.second));
            sel_pu_binned.push_back(std::move(pu_sel));
            std::unique_ptr<QGAnalysisHists> zpj(new QGAnalysisHists(ctx, TString::Format("ZPlusJets_QG_PU_%d_to_%d", puBin.first, puBin.second).Data(), NJETS_ZPJ, "zplusjets", "ZPlusJetsSelection", "ZPlusJetsGenSelection"));
            zplusjets_qg_hists_pu_binned.push_back(std::move(zpj));
            std::unique_ptr<QGAnalysisHists> dj(new QGAnalysisHists(ctx, TString::Format("Dijet_QG_PU_%d_to_%d", puBin.first, puBin.second).Data(), NJETS_DIJET, "dijet", "DijetSelection", "DijetGenSelection"));
            dijet_qg_hists_pu_binned.push_back(std::move(dj));
        }
    }

    // event_sel.reset(new EventNumberSelection({111}));

    nOverlapEvents = 0;
}


bool QGAnalysisMCModule::process(Event & event) {
    // if (!event_sel->passes(event)) return false;

    double orig_weight = 1.;
    event.set(gen_weight_handle, orig_weight); // need to set this at the start

    if (PRINTOUT) { cout << "-- Event: " << event.event << endl; }
    // cout << "-- Event: " << event.event << endl;

    if (!(njet_sel->passes(event) || ngenjet_sel->passes(event))) return false;

    if (PRINTOUT) printMuons(*event.muons, "Precleaning");
    if (PRINTOUT) printElectrons(*event.electrons, "Precleaning");
    if (PRINTOUT) printJets(*event.jets, "Precleaning");

    // Gen-level HT cut if necessary
    // -------------------------------------------------------------------------
    float genHT = 0;
    if (!event.isRealData) {
        genHT = calcGenHT(*(event.genparticles));
    }
    if ((htMax > 0) && (genHT > htMax)) { return false; }

    // Common things
    // -------------------------------------------------------------------------
    // Note that we only care about this for reco-specific bits,
    // not gen-specific (only false if fails MET filters)
    bool passCommonRecoSetup = common_setup->process(event);
    if (!(njet_sel->passes(event) || ngenjet_sel->passes(event))) return false;
    // MC-specific parts like reweighting for SF, for muR/F scale, etc
    mc_reweight->process(event);
    mc_scalevar->process(event);

    if (PRINTOUT) printMuons(*event.muons);
    if (PRINTOUT) printElectrons(*event.electrons);
    // if (PRINTOUT) printGenParticles(*event.genparticles);

    // Get Gen muons
    // -------------------------------------------------------------------------
    std::vector<GenParticle> goodGenMuons = getGenMuons(event.genparticles, 5., 2.4+(jetRadius/2.));
    event.set(genmuons_handle, std::move(goodGenMuons));

    // Get good GenJets, store in event
    // -------------------------------------------------------------------------
    double genjet_pt_cut = 15.;
    double genjet_eta_cut = 2.4+(jetRadius/2.);
    std::vector<GenJetWithParts> goodGenJets = getGenJets(event.genjets, &event.get(genmuons_handle), genjet_pt_cut, genjet_eta_cut, jetRadius);
    sort_by_pt(goodGenJets);
    event.set(genjets_handle, std::move(goodGenJets));

    // Need these as loosest possible requirement to run reco- or gen-specific bits
    bool hasRecoJets = njet_sel->passes(event) && passCommonRecoSetup; // commonReco bit here as common for all reco parts
    bool hasGenJets = ngenjet_good_sel->passes(event);

    // We need recojets and/or genjets (want both fakes and miss-recos)
    if (!(hasRecoJets || hasGenJets)) return false;

    // Cuts to throw away high-weight events from lower pT bins
    // (e.g. where leading jet actually PU jet)
    // -------------------------------------------------------------------------

    // 1. Cut on pt/genHt to avoid weird events
    if (genHT > 0 && (hasRecoJets && ((event.jets->at(0).pt() / genHT) > 2.5))) { return false; }
    if (genHT > 0 && (hasGenJets && ((event.get(genjets_handle)[0].pt() / genHT) > 2.5))) { return false; }

    // 2. Check event weight is sensible based on pthat - but isn't always available
    if (event.genInfo->binningValues().size() > 0) {
        double ptHat = event.genInfo->binningValues().at(0); // yes this is correct. no idea why
        if (hasRecoJets && (event.jets->at(0).pt() / ptHat > 2)) return false;
        if (hasGenJets && (event.get(genjets_handle)[0].pt() / ptHat > 2)) return false;
    }

    // Determine if good event if leading jet is a true jet or not (ie PU)
    // But really shouldn't need this?
    // bool goodEvent = false;
    // for (const auto genjtr: event.get(genjets_handle)) {
    //     if (deltaR(event.jets->at(0), genjtr) < jetRadius){
    //         goodEvent = true;
    //         break;
    //     }
    // }
    // if (!goodEvent) return false;

    // Theory-specific selection & hists
    // (Following selection etc in paper)
    // -------------------------------------------------------------------------
    // bool zpj_th = zplusjets_theory_sel->passes(event);
    // bool dj_th = dijet_theory_sel->passes(event);

    // if (zpj_th && dj_th) {
    //     cout << "Warning: event (runid, eventid) = ("  << event.run << ", "
    //          << event.event << ") passes both Z+jets and Dijet theory criteria" << endl;
    // }

    // if (zpj_th) {
    //     zplusjets_hists_theory->fill(event);
    // }

    // if (dj_th) {
    //     dijet_hists_theory->fill(event);
    // }

    // ptMin binned hists
    // for (uint i=0; i < theory_pt_bins.size(); i++) {
    //     if (zplusjets_theory_sel_pt_binned.at(i)->passes(event)) {
    //         zplusjets_hists_theory_pt_binned.at(i)->fill(event);
    //     }
    //     if (dijet_theory_sel_pt_binned.at(i)->passes(event)) {
    //         dijet_hists_theory_pt_binned.at(i)->fill(event);
    //     }
    // }

    // RECO PART

    if (PRINTOUT) printJets(*event.jets, "Original jets");

    // Get matching GenJets for reco jets
    // -------------------------------------------------------------------------
    // But we still use the event.jets as all interesting
    std::vector<Jet> goodJets = getMatchedJets(event.jets, &event.get(genjets_handle), jetRadius/2.);
    // std::swap(goodJets, *event.jets); // only save recojets with a match

    // if (PRINTOUT) printJets(*event.jets, "Matched Jets");
    // if (PRINTOUT) printGenJets(event.get(genjets_handle), "GoodGenJets");
    if (PRINTOUT) printJetsWithParts(*event.jets, event.pfparticles, "Matched Jets");
    if (PRINTOUT) printGenJetsWithParts(event.get(genjets_handle), event.genparticles, "GoodGenJets");

    // Save selection flags
    // At this point, all objects should have had all necessary corrections, filtering, etc
    // ------------------------------------------------------------------------------------
    bool pass_zpj_reco(false), pass_dj_reco(false);
    // incase they doesn't get set later
    event.set(pass_zpj_sel_handle, pass_zpj_reco);
    event.set(pass_dj_sel_handle, pass_dj_reco);

    // bool pass_dj_highPt(false);

    bool pass_zpj_gen = zplusjets_gen_sel->passes(event);
    event.set(pass_zpj_gen_sel_handle, pass_zpj_gen);

    bool pass_dj_gen = dijet_gen_sel->passes(event);
    event.set(pass_dj_gen_sel_handle, pass_dj_gen);

    // Calculate lambda vars for jets & genjets
    // These will be used in various histogram classes
    // -------------------------------------------------------------------------
    jetLambdaCreator->process(event);
    jetChargedLambdaCreator->process(event);
    genjetLambdaCreator->process(event);
    genjetChargedLambdaCreator->process(event);

    // Do reco-specific Z+Jet hists & selection
    // -------------------------------------------------------------------------
    if (hasRecoJets) {
        // flav-specific preselection hists, useful for optimising selection
        uint flav1 = event.jets->at(0).flavor();
        if (zFinder->process(event)) {
            zplusjets_hists_presel->fill(event);
            if (zplusjets_presel->passes(event)) {
                if (flav1 == PDGID::GLUON) {
                    zplusjets_hists_presel_g->fill(event);
                } else if (flav1 > PDGID::UNKNOWN && flav1 < PDGID::CHARM_QUARK) {
                    zplusjets_hists_presel_q->fill(event);
                } else if (flav1 == PDGID::UNKNOWN) {
                    zplusjets_hists_presel_unknown->fill(event);
                }

                pass_zpj_reco = zplusjets_sel->passes(event);
                event.set(pass_zpj_sel_handle, pass_zpj_reco);
                if (pass_zpj_reco) {
                    zplusjets_hists->fill(event);
                    zplusjets_qg_hists->fill(event);
                }
            }
        }
    }

    if (pass_zpj_reco || pass_zpj_gen) zplusjets_qg_unfold_hists->fill(event);

    // Do reco-specific DiJet hists & selection
    // -------------------------------------------------------------------------
    if (hasRecoJets) {
        // flav-specific preselection hists, useful for optimising selection
        uint flav1 = event.jets->at(0).flavor();
        uint flav2(99999999);
        if (event.jets->size() > 1) {
            dijet_hists_presel->fill(event);
            flav2 = event.jets->at(1).flavor();
            if (flav1 == PDGID::GLUON) {
                if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
                    dijet_hists_presel_gq->fill(event);
                } else if (flav2 == PDGID::GLUON) {
                    dijet_hists_presel_gg->fill(event);
                } else if (flav2 == PDGID::UNKNOWN) {
                    dijet_hists_presel_g_unknown->fill(event);
                }
            } else if (flav1 > PDGID::UNKNOWN && flav1 < PDGID::CHARM_QUARK) {
                if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
                    dijet_hists_presel_qq->fill(event);
                } else if (flav2 == PDGID::GLUON) {
                    dijet_hists_presel_qg->fill(event);
                } else if (flav2 == PDGID::UNKNOWN) {
                    dijet_hists_presel_q_unknown->fill(event);
                }
            } else if (flav1 == PDGID::UNKNOWN) {
                if (flav2 == PDGID::GLUON) {
                    dijet_hists_presel_unknown_g->fill(event);
                } else if (flav2 > PDGID::UNKNOWN && flav2 < PDGID::CHARM_QUARK) {
                    dijet_hists_presel_unknown_q->fill(event);
                } else if (flav2 == PDGID::UNKNOWN) {
                    dijet_hists_presel_unknown_unknown->fill(event);
                }
            }

            // pass_dj_reco = dijet_sel->passes(event);
            pass_dj_reco = dijet_sel_tighter->passes(event);
            event.set(pass_dj_sel_handle, pass_dj_reco);

            if (dijet_sel->passes(event)) {
                dijet_hists->fill(event);
                dijet_qg_hists->fill(event);
            }
            if (dijet_sel_tighter->passes(event)) {
                dijet_hists_tighter->fill(event);
                dijet_qg_hists_tighter->fill(event);
            }
        }
    }

    if (pass_dj_reco || pass_dj_gen) dijet_qg_unfold_hists_tighter->fill(event);

    // Do pu-binned hists
    // -------------------------------------------------------------------------
    if (DO_PU_BINNED_HISTS) {
        for (uint i=0; i<sel_pu_binned.size(); i++) {
            if (sel_pu_binned.at(i)->passes(event)) {
                if (pass_zpj_reco) zplusjets_qg_hists_pu_binned.at(i)->fill(event);
                if (pass_dj_reco) dijet_qg_hists_pu_binned.at(i)->fill(event);
            }
        }
    }

/*
    // Do high pt jet version
    // -------------------------------------------------------------------------
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

    return pass_zpj_reco || pass_dj_reco || pass_zpj_gen || pass_dj_gen;
}


/**
 * Get GenJets, ignoring those that are basically a lepton, and have some minimum pt and maximum eta.
 */
std::vector<GenJetWithParts> QGAnalysisMCModule::getGenJets(std::vector<GenJetWithParts> * genjets_in, std::vector<GenParticle> * genparticles, float pt_min, float eta_max, float lepton_overlap_dr) {
    std::vector<GenJetWithParts> genjets_out;
    for (const auto jet : *genjets_in) {
        bool found = (std::find(genjets_out.begin(), genjets_out.end(), jet) != genjets_out.end()); // avoid duplicate genjets
        // avoid jets that are just leptons + a few spurious gluons, i.e. check their pT fraction
        bool leptonOverlap = false;
        if (genparticles != nullptr) {
            for (const auto & ptr : *genparticles) {
                bool isLepton = ((abs(ptr.pdgId()) == PDGID::MUON) || (abs(ptr.pdgId()) == PDGID::ELECTRON) || (abs(ptr.pdgId()) == PDGID::TAU));
                if (!isLepton) continue;
                bool thisLeptonOverlap = isLepton && (deltaR(ptr.v4(), jet.v4()) < lepton_overlap_dr) && ((ptr.pt() / jet.pt()) > 0.5);
                leptonOverlap = leptonOverlap || thisLeptonOverlap;
            }
        }
        if ((jet.pt() > pt_min) && (fabs(jet.eta()) < eta_max) && !found && !leptonOverlap) genjets_out.push_back(jet);
    }
    sort_by_pt(genjets_out);
    return genjets_out;
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
            muons.push_back(*itr);
        }
    }
    // turn this off, because we want the latest muons, not similar ones
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
            auto thisDR = deltaR(jtr, genjtr);
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

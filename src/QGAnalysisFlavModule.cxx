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
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/QGAnalysis/include/QGAnalysisFlavCompHists.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

/** \brief Basic analysis preselection
 *
 */
class QGAnalysisFlavModule: public AnalysisModule {
public:

    explicit QGAnalysisFlavModule(Context & ctx);
    virtual bool process(Event & event) override;
    std::vector<GenJetWithParts> getGenJets(std::vector<GenJetWithParts> * jets, std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=1.5, float lepton_overlap_dr=0.2);
    std::vector<GenParticle> getGenMuons(std::vector<GenParticle> * genparticles, float pt_min=5., float eta_max=2.5);
    std::vector<Jet> getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax=0.8, bool uniqueMatch=true);

private:

    std::unique_ptr<GeneralEventSetup> common_setup;
    std::unique_ptr<MCReweighting> mc_reweight;

    // Reco selections/hists
    std::unique_ptr<Selection> njet_sel, zplusjets_sel, dijet_sel;

    Event::Handle<std::vector<GenJetWithParts>> genjets_handle;

    bool is_mc;
    float jetRadius;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;

    std::vector<std::pair<float, float>> pt_bins = {{30, 50}, {75, 100}, {100, 200}, {300, 500}, {500, 1000}};
    std::vector<std::pair<float, float>> eta_bins = {{0.0, 1.4}, {1.4, 2.6}, {2.6, 3.2}, {3.2, 4.7}};
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


    string jet_cone = ctx.get("JetCone", "AK4");
    string pu_removal = ctx.get("PURemoval", "CHS");
    if (pu_removal != "CHS" && pu_removal != "PUPPI") {
        throw runtime_error("Only PURemoval == CHS, PUPPI supported for now");
    }
    jetRadius = get_jet_radius(jet_cone);

    cout << "Running with jet cone: " << jet_cone << endl;
    cout << "Running with PUS: " << pu_removal << endl;

    common_setup.reset(new GeneralEventSetup(ctx, pu_removal, jet_cone, jetRadius));
    mc_reweight.reset(new MCReweighting(ctx));

    genjets_handle = ctx.declare_event_output< std::vector<GenJetWithParts> > ("GoodGenJets");

    // Event Selections
    njet_sel.reset(new NJetSelection(1));

    zplusjets_sel.reset(new ZplusJetsSelection(20, 20, 20, 2, 1));
    float deta = 1.2;
    float sumEta = 10.;
    dijet_sel.reset(new DijetSelection(2, 1, 1, false, 100, 100));

    // Hists
    for (auto ptBin : pt_bins ) {
        float ptMin = ptBin.first;
        float ptMax = ptBin.second;
        for (auto etaBin : eta_bins) {
            float etaMin = etaBin.first;
            float etaMax = etaBin.second;
            std::unique_ptr<QGAnalysisFlavCompHists> dj(new QGAnalysisFlavCompHists(ctx, TString::Format("Dijet_Pt%dto%d_Eta%gto%g", int(ptMin), int(ptMax), etaMin, etaMax).Data(), 2, "dijet", ptMin, ptMax, etaMin, etaMax));
            dijet_hists_binned.push_back(std::move(dj));
            std::unique_ptr<QGAnalysisFlavCompHists> zpj(new QGAnalysisFlavCompHists(ctx, TString::Format("ZPlusJets_Pt%dto%d_Eta%gto%g", int(ptMin), int(ptMax), etaMin, etaMax).Data(), 1, "zplusjets", ptMin, ptMax, etaMin, etaMax));
            zplusjets_hists_binned.push_back(std::move(zpj));
        }
    }
}


bool QGAnalysisFlavModule::process(Event & event) {
    // cout << "event" << endl;
    if (!common_setup->process(event)) {return false;}
    mc_reweight->process(event);

    if (!njet_sel->passes(event)) return false;

    std::vector<GenJetWithParts> goodGenJets = getGenJets(event.genjets, event.genparticles, 5., 5., jetRadius);
    if (goodGenJets.size() == 0) return false;
    if (event.jets->size() == 0) return false;
    sort_by_pt(goodGenJets);

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

    std::vector<Jet> goodJets = getMatchedJets(event.jets, &event.get(genjets_handle), jetRadius/2.);
    std::swap(goodJets, *event.jets);

    // RECO PART
    if (!njet_sel->passes(event)) return false;
    
    bool zpj = zplusjets_sel->passes(event);
    bool dj = dijet_sel->passes(event);
    if (!zpj && !dj) return false;

    for (uint i=0; i < pt_bins.size(); i++) {
        for (uint j=0; j < eta_bins.size(); j++) {
            uint ind = i*(pt_bins.size()-1) + j;
            if (dj) {
                dijet_hists_binned.at(ind)->fill(event);
            }
            if (zpj) {
                zplusjets_hists_binned.at(ind)->fill(event);
            }
        }
    }

    return true;
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
                leptonOverlap = leptonOverlap || (((abs(ptr.pdgId()) == PDGID::MUON) || (abs(ptr.pdgId()) == PDGID::ELECTRON) || (abs(ptr.pdgId()) == PDGID::TAU)) && (deltaR(ptr.v4(), itr.v4()) < lepton_overlap_dr) && (abs(itr.pdgId()) == 23));
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
std::vector<GenParticle> QGAnalysisFlavModule::getGenMuons(std::vector<GenParticle> * genparticles, float pt_min, float eta_max) {
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
 * Also stores index of matching GenJets in the passed Jet collection
 * Will take the closest matching GenJet as the match, provided it is within drMax.
 * uniqueMatch controls whether matching GenJets must be unique (i.e 2 reco jets can't match the same GenJet)
 *
 */
std::vector<Jet> QGAnalysisFlavModule::getMatchedJets(std::vector<Jet> * jets, std::vector<GenJetWithParts> * genjets, float drMax, bool uniqueMatch) {
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
        }
    }
    return goodJets;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the QGAnalysisFlavModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(QGAnalysisFlavModule)

}

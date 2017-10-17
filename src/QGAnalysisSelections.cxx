#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>
#include <algorithm>

using namespace uhh2examples;
using namespace uhh2;


ZplusJetsSelection::ZplusJetsSelection(float mu1_pt, float mu2_pt, float mZ_window, float dphi_jet_z_min, float second_jet_frac_max):
    mu1_pt_(mu1_pt),
    mu2_pt_(mu2_pt),
    mZ_window_(mZ_window),
    dphi_jet_z_min_(dphi_jet_z_min),
    second_jet_frac_max_(second_jet_frac_max)
    {}

bool ZplusJetsSelection::passes(const Event & event){
    if (event.muons->size() < 2) return false;
    std::vector<Muon> muons;
    for (const auto & itr : *event.muons) {
        if (itr.pt() > mu2_pt_) muons.push_back(itr);
    }
    if (muons.size() < 2) return false;

    const auto & muon1 = muons.at(0);
    const auto & muon2 = muons.at(1);

    if (muon1.pt() < mu1_pt_ || muon2.pt() < mu2_pt_) return false;

    // TODO: what if > 2 muons? how to pick the Z candidate?

    if (muon1.charge() == muon2.charge()) return false;

    const auto z_cand = muon1.v4() + muon2.v4();

    float m_mumu = z_cand.M();
    if (fabs(m_mumu - 90) > mZ_window_) return false;

    const auto & jet1 = event.jets->at(0);
    if (fabs(deltaPhi(jet1, z_cand)) < dphi_jet_z_min_) return false;

    if (event.jets->size() > 1) {
        const auto & jet2 = event.jets->at(1);
        auto jet_frac = jet2.pt() / z_cand.pt();
        if (jet_frac > second_jet_frac_max_) return false;
    }

    return true;
}


DijetSelection::DijetSelection(float dphi_min, float third_frac_max):
    dphi_min_(dphi_min),
    third_frac_max_(third_frac_max){}

bool DijetSelection::passes(const Event & event){
    if (event.jets->size() < 2) return false;
    const auto & jet1 = event.jets->at(0);
    const auto & jet2 = event.jets->at(1);

    auto dphi = fabs(deltaPhi(jet1, jet2));
    if (dphi < dphi_min_) return false;

    if (event.jets->size() == 2) return true;

    const auto & jet3 = event.jets->at(2);
    auto third_jet_frac = jet3.pt() / (0.5 * (jet1.pt() + jet2.pt()));
    return third_jet_frac < third_frac_max_;
}


JetFlavourSelection::JetFlavourSelection(std::vector<int> pdgids, uint jet_index): jet_index_(jet_index), flavours_(pdgids) {}

bool JetFlavourSelection::passes(const Event & event) {
    if (event.jets->size() <= jet_index_) return false;
    const auto & jet = event.jets->at(jet_index_);
    uint flav = abs(jet.genPartonFlavor());
    return std::find(flavours_.begin(), flavours_.end(), flav) != flavours_.end();
}


ZplusJetsTheorySelection::ZplusJetsTheorySelection(Context & ctx, float pt_min, float jet_frac_min, float jet_z_deta_max, float second_jet_frac_max, float mZ_window):
    pt_min_(pt_min),
    jet_frac_min_(jet_frac_min),
    jet_z_deta_max_(jet_z_deta_max),
    second_jet_frac_max_(second_jet_frac_max),
    mZ_window_(mZ_window) {
        genJets_handle = ctx.get_handle<std::vector<GenJetWithParts>>("GoodGenJets");
        genMuons_handle = ctx.get_handle<std::vector<GenParticle>>("GoodGenMuons");
    }

bool ZplusJetsTheorySelection::passes(const Event & event) {
    if(!event.is_valid(genJets_handle) || !event.is_valid(genMuons_handle)) return false;

    const auto & genjets = event.get(genJets_handle);
    const auto & muons = event.get(genMuons_handle);

    if ((genjets.size() < 1) || (muons.size() < 2)) return false;

    // Get Z candidate
    const auto & muon1 = muons.at(0);
    const auto & muon2 = muons.at(1);

    if (muon1.charge() == muon2.charge()) return false;

    const auto z_cand = muon1.v4() + muon2.v4();

    float m_mumu = z_cand.M();
    if (fabs(m_mumu - 90) > mZ_window_) return false;

    if (z_cand.pt() < pt_min_) return false;

    const auto & jet = genjets.at(0);

    if ((jet.pt() / z_cand.pt()) < jet_frac_min_) return false;

    if (fabs(jet.eta() - z_cand.eta()) > jet_z_deta_max_) return false;

    if (genjets.size() >=2) {
        const auto & jet2 = genjets.at(1);
        if ((jet2.pt() / z_cand.pt()) > second_jet_frac_max_) return false;
    }

    return true;
}


DijetTheorySelection::DijetTheorySelection(Context & ctx, float pt_min, float jet_frac_min, float jet_deta_max, float third_frac_max):
    pt_min_(pt_min),
    jet_frac_min_(jet_frac_min),
    jet_deta_max_(jet_deta_max),
    third_frac_max_(third_frac_max) {
        genJets_handle = ctx.get_handle<std::vector<GenJetWithParts>>("GoodGenJets");
    }

bool DijetTheorySelection::passes(const Event & event) {
    if(!event.is_valid(genJets_handle)) return false;

    const auto & genjets = event.get(genJets_handle);

    if (genjets.size() < 2) return false;

    const auto & jet1 = genjets.at(0);
    const auto & jet2 = genjets.at(1);

    if ((jet1.pt() + jet2.pt()) < (2. * pt_min_)) return false;

    if (jet2.pt() < (jet1.pt() * jet_frac_min_)) return false;

    if (fabs(jet1.eta() - jet2.eta()) > jet_deta_max_) return false;

    if (genjets.size()>=3) {
        const auto & jet3 = genjets.at(2);
        if ((jet3.pt() / (0.5*(jet1.pt()+jet2.pt()))) > third_frac_max_) return false;
    }
    return true;
}
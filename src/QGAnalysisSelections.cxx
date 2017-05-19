#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>

using namespace uhh2examples;
using namespace uhh2;


ZplusJetsSelection::ZplusJetsSelection(float mZ_window, float dphi_jet_z_min, float second_jet_frac_max):
    mZ_window_(mZ_window),
    dphi_jet_z_min_(dphi_jet_z_min),
    second_jet_frac_max_(second_jet_frac_max)
    {}

bool ZplusJetsSelection::passes(const Event & event){
    if (event.muons->size() < 2) return false;
    const auto & muon0 = event.muons->at(0);
    const auto & muon1 = event.muons->at(1);

    // TODO: what if > 2 muons? how to pick the Z candidate?

    if (muon0.charge() == muon1.charge()) return false;

    const auto z_cand = muon0.v4() + muon1.v4();

    float m_mumu = z_cand.M();
    if (fabs(m_mumu - 90) > mZ_window_) return false;

    const auto & jet0 = event.jets->at(0);
    if (deltaPhi(jet0, z_cand) < dphi_jet_z_min_) return false;

    if (event.jets->size() > 1) {
        const auto & jet2 = event.jets->at(1);
        auto jet_frac = jet2.pt() / z_cand.pt();
        if (jet_frac > second_jet_frac_max_) return false;
    }

    return true;
}


DijetSelection::DijetSelection(float dphi_min, float third_frac_max): dphi_min_(dphi_min), third_frac_max_(third_frac_max){}

bool DijetSelection::passes(const Event & event){
    if (event.jets->size() < 2) return false;
    const auto & jet0 = event.jets->at(0);
    const auto & jet1 = event.jets->at(1);

    auto dphi = deltaPhi(jet0, jet1);
    if (dphi < dphi_min_) return false;

    if (event.jets->size() == 2) return true;

    const auto & jet2 = event.jets->at(2);
    auto third_jet_frac = jet2.pt() / (0.5 * (jet0.pt() + jet1.pt()));
    return third_jet_frac < third_frac_max_;
}




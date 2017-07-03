#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"


namespace uhh2examples {

/*
 * Select Z->mumu + jet candidate events
 */
class ZplusJetsSelection: public uhh2::Selection {
public:
    ZplusJetsSelection(float mu1_pt=20., float mu2_pt=20., float mZ_window=20., float dphi_jet_z_min=2.0, float second_jet_frac_max=0.3);
    virtual bool passes(const uhh2::Event & event) override;
private:
    float mu1_pt_, mu2_pt_, mZ_window_, dphi_jet_z_min_, second_jet_frac_max_;
};


/*
 * Select dijet candidate events
 */
class DijetSelection: public uhh2::Selection {
public:
    DijetSelection(float dphi_min=2.0, float third_frac_max=0.3);
    virtual bool passes(const uhh2::Event & event) override;
private:
    float dphi_min_, third_frac_max_;
};


/*
 * Select event where the nth jet has suitable flavour
 */
class JetFlavourSelection: public uhh2::Selection {
public:
    JetFlavourSelection(std::vector<int> pdgids, uint jet_index=0);
    virtual bool passes(const uhh2::Event & event) override;
private:
    uint jet_index_;
    std::vector<int> flavours_;
};


}

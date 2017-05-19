#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"

namespace uhh2examples {

/*
 * Select Z->mumu + jet candidate events
 */
class ZplusJetsSelection: public uhh2::Selection {
public:
    ZplusJetsSelection();
    virtual bool passes(const uhh2::Event & event) override;

};

/*
 * Select dijet candidate events
 */
class DijetSelection: public uhh2::Selection {
public:
    DijetSelection(float dphi_min = 2.7f, float third_frac_max = 0.2f);
    virtual bool passes(const uhh2::Event & event) override;
private:
    float dphi_min, third_frac_max;
};

}

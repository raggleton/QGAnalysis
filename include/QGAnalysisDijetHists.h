#pragma once

#include "UHH2/core/include/Hists.h"

namespace uhh2examples {

/**  \brief Control histograms for dijet selection. Basically a sanity check.
 *
 */
class QGAnalysisDijetHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisDijetHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisDijetHists();

protected:
    TH1F *n_jets;
    TH1F *pt_jet1, *eta_jet1, *phi_jet1, *pt_jet2, *eta_jet2, *phi_jet2, *m_jj;
    TH2F *deta_dphi_jj;
};

}

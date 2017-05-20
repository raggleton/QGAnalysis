#pragma once

#include "UHH2/core/include/Hists.h"

namespace uhh2examples {

/**  \brief Control histograms for Z+jets selection. Basically a sanity check.
 *
 */
class QGAnalysisZPlusJetsHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisZPlusJetsHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisZPlusJetsHists();

protected:
    TH1F *n_jets;
    TH1F *n_mu, *pt_mu1, *eta_mu1, *reliso_mu1, *pt_mu2, *eta_mu2, *reliso_mu2, *m_mumu;
    TH2F *deta_dphi_mumu;
};

}

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
    TH1F *pt_jet1, *eta_jet1, *phi_jet1, *pt_jet2, *eta_jet2, *phi_jet2, *m_jj, *pt_jet1_jet2_ratio;
    TH1F *pt_jet3, *eta_jet3, *pt_jet3_frac;
    TH2F *deta_dphi_jj, *pt_jet1_jet2, *flav_jet1_jet2;
    TH1F *met_sig;
    TH2F *met_sig_pt_jet1;

    bool doHerwigReweighting;
    TH1F * reweightHist;
};

}

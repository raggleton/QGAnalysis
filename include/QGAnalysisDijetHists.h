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
    TH2F *n_jets_vs_pt_jet1;
    TH1F *pt_jet1;
    TH2F *eta_jet1_vs_pt_jet1, *phi_jet1_vs_pt_jet1;
    TH2F *pt_jet2_vs_pt_jet1, *eta_jet2_vs_pt_jet1, *phi_jet2_vs_pt_jet1;
    TH2F *m_jj_vs_pt_jet1, *pt_jet1_jet2_ratio_vs_pt_jet1, *deta_jj_vs_pt_jet1, *dphi_jj_vs_pt_jet1;
    TH2F *pt_jet3_vs_pt_jet1, *eta_jet3_vs_pt_jet1, *pt_jet3_frac_vs_pt_jet1;

    TH2F *deta_dphi_jj, *flav_jet1_jet2;

    TH2F *met_sig_vs_pt_jet1;

    bool doHerwigReweighting;
    TH1F * reweightHist;
};

}

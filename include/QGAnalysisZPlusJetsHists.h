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
    TH2F *n_jets_vs_pt_jet1;
    TH1F *pt_jet1;
    TH2F *eta_jet1_vs_pt_jet1, *pt_jet1_z_ratio_vs_pt_jet1;
    TH2F *pt_jet2_vs_pt_jet1, *eta_jet2_vs_pt_jet1, *pt_jet2_z_ratio_vs_pt_jet1;
    TH2F *n_mu_vs_pt_jet1, *pt_mu1_vs_pt_jet1, *eta_mu1_vs_pt_jet1, *reliso_mu1_vs_pt_jet1;
    TH2F *pt_mu2_vs_pt_jet1, *eta_mu2_vs_pt_jet1, *reliso_mu2_vs_pt_jet1, *m_mumu_vs_pt_jet1, *pt_mumu_vs_pt_jet1, *eta_mumu_vs_pt_jet1;
    TH2F *dphi_j_z_vs_pt_jet1;
    TH2F *deta_mumu_vs_pt_jet1, *dphi_mumu_vs_pt_jet1;
    TH2F *deta_mumu_jet1_vs_pt_jet1, *dphi_mumu_jet1_vs_pt_jet1;
    TH2F *pt_jet1_z_pt_jet2_z_ratio;

    bool doHerwigReweighting;
    TH1F * reweightHist;
};

}

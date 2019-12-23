#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"

namespace uhh2examples {

/**  \brief Control histograms for Z+jets selection. Basically a sanity check.
 *
 */
class QGAnalysisZPlusJetsHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisZPlusJetsHists(uhh2::Context & ctx, const std::string & dirname, const std::string & zLabel_, const std::string & genjets_name="GoodGenJets");

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisZPlusJetsHists();

protected:
    TH2F *n_jets_vs_pt;
    TH1F *pt_jet1, *pt_jet1_unweighted, *pt_jet_response_binning, *pt_genjet_response_binning;
    TH1F *pt_mumu;
    TH1F *gen_ht;
    TH2F *eta_jet1_vs_pt;
    TH2F *flav_jet1_vs_pt_jet1;
    TH2F *pt_jet1_z_ratio_vs_pt, *pt_jet1_z_ratio_vs_pt_jet1;
    TH2F *jet1_z_asym_vs_pt, *jet1_z_asym_vs_pt_jet1;
    TH2F *pt_jet2_vs_pt, *eta_jet2_vs_pt, *pt_jet2_z_ratio_vs_pt;
    TH2F *met_vs_pt, *met_sig_vs_pt;
    TH2F *n_mu_vs_pt, *pt_mu1_vs_pt, *eta_mu1_vs_pt, *reliso_mu1_vs_pt;
    TH2F *pt_mu2_vs_pt, *eta_mu2_vs_pt, *reliso_mu2_vs_pt, *m_mumu_vs_pt, *pt_jet1_vs_pt, *eta_mumu_vs_pt;
    TH2F *dphi_j_z_vs_pt;
    TH2F *deta_mumu_vs_pt, *dphi_mumu_vs_pt;
    TH2F *deta_mumu_jet1_vs_pt, *dphi_mumu_jet1_vs_pt;
    TH2F *pt_jet1_pt_jet2_ratio_vs_pt;
    TH2F *genjet1_ind_vs_pt_jet1, *pt_jet_response, *pt_jet_response_fine, *eta_jet_response;
    TH1F *pt_jet_genHT_ratio;
    std::vector<TH2F*> genjet_recojet_ind_binned;

    TH1F * n_pv;

    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;
    uhh2::Event::Handle<std::vector<Muon>> hndlZ;

    bool DO_MATCHING_INDS;
    float jetRadius;
};

}

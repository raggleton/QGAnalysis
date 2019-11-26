#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"

namespace uhh2examples {

/**  \brief Control histograms for dijet selection. Basically a sanity check.
 *
 */
class QGAnalysisDijetHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisDijetHists(uhh2::Context & ctx, const std::string & dirname, const std::string & binning_);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisDijetHists();

protected:
    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;

    TH2F *n_jets_vs_pt_jet;
    TH1F *pt_jet, *pt_jet_response_binning, *pt_genjet_response_binning;
    TH2F *eta_jet1_vs_pt_jet, *phi_jet1_vs_pt_jet;
    TH2F *eta_jet1_vs_eta_jet2;
    TH2F *pt_jet2_vs_pt_jet, *eta_jet2_vs_pt_jet, *phi_jet2_vs_pt_jet;
    TH2F *m_jj_vs_pt_jet, *pt_jet1_jet2_ratio_vs_pt_jet, *deta_jj_vs_pt_jet, *dphi_jj_vs_pt_jet, *sumeta_jj_vs_pt_jet, *jet1_jet2_asym_vs_pt_jet;
    TH2F *pt_jet3_vs_pt_jet, *eta_jet3_vs_pt_jet, *pt_jet3_frac_vs_pt_jet;

    TH2F *deta_dphi_jj, *flav_jet1_jet2;

    TH2F *met_vs_pt_jet, *met_sig_vs_pt_jet;

    TH2F *pt_jet_response_fine, *pt_jet_response, *eta_jet_response;
    TH1F *pt_jet_qScale_ratio;
    TH1F *pt_jet_genHT_ratio;
    TH2F *pt_jet_vs_pdf_scalePDF;
    TH2F *pt_jet_vs_genHT;
    TH2F *weight_vs_puHat_genHT_ratio;

    TH2F *genjet1_ind_vs_pt_jet1, *genjet2_ind_vs_pt_jet2;
    std::vector<TH2F*> genjet_recojet_ind_binned;

    bool is_mc_;
    TH1F * n_pv;

    std::vector<double> pt_bin_edges;
    int nbins_pt;

    std::string binning;
    bool DO_MATCHING_INDS;
    float jetRadius;
};

}

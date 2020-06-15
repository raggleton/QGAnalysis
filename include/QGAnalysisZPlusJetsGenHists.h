#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"

namespace uhh2examples {

/**  \brief Gen level histograms for Z+jets
 *
 */
class QGAnalysisZPlusJetsGenHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisZPlusJetsGenHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisZPlusJetsGenHists();
    GenParticle findGenZ(std::vector<GenParticle> & gps);

protected:
    TH1F *deta_mumu;
    TH1F *deta_mumu_jet1;
    TH1F *dphi_mumu;
    TH1F *dphi_mumu_jet1;
    TH1F *eta_jet1;
    TH1F *eta_jet2;
    TH1F *eta_mu1;
    TH1F *eta_mu2;
    TH1F *eta_mumu;
    TH1F *eta_z;
    TH1F *gen_ht;
    TH1F *m_mumu;
    TH1F *n_jets;
    TH1F *n_mu;
    TH1F *pt_hat;
    TH1F *pt_jet1;
    TH1F *pt_jet1_pt_jet2_ratio;
    TH1F *pt_jet1_z_ratio;
    TH1F *pt_jet2;
    TH1F *pt_jet2_z_ratio;
    TH1F *pt_jet_genHT_ratio;
    TH1F *pt_hat_pt_jet_ratio;
    TH1F *pt_mu1;
    TH1F *pt_mu2;
    TH1F *pt_mumu;
    TH1F *pt_z;
    TH1F *q_scale;
    TH1F *jet_kt;
    TH1F *jet_kt_pt_z_ratio;

    uhh2::Event::Handle<std::vector<GenJet> > genJets_handle;
    uhh2::Event::Handle<std::vector<GenParticle>> genMuons_handle;
    uhh2::Event::Handle<double> gen_weight_handle;
};

}

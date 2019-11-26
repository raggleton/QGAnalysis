#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"

namespace uhh2examples {

/**  \brief Gen level histograms for Dijets
 *
 */
class QGAnalysisDijetGenHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisDijetGenHists(uhh2::Context & ctx, const std::string & dirname, const std::string & genjets_name="GoodGenJets");

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisDijetGenHists();
    // const GenParticle & findGenZ(std::vector<GenParticle> & gps);

protected:
    TH1F *deta_jet;
    TH1F *dphi_jet;
    TH1F *dr_jet;
    TH1F *eta_jet1;
    TH1F *eta_jet2;
    TH1F *gen_ht;
    TH1F *n_jets;
    TH1F *pt_hat;
    TH1F *pt_jet1;
    TH1F *pt_jet1_pt_jet2_ratio;
    TH1F *pt_jet2;
    TH1F *pt_dijet_ave;
    TH1F *pt_jet_genHT_ratio;
    TH1F *pt_hat_pt_jet_ratio;
    TH1F *pt_asym;
    TH1F *q_scale;

    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;
    uhh2::Event::Handle<double> gen_weight_handle;
};

}

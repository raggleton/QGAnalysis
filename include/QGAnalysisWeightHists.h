#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"

#include "TH3.h"


namespace uhh2examples {

/**  \brief Control histograms for weights & correlations
 *
 */
class QGAnalysisWeightHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisWeightHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisWeightHists();

protected:
    uhh2::Event::Handle<std::vector<GenJet> > genJets_handle;
    TH3F *weight_vs_pt_vs_pt_jet_genHT_ratio,
         *weight_vs_pt_vs_pt_jet_genHT_ratio_unweighted,
         *weight_vs_pt_genjet_vs_pt_jet_genHT_ratio,
         *weight_vs_pt_genjet_vs_pt_jet_genHT_ratio_unweighted,
         *weight_vs_pt_vs_pt_genjet_genHT_ratio,
         *weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio,
         *weight_vs_pt_vs_pt_genjet_genHT_ratio_unweighted,
         *weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio_unweighted,
         *weight_vs_pt_vs_pt_jet_qScale_ratio,
         *weight_vs_pt_vs_pt_jet_qScale_ratio_unweighted,
         *weight_vs_pt_genjet_vs_pt_jet_qScale_ratio,
         *weight_vs_pt_genjet_vs_pt_jet_qScale_ratio_unweighted,
         *weight_vs_pt_vs_pt_genjet_qScale_ratio,
         *weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio,
         *weight_vs_pt_vs_pt_genjet_qScale_ratio_unweighted,
         *weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio_unweighted,
         *weight_vs_pt_vs_pt_jet_ptHat_ratio,
         *weight_vs_pt_vs_pt_jet_ptHat_ratio_unweighted,
         *weight_vs_pt_genjet_vs_pt_jet_ptHat_ratio,
         *weight_vs_pt_genjet_vs_pt_jet_ptHat_ratio_unweighted,
         *weight_vs_pt_vs_pt_genjet_ptHat_ratio,
         *weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio,
         *weight_vs_pt_vs_pt_genjet_ptHat_ratio_unweighted,
         *weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio_unweighted,
         *weight_vs_pt_vs_PU_ptHat_genHT_ratio,
         *weight_vs_pt_vs_PU_ptHat_genHT_ratio_unweighted,
         *weight_vs_pt_genjet_vs_PU_ptHat_genHT_ratio,
         *weight_vs_pt_genjet_vs_PU_ptHat_genHT_ratio_unweighted,
         *weight_vs_pt_vs_PU_ptHat_ptHat_ratio,
         *weight_vs_pt_vs_PU_ptHat_ptHat_ratio_unweighted,
         *weight_vs_pt_genjet_vs_PU_ptHat_ptHat_ratio,
         *weight_vs_pt_genjet_vs_PU_ptHat_ptHat_ratio_unweighted,
         *weight_vs_pt_vs_PU_ptHat_jetkT_ratio,
         *weight_vs_pt_vs_PU_ptHat_jetkT_ratio_unweighted,
         *weight_vs_pt_vs_pt_jet_jetkT_ratio,
         *weight_vs_pt_vs_pt_jet_jetkT_ratio_unweighted;

    bool is_mc_;
};

}

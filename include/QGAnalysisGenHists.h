#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/QGAnalysis/include/QGAddModules.h"

namespace uhh2examples {

/**  \brief Book and fill main analysis histograms for genjets, passing gen selection
 */
class QGAnalysisGenHists: public uhh2::Hists {
public:
    QGAnalysisGenHists(uhh2::Context & ctx,
                       const std::string & dirname,
                       int useNJets,
                       bool doGroomed,
                       bool useStatus23Flavour,
                       const std::string & selection,
                       const std::string & gen_sel_handle_name,
                       const std::string & gen_jetlambda_handle_name);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisGenHists();
protected:
    // std::vector<GenParticle*> get_genjet_genparticles(const GenJet &, std::vector<GenParticle>*);
    int get_jet_flavour(const Particle & obj, std::vector<GenParticle>* genparticles, float dr_max, bool pythiaMode);
    uint get_num_outgoing_partons(const std::vector<GenParticle> & genparticles);

    float jetRadius;
    TH1F *h_weights;
    TH2F *h_weights_vs_pt, *h_pthat_vs_weight, *h_pthat_vs_jet_pt;
    TH1F *h_jet_pt, *h_jet_pt_unweighted, *h_jet_eta, *h_jet_y, *h_jet_flavour;

    TH2F *h_jet_puppiMultiplicity_vs_pt, *h_jet_LHA_vs_pt, *h_jet_pTD_vs_pt, *h_jet_width_vs_pt, *h_jet_thrust_vs_pt;
    TH2F *h_jet_puppiMultiplicity_charged_vs_pt, *h_jet_LHA_charged_vs_pt, *h_jet_pTD_charged_vs_pt, *h_jet_width_charged_vs_pt, *h_jet_thrust_charged_vs_pt;

    TH2F *h_jet_flavour_vs_pt, *h_jet1_flavour_vs_pt, *h_jet2_flavour_vs_pt;
    std::vector<TH2F *> h_jet_flavour_vs_pt_nPartons, h_jet1_flavour_vs_pt_nPartons, h_jet2_flavour_vs_pt_nPartons;
    TH2F *h_jet_flavour_vs_eta_lowPt, *h_jet_flavour_vs_eta_midPt, *h_jet_flavour_vs_eta_highPt, *h_jet_flavour_vs_eta_highPt2;

    TH2F *h_gjet_puppiMultiplicity_vs_pt, *h_gjet_LHA_vs_pt, *h_gjet_pTD_vs_pt, *h_gjet_width_vs_pt, *h_gjet_thrust_vs_pt;
    TH2F *h_gjet_puppiMultiplicity_charged_vs_pt, *h_gjet_LHA_charged_vs_pt, *h_gjet_pTD_charged_vs_pt, *h_gjet_width_charged_vs_pt, *h_gjet_thrust_charged_vs_pt;

    TH2F *h_qjet_puppiMultiplicity_vs_pt, *h_qjet_LHA_vs_pt, *h_qjet_pTD_vs_pt, *h_qjet_width_vs_pt, *h_qjet_thrust_vs_pt;
    TH2F *h_qjet_puppiMultiplicity_charged_vs_pt, *h_qjet_LHA_charged_vs_pt, *h_qjet_pTD_charged_vs_pt, *h_qjet_width_charged_vs_pt, *h_qjet_thrust_charged_vs_pt;

    TH2F *h_bcjet_puppiMultiplicity_vs_pt, *h_bcjet_LHA_vs_pt, *h_bcjet_pTD_vs_pt, *h_bcjet_width_vs_pt, *h_bcjet_thrust_vs_pt;
    TH2F *h_bcjet_puppiMultiplicity_charged_vs_pt, *h_bcjet_LHA_charged_vs_pt, *h_bcjet_pTD_charged_vs_pt, *h_bcjet_width_charged_vs_pt, *h_bcjet_thrust_charged_vs_pt;

    std::string dirname_;
    int useNJets_;
    bool doGroomed_;

    uhh2::Event::Handle<std::vector<GenJetLambdaBundle> > genJetsLambda_handle;
    uhh2::Event::Handle<bool> pass_gen_handle;

    bool useStatus23Flavour_;
    uint N_PARTONS_MAX;
    MC::Name dataset_;
};

}

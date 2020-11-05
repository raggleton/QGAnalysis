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

/**  \brief Book and fill main analysis histograms (reco)
 */
class QGAnalysisHists: public uhh2::Hists {
public:
    QGAnalysisHists(uhh2::Context & ctx,
                    const std::string & dirname,
                    int useNJets,
                    bool doGroomed,
                    bool useStatus23Flavour,
                    const std::string & selection,
                    const std::string & reco_sel_handle_name,
                    const std::string & gen_sel_handle_name,
                    const std::string & reco_jetlambda_handle_name,
                    const std::string & gen_jetlambda_handle_name);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisHists();
protected:
    // std::vector<GenParticle*> get_genjet_genparticles(const GenJet &, std::vector<GenParticle>*);
    int get_jet_flavour(const Particle & obj, std::vector<GenParticle>* genparticles, float dr_max, bool pythiaMode);
    std::vector<PFParticle*> get_jet_pfparticles(const Jet &, std::vector<PFParticle>*);
    void fill_lambda_rsp_hists(float reco_val,
                               float gen_val,
                               float weight,
                               TH2F * response,
                               TH2F * rel_response,
                               float jet_pt,
                               TH2F * response_lowPt,
                               TH2F * response_midPt,
                               TH2F * response_highPt,
                               TH2F * rel_response_lowPt,
                               TH2F * rel_response_midPt,
                               TH2F * rel_response_highPt);
    uint get_num_outgoing_partons(const std::vector<GenParticle> & genparticles);
    // reco jet hists
    float jetRadius;
    TH1F * h_weights;
    TH2F * h_weights_vs_pt, *h_pthat_vs_weight, *h_pthat_vs_jet_pt;

    TH1F *h_jet_pt, *h_jet_pt_unweighted, *h_jet_eta, *h_jet_y, *h_jet_flavour;

    // TH2F *h_jet_flavour_vs_pt, *h_jet1_flavour_vs_pt, *h_jet2_flavour_vs_pt, *h_jet_flavour_vs_eta;
    TH2F *h_genjet_flavour_vs_pt, *h_genjet1_flavour_vs_pt, *h_genjet2_flavour_vs_pt;
    std::vector<TH2F *> h_genjet_flavour_vs_pt_nPartons, h_genjet1_flavour_vs_pt_nPartons, h_genjet2_flavour_vs_pt_nPartons;
    TH2F *h_genjet_flavour_vs_eta_lowPt, *h_genjet_flavour_vs_eta_midPt, *h_genjet_flavour_vs_eta_highPt, *h_genjet_flavour_vs_eta_highPt2;
    TH2F *h_jet_response_vs_genjet_pt, *h_jet_pt_vs_genjet_pt;

    TH2F *h_jet_puppiMultiplicity_vs_pt, *h_jet_LHA_vs_pt, *h_jet_pTD_vs_pt, *h_jet_width_vs_pt, *h_jet_thrust_vs_pt;
    TH2F *h_jet_puppiMultiplicity_charged_vs_pt, *h_jet_LHA_charged_vs_pt, *h_jet_pTD_charged_vs_pt, *h_jet_width_charged_vs_pt, *h_jet_thrust_charged_vs_pt;

    TH2F *h_gjet_puppiMultiplicity_vs_pt, *h_gjet_LHA_vs_pt, *h_gjet_pTD_vs_pt, *h_gjet_width_vs_pt, *h_gjet_thrust_vs_pt;
    TH2F *h_gjet_puppiMultiplicity_charged_vs_pt, *h_gjet_LHA_charged_vs_pt, *h_gjet_pTD_charged_vs_pt, *h_gjet_width_charged_vs_pt, *h_gjet_thrust_charged_vs_pt;
    TH2F *h_gjet_response_vs_genjet_pt;

    TH2F *h_qjet_puppiMultiplicity_vs_pt, *h_qjet_LHA_vs_pt, *h_qjet_pTD_vs_pt, *h_qjet_width_vs_pt, *h_qjet_thrust_vs_pt;
    TH2F *h_qjet_puppiMultiplicity_charged_vs_pt, *h_qjet_LHA_charged_vs_pt, *h_qjet_pTD_charged_vs_pt, *h_qjet_width_charged_vs_pt, *h_qjet_thrust_charged_vs_pt;
    TH2F *h_qjet_response_vs_genjet_pt;

    TH2F *h_bcjet_puppiMultiplicity_vs_pt, *h_bcjet_LHA_vs_pt, *h_bcjet_pTD_vs_pt, *h_bcjet_width_vs_pt, *h_bcjet_thrust_vs_pt;
    TH2F *h_bcjet_puppiMultiplicity_charged_vs_pt, *h_bcjet_LHA_charged_vs_pt, *h_bcjet_pTD_charged_vs_pt, *h_bcjet_width_charged_vs_pt, *h_bcjet_thrust_charged_vs_pt;
    TH2F *h_bcjet_response_vs_genjet_pt;

    // gen-reco response hists
    TH2F *h_jet_puppiMultiplicity_response, *h_jet_LHA_response, *h_jet_pTD_response, *h_jet_width_response, *h_jet_thrust_response;
    TH2F *h_jet_puppiMultiplicity_charged_response, *h_jet_LHA_charged_response, *h_jet_pTD_charged_response, *h_jet_width_charged_response, *h_jet_thrust_charged_response;
    TH2F *h_jet_puppiMultiplicity_rel_response, *h_jet_LHA_rel_response, *h_jet_pTD_rel_response, *h_jet_width_rel_response, *h_jet_thrust_rel_response;
    TH2F *h_jet_puppiMultiplicity_charged_rel_response, *h_jet_LHA_charged_rel_response, *h_jet_pTD_charged_rel_response, *h_jet_width_charged_rel_response, *h_jet_thrust_charged_rel_response;

    // and split by low, mid & high pt
    TH2F *h_jet_puppiMultiplicity_lowPt_response, *h_jet_LHA_lowPt_response, *h_jet_pTD_lowPt_response, *h_jet_width_lowPt_response, *h_jet_thrust_lowPt_response;
    TH2F *h_jet_puppiMultiplicity_charged_lowPt_response, *h_jet_LHA_charged_lowPt_response, *h_jet_pTD_charged_lowPt_response, *h_jet_width_charged_lowPt_response, *h_jet_thrust_charged_lowPt_response;
    TH2F *h_jet_puppiMultiplicity_lowPt_rel_response, *h_jet_LHA_lowPt_rel_response, *h_jet_pTD_lowPt_rel_response, *h_jet_width_lowPt_rel_response, *h_jet_thrust_lowPt_rel_response;
    TH2F *h_jet_puppiMultiplicity_charged_lowPt_rel_response, *h_jet_LHA_charged_lowPt_rel_response, *h_jet_pTD_charged_lowPt_rel_response, *h_jet_width_charged_lowPt_rel_response, *h_jet_thrust_charged_lowPt_rel_response;

    TH2F *h_jet_puppiMultiplicity_midPt_response, *h_jet_LHA_midPt_response, *h_jet_pTD_midPt_response, *h_jet_width_midPt_response, *h_jet_thrust_midPt_response;
    TH2F *h_jet_puppiMultiplicity_charged_midPt_response, *h_jet_LHA_charged_midPt_response, *h_jet_pTD_charged_midPt_response, *h_jet_width_charged_midPt_response, *h_jet_thrust_charged_midPt_response;
    TH2F *h_jet_puppiMultiplicity_midPt_rel_response, *h_jet_LHA_midPt_rel_response, *h_jet_pTD_midPt_rel_response, *h_jet_width_midPt_rel_response, *h_jet_thrust_midPt_rel_response;
    TH2F *h_jet_puppiMultiplicity_charged_midPt_rel_response, *h_jet_LHA_charged_midPt_rel_response, *h_jet_pTD_charged_midPt_rel_response, *h_jet_width_charged_midPt_rel_response, *h_jet_thrust_charged_midPt_rel_response;

    TH2F *h_jet_puppiMultiplicity_highPt_response, *h_jet_LHA_highPt_response, *h_jet_pTD_highPt_response, *h_jet_width_highPt_response, *h_jet_thrust_highPt_response;
    TH2F *h_jet_puppiMultiplicity_charged_highPt_response, *h_jet_LHA_charged_highPt_response, *h_jet_pTD_charged_highPt_response, *h_jet_width_charged_highPt_response, *h_jet_thrust_charged_highPt_response;
    TH2F *h_jet_puppiMultiplicity_highPt_rel_response, *h_jet_LHA_highPt_rel_response, *h_jet_pTD_highPt_rel_response, *h_jet_width_highPt_rel_response, *h_jet_thrust_highPt_rel_response;
    TH2F *h_jet_puppiMultiplicity_charged_highPt_rel_response, *h_jet_LHA_charged_highPt_rel_response, *h_jet_pTD_charged_highPt_rel_response, *h_jet_width_charged_highPt_rel_response, *h_jet_thrust_charged_highPt_rel_response;

    std::string dirname_;
    int useNJets_;
    bool doPuppi_;
    bool doGroomed_;

    uhh2::Event::Handle<std::vector<GenJetLambdaBundle> > genJetsLambda_handle;
    uhh2::Event::Handle<std::vector<JetLambdaBundle> > jetsLambda_handle;
    uhh2::Event::Handle<bool> pass_reco_handle, pass_gen_handle;

    // uhh2::Event::Handle<std::vector<GenJet> > genJets_handle;

    bool is_mc_;
    float rsp_lowPt_cut_, rsp_midPt_cut_, rsp_highPt_cut_, rsp_highPt2_cut_;
    bool useStatus23Flavour_;
    uint N_PARTONS_MAX;
};


class QGJetTrigHists: public uhh2::Hists {
public:
    QGJetTrigHists(uhh2::Context & ctx, const std::string & dirname, const std::vector<std::string> & trigNames_);
    virtual void fill(const uhh2::Event & event) override;
    virtual ~QGJetTrigHists();

protected:
    std::vector<TH2F*> hTrigs;
    std::vector<TriggerSelection> trigSels;
    std::vector<std::string> trigNames;
    TH2F * hAll;
};

}

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

/**  \brief Book and fill main analysis histograms
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class QGAnalysisHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisHists(uhh2::Context & ctx, const std::string & dirname, int useNJets, const std::string & selection, const std::string & reco_sel_handle_name, const std::string & gen_sel_handle_name);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisHists();
protected:
    // std::vector<GenParticle*> get_genjet_genparticles(const GenJetWithParts &, std::vector<GenParticle>*);
    std::vector<PFParticle*> get_jet_pfparticles(const Jet &, std::vector<PFParticle>*);
    void fill_lambda_rsp_hists(float reco_val, float gen_val, float weight,
                               TH2F * response, TH2F * rel_response,
                               float jet_pt,
                               TH2F * response_lowPt, TH2F * response_midPt, TH2F * response_highPt,
                               TH2F * rel_response_lowPt, TH2F * rel_response_midPt, TH2F * rel_response_highPt);

    // reco jet hists
    float jetRadius, LHA_rescale, width_rescale, thrust_rescale;
    TH1F * h_weights;
    TH2F * h_weights_vs_pt, *h_pthat_vs_weight, *h_pthat_vs_jet_pt;

    TH1F *h_jet_pt, *h_jet_pt_unweighted, *h_jet_eta, *h_jet_flavour;
    TH1F *h_jet_multiplicity, *h_jet_puppiMultiplicity, *h_jet_LHA, *h_jet_pTD, *h_jet_width, *h_jet_thrust;
    TH1F *h_jet_multiplicity_charged, *h_jet_puppiMultiplicity_charged, *h_jet_LHA_charged, *h_jet_pTD_charged, *h_jet_width_charged, *h_jet_thrust_charged;
    TH1F *h_qjet_multiplicity, *h_qjet_LHA, *h_qjet_pTD, *h_qjet_width, *h_qjet_thrust;
    TH1F *h_gjet_multiplicity, *h_gjet_LHA, *h_gjet_pTD, *h_gjet_width, *h_gjet_thrust;

    TH2F *h_jet_multiplicity_vs_pt, *h_jet_LHA_vs_pt, *h_jet_pTD_vs_pt, *h_jet_width_vs_pt, *h_jet_thrust_vs_pt;
    TH2F *h_jet_multiplicity_charged_vs_pt, *h_jet_LHA_charged_vs_pt, *h_jet_pTD_charged_vs_pt, *h_jet_width_charged_vs_pt, *h_jet_thrust_charged_vs_pt;
    TH2F *h_jet_flavour_vs_pt, *h_jet1_flavour_vs_pt, *h_jet2_flavour_vs_pt;
    TH2F *h_jet_flavour_vs_eta, *h_jet_response_vs_genjet_pt;
    TH2F *h_qjet_multiplicity_vs_pt, *h_qjet_LHA_vs_pt, *h_qjet_pTD_vs_pt, *h_qjet_width_vs_pt, *h_qjet_thrust_vs_pt, *h_qjet_response_vs_genjet_pt;
    TH2F *h_gjet_multiplicity_vs_pt, *h_gjet_LHA_vs_pt, *h_gjet_pTD_vs_pt, *h_gjet_width_vs_pt, *h_gjet_thrust_vs_pt, *h_gjet_response_vs_genjet_pt;
    TH2F *h_jet_puppiMultiplicity_vs_pt, *h_jet_puppiMultiplicity_charged_vs_pt, *h_qjet_puppiMultiplicity_vs_pt, *h_gjet_puppiMultiplicity_vs_pt;

    // gen-reco response hists
    TH2F * h_jet_multiplicity_response, *h_jet_puppiMultiplicity_response, *h_jet_LHA_response, *h_jet_pTD_response, *h_jet_width_response, *h_jet_thrust_response;
    TH2F * h_jet_multiplicity_charged_response, *h_jet_puppiMultiplicity_charged_response, *h_jet_LHA_charged_response, *h_jet_pTD_charged_response, *h_jet_width_charged_response, *h_jet_thrust_charged_response;
    TH2F * h_jet_multiplicity_rel_response, *h_jet_puppiMultiplicity_rel_response, *h_jet_LHA_rel_response, *h_jet_pTD_rel_response, *h_jet_width_rel_response, *h_jet_thrust_rel_response;
    TH2F * h_jet_multiplicity_charged_rel_response, *h_jet_puppiMultiplicity_charged_rel_response, *h_jet_LHA_charged_rel_response, *h_jet_pTD_charged_rel_response, *h_jet_width_charged_rel_response, *h_jet_thrust_charged_rel_response;

    // and split by low, mid & high pt
    TH2F * h_jet_multiplicity_lowPt_response, *h_jet_puppiMultiplicity_lowPt_response, *h_jet_LHA_lowPt_response, *h_jet_pTD_lowPt_response, *h_jet_width_lowPt_response, *h_jet_thrust_lowPt_response;
    TH2F * h_jet_multiplicity_charged_lowPt_response, *h_jet_puppiMultiplicity_charged_lowPt_response, *h_jet_LHA_charged_lowPt_response, *h_jet_pTD_charged_lowPt_response, *h_jet_width_charged_lowPt_response, *h_jet_thrust_charged_lowPt_response;
    TH2F * h_jet_multiplicity_lowPt_rel_response, *h_jet_puppiMultiplicity_lowPt_rel_response, *h_jet_LHA_lowPt_rel_response, *h_jet_pTD_lowPt_rel_response, *h_jet_width_lowPt_rel_response, *h_jet_thrust_lowPt_rel_response;
    TH2F * h_jet_multiplicity_charged_lowPt_rel_response, *h_jet_puppiMultiplicity_charged_lowPt_rel_response, *h_jet_LHA_charged_lowPt_rel_response, *h_jet_pTD_charged_lowPt_rel_response, *h_jet_width_charged_lowPt_rel_response, *h_jet_thrust_charged_lowPt_rel_response;

    TH2F * h_jet_multiplicity_midPt_response, *h_jet_puppiMultiplicity_midPt_response, *h_jet_LHA_midPt_response, *h_jet_pTD_midPt_response, *h_jet_width_midPt_response, *h_jet_thrust_midPt_response;
    TH2F * h_jet_multiplicity_charged_midPt_response, *h_jet_puppiMultiplicity_charged_midPt_response, *h_jet_LHA_charged_midPt_response, *h_jet_pTD_charged_midPt_response, *h_jet_width_charged_midPt_response, *h_jet_thrust_charged_midPt_response;
    TH2F * h_jet_multiplicity_midPt_rel_response, *h_jet_puppiMultiplicity_midPt_rel_response, *h_jet_LHA_midPt_rel_response, *h_jet_pTD_midPt_rel_response, *h_jet_width_midPt_rel_response, *h_jet_thrust_midPt_rel_response;
    TH2F * h_jet_multiplicity_charged_midPt_rel_response, *h_jet_puppiMultiplicity_charged_midPt_rel_response, *h_jet_LHA_charged_midPt_rel_response, *h_jet_pTD_charged_midPt_rel_response, *h_jet_width_charged_midPt_rel_response, *h_jet_thrust_charged_midPt_rel_response;

    TH2F * h_jet_multiplicity_highPt_response, *h_jet_puppiMultiplicity_highPt_response, *h_jet_LHA_highPt_response, *h_jet_pTD_highPt_response, *h_jet_width_highPt_response, *h_jet_thrust_highPt_response;
    TH2F * h_jet_multiplicity_charged_highPt_response, *h_jet_puppiMultiplicity_charged_highPt_response, *h_jet_LHA_charged_highPt_response, *h_jet_pTD_charged_highPt_response, *h_jet_width_charged_highPt_response, *h_jet_thrust_charged_highPt_response;
    TH2F * h_jet_multiplicity_highPt_rel_response, *h_jet_puppiMultiplicity_highPt_rel_response, *h_jet_LHA_highPt_rel_response, *h_jet_pTD_highPt_rel_response, *h_jet_width_highPt_rel_response, *h_jet_thrust_highPt_rel_response;
    TH2F * h_jet_multiplicity_charged_highPt_rel_response, *h_jet_puppiMultiplicity_charged_highPt_rel_response, *h_jet_LHA_charged_highPt_rel_response, *h_jet_pTD_charged_highPt_rel_response, *h_jet_width_charged_highPt_rel_response, *h_jet_thrust_charged_highPt_rel_response;

    // TH2F *h_jet_multiplicity_vs_LHA, *h_jet_multiplicity_vs_pTD, *h_jet_multiplicity_vs_width, *h_jet_multiplicity_vs_thrust;
    // TH2F *h_jet_LHA_vs_pTD, *h_jet_LHA_vs_width, *h_jet_LHA_vs_thrust;
    // TH2F *h_jet_pTD_vs_width, *h_jet_pTD_vs_thrust;
    // TH2F *h_jet_width_vs_thrust;

    // TH2F *h_qjet_multiplicity_vs_LHA, *h_qjet_multiplicity_vs_pTD, *h_qjet_multiplicity_vs_width, *h_qjet_multiplicity_vs_thrust;
    // TH2F *h_qjet_LHA_vs_pTD, *h_qjet_LHA_vs_width, *h_qjet_LHA_vs_thrust;
    // TH2F *h_qjet_pTD_vs_width, *h_qjet_pTD_vs_thrust;
    // TH2F *h_qjet_width_vs_thrust;

    // TH2F *h_gjet_multiplicity_vs_LHA, *h_gjet_multiplicity_vs_pTD, *h_gjet_multiplicity_vs_width, *h_gjet_multiplicity_vs_thrust;
    // TH2F *h_gjet_LHA_vs_pTD, *h_gjet_LHA_vs_width, *h_gjet_LHA_vs_thrust;
    // TH2F *h_gjet_pTD_vs_width, *h_gjet_pTD_vs_thrust;
    // TH2F *h_gjet_width_vs_thrust;

    // genjet hists
    TH1F *h_genjet_pt, *h_genjet_eta, *h_genjet_flavour;
    TH1F *h_genjet_multiplicity, *h_genjet_LHA, *h_genjet_pTD, *h_genjet_width, *h_genjet_thrust;
    TH1F *h_qgenjet_multiplicity, *h_qgenjet_LHA, *h_qgenjet_pTD, *h_qgenjet_width, *h_qgenjet_thrust;
    TH1F *h_ggenjet_multiplicity, *h_ggenjet_LHA, *h_ggenjet_pTD, *h_ggenjet_width, *h_ggenjet_thrust;

    TH2F *h_genjet_multiplicity_vs_pt, *h_genjet_LHA_vs_pt, *h_genjet_pTD_vs_pt, *h_genjet_width_vs_pt, *h_genjet_thrust_vs_pt;

    int useNJets_;
    bool doPuppi_;
    bool doHerwigReweighting;
    TH1F * reweightHist;

    uhh2::Event::Handle<std::vector<GenJetLambdaBundle> > genJetsLambda_handle, genJetsChargedLambda_handle;
    uhh2::Event::Handle<std::vector<JetLambdaBundle> > jetsLambda_handle, jetsChargedLambda_handle;
    uhh2::Event::Handle<double> gen_weight_handle;
    uhh2::Event::Handle<bool> pass_reco_handle;
    uhh2::Event::Handle<bool> pass_gen_handle;

    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;

    bool is_mc_;
    float rsp_midPt_cut_, rsp_highPt_cut_;
    float recoDauPtCut_;
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

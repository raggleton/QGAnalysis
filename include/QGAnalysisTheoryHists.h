#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"

namespace uhh2examples {

struct LambdaVariable {
    std::string shortName;
    std::string label;
    int nBins;
    double min;
    double max;
};

/**  \brief Book and fill main analysis histograms for theory selections
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class QGAnalysisTheoryHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisTheoryHists(uhh2::Context & ctx, const std::string & dirname, int useNJets, const std::string & selection);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisTheoryHists();
protected:
    std::vector<GenParticle*> get_genjet_genparticles(const GenJetWithParts &, std::vector<GenParticle>*);
    int get_jet_flavour(const GenJetWithParts & jet, std::vector<GenParticle>* genparticles, float dr_max=0.2, bool PythiaMode=true);

    float jetRadius;

    TH1F *h_weights;
    TH2F *h_weights_vs_pt;
    TH2F *h_pthat_vs_weight, *h_pthat_vs_genjet_pt;

    // genjet hists
    TH1F *h_genjet_pt, *h_genjet_pt_unweighted, *h_genjet_pt_all, *h_genjet_ht, *h_genjet_eta, *h_genjet_flavour;
    TH2F *h_genjet_flavour_vs_pt, *h_genjet_flavour_vs_eta;
    TH1F *h_genjet_multiplicity, *h_genjet_LHA, *h_genjet_pTD, *h_genjet_width, *h_genjet_thrust;
    TH1F *h_qgenjet_multiplicity, *h_qgenjet_LHA, *h_qgenjet_pTD, *h_qgenjet_width, *h_qgenjet_thrust;
    TH1F *h_ggenjet_multiplicity, *h_ggenjet_LHA, *h_ggenjet_pTD, *h_ggenjet_width, *h_ggenjet_thrust;

    TH2F *h_genjet_pt_vs_const_pt, *h_genjet_pt_vs_const_zi, *h_genjet_pt_vs_const_deta, *h_genjet_pt_vs_const_dphi, *h_genjet_pt_vs_const_thetai, *h_genjet_const_zi_vs_const_thetai_pt100to200, *h_genjet_const_zi_vs_const_thetai_pt800to1000;
    TH2F *h_qgenjet_pt_vs_const_pt, *h_qgenjet_pt_vs_const_zi, *h_qgenjet_pt_vs_const_deta, *h_qgenjet_pt_vs_const_dphi, *h_qgenjet_pt_vs_const_thetai, *h_qgenjet_const_zi_vs_const_thetai_pt100to200, *h_qgenjet_const_zi_vs_const_thetai_pt800to1000;
    TH2F *h_ggenjet_pt_vs_const_pt, *h_ggenjet_pt_vs_const_zi, *h_ggenjet_pt_vs_const_deta, *h_ggenjet_pt_vs_const_dphi, *h_ggenjet_pt_vs_const_thetai, *h_ggenjet_const_zi_vs_const_thetai_pt100to200, *h_ggenjet_const_zi_vs_const_thetai_pt800to1000;

    TH2F *h_genjet_multiplicity_vs_pt, *h_genjet_LHA_vs_pt, *h_genjet_pTD_vs_pt, *h_genjet_width_vs_pt, *h_genjet_thrust_vs_pt;
    TH2F *h_qgenjet_multiplicity_vs_pt, *h_qgenjet_LHA_vs_pt, *h_qgenjet_pTD_vs_pt, *h_qgenjet_width_vs_pt, *h_qgenjet_thrust_vs_pt;
    TH2F *h_ggenjet_multiplicity_vs_pt, *h_ggenjet_LHA_vs_pt, *h_ggenjet_pTD_vs_pt, *h_ggenjet_width_vs_pt, *h_ggenjet_thrust_vs_pt;

    TH2F *h_genjet_LHA_vs_zi, *h_genjet_LHA_vs_thetai;
    TH2F *h_qgenjet_LHA_vs_zi, *h_qgenjet_LHA_vs_thetai;
    TH2F *h_ggenjet_LHA_vs_zi, *h_ggenjet_LHA_vs_thetai;

    // lambda correlation hists
    TH2F *h_genjet_multiplicity_vs_LHA, *h_genjet_multiplicity_vs_pTD, *h_genjet_multiplicity_vs_width, *h_genjet_multiplicity_vs_thrust;
    TH2F *h_genjet_LHA_vs_pTD, *h_genjet_LHA_vs_width, *h_genjet_LHA_vs_thrust;
    TH2F *h_genjet_pTD_vs_width, *h_genjet_pTD_vs_thrust;
    TH2F *h_genjet_width_vs_thrust;
    
    TH2F *h_qgenjet_multiplicity_vs_LHA, *h_qgenjet_multiplicity_vs_pTD, *h_qgenjet_multiplicity_vs_width, *h_qgenjet_multiplicity_vs_thrust;
    TH2F *h_qgenjet_LHA_vs_pTD, *h_qgenjet_LHA_vs_width, *h_qgenjet_LHA_vs_thrust;
    TH2F *h_qgenjet_pTD_vs_width, *h_qgenjet_pTD_vs_thrust;
    TH2F *h_qgenjet_width_vs_thrust;

    TH2F *h_ggenjet_multiplicity_vs_LHA, *h_ggenjet_multiplicity_vs_pTD, *h_ggenjet_multiplicity_vs_width, *h_ggenjet_multiplicity_vs_thrust;
    TH2F *h_ggenjet_LHA_vs_pTD, *h_ggenjet_LHA_vs_width, *h_ggenjet_LHA_vs_thrust;
    TH2F *h_ggenjet_pTD_vs_width, *h_ggenjet_pTD_vs_thrust;
    TH2F *h_ggenjet_width_vs_thrust;

    std::map<std::string, LambdaVariable> lambda_variables;
    std::map<std::string, TH2F*> lambda_correlation_hists;

    int useNJets_;

    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;

    bool doHerwigReweighting;
    TH1F * reweightHist;
};

}

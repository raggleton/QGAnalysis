#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"


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
    QGAnalysisHists(uhh2::Context & ctx, const std::string & dirname, int useNJets, const std::string & selection);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisHists();
protected:
    std::vector<GenParticle*> get_genjet_genparticles(const GenJetWithParts &, std::vector<GenParticle>*);
    int get_jet_flavour(const Jet & jet, std::vector<GenParticle> * genparticles);

    // reco jet hists
    float jetRadius, LHA_rescale, width_rescale, thrust_rescale;
    TH1F * h_weights;
    TH2F * h_weights_vs_pt, *h_pthat_vs_weight, *h_pthat_vs_jet_pt;

    TH1F *h_jet_pt, *h_jet_pt_unweighted, *h_jet_eta, *h_jet_flavour, *h_jet_genParton_flavour;
    TH1F *h_jet_multiplicity, *h_jet_LHA, *h_jet_pTD, *h_jet_width, *h_jet_thrust;
    TH1F *h_qjet_multiplicity, *h_qjet_LHA, *h_qjet_pTD, *h_qjet_width, *h_qjet_thrust;
    TH1F *h_gjet_multiplicity, *h_gjet_LHA, *h_gjet_pTD, *h_gjet_width, *h_gjet_thrust;

    TH2F *h_jet_multiplicity_vs_pt, *h_jet_LHA_vs_pt, *h_jet_pTD_vs_pt, *h_jet_width_vs_pt, *h_jet_thrust_vs_pt;
    TH2F *h_jet_flavour_vs_pt, *h_jet1_flavour_vs_pt, *h_jet2_flavour_vs_pt, *h_jet_genParton_flavour_vs_pt, *h_jet1_genParton_flavour_vs_pt, *h_jet2_genParton_flavour_vs_pt;
    TH2F *h_jet_flavour_vs_eta, *h_jet_genParton_flavour_vs_eta, *h_jet_response_vs_genjet_pt;
    TH2F *h_qjet_multiplicity_vs_pt, *h_qjet_LHA_vs_pt, *h_qjet_pTD_vs_pt, *h_qjet_width_vs_pt, *h_qjet_thrust_vs_pt, *h_qjet_response_vs_genjet_pt;
    TH2F *h_gjet_multiplicity_vs_pt, *h_gjet_LHA_vs_pt, *h_gjet_pTD_vs_pt, *h_gjet_width_vs_pt, *h_gjet_thrust_vs_pt, *h_gjet_response_vs_genjet_pt;
    TH2F *h_jet_puppiMultiplicity_vs_pt, *h_qjet_puppiMultiplicity_vs_pt, *h_gjet_puppiMultiplicity_vs_pt;

    // lambda correlation hists
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
    bool useGenPartonFlav_;
    bool doHerwigReweighting;
    TH1F * reweightHist;

    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;
    bool is_mc_;
};

}

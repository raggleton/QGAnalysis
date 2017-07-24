#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"

namespace uhh2examples {

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

    // genjet hists
    TH1F *h_genjet_pt, *h_genjet_pt_unweighted, *h_genjet_pt_all, *h_genjet_ht, *h_genjet_eta, *h_genjet_flavour;
    TH2F *h_genjet_flavour_vs_pt;
    TH1F *h_genjet_multiplicity, *h_genjet_LHA, *h_genjet_pTD, *h_genjet_width, *h_genjet_thrust;
    TH1F *h_qgenjet_multiplicity, *h_qgenjet_LHA, *h_qgenjet_pTD, *h_qgenjet_width, *h_qgenjet_thrust;
    TH1F *h_ggenjet_multiplicity, *h_ggenjet_LHA, *h_ggenjet_pTD, *h_ggenjet_width, *h_ggenjet_thrust;

    TH2F *h_genjet_multiplicity_vs_pt, *h_genjet_LHA_vs_pt, *h_genjet_pTD_vs_pt, *h_genjet_width_vs_pt, *h_genjet_thrust_vs_pt;
    TH2F *h_qgenjet_multiplicity_vs_pt, *h_qgenjet_LHA_vs_pt, *h_qgenjet_pTD_vs_pt, *h_qgenjet_width_vs_pt, *h_qgenjet_thrust_vs_pt;
    TH2F *h_ggenjet_multiplicity_vs_pt, *h_ggenjet_LHA_vs_pt, *h_ggenjet_pTD_vs_pt, *h_ggenjet_width_vs_pt, *h_ggenjet_thrust_vs_pt;

    int useNJets_;

    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;

    bool doHerwigReweighting;
    TH1F * reweightHist;
};

}

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
class QGAnalysisFlavCompHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    QGAnalysisFlavCompHists(uhh2::Context & ctx, const std::string & dirname, int useNJets, const std::string & selection, float jetPtMin);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~QGAnalysisFlavCompHists();
protected:
    TH2F * h_jet_flavour_vs_genparton_flavour;

    int useNJets_;
    float jetPtMin_;
    bool doHerwigReweighting;
    TH1F * reweightHist;

};

}

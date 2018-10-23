#pragma once

#include <vector>

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/TriggerSelection.h"


namespace uhh2examples {

/*
 * Select Z->mumu + jet candidate events
 */
class ZplusJetsSelection: public uhh2::Selection {
public:
    ZplusJetsSelection(uhh2::Context & ctx, const std::string & zLabel_, float mu1_pt=20., float mu2_pt=20., float mZ_window=20., float dphi_jet_z_min=2.0, float second_jet_frac_max=0.3);
    virtual bool passes(const uhh2::Event & event) override;
private:
    uhh2::Event::Handle<std::vector<Muon>> hndlZ;
    float mu1_pt_, mu2_pt_, mZ_window_, dphi_jet_z_min_, second_jet_frac_max_;
};


/*
 * Select dijet candidate events
 */
class DijetSelection: public uhh2::Selection {
public:
    DijetSelection(float dphi_min=2.0, float second_jet_frac_max=0.94, float jet_asym_max=0.3, bool ss_eta=false, float deta_max=10, float sum_eta=10);
    virtual bool passes(const uhh2::Event & event) override;
private:
    float dphi_min_, second_jet_frac_max_, jet_asym_max_, ss_eta_, deta_max_, sum_eta_;
};


/*
 * Select event where the nth jet has suitable flavour
 */
class JetFlavourSelection: public uhh2::Selection {
public:
    JetFlavourSelection(std::vector<int> pdgids, uint jet_index=0);
    virtual bool passes(const uhh2::Event & event) override;
private:
    uint jet_index_;
    std::vector<int> flavours_;
};


/*
 * Select Z+jets events using the theory paper method
 */
class ZplusJetsTheorySelection: public uhh2::Selection {
public:
    ZplusJetsTheorySelection(uhh2::Context & ctx, float pt_min=100., float jet_frac_min=0.8, float jet_z_deta_max_=1.0, float second_jet_frac_max=999., float mZ_window=20);
    virtual bool passes(const uhh2::Event & event) override;
private:
    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;
    uhh2::Event::Handle<std::vector<GenParticle> > genMuons_handle;
    float pt_min_;
    float jet_frac_min_;
    float jet_z_deta_max_;
    float second_jet_frac_max_;
    float mZ_window_;
};


/*
 * Select dijet events using the theory paper method
 */
class DijetTheorySelection: public uhh2::Selection {
public:
    DijetTheorySelection(uhh2::Context & ctx, float pt_min=100., float jet_frac_min=0.8, float jet_deta_max=1.0, float third_frac_max=999.);
    virtual bool passes(const uhh2::Event & event) override;
private:
    uhh2::Event::Handle<std::vector<GenJetWithParts> > genJets_handle;
    float pt_min_;
    float jet_frac_min_;
    float jet_deta_max_;
    float third_frac_max_;
};


/**
 * Allows selection of events by Event number
 */
class EventNumberSelection: public uhh2::Selection {
public:
    EventNumberSelection(std::vector<unsigned long> eventNums);
    virtual bool passes(const uhh2::Event & event) override;
private:
    std::vector<unsigned long> eventNums_;
};



/**
 * Select event based on trigger and leading jet pt range
 * Also return bin index to say which passed.
 */
class DataJetSelection: public uhh2::Selection {
public:
    DataJetSelection(const std::vector<std::string> & triggers, const std::vector<std::pair<float, float>> & ptBins);
    virtual bool passes(const uhh2::Event & event) override;
    virtual int passIndex();
private:
    std::vector<TriggerSelection> trigSel_;
    std::vector<PtEtaCut> ptSel_;
    int passBinInd_;
};
}

#include "UHH2/QGAnalysis/include/QGAnalysisSelections.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>
#include <algorithm>

using namespace std;
using namespace uhh2examples;
using namespace uhh2;


ZplusJetsSelection::ZplusJetsSelection(uhh2::Context & ctx, const std::string & zLabel_, float mu1_pt, float mu2_pt, float mZ_window, float dphi_jet_z_min, float second_jet_frac_max, float z_pt_min, float z_jet_asym_max, const std::string & cutflow_hname):
    hndlZ(ctx.get_handle<std::vector<Muon>>(zLabel_)),
    mu1_pt_(mu1_pt),
    mu2_pt_(mu2_pt),
    mZ_window_(mZ_window),
    dphi_jet_z_min_(dphi_jet_z_min),
    second_jet_frac_max_(second_jet_frac_max),
    z_pt_min_(z_pt_min),
    z_jet_asym_max_(z_jet_asym_max)
    {
        // Remember to update when you update passes() ! Yes is horrible
        std::vector<std::string> descriptions = {
            "nJets>=1",
            "nZMuons>=2",
            "muon1_pt",
            "muon2_pt",
            "muon_os",
            "z_pt_min",
            "m_mumu",
            "dphi_jet_z_min",
            "second_jet_frac_max",
            "z_jet_asym_max"
        };
        cutflow_raw = new TH1D(("cf_" + cutflow_hname + "_raw").c_str(), ("Cutflow '" + cutflow_hname + "' unweighted").c_str(), descriptions.size()+1, -1, descriptions.size());
        cutflow_weighted = new TH1D(("cf_" + cutflow_hname).c_str(), ("Cutflow '" + cutflow_hname + "' using weights").c_str(), descriptions.size()+1, -1, descriptions.size());
        for(TAxis * ax : {cutflow_raw->GetXaxis(), cutflow_weighted->GetXaxis()}){
            ax->SetBinLabel(1, "all");
            for(size_t i=0; i<descriptions.size(); ++i){
                ax->SetBinLabel(i+2, descriptions[i].c_str());
            }
        }
        ctx.put("", cutflow_raw);
        ctx.put("", cutflow_weighted);

    }

bool ZplusJetsSelection::passes(const Event & event){
    int i = -1;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (event.jets->size() < 1) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    auto & zMuons = event.get(hndlZ);
    if (zMuons.size() < 2) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto & muon1 = zMuons.at(0);
    const auto & muon2 = zMuons.at(1);

    if (muon1.pt() < mu1_pt_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (muon2.pt() < mu2_pt_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    // TODO: what if > 2 muons? how to pick the Z candidate?

    if (muon1.charge() == muon2.charge()) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto z_cand = muon1.v4() + muon2.v4();

    if (z_cand.pt() < z_pt_min_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    float m_mumu = z_cand.M();
    if (fabs(m_mumu - 90) > mZ_window_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto & jet1 = event.jets->at(0);
    if (fabs(deltaPhi(jet1, z_cand)) < dphi_jet_z_min_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (event.jets->size() > 1) {
        const auto & jet2 = event.jets->at(1);
        auto jet_frac = jet2.pt() / z_cand.pt();
        if (jet_frac > second_jet_frac_max_) return false;
    }
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    float asym = (jet1.pt() - z_cand.pt()) / (jet1.pt() + z_cand.pt());
    if (asym > z_jet_asym_max_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    return true;
}


ZplusJetsGenSelection::ZplusJetsGenSelection(uhh2::Context & ctx,
                                             float mu1_pt,
                                             float mu2_pt,
                                             float mZ_window,
                                             float dphi_jet_z_min,
                                             float second_jet_frac_max,
                                             float z_pt_min,
                                             float z_jet_asym_max,
                                             const std::string & cutflow_hname,
                                             const std::string & genjet_name,
                                             const std::string & genmuon_name):
    genJets_handle(ctx.get_handle<std::vector<GenJetWithParts>>(genjet_name)),
    zMuons_handle(ctx.get_handle<std::vector<GenParticle>>(genmuon_name)),
    mu1_pt_(mu1_pt),
    mu2_pt_(mu2_pt),
    mZ_window_(mZ_window),
    dphi_jet_z_min_(dphi_jet_z_min),
    second_jet_frac_max_(second_jet_frac_max),
    z_pt_min_(z_pt_min),
    z_jet_asym_max_(z_jet_asym_max)
    {
        // Remember to update when you update passes() ! Yes is horrible
        std::vector<std::string> descriptions = {
            "nGenJets>=1",
            "nZMuons>=2",
            "muon1_pt",
            "muon2_pt",
            "muon_os",
            "z_pt_min",
            "m_mumu",
            "dphi_jet_z_min",
            "second_jet_frac_max",
            "z_jet_asym_max"
        };
        cutflow_raw = new TH1D(("cf_" + cutflow_hname + "_raw").c_str(), ("Cutflow '" + cutflow_hname + "' unweighted").c_str(), descriptions.size()+1, -1, descriptions.size());
        cutflow_weighted = new TH1D(("cf_" + cutflow_hname).c_str(), ("Cutflow '" + cutflow_hname + "' using weights").c_str(), descriptions.size()+1, -1, descriptions.size());
        for(TAxis * ax : {cutflow_raw->GetXaxis(), cutflow_weighted->GetXaxis()}){
            ax->SetBinLabel(1, "all");
            for(size_t i=0; i<descriptions.size(); ++i){
                ax->SetBinLabel(i+2, descriptions[i].c_str());
            }
        }
        ctx.put("", cutflow_raw);
        ctx.put("", cutflow_weighted);
    }

bool ZplusJetsGenSelection::passes(const Event & event){
    const auto & jets = event.get(genJets_handle);
    int i = -1;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (jets.size() < 1) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto & zMuons = event.get(zMuons_handle);
    if (zMuons.size() < 2) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto & muon1 = zMuons.at(0);
    const auto & muon2 = zMuons.at(1);

    if (muon1.pt() < mu1_pt_ ) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (muon2.pt() < mu2_pt_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    // TODO: what if > 2 muons? how to pick the Z candidate?

    if (muon1.charge() == muon2.charge()) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto z_cand = muon1.v4() + muon2.v4();

    if (z_cand.pt() < z_pt_min_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    float m_mumu = z_cand.M();
    if (fabs(m_mumu - 90) > mZ_window_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto & jet1 = jets.at(0);
    if (fabs(deltaPhi(jet1, z_cand)) < dphi_jet_z_min_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (jets.size() > 1) {
        const auto & jet2 = jets.at(1);
        auto jet_frac = jet2.pt() / z_cand.pt();
        if (jet_frac > second_jet_frac_max_) return false;
    }
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    float asym = (jet1.pt() - z_cand.pt()) / (jet1.pt() + z_cand.pt());
    if (asym > z_jet_asym_max_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    return true;
}


DijetSelection::DijetSelection(uhh2::Context & ctx,
                               float jet_pt_min,
                               float jet_y_max,
                               float dphi_min,
                               float second_jet_frac_max,
                               float jet_asym_max,
                               bool ss_eta,
                               float deta_max,
                               float sum_eta,
                               const std::string & cutflow_hname):
    jet_pt_min_(jet_pt_min),
    jet_y_max_(jet_y_max),
    dphi_min_(dphi_min),
    second_jet_frac_max_(second_jet_frac_max),
    jet_asym_max_(jet_asym_max),
    ss_eta_(ss_eta),
    deta_max_(deta_max),
    sum_eta_(sum_eta)
    {
        // Remember to update when you update passes() ! Yes is horrible
        std::vector<std::string> descriptions = {
            "nJets>=2",
            "jet_pt_min",
            "jet_y_max",
            "second_jet_frac_max",
            "dphi_min",
            "ss_eta",
            "deta_max",
            "jet_asym_max"
        };
        cutflow_raw = new TH1D(("cf_" + cutflow_hname + "_raw").c_str(), ("Cutflow '" + cutflow_hname + "' unweighted").c_str(), descriptions.size()+1, -1, descriptions.size());
        cutflow_weighted = new TH1D(("cf_" + cutflow_hname).c_str(), ("Cutflow '" + cutflow_hname + "' using weights").c_str(), descriptions.size()+1, -1, descriptions.size());
        for(TAxis * ax : {cutflow_raw->GetXaxis(), cutflow_weighted->GetXaxis()}){
            ax->SetBinLabel(1, "all");
            for(size_t i=0; i<descriptions.size(); ++i){
                ax->SetBinLabel(i+2, descriptions[i].c_str());
            }
        }
        ctx.put("", cutflow_raw);
        ctx.put("", cutflow_weighted);
    }

bool DijetSelection::passes(const Event & event){
    int i = -1;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (event.jets->size() < 2) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto & jet1 = event.jets->at(0);
    const auto & jet2 = event.jets->at(1);

    // Should be pt-ordered
    if (jet1.pt() < jet2.pt())
        throw runtime_error("jet1.pt < jet2.pt, expecting descending pt order");
    // so only need to check jet2
    if (jet2.pt() < jet_pt_min_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if ((fabs(jet1.Rapidity()) > jet_y_max_) || (fabs(jet2.Rapidity()) > jet_y_max_)) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    float ratio = (jet1.pt() > jet2.pt()) ? jet2.pt() / jet1.pt() : jet1.pt() / jet2.pt();
    if (ratio > second_jet_frac_max_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    auto dphi = fabs(deltaPhi(jet1, jet2));
    if (dphi < dphi_min_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    auto eta1 = jet1.Rapidity();
    auto eta2 = jet2.Rapidity();

    if (ss_eta_ && ((eta1 * eta2) < 0)) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (fabs(eta1 - eta2) > deta_max_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if ((fabs(jet1.pt() - jet2.pt())/(jet1.pt() + jet2.pt())) > jet_asym_max_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    return true;
}

DijetGenSelection::DijetGenSelection(uhh2::Context & ctx,
                                     float jet_pt_min,
                                     float jet_y_max,
                                     float dphi_min,
                                     float second_jet_frac_max,
                                     float jet_asym_max,
                                     bool ss_eta,
                                     float deta_max,
                                     float sum_eta,
                                     const std::string & cutflow_hname,
                                     const std::string & genjet_name):
    genJets_handle(ctx.get_handle<std::vector<GenJetWithParts>>(genjet_name)),
    jet_pt_min_(jet_pt_min),
    jet_y_max_(jet_y_max),
    dphi_min_(dphi_min),
    second_jet_frac_max_(second_jet_frac_max),
    jet_asym_max_(jet_asym_max),
    ss_eta_(ss_eta),
    deta_max_(deta_max),
    sum_eta_(sum_eta)
    {
        // Remember to update when you update passes() ! Yes is horrible
        std::vector<std::string> descriptions = {
            "nGenJets>=2",
            "jet_pt_min",
            "jet_y_max",
            "second_jet_frac_max",
            "dphi_min",
            "ss_eta",
            "deta_max",
            "jet_asym_max"
        };
        cutflow_raw = new TH1D(("cf_" + cutflow_hname + "_raw").c_str(), ("Cutflow '" + cutflow_hname + "' unweighted").c_str(), descriptions.size()+1, -1, descriptions.size());
        cutflow_weighted = new TH1D(("cf_" + cutflow_hname).c_str(), ("Cutflow '" + cutflow_hname + "' using weights").c_str(), descriptions.size()+1, -1, descriptions.size());
        for(TAxis * ax : {cutflow_raw->GetXaxis(), cutflow_weighted->GetXaxis()}){
            ax->SetBinLabel(1, "all");
            for(size_t i=0; i<descriptions.size(); ++i){
                ax->SetBinLabel(i+2, descriptions[i].c_str());
            }
        }
        ctx.put("", cutflow_raw);
        ctx.put("", cutflow_weighted);
    }

bool DijetGenSelection::passes(const Event & event){
    const auto & jets = event.get(genJets_handle);
    int i = -1;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (jets.size() < 2) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    const auto & jet1 = jets.at(0);
    const auto & jet2 = jets.at(1);

    // Should be pt-ordered
    if (jet1.pt() < jet2.pt())
        throw runtime_error("jet1.pt < jet2.pt, expecting descending pt order");
    // so only need to check jet2
    if (jet2.pt() < jet_pt_min_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if ((fabs(jet1.Rapidity()) > jet_y_max_) || (fabs(jet2.Rapidity()) > jet_y_max_)) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    float ratio = (jet1.pt() > jet2.pt()) ? jet2.pt() / jet1.pt() : jet1.pt() / jet2.pt();
    if (ratio > second_jet_frac_max_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    auto dphi = fabs(deltaPhi(jet1, jet2));
    if (dphi < dphi_min_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    auto eta1 = jet1.Rapidity();
    auto eta2 = jet2.Rapidity();

    if (ss_eta_ && ((eta1 * eta2) < 0)) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if (fabs(eta1 - eta2) > deta_max_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    if ((fabs(jet1.pt() - jet2.pt())/(jet1.pt() + jet2.pt())) > jet_asym_max_) return false;
    i++;
    cutflow_raw->Fill(i);
    cutflow_weighted->Fill(i, event.weight);

    return true;
}


JetFlavourSelection::JetFlavourSelection(std::vector<int> pdgids, uint jet_index): jet_index_(jet_index), flavours_(pdgids) {}

bool JetFlavourSelection::passes(const Event & event) {
    if (event.jets->size() <= jet_index_) return false;
    const auto & jet = event.jets->at(jet_index_);
    uint flav = abs(jet.genPartonFlavor());
    return std::find(flavours_.begin(), flavours_.end(), flav) != flavours_.end();
}


ZplusJetsTheorySelection::ZplusJetsTheorySelection(Context & ctx, float pt_min, float jet_frac_min, float jet_z_deta_max, float second_jet_frac_max, float mZ_window):
    pt_min_(pt_min),
    jet_frac_min_(jet_frac_min),
    jet_z_deta_max_(jet_z_deta_max),
    second_jet_frac_max_(second_jet_frac_max),
    mZ_window_(mZ_window) {
        genJets_handle = ctx.get_handle<std::vector<GenJetWithParts>>("GoodGenJets");
        genMuons_handle = ctx.get_handle<std::vector<GenParticle>>("GoodGenMuons");
    }

bool ZplusJetsTheorySelection::passes(const Event & event) {
    if(!event.is_valid(genJets_handle) || !event.is_valid(genMuons_handle)) return false;

    const auto & genjets = event.get(genJets_handle);
    const auto & muons = event.get(genMuons_handle);

    if ((genjets.size() < 1) || (muons.size() < 2)) return false;

    // Get Z candidate
    const auto & muon1 = muons.at(0);
    const auto & muon2 = muons.at(1);

    if (muon1.charge() == muon2.charge()) return false;

    const auto z_cand = muon1.v4() + muon2.v4();

    float m_mumu = z_cand.M();
    if (fabs(m_mumu - 90) > mZ_window_) return false;

    if (z_cand.pt() < pt_min_) return false;

    const auto & jet = genjets.at(0);

    if ((jet.pt() / z_cand.pt()) < jet_frac_min_) return false;

    if (fabs(jet.Rapidity() - z_cand.Rapidity()) > jet_z_deta_max_) return false;

    if (genjets.size() >=2) {
        const auto & jet2 = genjets.at(1);
        if ((jet2.pt() / z_cand.pt()) > second_jet_frac_max_) return false;
    }

    return true;
}


DijetTheorySelection::DijetTheorySelection(Context & ctx, float pt_min, float jet_frac_min, float jet_deta_max, float third_frac_max):
    pt_min_(pt_min),
    jet_frac_min_(jet_frac_min),
    jet_deta_max_(jet_deta_max),
    third_frac_max_(third_frac_max) {
        genJets_handle = ctx.get_handle<std::vector<GenJetWithParts>>("GoodGenJets");
    }

bool DijetTheorySelection::passes(const Event & event) {
    if(!event.is_valid(genJets_handle)) return false;

    const auto & genjets = event.get(genJets_handle);

    if (genjets.size() < 2) return false;

    const auto & jet1 = genjets.at(0);
    const auto & jet2 = genjets.at(1);

    if ((jet1.pt() + jet2.pt()) < (2. * pt_min_)) return false;

    if (jet2.pt() < (jet1.pt() * jet_frac_min_)) return false;

    if (fabs(jet1.Rapidity() - jet2.Rapidity()) > jet_deta_max_) return false;

    if (genjets.size()>=3) {
        const auto & jet3 = genjets.at(2);
        if ((jet3.pt() / (0.5*(jet1.pt()+jet2.pt()))) > third_frac_max_) return false;
    }
    return true;
}


EventNumberSelection::EventNumberSelection(std::vector<unsigned long> eventNums) {
    eventNums_ = eventNums;
}

bool EventNumberSelection::passes(const Event & event) {
    return std::find(eventNums_.begin(), eventNums_.end(), event.event) != eventNums_.end();
}


DataJetSelection::DataJetSelection(const std::vector<std::string> & triggers, const std::vector<std::pair<float, float>> & ptBins):
    passBinInd_(-1)
{
    for (const auto & itr : triggers) {
        TriggerSelection trig(itr);
        std::cout << "Adding trigger " << itr << std::endl;
        trigSel_.push_back(trig);
    }
    for (const auto & ptValues : ptBins) {
        if (ptValues.second < ptValues.first) {
            throw std::runtime_error("trigger thresholds must be ordered ascending");
        }
        std::cout << "Adding pt cut " << ptValues.first << " - " << ptValues.second << std::endl;
        PtEtaCut ptcut(ptValues.first, 999, ptValues.second, -999);
        ptSel_.push_back(ptcut);
    }
}

bool DataJetSelection::passes(const Event & event) {
    int ind = -1;
    // Do basic check to see if leading jet exists
    if (event.jets->size() == 0) return false;
    // iterate in reverse so as to get the highest trigger threshold that fires
    // (assumes you loaded the triggers in ascending order!)
    // std::cout << "Leading jet: "<< event.jets->at(0).pt() << " : " << event.jets->at(0).Rapidity() << std::endl;
    for (ind = trigSel_.size()-1; ind >= 0; ind--) {
        // if (trigSel_[ind].passes(event)) std::cout << "Fired " << ind << std::endl;
        if (trigSel_[ind].passes(event) && ptSel_[ind](event.jets->at(0), event)) {
            passBinInd_ = ind;
            return true;
        }
    }
    passBinInd_ = ind;
    return false;
}

int DataJetSelection::passIndex() {
    return passBinInd_;
}

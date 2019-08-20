#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"


using namespace std;
using namespace uhh2;
using namespace fastjet;


DataJetMetCorrector::DataJetMetCorrector(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone){
  std::vector<std::string> JEC_BCD, JEC_EFearly, JEC_GH;
  if (pu_removal == "CHS") {
    if (jet_cone == "AK4") {
      JEC_BCD = JERFiles::Summer16_07Aug2017_V11_BCD_L123_AK4PFchs_DATA;
      JEC_EFearly = JERFiles::Summer16_07Aug2017_V11_EF_L123_AK4PFchs_DATA;
      JEC_GH = JERFiles::Summer16_07Aug2017_V11_GH_L123_AK4PFchs_DATA;
    } else if (jet_cone == "AK8") {
      JEC_BCD = JERFiles::Summer16_07Aug2017_V11_BCD_L123_AK8PFchs_DATA;
      JEC_EFearly = JERFiles::Summer16_07Aug2017_V11_EF_L123_AK8PFchs_DATA;
      JEC_GH = JERFiles::Summer16_07Aug2017_V11_GH_L123_AK8PFchs_DATA;
    } else {
      throw runtime_error("CHS must have jet_cone of AK4 or AK8");
    }
  } else if (pu_removal == "PUPPI") {
    if (jet_cone == "AK4") {
      JEC_BCD = JERFiles::Summer16_07Aug2017_V11_BCD_L123_AK4PFPuppi_DATA;
      JEC_EFearly = JERFiles::Summer16_07Aug2017_V11_EF_L123_AK4PFPuppi_DATA;
      JEC_GH = JERFiles::Summer16_07Aug2017_V11_GH_L123_AK4PFPuppi_DATA;
    } else if (jet_cone == "AK8") {
      JEC_BCD = JERFiles::Summer16_07Aug2017_V11_BCD_L123_AK8PFPuppi_DATA;
      JEC_EFearly = JERFiles::Summer16_07Aug2017_V11_EF_L123_AK8PFPuppi_DATA;
      JEC_GH = JERFiles::Summer16_07Aug2017_V11_GH_L123_AK8PFPuppi_DATA;
    } else {
      throw runtime_error("PUPPI must have jet_cone of AK4 or AK8");
    }
  }

  jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_BCD));
  jet_corrector_EFearly.reset(new JetCorrector(ctx, JEC_EFearly));
  jet_corrector_GH.reset(new JetCorrector(ctx, JEC_GH));
}

bool DataJetMetCorrector::process(uhh2::Event & event) {
  if (event.run <= runnr_BCD) {
    jet_corrector_BCD->process(event);
    jet_corrector_BCD->correct_met(event);
  } else if (event.run <= runnr_EFearly) { //< is correct, not <=
    jet_corrector_EFearly->process(event);
    jet_corrector_EFearly->correct_met(event);
  } else if (event.run > runnr_EFearly) {
    jet_corrector_GH->process(event);
    jet_corrector_GH->correct_met(event);
  } else {
    throw runtime_error("Run number not covered by if-statements in DataJetMetCorrector.");
  }
  return true;
}


MCJetMetCorrector::MCJetMetCorrector(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, const std::string & jet_coll_name, const std::string & genjet_coll_name) {
  std::vector<std::string> JEC_MC;
  std::string resolutionFilename;
  if (pu_removal == "CHS") {
    if (jet_cone == "AK4") {
      JEC_MC = JERFiles::Summer16_07Aug2017_V11_L123_AK4PFchs_MC;
      resolutionFilename = "Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
    } else if (jet_cone == "AK8") {
      JEC_MC = JERFiles::Summer16_07Aug2017_V11_L123_AK8PFchs_MC;
      resolutionFilename = "Summer16_25nsV1_MC_PtResolution_AK8PFchs.txt";
    } else {
      throw runtime_error("CHS must have jet_cone of AK4 or AK8");
    }
  } else if (pu_removal == "PUPPI") {
    if (jet_cone == "AK4") {
      JEC_MC = JERFiles::Summer16_07Aug2017_V11_L123_AK4PFPuppi_MC;
      resolutionFilename = "Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
    } else if (jet_cone == "AK8") {
      JEC_MC = JERFiles::Summer16_07Aug2017_V11_L123_AK8PFPuppi_MC;
      resolutionFilename = "Summer16_25nsV1_MC_PtResolution_AK8PFchs.txt";
    } else {
      throw runtime_error("PUPPI must have jet_cone of AK4 or AK8");
    }
  }
  jet_corrector.reset(new JetCorrector(ctx, JEC_MC));
  jet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, jet_coll_name, genjet_coll_name, true, JERSmearing::SF_13Tev_Summer16_25nsV1, resolutionFilename));
}

bool MCJetMetCorrector::process(uhh2::Event & event) {
  jet_corrector->process(event);
  jet_corrector->correct_met(event);
  jet_resolution_smearer->process(event);
  return true;
}


GeneralEventSetup::GeneralEventSetup(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, float jet_radius, float jet_pt_min) {
  is_mc = ctx.get("dataset_type") == "MC";

  if (!is_mc) lumi_selection.reset(new LumiSelection(ctx));

  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("globalTightHalo2016Filter", "Flag_globalTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
  PrimaryVertexId pvid = StandardPrimaryVertexId();
  metfilters_selection->add<NPVSelection>("1 good PV", 1, -1, pvid);

  pv_cleaner.reset(new PrimaryVertexCleaner(pvid));

  electron_cleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_Spring16_medium, PtEtaCut(20.0, 2.5))));

  muon_cleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(26.0, 2.4), MuonIso(0.25))));
  // muon_cleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(26.0, 2.4))));

  jet_pf_id.reset(new JetCleaner(ctx, JetPFID(JetPFID::wp::WP_LOOSE)));

  if (is_mc) {
    jet_met_corrector.reset(new MCJetMetCorrector(ctx, pu_removal, jet_cone));
  } else {
    jet_met_corrector.reset(new DataJetMetCorrector(ctx, pu_removal, jet_cone));
  }

  // eta cut to allow for all of the jet to fall inside tracker (<2.5)
  jet_cleaner.reset(new JetCleaner(ctx, PtEtaCut(jet_pt_min, 2.5-jet_radius)));

  jet_ele_cleaner.reset(new JetElectronOverlapRemoval(jet_radius));

  jet_mu_cleaner.reset(new JetMuonOverlapRemoval(jet_radius));
}

bool GeneralEventSetup::process(uhh2::Event & event) {

  if(event.isRealData && !lumi_selection->passes(event)) return false;

  pv_cleaner->process(event);

  electron_cleaner->process(event);

  muon_cleaner->process(event);

  jet_pf_id->process(event);

  jet_met_corrector->process(event);

  jet_cleaner->process(event);

  jet_ele_cleaner->process(event);

  jet_mu_cleaner->process(event);

  sort_by_pt(*event.jets);

  // put this last, so objects are correctly cleaned, etc, for MCWeight afterwards
  if(!metfilters_selection->passes(event)) return false;

  return true;
}


MCReweighting::MCReweighting(uhh2::Context & ctx) {
  gen_weight_handle = ctx.get_handle<double>("gen_weight");

  lumi_weighter.reset(new MCLumiWeight(ctx));

  pileup_reweighter.reset(new MCPileupReweight(ctx, ctx.get("pileup_direction", "central")));

  std::string sf_path_name = locate_file("common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root");
  std::string sf_name = "MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta";
  muon_id_reweighter_pt_eta.reset(new MCMuonScaleFactor(ctx, sf_path_name, sf_name, 100));

  sf_name = "MC_NUM_MediumID2016_DEN_genTracks_PAR_vtx";
  // Doesn't work
  // muon_id_reweighter_vtx.reset(new MCMuonScaleFactor(ctx, sf_path_name, sf_name, 100));

  sf_path_name = locate_file("common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root");
  sf_name = "IsoMu24_OR_IsoTkMu24_PtEtaBins";
  muon_trg_reweighter.reset(new MCMuonScaleFactor(ctx, sf_path_name, sf_name, 100));

  std::string trk_path_name = locate_file("common/data/general_eff_aeta_dr030e030_corr_ratio.txt");
  muon_trk_reweighter.reset(new MCMuonTrkScaleFactor(ctx, trk_path_name, 100));
}


bool MCReweighting::process(uhh2::Event & event) {
  // store gen-only weights in separate variable
  // event.weight is then the product of reco weights & gen weights
  double old_gen_weight = event.get(gen_weight_handle);
  double old_weight = event.weight;
  if (!event.isRealData){
    lumi_weighter->process(event);
    old_gen_weight *= (event.weight / old_weight);

    pileup_reweighter->process(event);

    muon_id_reweighter_pt_eta->process(event);

    // muon_id_reweighter_vtx->process(event);

    muon_trg_reweighter->process(event);

    muon_trk_reweighter->process(event);
  }
  event.set(gen_weight_handle, old_gen_weight);
  return true;
}


float get_jet_radius(const std::string & jet_cone) {
  if (jet_cone.find("AK4") != std::string::npos)
    return 0.4;
  else if (jet_cone.find("AK8") != std::string::npos)
    return 0.8;
  else if (jet_cone.find("ca15") != std::string::npos)
    return 1.5;
  else
    throw std::runtime_error("Cannot determine jet radius from jet_cone string");
}


ZFinder::ZFinder(uhh2::Context & ctx, const std::string & inputLabel_, const std::string & outputLabel_, const std::string & weightFilename_):
  // Would like more generic FlavourParticle handle, but may need to do additional declare_event_input?
  hndlInput(ctx.get_handle<vector<Muon>>(inputLabel_)),
  hndlZ(ctx.get_handle<vector<Muon>>(outputLabel_)),
  z_weight_handle(ctx.get_handle<double>("z_weight"))

{
  if (weightFilename_ != "")
    zReweight.reset(new ZllKFactor(weightFilename_));
}

bool ZFinder::process(uhh2::Event & event) {
  // Reweight to higher order cross-section using gen level Z pT
  if (zReweight) {
    // double realZPt = zCand.pt();
    double realZPt = 0;
    for (const auto & itr : *event.genparticles) {
      if (itr.pdgId() == 23) {
        realZPt = itr.pt();
        // cout << "Found Gen Z with pT " << realZPt << " with status " << itr.status() << endl;
        break;
        // use 1st Z occurence in GenParticles (should be status 22)
        // later Z is evolved (status 44) and then finalised (status 62),
        // but I guess we want to use the 1st occurence?
      }
    }
    double zWeight = 1;
    if (realZPt > 0) zWeight = zReweight->getKFactor(realZPt);
    event.set(z_weight_handle, zWeight);
  }

  // Now look for a reconstructed Z using inputs
  auto inputs = event.get(hndlInput);
  if (inputs.size() < 2) return false;
  // Do we also want to consider more than leading & subleading?
  // Prob v.v.litle diff as not often > 2 leptons that pass selection
  auto zCand = inputs[0].v4() + inputs[1].v4();
  if ((fabs(zCand.M() - 90) < 20) && (inputs[0].charge() * inputs[1].charge() < 0)) {
    std::vector<Muon> cands = {inputs[0], inputs[1]};
    event.set(hndlZ, cands);
    return true;
  }
  return false;
}


float calcGenHT(const std::vector<GenParticle> & gps) {
  float ht = 0.;
  for (const auto & itr: gps) {
    if (abs(itr.status()) != 23) continue;
    uint pdg = abs(itr.pdgId());
    if (( pdg <= PDGID::TOP_QUARK && pdg >= PDGID::DOWN_QUARK) || pdg == PDGID::GLUON) {
      ht += itr.pt();
    }
  }
  return ht;
}

ZllKFactor::ZllKFactor(const std::string & weightFilename_)
{
  file.reset(TFile::Open(locate_file(weightFilename_).c_str()));
  grNNLO.reset((TGraph*) file->Get("kNNLO"));
}

float ZllKFactor::getKFactor(float zPt) {
  float factor = 1.;
  if (zPt > grNNLO->GetXaxis()->GetXmin() && zPt < grNNLO->GetXaxis()->GetXmax()) {
    factor = grNNLO->Eval(zPt);
  } else if (zPt < grNNLO->GetXaxis()->GetXmin()) {
    factor = grNNLO->Eval(grNNLO->GetXaxis()->GetXmin());
  } else {
    factor = grNNLO->Eval(grNNLO->GetXaxis()->GetXmax());
  }
  return factor;
}



PtReweight::PtReweight(uhh2::Context & ctx, const std::string & selection, const std::string & weightFilename):
  gen_weight_handle(ctx.get_handle<double>("gen_weight"))
{
  f_weight.reset(TFile::Open(weightFilename.c_str()));
  if (selection == "dijet")
    reweightHist.reset((TH1F*) f_weight->Get("dijet_reco"));
  else if (selection == "zplusjets")
    reweightHist.reset((TH1F*) f_weight->Get("zpj_reco"));
  if (reweightHist) reweightHist->SetDirectory(0);
}


bool PtReweight::process(uhh2::Event & event, float value) {
  (void) event;
  if (value >= reweightHist->GetXaxis()->GetXmax()) {
    value = reweightHist->GetXaxis()->GetXmax() - 0.1;
  }
  int bin_num = reweightHist->GetXaxis()->FindBin(value);
  double new_weight = reweightHist->GetBinContent(bin_num);

  // Update the event weight & gen_weight stored in the event
  double gen_weight = event.get(gen_weight_handle);
  event.weight *= new_weight;
  event.set(gen_weight_handle, gen_weight*new_weight);

  return true;
}

// namespace uhh2examples {


template<class T>
LambdaCalculator<T>::LambdaCalculator(std::vector<T> & daughters, float jet_radius, const LorentzVector & jet_vector, bool usePuppiWeight):
  daughters_(daughters),
  jetRadius_(jet_radius),
  ptSum_(0),
  jetVector_(jet_vector),
  usePuppiWeight_(usePuppiWeight)
{
  // cache pt_sum
  for (auto & dtr : daughters_) {
    float weight = usePuppiWeight_ ? dtr.puppiWeight() : 1.;
    ptSum_ += weight*dtr.pt();
  }
  // cout << "calc ptSum: " << ptSum_ << endl;
}

template<>
LambdaCalculator<GenParticle>::LambdaCalculator(std::vector<GenParticle> & daughters, float jet_radius, const LorentzVector & jet_vector, bool usePuppiWeight):
  daughters_(daughters),
  jetRadius_(jet_radius),
  ptSum_(0),
  jetVector_(jet_vector),
  usePuppiWeight_(usePuppiWeight)
{
  // cache pt_sum
  for (auto & dtr : daughters_) {
    ptSum_ += dtr.pt();
  }
}

// TODO: unify this into one calculation, otherwise kinda useless
template<class T>
float LambdaCalculator<T>::getLambda(float kappa, float beta)
{
  // Check if result already exists in cache
  auto thisArgs = std::make_pair(kappa, beta);
  auto findRes = resultsCache_.find(thisArgs);
  if (findRes != resultsCache_.end()) {
    return findRes->second;
  }

  // If not, calculate it and store in cache
  // Special case if both 0 ie multiplicity
  // Do it this way to ensure puppi weights correctly accounted for
  float result = 0.;
  if (kappa == 0 && beta == 0) {
    if (usePuppiWeight_) {
      for (auto dtr : daughters_) {
        result += dtr.puppiWeight();
      }
    } else {
      result = daughters_.size();
    }
    resultsCache_[thisArgs] = result;
    return result;
  }

  for (auto dtr : daughters_) {
    float weight = usePuppiWeight_ ? dtr.puppiWeight() : 1.;
    float z = (kappa != 0) ? dtr.pt() / ptSum_ : 1.;
    z *= weight;
    float theta = (beta != 0) ? deltaR(dtr.v4(), jetVector_) / jetRadius_ : 1.; // 1 as puppi doesn't change direction
    result += (pow(z, kappa) * pow(theta, beta));
  }
  resultsCache_[thisArgs] = result;
  return result;
}

template<>
float LambdaCalculator<GenParticle>::getLambda(float kappa, float beta)
{
  // Check if result already exists in cache
  auto thisArgs = std::make_pair(kappa, beta);
  auto findRes = resultsCache_.find(thisArgs);
  if (findRes != resultsCache_.end()) {
    return findRes->second;
  }

  // If not, calculate it and store in cache
  float result = 0.;
  for (auto dtr : daughters_) {
    float z = (kappa != 0) ? dtr.pt() / ptSum_ : 1.;
    float theta = (beta != 0) ? deltaR(dtr.v4(), jetVector_) / jetRadius_ : 1.;
    result += (pow(z, kappa) * pow(theta, beta));
  }
  resultsCache_[thisArgs] = result;
  return result;
}

template<class T>
void LambdaCalculator<T>::clearCache()
{
  resultsCache_.clear();
}



QGAnalysisJetLambda::QGAnalysisJetLambda(uhh2::Context & ctx, float jetRadius, int nJetsMax, bool doPuppi, const PFParticleId & pfId, const std::string & jet_coll_name, const std::string & output_coll_name):
  ca_wta_cluster_(JetDefinition(cambridge_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme)),
  jetRadius_(jetRadius),
  nJetsMax_(nJetsMax),
  doPuppi_(doPuppi),
  pfId_(pfId),
  neutralHadronShift_(0), // these are the fractional shift, e.g. 0.2 for 20%
  photonShift_(0), // must init to 0 otherwise bad things will happen, will go haywire!
  jet_handle(ctx.get_handle<std::vector<Jet>>(jet_coll_name)),
  output_handle(ctx.get_handle<std::vector<JetLambdaBundle>>(output_coll_name))
{}


PseudoJet QGAnalysisJetLambda::convert_uhh_pfparticle_to_pseudojet(const PFParticle & particle, bool applyPuppiWeight) {
  float weight = (applyPuppiWeight) ? particle.puppiWeight() : 1.;
  LorentzVectorXYZE lv = toXYZ(particle.v4()) * weight;
  PseudoJet outputParticle(lv.px(), lv.py(), lv.pz(), lv.E());
  return outputParticle;
}

bool QGAnalysisJetLambda::process(uhh2::Event & event) {
  std::vector<Jet> jets = event.get(jet_handle);
  std::vector<JetLambdaBundle> outputs;
  int nJetCounter = 0;
  for (auto & jet : jets) {
    if (nJetsMax_ > 0 && nJetCounter == nJetsMax_) break;
    // Get constituents
    // don't apply puppi weight here, since we already account for it in the LambdaCalculator
    std::vector<PFParticle> daughters = get_jet_pfparticles(jet, event, false);
    // Shift energies if appropriate
    if (fabs(neutralHadronShift_) > 1E-6) { shift_neutral_hadron_pfparticles(daughters, neutralHadronShift_); }
    if (fabs(photonShift_) > 1E-6) { shift_photon_pfparticles(daughters, photonShift_); }

    // Calculate the WTA axis
    // First convert Jet to PseudoJet, then convert it
    // PseudoJet origJet = convert_uhh_jet_to_pseudojet(jet);
    vector<PseudoJet> constits = {};
    int dauCounter = 0;
    for (const auto & dau : daughters) {
      PseudoJet thisDau = convert_uhh_pfparticle_to_pseudojet(dau, doPuppi_); // apply puppi weight here to recluster
      thisDau.set_user_index(dauCounter);
      constits.push_back(thisDau);
      dauCounter++;
    }
    ClusterSequence cs(constits, ca_wta_cluster_);
    vector<PseudoJet> wtaJets = sorted_by_pt(cs.inclusive_jets());
    if (wtaJets.size() > 1) {
      cout << " >1 WTA jets" << endl;
    } else if (wtaJets.size() == 0) {
      cout << constits.size() << " daughters" << endl;
      throw std::runtime_error("WTA reclustering failed - no jets");
    }
    LorentzVectorXYZE wtaJetAxis(wtaJets[0].px(), wtaJets[0].py(), wtaJets[0].pz(), wtaJets[0].E());

    // cout << "Original jet axis: " << jet.v4().px() << " : " << jet.v4().py() << " : " << jet.v4().pz() << endl;
    // cout << "Original jet axis: " << jet.eta() << " : " << jet.phi() << endl;
    // cout << "WTA jet axis: " << wtaJetAxis.px() << " : " << wtaJetAxis.py() << " : " << wtaJetAxis.pz() << endl;
    // cout << "WTA jet axis: " << wtaJetAxis.eta() << " : " << wtaJetAxis.phi() << endl;

    // Finally apply any pt cuts etc
    if (pfId_) clean_collection<PFParticle>(daughters, event, pfId_);
    LambdaCalculator<PFParticle> recoJetCalc(daughters, jetRadius_, toPtEtaPhi(wtaJetAxis), doPuppi_);
    JetLambdaBundle thisBundle{jet, recoJetCalc};
    outputs.push_back(thisBundle);
    nJetCounter++;
  }
  event.set(output_handle, std::move(outputs));
  return true;
}

std::vector<PFParticle> QGAnalysisJetLambda::get_jet_pfparticles(const Jet & jet, uhh2::Event & event, bool applyPuppiWeight) {
  std::vector<PFParticle> * pfparticles = event.pfparticles;
  std::vector<PFParticle> pfCopy;
  // Create a copy, since we might modify it later and we want to keep the original event.pfparticles intact
  for (const uint i : jet.daughterIndices()) {
    PFParticle newPf;
    newPf.set_particleID(pfparticles->at(i).particleID());
    newPf.set_puppiWeight(pfparticles->at(i).puppiWeight());
    if (applyPuppiWeight) {
      LorentzVectorXYZE v4XYZ = toXYZ(pfparticles->at(i).v4());
      v4XYZ * pfparticles->at(i).puppiWeight();
      newPf.set_v4(toPtEtaPhi(v4XYZ));
    } else {
      newPf.set_v4(pfparticles->at(i).v4());
    }
    newPf.set_charge(pfparticles->at(i).charge());
    pfCopy.push_back(newPf);
  }
  return pfCopy;
}

void QGAnalysisJetLambda::set_neutral_hadron_shift(int direction, float rel_shift) {
  neutralHadronShift_ = (direction * rel_shift);
}

void QGAnalysisJetLambda::shift_neutral_hadron_pfparticles(std::vector<PFParticle> & pfparticles, float shift) {
  for (auto & itr : pfparticles) {
    if (itr.particleID() == PFParticle::eH0) {
      LorentzVectorXYZE lv = toXYZ(itr.v4());
      lv *= (1+shift);
      LorentzVector newLv = toPtEtaPhi(lv);
      itr.set_pt(newLv.pt());
      itr.set_eta(newLv.eta());
      itr.set_phi(newLv.phi());
      itr.set_energy(newLv.energy());
    }
  }
}

void QGAnalysisJetLambda::set_photon_shift(int direction, float rel_shift) {
  photonShift_ = (direction * rel_shift);
}

void QGAnalysisJetLambda::shift_photon_pfparticles(std::vector<PFParticle> & pfparticles, float shift) {
  for (auto & itr : pfparticles) {
    if (itr.particleID() == PFParticle::eGamma) {
      LorentzVectorXYZE lv = toXYZ(itr.v4());
      lv *= (1+shift);
      LorentzVector newLv = toPtEtaPhi(lv);
      itr.set_pt(newLv.pt());
      itr.set_eta(newLv.eta());
      itr.set_phi(newLv.phi());
      itr.set_energy(newLv.energy());
    }
  }
}


QGAnalysisGenJetLambda::QGAnalysisGenJetLambda(uhh2::Context & ctx, float jetRadius, int nJetsMax, const GenParticleId & genId, const std::string & jet_coll_name, const std::string & output_coll_name):
  ca_wta_cluster_(JetDefinition(cambridge_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme)),
  jetRadius_(jetRadius),
  nJetsMax_(nJetsMax),
  genId_(genId),
  neutralHadronShift_(0), // these are the fractional shift, e.g. 0.2 for 20%
  photonShift_(0),
  genjet_handle(ctx.get_handle<std::vector<GenJetWithParts>>(jet_coll_name)),
  output_handle(ctx.get_handle<std::vector<GenJetLambdaBundle>>(output_coll_name))
{}


PseudoJet QGAnalysisGenJetLambda::convert_uhh_genparticle_to_pseudojet(const GenParticle & particle) {
  LorentzVectorXYZE lv = toXYZ(particle.v4());
  PseudoJet outputParticle(lv.px(), lv.py(), lv.pz(), lv.E());
  return outputParticle;
}


bool QGAnalysisGenJetLambda::process(uhh2::Event & event) {
  std::vector<GenJetWithParts> jets = event.get(genjet_handle);
  std::vector<GenJetLambdaBundle> outputs;
  int nJetCounter = 0;
  for (auto & jet : jets) {
    if (nJetsMax_ > 0 && nJetCounter == nJetsMax_) break;
    // Get constituents
    std::vector<GenParticle> daughters = get_jet_genparticles(jet, event);
    // Shift energies if appropriate
    // if (fabs(neutralHadronShift_) > 1E-6) { shift_neutral_hadron_genparticles(daughters, neutralHadronShift_); }
    // if (fabs(photonShift_) > 1E-6) { shift_photon_genparticles(daughters, photonShift_); }

    // FIXME: handle when 0 leftover constituents

    // Calculate the WTA axis
    // First convert Jet to PseudoJet, then convert it
    // PseudoJet origJet = convert_uhh_jet_to_pseudojet(jet);
    vector<PseudoJet> constits = {};
    int dauCounter = 0;
    for (const auto & dau : daughters) {
      PseudoJet thisDau = convert_uhh_genparticle_to_pseudojet(dau);
      thisDau.set_user_index(dauCounter);
      constits.push_back(thisDau);
      dauCounter++;
    }
    ClusterSequence cs(constits, ca_wta_cluster_);
    vector<PseudoJet> wtaJets = sorted_by_pt(cs.inclusive_jets());
    if (wtaJets.size() > 1) {
      cout << " >1 WTA jets" << endl;
    } else if (wtaJets.size() == 0) {
      cout << jet.genparticles_indices().size() << " genjet gp indices" << endl;
      cout << constits.size() << " gen pseudo daughters" << endl;
      std::vector<GenParticle> * genparticles = event.genparticles;
      for (const uint i : jet.genparticles_indices()) {
        auto gp = genparticles->at(i);
        cout << gp.pt() << " : " << gp.eta() << " : " << gp.phi() << " : " << gp.charge() << " : " << gp.pdgId() << " : " << gp.status() << endl;
        if (genId_ && !genId_(genparticles->at(i), event)) {
          cout << "FAIL GEN ID" << endl;
        }
      }
      throw std::runtime_error("WTA reclustering failed - no genjets");
    }
    LorentzVectorXYZE wtaJetAxis(wtaJets[0].px(), wtaJets[0].py(), wtaJets[0].pz(), wtaJets[0].E());

    // cout << "Original jet axis: " << jet.v4().px() << " : " << jet.v4().py() << " : " << jet.v4().pz() << endl;
    // cout << "Original jet axis: " << jet.eta() << " : " << jet.phi() << endl;
    // cout << "WTA jet axis: " << wtaJetAxis.px() << " : " << wtaJetAxis.py() << " : " << wtaJetAxis.pz() << endl;
    // cout << "WTA jet axis: " << wtaJetAxis.eta() << " : " << wtaJetAxis.phi() << endl;

    // Finally apply any pt charge cuts etc
    if (genId_) clean_collection<GenParticle>(daughters, event, genId_);

    LambdaCalculator<GenParticle> genJetCalc(daughters, jetRadius_, toPtEtaPhi(wtaJetAxis), false);
    GenJetLambdaBundle thisBundle{jet, genJetCalc};
    outputs.push_back(thisBundle);
    nJetCounter++;
  }
  event.set(output_handle, std::move(outputs));
  return true;
}

std::vector<GenParticle> QGAnalysisGenJetLambda::get_jet_genparticles(const GenJetWithParts & genjet, uhh2::Event & event) {
  std::vector<GenParticle> * genparticles = event.genparticles;
  std::vector<GenParticle> gp;
  for (const uint i : genjet.genparticles_indices()) {
    gp.push_back(genparticles->at(i)); // TODO store copy incase we shift it?
  }
  return gp;
}

// void QGAnalysisGenJetLambda::set_neutral_hadron_shift(float direction, float rel_shift) {
//   neutralHadronShift_ = (direction * rel_shift);
// }

// void QGAnalysisGenJetLambda::shift_neutral_hadron_genparticles(std::vector<GenParticle*> genparticles, float shift) {
//   for (auto & itr : genparticles) {
//     if (itr->pdgId() == PFParticle::eH0) {
//       itr->set_pt(itr->pt() * (1 + shift));
//     }
//   }
// }

// void QGAnalysisGenJetLambda::set_photon_shift(float direction, float rel_shift) {
//   photonShift_ = (direction * rel_shift);
// }

// void QGAnalysisGenJetLambda::shift_photon_genparticles(std::vector<GenParticle*> genparticles, float shift) {
//   for (auto & itr : genparticles) {
//     if (itr->pdgId() == PFParticle::eGamma) {
//       itr->set_pt(itr->pt() * (1 + shift));
//     }
//   }
// }


std::vector<double> Binning::calculate_fine_binning(const std::vector<double> & coarse_bin_edges)
{
  std::vector<double> fine_bin_edges;
  for (uint i=0; i<coarse_bin_edges.size()-1; i++) {
    fine_bin_edges.push_back(coarse_bin_edges[i]);
    fine_bin_edges.push_back(0.5*(coarse_bin_edges[i+1] + coarse_bin_edges[i]));
  }
  fine_bin_edges.push_back(coarse_bin_edges[coarse_bin_edges.size()-1]);
  return fine_bin_edges;
}

std::vector<double> Binning::sum_vectors(const std::vector<double> & vec1, const std::vector<double> & vec2)
{
  std::vector<double> result(vec1);
  result.insert(result.end(), vec2.begin(), vec2.end());
  return result;
}

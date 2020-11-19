#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"

#include "UHH2/core/include/Utils.h"

#include "TSystem.h"

using namespace std;
using namespace uhh2;
using namespace fastjet;
// using namespace fastjet::contrib;


MC::Name matchDatasetName(const std::string & name) {
    if (name.find("MGPYTHIA_QCD") != string::npos) {
        return MC::MGPYTHIA_QCD;
    } else if (name.find("PYTHIA-QCD-Pt") != string::npos) {
        return MC::PYTHIA_QCD_BINNED;
    } else if (name.find("HERWIG_QCD") != string::npos) {
        return MC::HERWIG_QCD;
    } else if (name.find("MGPYTHIA_DYJetsToLL") != string::npos) {
        return MC::MGPYTHIA_DY;
    } else if (name.find("HERWIG_DYJetsToLL") != string::npos) {
        return MC::HERWIG_DY;
    } else {
        throw std::runtime_error("Cannot understand MC dataset with name " + name);
    }
}

DataJetMetCorrector::DataJetMetCorrector(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone){
  if (jet_cone != "AK4" && jet_cone != "AK8") {
    throw runtime_error("jet_cone must be AK4 or AK8");
  }

  std::string puName = "chs";
  if (pu_removal == "PUPPI") {
    puName = "Puppi";
  } else if (pu_removal != "CHS") {
    throw runtime_error("pu_removal must be CHS or PUPPI");
  }

  std::string jec_jet_coll = jet_cone + "PF" + puName;

  jec_switcher_16.reset(new RunSwitcher(ctx, "2016"));
  for (const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
    jec_switcher_16->setupRun(runItr, std::make_shared<JetCorrector>(ctx, JERFiles::JECFilesDATA(Cuts::jec_tag_2016, Cuts::jec_ver_2016, jec_jet_coll, runItr)));
  }
}

bool DataJetMetCorrector::process(uhh2::Event & event) {
  jec_switcher_16->process(event);
  return true;
}


MCJetMetCorrector::MCJetMetCorrector(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, const std::string & jet_coll_name, const std::string & genjet_coll_name) {
  if (jet_cone != "AK4" && jet_cone != "AK8") {
    throw runtime_error("jet_cone must be AK4 or AK8");
  }

  std::string puName = "chs";
  if (pu_removal == "PUPPI") {
    puName = "Puppi";
  } else if (pu_removal != "CHS") {
    throw runtime_error("pu_removal must be CHS or PUPPI");
  }
  // TODO: move to YearSwitcher...?
  std::string jec_jet_coll = jet_cone + "PF" + puName;
  std::vector<std::string> JEC_MC = JERFiles::JECFilesMC(Cuts::jec_tag_2016, Cuts::jec_ver_2016, jec_jet_coll);
  jet_corrector.reset(new JetCorrector(ctx, JEC_MC));
  std::string resolutionFilename = "common/data/2016/Summer16_25nsV1_MC_PtResolution_" + jet_cone + "PF" + puName + ".txt";
  // cout << resolutionFilename << endl;
  jet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, jet_coll_name, genjet_coll_name, JERSmearing::SF_13TeV_Summer16_25nsV1, resolutionFilename));
}

bool MCJetMetCorrector::process(uhh2::Event & event) {
  jet_corrector->process(event);
  jet_corrector->correct_met(event);
  jet_resolution_smearer->process(event);
  return true;
}


GeneralEventSetup::GeneralEventSetup(uhh2::Context & ctx) {
  bool is_mc = ctx.get("dataset_type") == "MC";

  if (!is_mc) lumi_selection.reset(new LumiSelection(ctx));

  auto year = extract_year(ctx);

  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  // metfilters_selection->add<TriggerSelection>("globalTightHalo2016Filter", "Flag_globalTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  if (!is_mc) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
  metfilters_selection->add<TriggerSelection>("goodVertices", "Flag_goodVertices");
  metfilters_selection->add<EcalBadCalibSelection>("EcalBadCalibSelection"); // Use this instead of Flag_ecalBadCalibFilter, uses ecalBadCalibReducedMINIAODFilter in ntuple_generator
  if (year != Year::is2016v2) {
    metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter");
  } else {
    metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Extra_BadPFMuonFilter");
  }
  PrimaryVertexId pvid = StandardPrimaryVertexId();
  metfilters_selection->add<NPVSelection>("1 good PV", 1, -1, pvid);

  pv_cleaner.reset(new PrimaryVertexCleaner(pvid));

  electron_cleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronTagID(Electron::cutBasedElectronID_Summer16_80X_V1_medium), PtEtaCut(Cuts::reco_electron_pt_min, Cuts::electron_eta_max))));

  if (year == Year::is2016v2) {
    muon_cleaner.reset(new MuonCleaner(AndId<Muon>(MuonID(Muon::Medium), PtEtaCut(Cuts::reco_muon_pt_min, Cuts::muon_eta_max), MuonIso(0.2))));
  } else {
    muon_cleaner.reset(new MuonCleaner(AndId<Muon>(MuonID(Muon::CutBasedIdMedium), PtEtaCut(Cuts::reco_muon_pt_min, Cuts::muon_eta_max), MuonID(Muon::PFIsoMedium))));
  }
}

bool GeneralEventSetup::process(uhh2::Event & event) {

  if(event.isRealData && !lumi_selection->passes(event)) return false;

  pv_cleaner->process(event);

  electron_cleaner->process(event);

  muon_cleaner->process(event);

  // put this last, so objects are correctly cleaned, etc, for MCWeight afterwards
  if(!metfilters_selection->passes(event)) return false;

  return true;
}


RecoJetSetup::RecoJetSetup(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, float jet_radius, float jet_pt_min, float jet_y_max, bool doJetID, const std::string & muonHandleName) {
  bool is_mc = ctx.get("dataset_type") == "MC";

  if (doJetID) jet_pf_id.reset(new JetCleaner(ctx, JetPFID(Cuts::RECO_JET_ID)));

  if (is_mc) {
    jet_met_corrector.reset(new MCJetMetCorrector(ctx, pu_removal, jet_cone));
  } else {
    jet_met_corrector.reset(new DataJetMetCorrector(ctx, pu_removal, jet_cone));
  }

  // note we cut on y not eta, since our jets can be massive
  jet_cleaner.reset(new JetCleaner(ctx, PtYCut(jet_pt_min, jet_y_max)));

  // jet_ele_cleaner.reset(new JetElectronOverlapRemoval(ctx, jet_radius));

  if (string2bool(ctx.get("isZPlusJetsMuons")) && (muonHandleName != "")) {
    jet_mu_cleaner.reset(new JetMuonOverlapRemoval(ctx, jet_radius, muonHandleName));
  }
}

bool RecoJetSetup::process(uhh2::Event & event) {

  if (jet_pf_id) jet_pf_id->process(event);

  jet_met_corrector->process(event);

  jet_cleaner->process(event);

  // jet_ele_cleaner->process(event);

  if (jet_mu_cleaner) jet_mu_cleaner->process(event);

  sort_by_pt(*event.jets);

  return true;
}


ZFinder::ZFinder(uhh2::Context & ctx, const std::string & inputLabel, const std::string & zLeptonLabel):
  input_handle(ctx.get_handle<vector<Muon>>(inputLabel)),
  z_leptons_handle(ctx.get_handle<vector<Muon>>(zLeptonLabel))
{}

bool ZFinder::process(uhh2::Event & event) {
  event.set(z_leptons_handle, std::vector<Muon>()); // set handle to empty vector if fails later

  // Look for a reconstructed Z using inputs
  auto inputs = event.get(input_handle);
  if (inputs.size() < 2) return false;
  // Do we also want to consider more than leading & subleading?
  // Prob v.v.litle diff as not often > 2 leptons that pass selection
  auto zCand = inputs[0].v4() + inputs[1].v4();
  if ((fabs(zCand.M() - 90) < Cuts::mZ_window) && (inputs[0].charge() * inputs[1].charge() < 0)) {
    std::vector<Muon> cands = {inputs[0], inputs[1]};
    event.set(z_leptons_handle, cands);
    return true;
  }
  return false;
}


GenZFinder::GenZFinder(uhh2::Context & ctx, const std::string & inputLabel, const std::string & genZLabel, const std::string & genZLeptonLabel):
  input_handle(ctx.get_handle<std::vector<GenParticle>>(inputLabel)),
  z_handle(ctx.get_handle<GenParticle>(genZLabel)),
  z_leptons_handle(ctx.get_handle<std::vector<GenParticle>>(genZLeptonLabel))
{}


bool GenZFinder::process(uhh2::Event & event) {
  event.set(z_handle, GenParticle()); // set handle to empty if fails later
  event.set(z_leptons_handle, std::vector<GenParticle>()); // set handle to empty vector if fails later

  // Look for a Z using inputs
  auto inputs = event.get(input_handle);
  if (inputs.size() < 2) return false;
  // Do we also want to consider more than leading & subleading?
  // Prob v.v.litle diff as not often > 2 leptons that pass selection
  auto zCand = inputs[0].v4() + inputs[1].v4();
  if ((fabs(zCand.M() - 90) < Cuts::mZ_window) && (inputs[0].charge() * inputs[1].charge() < 0)) {
    std::vector<GenParticle> cands = {inputs[0], inputs[1]};
    GenParticle genZ;
    genZ.set_v4(zCand);
    genZ.set_pdgId(PDGID::Z); // in reality could be a photon
    event.set(z_handle, genZ);
    event.set(z_leptons_handle, cands);
    return true;
  }
  return false;
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


ZkFactorReweight::ZkFactorReweight(uhh2::Context & ctx, const std::string & weightFilename_, const std::string & genZLabel):
  z_weight_handle(ctx.get_handle<double>("z_weight"))
{
  if (weightFilename_ != "") {
    zReweight.reset(new ZllKFactor(weightFilename_));
  }
  if (genZLabel != "") {
    gen_z_handle = ctx.get_handle<GenParticle>(genZLabel);
  }
}

bool ZkFactorReweight::process(uhh2::Event & event) {
  if (zReweight) {
    // the k factor uses the Z as reconstructed by muons
    // TODO: dress them with photons?
    GenParticle genZ = event.get(gen_z_handle);
    double zPt = genZ.pt();
    double zWeight = 1;
    if (zPt > 0) zWeight = zReweight->getKFactor(zPt);
    event.set(z_weight_handle, zWeight);
    // cout << "z k factor: " << zWeight << endl;
    event.weight *= zWeight;
    return true;
  }
  return false;
}


PtReweight::PtReweight(uhh2::Context & ctx, const std::string & genjets_name, const std::string & weightFilename_, const std::string & region_):
genjets_handle(ctx.get_handle<std::vector<GenJet>>(genjets_name))
{
  if (weightFilename_ != "" && region_ != "") {
    TFile f_weight(weightFilename_.c_str());
    if (f_weight.IsZombie()) {
      throw std::runtime_error("Cannot open " + weightFilename_);
    }
    if (region_ == "dijet") {
      reweightHist.reset((TH1F*) f_weight.Get("dijet"));
      method = "dijet";
    } else if (region_ == "zplusjets") {
      reweightHist.reset((TH1F*) f_weight.Get("zplusjets"));
      method = "jet";
    } else {
      throw std::runtime_error("PtReweight: region_ argument not valid");
    }
    if (reweightHist) {
      reweightHist->SetDirectory(0);
    }
  }
}

bool PtReweight::process(uhh2::Event & event) {
  if (reweightHist) {
    double val = -1;
    std::vector<GenJet> genjets = event.get(genjets_handle);
    if (method == "jet" && genjets.size() >= 1) {
      val = genjets.at(0).pt();
    } else if (method == "dijet" && genjets.size() >=2 ) {
      val = 0.5*(genjets.at(0).pt() + genjets.at(1).pt());
    }
    if (val > reweightHist->GetXaxis()->GetXmax()) {
      // overflow protection
      val = reweightHist->GetXaxis()->GetXmax() - 0.1;
    }
    double weight = 1.;
    if (val > 0) {
      int bin_num = reweightHist->GetXaxis()->FindBin(val);
      weight = reweightHist->GetBinContent(bin_num);
    }
    event.weight *= weight;
    return true;
  }
  return false;
}

MCTrackScaleFactor::MCTrackScaleFactor(uhh2::Context & ctx, const std::string & direction):
  eta_regions{0.00, 0.80, 1.50, 2.5},
  SF{
    {1.01, 1.08, 0.93}, // Run B to F
    {1.04, 1.07, 1.12}  // Run G, H
  },
  SF_uncert{
    {0.03, 0.04, 0.04}, // Run B to F
    {0.03, 0.06, 0.05}  // Run G, H
  },
  run_BtoF_lumi(19284.990943025),  // taken from IsoMu24 trigger lumis
  run_GtoH_lumi(16633.228549923), // accounts for the HIP fix partway through Run F
  // run_BtoF_lumi(19691.766856822),  // taken from IsoMu24 trigger lumis
  // run_GtoH_lumi(16226.452636126),
  direction_(direction),
  drMax_(0.05),
  dropped_pf_handle(ctx.get_handle<std::vector<PFParticle>>("dropped_pfparticles")),
  promoted_pf_handle(ctx.get_handle<std::vector<PFParticle>>("promoted_genparticles"))
{
  total_lumi = run_BtoF_lumi + run_GtoH_lumi;
  random_.reset(new TRandom3(0));
  if (!(direction == "nominal" || direction == "up" || direction == "down")) {
    throw std::runtime_error("MCTrackScaleFactor direction should be one of: nominal, up, down");
  }
  int nbins_pt = 100;
  float pt_min = 0;
  float pt_max = 20;
  int nbins_dr = 20;
  for (int i=PFParticle::eX; i<=PFParticle::eH0; i++) {
    matching_pf_hists[i] = new TH3D(TString::Format("MCTrackScaleFactor_PF_%d", i), TString::Format("GenParticle-PF matching for PF ID %d;GenParticle p_{T} [GeV];PF particle p_{T} [GeV];#DeltaR", i), nbins_pt, pt_min, pt_max, nbins_pt, pt_min, pt_max, nbins_dr, 0, drMax_);;
    ctx.put("", matching_pf_hists[i]);
  }
}

PFParticle::EParticleID MCTrackScaleFactor::pdgid_to_pfid(int pdgid, int charge){
  int abspdgid = abs(pdgid);
  if (abspdgid == 11) {
    return PFParticle::eE;
  } else if (abspdgid == 13) {
    return PFParticle::eMu;
  } else if (abspdgid == 22) {
    return PFParticle::eGamma;
  } else if (abspdgid == 15) {
    return PFParticle::eH;
  } else if (charge != 0) {
    return PFParticle::eH;
  } else if (abspdgid != 12 && abspdgid != 14 && abspdgid != 16) {
    return PFParticle::eH0;
  } else {
    return PFParticle::eX;
  }
}

bool MCTrackScaleFactor::process(uhh2::Event & event) {
  // Decide randomly which run period we're in for this event
  // Done proportional to total lumi in each run period
  bool isRunGtoH = (random_->Rndm() < (run_GtoH_lumi / total_lumi));
  // Go through each eta region, and apply the SF
  // If SF < 1, then need to drop a track with a certain probability
  // If SF > 1, then need to promote a gen particle (that doesn't have a matching track)
  // to become a track, and somehow add it to a jet
  std::vector<PFParticle> dropped_pfparticles;
  std::vector<GenParticle> promoted_genparticles;
  for (uint eta_bin=0; eta_bin<eta_regions.size()-1; eta_bin++) {
    float eta_min = eta_regions[eta_bin];
    float eta_max = eta_regions[eta_bin+1];
    float sf = SF[isRunGtoH][eta_bin];
    if (direction_ == "up") {
      sf += SF_uncert[isRunGtoH][eta_bin];
    } else if (direction_ == "down") {
      sf -= SF_uncert[isRunGtoH][eta_bin];
    }

    if (sf < 1) {
      // Drop track with probability proportional to SF
      // FIXME: should this only apply to matched PF particles?
      for (uint pf_ind=0; pf_ind<event.pfparticles->size(); pf_ind++) {
        PFParticle pf = event.pfparticles->at(pf_ind);
        if (fabs(pf.eta()) > eta_max || fabs(pf.eta()) < eta_min) continue; // ignore if outside this eta bin
        if (pf.particleID() != PFParticle::eH) continue; //ignore if not charged hadron
        if (random_->Rndm() > sf) {
          // reject this particle
          dropped_pfparticles.push_back(pf);
        }
      }
    } else if (sf > 1) {
      // Promote charged genparticle to PF particle with probability proportional to SF,
      // but only if there's no match to a PF particle already
      int genCounter = 0;
      int matchedGenCounter = 0;
      std::vector<GenParticle> consideredGP; // needed since there are duplicate genparticles :((
      std::vector<GenParticle> unmatchedGP;
      for (uint gp_ind=0; gp_ind<event.genparticles->size(); gp_ind++) {
        GenParticle gp = event.genparticles->at(gp_ind);
        if (gp.status() != 1) continue; // final state only
        if (fabs(gp.eta()) > eta_max || fabs(gp.eta()) < eta_min) continue; // ignore if outside this eta bin
        if ((abs(gp.pdgId()) < 100) || (gp.charge() == 0)) continue; // only care about charged hadrons
        if (std::find(consideredGP.begin(), consideredGP.end(), gp) == consideredGP.end()) {
          // only count if non-duplicate
          consideredGP.push_back(gp);
        } else {
          continue;
        }

        genCounter++;
        // Look for matching pf particle
        bool matched = false;
        float minDR = 99999;
        int pfMatch = -1;
        for (uint pf_ind=0; pf_ind<event.pfparticles->size(); pf_ind++) {
          PFParticle pf = event.pfparticles->at(pf_ind);
          float dr = deltaR(pf.v4(), gp.v4());
          if (dr < drMax_ && dr < minDR) {
            minDR = dr;
            matched = true;
            pfMatch = pf_ind;
            // cout << "Found a gp-pf match: " << dr << endl;
            // cout << "GP: " << gp.pt() << " : " << gp.eta() << " : " << gp.phi() << " : " << gp.pdgId() << endl;
            // cout << "PF: " << pf.pt() << " : " << pf.eta() << " : " << pf.phi() << " : " << pf.particleID() << endl;
            // if (pf.particleID() == 5) { cout << " ************************************************************************************************" << endl;}
          }
        }

        if (matched) {
          matchedGenCounter++;
          matching_pf_hists[event.pfparticles->at(pfMatch).particleID()]->Fill(gp.pt(), event.pfparticles->at(pfMatch).pt(), minDR, event.weight);
          continue;
        }
        bool alreadyAdded = (std::find(unmatchedGP.begin(), unmatchedGP.end(), gp) != unmatchedGP.end());
        if (!alreadyAdded) { unmatchedGP.push_back(gp); }
        // bool alreadyAdded = (std::find(promoted_genparticles.begin(), promoted_genparticles.end(), gp) != promoted_genparticles.end());
        // if ((random_->Rndm() < (sf-1)) && !alreadyAdded) { promoted_genparticles.push_back(gp); }
      }
      // cout << "genCounter: " << genCounter << endl;
      // cout << "matchedGenCounter: " << matchedGenCounter << endl;
      // now decide how many to promote.
      double probability = TMath::Min(1., (sf-1.) * (consideredGP.size() - unmatchedGP.size())/ unmatchedGP.size());
      for (uint iu=0; iu<unmatchedGP.size(); iu++) {
        if (random_->Rndm() < probability) { promoted_genparticles.push_back(unmatchedGP.at(iu)); }
      }
    }
  }
  // convert those GPs to PFs
  std::vector<PFParticle> promoted_genparticles_as_pf;
  for (const auto & gpItr : promoted_genparticles) {
    PFParticle pf;
    pf.set_v4(gpItr.v4());
    pf.set_particleID(pdgid_to_pfid(gpItr.pdgId(), gpItr.charge()));
    pf.set_charge(gpItr.charge());
  }
  event.set(dropped_pf_handle, dropped_pfparticles);
  event.set(promoted_pf_handle, promoted_genparticles_as_pf);
  return true;
}


JetPFUpdater::JetPFUpdater(uhh2::Context & ctx, const std::string & jet_coll_name, bool update4vec):
  jet_handle(ctx.get_handle<std::vector<Jet>>(jet_coll_name)),
  dropped_pf_handle(ctx.get_handle<std::vector<PFParticle>>("dropped_pfparticles")),
  promoted_pf_handle(ctx.get_handle<std::vector<PFParticle>>("promoted_genparticles")),
  update4vec_(update4vec)
{
}

bool JetPFUpdater::process(uhh2::Event & event) {
  std::vector<PFParticle> & promoted_pf_particles = event.get(promoted_pf_handle);
  std::vector<PFParticle> & dropped_pf_particles = event.get(dropped_pf_handle);
  if (promoted_pf_particles.size() == 0 && dropped_pf_particles.size() == 0) { return true; }

  // add promoted pfparticles to event.pfparticles
  uint old_pf_size = event.pfparticles->size();
  event.pfparticles->insert(event.pfparticles->end(), promoted_pf_particles.begin(), promoted_pf_particles.end());

  // cout << "++++++++++++++++++++++++++++= Promoting: " << promoted_pf_particles.size() << endl;
  // cout << "----------------------------= Dropping: " << dropped_pf_particles.size() << endl;
  // go through jets, for each remove constituent if in dropped collection, add if in promoted

  for (auto & jet: event.get(jet_handle)) {
    std::vector<long int> newDauIndices;

    // get un-corrected 4-vector, will be updated with dropped/promoted particles
    LorentzVectorXYZE jetv4xyz = toXYZ(jet.v4()) * jet.JEC_factor_raw();
    LorentzVector jetv4 = toPtEtaPhi(jetv4xyz);
    // cout << "Old v4: " << jetv4.pt() << " : " << jetv4.eta() << " : " << jetv4.phi() << endl;

    // remove dropped
    if (dropped_pf_particles.size() > 0) {
      for (const auto dInd : jet.pfcand_indexs()) {
        PFParticle dau = event.pfparticles->at(dInd);
        if (std::find(dropped_pf_particles.begin(), dropped_pf_particles.end(), dau) == dropped_pf_particles.end()) {
          newDauIndices.push_back(dInd);
        }
        else {
          jetv4 -= dau.v4();
          // cout << "*************************** removed dau" << endl;
        }
      }
    } else {
      newDauIndices = jet.pfcand_indexs();
    }

    if (promoted_pf_particles.size() > 0) {
      // add promoted particles, if they fall within jet radius
      float jetRadius = TMath::Sqrt(jet.jetArea() / TMath::Pi()); // take from user instead?
      for (uint promInd=0; promInd < promoted_pf_particles.size(); promInd++) {
        const PFParticle & promoted = promoted_pf_particles.at(promInd);
        if (deltaRUsingY(promoted.v4(), jet.v4()) < jetRadius) {
          // cout << "*************************** Adding new dau ind " << old_pf_size+promInd << " to jet " << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << endl;
          newDauIndices.push_back(old_pf_size+promInd);
          jetv4 += promoted.v4();
        }
      }
    }
    // cout << "New v4: " << jetv4.pt() << " : " << jetv4.eta() << " : " << jetv4.phi() << endl;
    // TODO: recalc energy fractions as well?
    if (update4vec_) jet.set_v4(jetv4); // update with uncorrected 4-vector
    jet.set_pfcand_indexs(newDauIndices);
    jet.set_numberOfDaughters(newDauIndices.size());
  }

  return true;
}


TrackingEfficiency::TrackingEfficiency(uhh2::Context & ctx, bool update4vec) {
  track_sf.reset(new MCTrackScaleFactor(ctx, ctx.get("track_direction", "nominal")));
  jet_updater.reset(new JetPFUpdater(ctx, "jets", update4vec));
}

bool TrackingEfficiency::process(uhh2::Event & event) {
  if (event.isRealData) return true;
  track_sf->process(event);
  jet_updater->process(event);
  return true;
}


PDFReweight::PDFReweight(uhh2::Context & ctx) {
  if (ctx.get("pdf_reweight", "") == "") {
    skip_ = true;
    return;
  }

  if (gSystem->Load("libLHAPDF") == -1) {
    skip_ = true;
    std::cerr << "libLHAPDF not found, no pdf weights will be applied. To apply pdf re-weighting, add path to libLHAPDF.so to LD_LIBRARY_PATH" << std::endl;
    return;
  }
  skip_ = false;

  const std::string & pdfname = ctx.get("pdf_reweight");
  cout << "Using " << pdfname << " as new PDF set" << endl;
  pdf_ = LHAPDF::mkPDF(pdfname, 0);

  const std::string & original_pdfname = ctx.get("original_pdf", "NNPDF30_lo_as_0130");
  cout << "Using " << original_pdfname << " as original PDF" << endl;
  original_pdf_ = LHAPDF::mkPDF(original_pdfname, 0);
}

bool PDFReweight::process(uhh2::Event & event) {
  if (skip_) return true;

  double x1 = event.genInfo->pdf_x1();
  double x2 = event.genInfo->pdf_x2();

  int id1 = event.genInfo->pdf_id1();
  int id2 = event.genInfo->pdf_id2();

  double q = event.genInfo->pdf_scalePDF();

  // Get new PDF weights
  double xpdf1 = pdf_->xfxQ(id1, x1, q);
  double xpdf2 = pdf_->xfxQ(id2, x2, q);

  // Get original PDF weights
  // Have to recalculate ourselves, since the ones in the ntuple are 0 annoyingly
  double xpdf_orig1 = original_pdf_->xfxQ(id1, x1, q);
  double xpdf_orig2 = original_pdf_->xfxQ(id2, x2, q);
  // std::cout << "xpdf1 " << xpdf1 << " xpdf2 " << xpdf2 << std::endl;
  // std::cout << "xpdf_orig1 " << xpdf_orig1 << " xpdf_orig2 " << xpdf_orig2 << std::endl;

  event.weight *= (xpdf1 * xpdf2) / (xpdf_orig1 * xpdf_orig2);

  return true;
}


MCReweighting::MCReweighting(uhh2::Context & ctx, const std::string & genjet_name, const std::string & genZ_name) {
  gen_weight_handle = ctx.get_handle<double>("gen_weight");

  lumi_weighter.reset(new MCLumiWeight(ctx));

  pileup_reweighter.reset(new MCPileupReweight(ctx, ctx.get("pileup_direction", "central")));

  std::string datasetVersion = ctx.get("dataset_version");

  is_DY = isSubstring(datasetVersion, "DYJetsToLL", true);
  doMuons = string2bool(ctx.get("isZPlusJets")); // this is set in MCModule
  cout << "doMuons: " << doMuons << endl;

  if (doMuons) {
    std::string sf_path_name = locate_file("common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root");
    std::string sf_name = "MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta";
    muon_id_reweighter_pt_eta.reset(new MCMuonScaleFactor(ctx, sf_path_name, sf_name, 100));

    sf_name = "MC_NUM_MediumID2016_DEN_genTracks_PAR_vtx";
    // Doesn't work
    // muon_id_reweighter_vtx.reset(new MCMuonScaleFactor(ctx, sf_path_name, sf_name, 100));

    sf_path_name = locate_file("common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root");
    sf_name = "IsoMu24_OR_IsoTkMu24_PtEtaBins";
    muon_trg_reweighter.reset(new MCMuonScaleFactor(ctx, sf_path_name, sf_name, 100));

    std::string trk_path_name = locate_file("common/data/2016/general_eff_aeta_dr030e030_corr_ratio.txt");
    muon_trk_reweighter.reset(new MCMuonTrkScaleFactor(ctx, trk_path_name, 100));
  }

  if (is_DY) z_reweighter.reset(new ZkFactorReweight(ctx, ctx.get("z_reweight_file", ""), genZ_name));

  std::string pt_filename = ctx.get("pt_reweight_file", "");
  if (doMuons && ! is_DY) {
    // don't apply to non-DY Z+Jets MC (e.g. ttbar)
    pt_filename = "";
  }
  std::string region = is_DY ? "zplusjets" : "dijet";
  pt_reweighter.reset(new PtReweight(ctx, genjet_name, pt_filename, region));

  mc_scalevar.reset(new MCScaleVariation(ctx));

  pdf_reweighter.reset(new PDFReweight(ctx));
}


bool MCReweighting::process(uhh2::Event & event) {
  // store gen-only weights in separate variable
  // event.weight is then the product of reco weights & gen weights
  double old_gen_weight = event.get(gen_weight_handle);
  double old_weight = event.weight;
  if (!event.isRealData){
    lumi_weighter->process(event);
    if (is_DY) {
      z_reweighter->process(event);
    }
    pt_reweighter->process(event);
    mc_scalevar->process(event);
    pdf_reweighter->process(event);

    // ONLY DO THIS AFTER ALL GEN-SPECIFIC REWEIGHTING
    // Take the existing gen_weight, and update it with the gen bits we just did
    old_gen_weight *= (event.weight / old_weight);

    pileup_reweighter->process(event);

    if (doMuons) {
      muon_id_reweighter_pt_eta->process(event);

      // muon_id_reweighter_vtx->process(event);

      muon_trg_reweighter->process(event);

      muon_trk_reweighter->process(event);
    }
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



float calcGenHT(const std::vector<GenParticle> & genparticles) {
  // Find scalar sum of all matrix-element partons
  // Only works for Pythia8 samples due to status check
  // Returns -1 if no partons found
  float ht = 0.;
  int nFound = 0;
  for (const auto & itr: genparticles) {
    if (abs(itr.status()) != 23) continue;
    uint pdg = abs(itr.pdgId());
    if (( pdg <= PDGID::TOP_QUARK && pdg >= PDGID::DOWN_QUARK) || pdg == PDGID::GLUON) {
      nFound++;
      ht += itr.pt();
    }
  }
  return (nFound > 0) ? ht : -1;
}

float calcJetKt(const std::vector<GenParticle> & genparticles) {
  // find hardest parton with status 11
  float highestPt  = -1;
  for (const auto & gp : genparticles) {
      if (gp.status() != 11) continue; // can't assume all the 11s first - some status 4 in there as well! so don't break
      int absId = abs(gp.pdgId());
      if (absId < PDGID::TOP_QUARK || absId == PDGID::GLUON) {
        highestPt = max(highestPt, gp.pt());
      }
      // cout << "jet kt max: " << highestPt << " : " << gp.pt() << " : " << absId << " : " << gp.status() << endl;
  }
  return highestPt;
}


// PtReweight::PtReweight(uhh2::Context & ctx, const std::string & selection, const std::string & weightFilename):
//   gen_weight_handle(ctx.get_handle<double>("gen_weight"))
// {
//   f_weight.reset(TFile::Open(weightFilename.c_str()));
//   if (selection == "dijet")
//     reweightHist.reset((TH1F*) f_weight->Get("dijet_reco"));
//   else if (selection == "zplusjets")
//     reweightHist.reset((TH1F*) f_weight->Get("zpj_reco"));
//   if (reweightHist) reweightHist->SetDirectory(0);
// }


// bool PtReweight::process(uhh2::Event & event, float value) {
//   (void) event;
//   if (value >= reweightHist->GetXaxis()->GetXmax()) {
//     value = reweightHist->GetXaxis()->GetXmax() - 0.1;
//   }
//   int bin_num = reweightHist->GetXaxis()->FindBin(value);
//   double new_weight = reweightHist->GetBinContent(bin_num);

//   // Update the event weight & gen_weight stored in the event
//   double gen_weight = event.get(gen_weight_handle);
//   event.weight *= new_weight;
//   event.set(gen_weight_handle, gen_weight*new_weight);

//   return true;
// }

// namespace uhh2examples {


template<class T>
LambdaCalculator<T>::LambdaCalculator(std::vector<T> & constits, float jet_radius, const LorentzVector & jet_vector, const LorentzVector & wta_vector, bool usePuppiWeight):
  constits_(constits),
  jetRadius_(jet_radius),
  jetVector_(jet_vector),
  wtaVector_(wta_vector),
  usePuppiWeight_(usePuppiWeight)
{}

template<>
LambdaCalculator<GenParticle>::LambdaCalculator(std::vector<GenParticle> & constits, float jet_radius, const LorentzVector & jet_vector, const LorentzVector & wta_vector, bool usePuppiWeight):
  constits_(constits),
  jetRadius_(jet_radius),
  jetVector_(jet_vector),
  wtaVector_(wta_vector),
  usePuppiWeight_(usePuppiWeight) // Not used for GenParticles
{}

template<>
double LambdaCalculator<PFParticle>::getLambda(LambdaArgs lambdaArgs) const
{
  // special case if no constits
  if (constits_.size() == 0) {
    return -1.;
  }

  float kappa = lambdaArgs.kappa;
  float beta = lambdaArgs.beta;
  ParticleId constitId = lambdaArgs.id;

  double result = 0.;
  // Special case if both 0 ie multiplicity
  // Do it this way to ensure puppi weights correctly accounted for
  if (kappa == 0 && beta == 0) {
    for (const auto & dtr : constits_) {
      if (!constitId(dtr)) continue;
      if (usePuppiWeight_) {
        result += dtr.puppiWeight();
      } else {
        result += 1;
      }
    }
    return result;
  }

  double numerator(0.), ptSum(0.);
  for (const auto & dtr : constits_) {
    if (!constitId(dtr))
      continue;

    double weight = usePuppiWeight_ ? dtr.puppiWeight() : 1.;
    double thisPt = (kappa != 0) ? dtr.pt() * weight : 1.;
    // for the reference jet vector for deltaR, we use the WTA vector
    // if beta <=1, and the normal jet vector otherwise
    // Note: better compute (dist^2)^(beta/2) to avoid an extra square root
    // Taken from Gregory
    auto referenceAxis = (beta <= 1) ? wtaVector_ : jetVector_;
    double theta2 = (beta != 0) ? deltaR2UsingY(dtr.v4(), referenceAxis) : 1.;
    // 1 as puppi shouldn't change direction
    numerator += (pow(thisPt, kappa) * pow(theta2, 0.5*beta));
    ptSum += thisPt;
  }
  result = numerator / (pow(ptSum, kappa) * pow(jetRadius_, beta));
  return result;
}

template<>
double LambdaCalculator<GenParticle>::getLambda(LambdaArgs lambdaArgs) const
{
  // special case if no constits
  if (constits_.size() == 0) {
    return -1.;
  }

  float kappa = lambdaArgs.kappa;
  float beta = lambdaArgs.beta;
  ParticleId constitId = lambdaArgs.id;

  double result = 0.;
  // Special case if both 0 ie multiplicity
  // Do it this way to save time
  if (kappa == 0 && beta == 0) {
    for (const auto & dtr : constits_) {
      if (!constitId(dtr)) continue;
      result++;
    }
    return result;
  }

  // If not, calculate it and store in cache
  double numerator(0.), ptSum(0.);
  for (const auto & dtr : constits_) {
    if (!constitId(dtr)) continue;
    double thisPt = dtr.pt();
    // for the reference jet vector for deltaR, we use the WTA vector if beta <=1,
    // and the normal jet vector otherwise
    // Note: better compute (dist^2)^(beta/2) to avoid an extra square root
    // Taken from Gregory
    auto referenceAxis = (beta <= 1) ? wtaVector_ : jetVector_;
    double theta2 = (beta != 0) ? deltaR2UsingY(dtr.v4(), referenceAxis) : 1.;
    numerator += (pow(thisPt, kappa) * pow(theta2, 0.5*beta));
    ptSum += thisPt;
  }
  result = numerator / (pow(ptSum, kappa) * pow(jetRadius_, beta));
  return result;
}


QGAnalysisJetLambda::QGAnalysisJetLambda(uhh2::Context & ctx,
                                         float jetRadius,
                                         int nJetsMax,
                                         bool doPuppi,
                                         const PFParticleId & pfId,
                                         const std::string & jet_coll_name,
                                         const std::string & output_coll_name):
  wta_cluster_(Recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme), false, Recluster::keep_only_hardest)),
  // note that if MMDT/SD does reclustering, then it will call
  // Recluster(cambridge_algorithm, JetDefinition::max_allowable_R)(jet);
  // in RecursiveSymmetryCutBase::_recluster_if_needed. And inside Recluster ctor,
  // it will do:
  // JetDefinition(new_jet_alg, new_jet_radius)
  // which by default has E-scheme recombination
  jetRadius_(jetRadius),
  nJetsMax_(nJetsMax),
  doPuppi_(doPuppi),
  pfId_(pfId),
  chargedHadronShift_(0), // these are the fractional shift, e.g. 0.2 for 20%
  neutralHadronShift_(0), // these are the fractional shift, e.g. 0.2 for 20%
  photonShift_(0), // must init to 0 otherwise bad things will happen, will go haywire!
  jet_handle_(ctx.get_handle<std::vector<Jet>>(jet_coll_name)),
  output_handle_(ctx.get_handle<std::vector<JetLambdaBundle>>(output_coll_name))
{}


PseudoJet QGAnalysisJetLambda::convert_uhh_pfparticle_to_pseudojet(const PFParticle & particle, bool applyPuppiWeight) {
  float weight = (applyPuppiWeight) ? particle.puppiWeight() : 1.;
  LorentzVectorXYZE lv = toXYZ(particle.v4()) * weight;
  PseudoJet outputParticle(lv.px(), lv.py(), lv.pz(), lv.E());
  return outputParticle;
}

bool QGAnalysisJetLambda::process(uhh2::Event & event) {
  std::vector<Jet> jets = event.get(jet_handle_);
  std::vector<JetLambdaBundle> outputs;
  int nJetCounter = 0;
  for (auto & jet : jets) {
    if (nJetsMax_ > 0 && nJetCounter == nJetsMax_) break;
    // For each jet, we:
    // - Get constits, apply energy shifts if need be
    // - Calculate the WTA axis
    // - Get charged constituents, recluster into jet, get WTA axis
    // - Apply grooming to create subset of groomed constits, get WTA axis
    // - Apply grooming to charged-only jet, get WTA axis
    // - Apply any other cuts (e.g. pt)
    // - Create LambdaCalculator for each constit collection, save to struct
    // {ungroomed, groomed} x {neutral+charged, charged-only}

    // Get constituents
    // don't apply puppi weight here, since we already account for it in the LambdaCalculator
    std::vector<PFParticle> constits = get_jet_pfparticles(jet, event, false);
    clean_collection<PFParticle>(constits, event, PtEtaCut(1E-8, 10.)); // basic cut to remove weird 0 pt constits

    LorentzVector sumv4 = jet_constit_4vec(jet, *event.pfparticles, doPuppi_);
    if (fabs(sumv4.pt() - jet.pt()*jet.JEC_factor_raw())/(jet.pt()*jet.JEC_factor_raw()) >  0.05) {
      cout << "Warning sumv4 - jet pt mismatch" << endl;
      cout << "jet: " << jet.pt()*jet.JEC_factor_raw() << " : " << jet.Rapidity() << " : " << jet.phi() << endl;
      cout << "sumv4: " << sumv4.pt() << " : " << sumv4.Rapidity() << " : " << sumv4.phi() << endl;
    }

    // Shift energies if appropriate
    if (fabs(chargedHadronShift_) > 1E-6) { shift_charged_hadron_pfparticles(constits, chargedHadronShift_); }
    if (fabs(neutralHadronShift_) > 1E-6) { shift_neutral_hadron_pfparticles(constits, neutralHadronShift_); }
    if (fabs(photonShift_) > 1E-6) { shift_photon_pfparticles(constits, photonShift_); }

    // Calculate the WTA axis
    // First convert constits to PseudoJets
    vector<PseudoJet> pjconstits = {};
    int dauCounter = 0;
    for (const auto & dau : constits) {
      PseudoJet thisDau = convert_uhh_pfparticle_to_pseudojet(dau, doPuppi_); // do apply puppi weight here to recluster
      thisDau.set_user_index(dauCounter);
      pjconstits.push_back(thisDau);
      dauCounter++;
    }

    // Cluster using AKx to make pseudojet with right history
    // don't use jet radius, as it can lead to several jets instead of the original
    // Uses E-scheme by default
    JetDefinition jet_def(antikt_algorithm, JetDefinition::max_allowable_R);
    fastjet::ClusterSequence clust_seq(pjconstits, jet_def);
    std::vector<PseudoJet> akJets = sorted_by_pt(clust_seq.inclusive_jets(1.));

    if (akJets.size() > 1) {
      cout << "QGAnalysisJetLambda : >1 ak jets" << endl;
      cout << "Original jet: " << jet.pt()*jet.JEC_factor_raw() << " : " << jet.Rapidity() << " : " << jet.phi() << endl;
      for (const auto & newjet : akJets) {
        cout << "ak jet: " << newjet.pt() << " : " << newjet.rap() << " : " << newjet.phi() << endl;
        for (const auto & d : newjet.constituents()) {
          cout << "   constit: " << d.pt() << " : " << d.rap() << " : " << d.phi() << endl;
        }
      }
    } else if (akJets.size() == 0) {
      cout << "uhh jet constits:" << endl;
      for (const auto & dau : constits) {
        cout << "    " << dau.pt() << " : " << dau.Rapidity() << " :" << dau.phi() << endl;
      }
      cout << "pjconstits:" << endl;
      for (const auto & dau : pjconstits) {
        cout << "    " << dau.pt() << " : " << dau.rap() << " :" << dau.phi() << endl;
      }
      cout << "Event number : run : lumi: " << event.event << " : " << event.run << " : " << event.luminosityBlock << endl;
      throw std::runtime_error("AKx reclustering failed QGAnalysisJetLambda - no jets");
    }

    // Check AK jet is same as original
    if ((fabs(akJets[0].pt() - sumv4.pt()) / akJets[0].pt() > 0.1) || (fabs(akJets[0].rap() - sumv4.Rapidity()) > 0.1)){
      cout << "Original jet: " << sumv4.pt() << " : " << sumv4.Rapidity() << " : " << sumv4.phi() << endl;
      for (const auto & dau : constits) {
        cout << "   dau: " << dau.pt() << " : " << dau.Rapidity() << " : " << dau.phi() << endl;
      }
      cout << "ak jet: " << akJets[0].pt() << " : " << akJets[0].rap() << " : " << akJets[0].phi() << endl;
      for (const auto & d : akJets[0].constituents()) {
        cout << "   constit: " << d.pt() << " : " << d.rap() << " : " << d.phi() << endl;
      }
      cout << "Event number : run : lumi: " << event.event << " : " << event.run << " : " << event.luminosityBlock << endl;
      throw std::runtime_error("Reclustered AK pT/y mismatch");
    }

    // Now recluster using WTA
    PseudoJet wtaJet = wta_cluster_(akJets[0]);

    // Get charged constits, recluster into jet, calculate WTA axis
    vector<PseudoJet> chargedConstitsPJ;
    for (auto pj : pjconstits){
      // have to get charge from original PFParticle, not stored in PseudoJet
      if (constits.at(pj.user_index()).charge() != 0) {
        chargedConstitsPJ.push_back(pj);
      }
    }
    JetDefinition jet_def_R(antikt_algorithm, jetRadius_);
    vector<PseudoJet> chargedJets = jet_def_R(chargedConstitsPJ);
    PseudoJet chargedJet;
    PseudoJet wtaJetCharged;
    std::vector<PFParticle> chargedConstits;
    if (chargedJets.size() > 0) {
      chargedJet = chargedJets.at(0);
      wtaJetCharged = wta_cluster_(chargedJet);
      // Store the corresponding PFParticles
      for (const auto & dItr : chargedJet.constituents()) {
        chargedConstits.push_back(constits.at(dItr.user_index()));
      }
    }

    // Apply grooming to jet (does CA reclustering internally with E-scheme)
    // Also get WTA axis, and constits left after grooming
    fastjet::contrib::SoftDrop sd(0, 0.1, jetRadius_);
    PseudoJet groomedJet = sd(akJets[0]);
    PseudoJet groomedJetWTA = wta_cluster_(groomedJet); // want the AK WTA axis for lambda variables
    // create collection with only those in groomed jet
    std::vector<PFParticle> groomedConstits;
    for (const auto & dItr : groomedJet.constituents()) {
      groomedConstits.push_back(constits.at(dItr.user_index()));
    }

    // Do charged-only verions
    PseudoJet groomedChargedJet;
    PseudoJet groomedChargedJetWTA;
    std::vector<PFParticle> groomedChargedConstits;
    if (chargedJets.size() > 0) {
      groomedChargedJet = sd(chargedJets[0]);
      groomedChargedJetWTA = wta_cluster_(groomedChargedJet);
      for (const auto & dItr : groomedChargedJet.constituents()) {
        groomedChargedConstits.push_back(constits.at(dItr.user_index()));
      }
    }

    // Finally apply any pt cuts etc
    if (pfId_) {
      clean_collection<PFParticle>(constits, event, pfId_);
      clean_collection<PFParticle>(chargedConstits, event, pfId_);
      clean_collection<PFParticle>(groomedConstits, event, pfId_);
      clean_collection<PFParticle>(groomedChargedConstits, event, pfId_);
    }

    // Note that constits/chargedConstits/groomedConstits/groomedChargedConstits are PFParticles WITHOUT Puppi weights applied to 4-vector
    auto wtaJetAxis = toPtEtaPhi(LorentzVectorXYZE(wtaJet.px(), wtaJet.py(), wtaJet.pz(), wtaJet.E()));
    LambdaCalculator<PFParticle> recoJetCalc(constits, jetRadius_, jet.v4(), wtaJetAxis, doPuppi_);

    auto jetAxisCharged = toPtEtaPhi(LorentzVectorXYZE(chargedJet.px(), chargedJet.py(), chargedJet.pz(), chargedJet.E()));
    auto wtaJetAxisCharged = toPtEtaPhi(LorentzVectorXYZE(wtaJetCharged.px(), wtaJetCharged.py(), wtaJetCharged.pz(), wtaJetCharged.E()));
    LambdaCalculator<PFParticle> recoJetCalcCharged(chargedConstits, jetRadius_, jetAxisCharged, wtaJetAxisCharged, doPuppi_);

    auto jetAxisGroomed = toPtEtaPhi(LorentzVectorXYZE(groomedJet.px(), groomedJet.py(), groomedJet.pz(), groomedJet.E()));
    auto wtaJetAxisGroomed = toPtEtaPhi(LorentzVectorXYZE(groomedJetWTA.px(), groomedJetWTA.py(), groomedJetWTA.pz(), groomedJetWTA.E()));
    LambdaCalculator<PFParticle> recoJetCalcGroomed(groomedConstits, jetRadius_, jetAxisGroomed, wtaJetAxisGroomed, doPuppi_);

    auto jetAxisGroomedCharged = toPtEtaPhi(LorentzVectorXYZE(groomedChargedJet.px(), groomedChargedJet.py(), groomedChargedJet.pz(), groomedChargedJet.E()));
    auto wtaJetAxisGroomedCharged = toPtEtaPhi(LorentzVectorXYZE(groomedChargedJetWTA.px(), groomedChargedJetWTA.py(), groomedChargedJetWTA.pz(), groomedChargedJetWTA.E()));
    LambdaCalculator<PFParticle> recoJetCalcGroomedCharged(groomedChargedConstits, jetRadius_, jetAxisGroomedCharged, wtaJetAxisGroomedCharged, doPuppi_);

    JetLambdaBundle thisBundle{jet, recoJetCalc, recoJetCalcCharged, recoJetCalcGroomed, recoJetCalcGroomedCharged};
    outputs.push_back(thisBundle);
    nJetCounter++;
  }
  event.set(output_handle_, std::move(outputs));
  return true;
}

std::vector<PFParticle> QGAnalysisJetLambda::get_jet_pfparticles(const Jet & jet, uhh2::Event & event, bool applyPuppiWeight) {
  std::vector<PFParticle> * pfparticles = event.pfparticles;
  std::vector<PFParticle> pfCopy;
  // Create a copy, since we might modify it later and we want to keep the original event.pfparticles intact
  for (const uint i : jet.pfcand_indexs()) {
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

void QGAnalysisJetLambda::set_charged_hadron_shift(int direction, float rel_shift) {
  chargedHadronShift_ = (direction * rel_shift);
}

void QGAnalysisJetLambda::shift_charged_hadron_pfparticles(std::vector<PFParticle> & pfparticles, float shift) {
  for (auto & itr : pfparticles) {
    if (itr.particleID() == PFParticle::eH) {
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


JetLambdaCopier::JetLambdaCopier(uhh2::Context & ctx,
                                 const std::string & jet_coll_name,
                                 const std::string & lambda_coll_name,
                                 const std::string & output_coll_name):
  jet_handle_(ctx.get_handle<std::vector<Jet>>(jet_coll_name)),
  lambda_handle_(ctx.get_handle<std::vector<JetLambdaBundle>>(lambda_coll_name)),
  output_handle_(ctx.get_handle<std::vector<JetLambdaBundle>>(output_coll_name))
{}

bool JetLambdaCopier::process(uhh2::Event & event) {
  std::vector<JetLambdaBundle> output = {};
  auto & bundles = event.get(lambda_handle_);
  auto & jets = event.get(jet_handle_);
  bool foundMatch = false;
  for (auto & bundle : bundles) {
    for (auto & jet : jets) {
      // check if any of the jets match the jet in this bundle
      if (deltaRUsingY(jet, bundle.jet) < 0.01) {
        JetLambdaBundle newBundle = bundle;
        newBundle.jet = jet;
        output.push_back(newBundle);
        foundMatch = true;
      }
    }
  }
  event.set(output_handle_, output);
  return foundMatch;
}


QGAnalysisGenJetLambda::QGAnalysisGenJetLambda(uhh2::Context & ctx,
                                               float jetRadius,
                                               int nJetsMax,
                                               const GenParticleId & genId,
                                               const std::string & jet_coll_name,
                                               const std::string & output_coll_name):
  wta_cluster_(Recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme), false, Recluster::keep_only_hardest)),
  jetRadius_(jetRadius),
  nJetsMax_(nJetsMax),
  genId_(genId),
  neutralHadronShift_(0), // these are the fractional shift, e.g. 0.2 for 20%
  photonShift_(0),
  genjet_handle_(ctx.get_handle<std::vector<GenJet>>(jet_coll_name)),
  output_handle_(ctx.get_handle<std::vector<GenJetLambdaBundle>>(output_coll_name))
{}


PseudoJet QGAnalysisGenJetLambda::convert_uhh_genparticle_to_pseudojet(const GenParticle & particle) {
  LorentzVectorXYZE lv = toXYZ(particle.v4());
  PseudoJet outputParticle(lv.px(), lv.py(), lv.pz(), lv.E());
  return outputParticle;
}


bool QGAnalysisGenJetLambda::process(uhh2::Event & event) {
  std::vector<GenJet> jets = event.get(genjet_handle_);
  std::vector<GenJetLambdaBundle> outputs = {};
  int nJetCounter = 0;
  // cout << "---- Doing QGAnalysisGenJetLambda::process ----" << endl;
  for (auto & jet : jets) {
    if (nJetsMax_ > 0 && nJetCounter == nJetsMax_) break;
    // For each jet, we:
    // - Get constits, apply energy shifts if need be
    // - Calculate the WTA axis
    // - Get charged constituents, recluster into jet, get WTA axis
    // - Apply grooming to create subset of groomed constits, get WTA axis
    // - Apply grooming to charged-only jet, get WTA axis
    // - Apply any other cuts (e.g. pt)
    // - Create LambdaCalculator for each constit collection, save to struct
    // {ungroomed, groomed} x {neutral+charged, charged-only}

    // Get constituents
    std::vector<GenParticle> constits = get_jet_genparticles(jet, event);
    clean_collection<GenParticle>(constits, event, PtEtaCut(1E-8, 10.)); // basic cut to remove weird 0 pt constits

    LorentzVector sumv4 = genjet_constit_4vec(jet, *event.genparticles);
    if (fabs(sumv4.pt() - jet.pt())/jet.pt() >  0.05) {
      cout << "Warning sumv4 - genjet pt mismatch" << endl;
      cout << "jet: " << jet.pt() << " : " << jet.Rapidity() << " : " << jet.phi() << endl;
      cout << "sumv4: " << sumv4.pt() << " : " << sumv4.Rapidity() << " : " << sumv4.phi() << endl;
    }

    // Calculate the WTA axis
    // First convert constits to PseudoJets
    vector<PseudoJet> pjconstits = {};
    int dauCounter = 0;
    for (const auto & dau : constits) {
      PseudoJet thisDau = convert_uhh_genparticle_to_pseudojet(dau);
      thisDau.set_user_index(dauCounter);
      pjconstits.push_back(thisDau);
      dauCounter++;
    }

    // Cluster using AKx to make pseudojet with right history
    JetDefinition jet_def(antikt_algorithm, JetDefinition::max_allowable_R); // don't use jet radius, as it can lead to several jets instead of the original
    fastjet::ClusterSequence clust_seq(pjconstits, jet_def);
    std::vector<PseudoJet> akJets = sorted_by_pt(clust_seq.inclusive_jets(1.));

    if (akJets.size() > 1) {
      cout << "QGAnalysisGenJetLambda : >1 ak jets" << endl;
      cout << "Original jet: " << jet.pt() << " : " << jet.Rapidity() << " : " << jet.phi() << endl;
      for (const auto & newjet : akJets) {
        cout << "ak jet: " << newjet.pt() << " : " << newjet.rap() << " : " << newjet.phi() << endl;
        for (const auto & d : newjet.constituents()) {
          cout << "   constit: " << d.pt() << " : " << d.rap() << " : " << d.phi() << endl;
        }
      }
    } else if (akJets.size() == 0) {
      cout << "uhh jet constits:" << endl;
      for (const auto & dau : constits) {
        cout << "    " << dau.pt() << " : " << dau.Rapidity() << " : " << dau.phi() << endl;
      }
      cout << "pjconstits:" << endl;
      for (const auto & dau : pjconstits) {
        cout << "    " << dau.pt() << " : " << dau.rap() << " : " << dau.phi() << endl;
      }
      cout << "Event number : run : lumi: " << event.event << " : " << event.run << " : " << event.luminosityBlock << endl;
      throw std::runtime_error("AKx reclustering failed QGAnalysisGenJetLambda - no jets");
    }

    // Check AK jet is same as original
    // Use summed 4vec as sometimes differs to original jet pt, better comparison
    if ((fabs(akJets[0].pt() - sumv4.pt()) / sumv4.pt() > 0.1) || (fabs(akJets[0].rap() - sumv4.Rapidity()) > 0.1))  {
      cout << "Original genjet: " << sumv4.pt() << " : " << sumv4.Rapidity() << " : " << sumv4.phi() << endl;
      for (const auto & dau : constits) {
        cout << "   dau: " << dau.pt() << " : " << dau.Rapidity() << " : " << dau.phi() << endl;
      }
      cout << "ak jet: " << akJets[0].pt() << " : " << akJets[0].rap() << " : " << akJets[0].phi() << endl;
      for (const auto & d : akJets[0].constituents()) {
        cout << "   constit: " << d.pt() << " : " << d.rap() << " : " << d.phi() << endl;
      }
      cout << "Event number : run : lumi: " << event.event << " : " << event.run << " : " << event.luminosityBlock << endl;
      throw std::runtime_error("Reclustered gen AK pT/y mismatch");
    }

    // Now recluster using WTA
    PseudoJet wtaJet = wta_cluster_(akJets[0]);

    // Get charged constits, recluster into jet, calculate WTA axis
    vector<PseudoJet> chargedConstitsPJ;
    for (auto pj : pjconstits){
      // have to get charge from original GenParticle, not stored in PseudoJet
      if (constits.at(pj.user_index()).charge() != 0) {
        chargedConstitsPJ.push_back(pj);
      }
    }
    JetDefinition jet_def_R(antikt_algorithm, jetRadius_);
    vector<PseudoJet> chargedJets = jet_def_R(chargedConstitsPJ);
    PseudoJet chargedJet;
    PseudoJet wtaJetCharged;
    std::vector<GenParticle> chargedConstits;
    if (chargedJets.size() > 0) {
      chargedJet = chargedJets.at(0);
      wtaJetCharged = wta_cluster_(chargedJet);
      // Store the corresponding GenParticles
      for (const auto & dItr : chargedJet.constituents()) {
        chargedConstits.push_back(constits.at(dItr.user_index()));
      }
    }

    // Apply grooming to jet (does CA reclustering internally with E-scheme)
    // Also get WTA axis, and constits left after grooming
    fastjet::contrib::SoftDrop sd(0, 0.1, jetRadius_);
    PseudoJet groomedJet = sd(akJets[0]);
    PseudoJet groomedJetWTA = wta_cluster_(groomedJet); // want the AK WTA axis for lambda variables
    // create collection with only those in groomed jet
    std::vector<GenParticle> groomedConstits;
    for (const auto & dItr : groomedJet.constituents()) {
      groomedConstits.push_back(constits.at(dItr.user_index()));
    }

    // Do charged-only verions
    PseudoJet groomedChargedJet;
    PseudoJet groomedChargedJetWTA;
    std::vector<GenParticle> groomedChargedConstits;
    if (chargedJets.size() > 0) {
      groomedChargedJet = sd(chargedJets[0]);
      groomedChargedJetWTA = wta_cluster_(groomedChargedJet);
      for (const auto & dItr : groomedChargedJet.constituents()) {
        groomedChargedConstits.push_back(constits.at(dItr.user_index()));
      }
    }

    // Finally apply any pt cuts etc
    if (genId_) {
      clean_collection<GenParticle>(constits, event, genId_);
      clean_collection<GenParticle>(chargedConstits, event, genId_);
      clean_collection<GenParticle>(groomedConstits, event, genId_);
      clean_collection<GenParticle>(groomedChargedConstits, event, genId_);
    }

    auto wtaJetAxis = toPtEtaPhi(LorentzVectorXYZE(wtaJet.px(), wtaJet.py(), wtaJet.pz(), wtaJet.E()));
    LambdaCalculator<GenParticle> genJetCalc(constits, jetRadius_, jet.v4(), wtaJetAxis, false);

    auto jetAxisCharged = toPtEtaPhi(LorentzVectorXYZE(chargedJet.px(), chargedJet.py(), chargedJet.pz(), chargedJet.E()));
    auto wtaJetAxisCharged = toPtEtaPhi(LorentzVectorXYZE(wtaJetCharged.px(), wtaJetCharged.py(), wtaJetCharged.pz(), wtaJetCharged.E()));
    LambdaCalculator<GenParticle> genJetCalcCharged(chargedConstits, jetRadius_, jetAxisCharged, wtaJetAxisCharged, false);

    auto jetAxisGroomed = toPtEtaPhi(LorentzVectorXYZE(groomedJet.px(), groomedJet.py(), groomedJet.pz(), groomedJet.E()));
    auto wtaJetAxisGroomed = toPtEtaPhi(LorentzVectorXYZE(groomedJetWTA.px(), groomedJetWTA.py(), groomedJetWTA.pz(), groomedJetWTA.E()));
    LambdaCalculator<GenParticle> genJetCalcGroomed(groomedConstits, jetRadius_, jetAxisGroomed, wtaJetAxisGroomed, false);

    auto jetAxisGroomedCharged = toPtEtaPhi(LorentzVectorXYZE(groomedChargedJet.px(), groomedChargedJet.py(), groomedChargedJet.pz(), groomedChargedJet.E()));
    auto wtaJetAxisGroomedCharged = toPtEtaPhi(LorentzVectorXYZE(groomedChargedJetWTA.px(), groomedChargedJetWTA.py(), groomedChargedJetWTA.pz(), groomedChargedJetWTA.E()));
    LambdaCalculator<GenParticle> genJetCalcGroomedCharged(groomedChargedConstits, jetRadius_, jetAxisGroomedCharged, wtaJetAxisGroomedCharged, false);

    GenJetLambdaBundle thisBundle{jet, genJetCalc, genJetCalcCharged, genJetCalcGroomed, genJetCalcGroomedCharged};
    outputs.push_back(thisBundle);
    nJetCounter++;
  }

  // cout << "---- End QGAnalysisGenJetLambda::process ----" << endl;
  event.set(output_handle_, std::move(outputs));
  return true;
}

std::vector<GenParticle> QGAnalysisGenJetLambda::get_jet_genparticles(const GenJet & genjet, uhh2::Event & event) {
  std::vector<GenParticle> * genparticles = event.genparticles;
  std::vector<GenParticle> gps;
  bool firstTime = false;
  for (const uint i : genjet.genparticles_indices()) {
    auto gp = genparticles->at(i);
    if (gp.pt() <  1E-7) {
      if (!firstTime) {
        cout << "Low pT GP in this GenJet: " << genjet.pt() << " : " << genjet.Rapidity() << " : " << genjet.phi() << " : " << genjet.genparticles_indices().size() << " : " << genjet.partonFlavour() << endl;
        firstTime = true;
      }
      cout << "Low pT GP: " << gp.pt() << " : " << gp.Rapidity() << " : " << gp.phi() << " : " << gp.status() << " : " << gp.pdgId() << endl;
    }
    gps.push_back(genparticles->at(i)); // TODO store copy incase we shift it?
  }
  return gps;
}


GenJetLambdaCopier::GenJetLambdaCopier(uhh2::Context & ctx,
                                       const std::string & jet_coll_name,
                                       const std::string & lambda_coll_name,
                                       const std::string & output_coll_name):
  genjet_handle_(ctx.get_handle<std::vector<GenJet>>(jet_coll_name)),
  lambda_handle_(ctx.get_handle<std::vector<GenJetLambdaBundle>>(lambda_coll_name)),
  output_handle_(ctx.get_handle<std::vector<GenJetLambdaBundle>>(output_coll_name))
{}

bool GenJetLambdaCopier::process(uhh2::Event & event) {
  std::vector<GenJetLambdaBundle> output = {};
  auto & bundles = event.get(lambda_handle_);
  auto & jets = event.get(genjet_handle_);
  bool foundMatch = false;
  for (auto & bundle : bundles) {
    for (auto & jet : jets) {
      // check if any of the jets match the jet in this bundle
      if (deltaRUsingY(jet, bundle.jet) < 0.01) {
        GenJetLambdaBundle newBundle = bundle;
        newBundle.jet = jet;
        output.push_back(newBundle);
        foundMatch = true;
      }
    }
  }
  event.set(output_handle_, output);
  return foundMatch;
}


JetMatcher::JetMatcher(uhh2::Context & ctx, const std::string & jet_coll_name, const std::string & genjet_coll_name, float matchRadius, bool uniqueMatch):
  recojet_handle_(ctx.get_handle<std::vector<Jet>>(jet_coll_name)),
  genjet_handle_(ctx.get_handle<std::vector<GenJet>>(genjet_coll_name)),
  matchRadius_(matchRadius),
  uniqueMatch_(uniqueMatch)
{}

bool JetMatcher::process(uhh2::Event & event) {
  std::vector<uint> matchedIndices;
  std::vector<Jet> & jets = event.get(recojet_handle_);
  std::vector<GenJet> & genjets = event.get(genjet_handle_);
  // For each reco jet, loop through genjets and find closest that is within matchRadius_
  // If uniqueMatch is true, then each genjet can only be matched to at most 1 reco jet
  for (auto & jtr: jets) {
    float minDR = 9999.;
    int matchInd = -1; // sensible default - not 0!

    for (uint gjInd=0; gjInd < genjets.size(); gjInd++) {
      // If we want unique matches and we've already matched then skip this genjet
      if (uniqueMatch_ && std::find(matchedIndices.begin(), matchedIndices.end(), gjInd) != matchedIndices.end())
        continue;

      const auto genjtr = genjets.at(gjInd);
      auto thisDR = deltaRUsingY(jtr, genjtr);
      if (thisDR < matchRadius_ && thisDR < minDR) {
        matchInd = gjInd;
        minDR = thisDR;
      }
    }

    jtr.set_genjet_index(matchInd);
  }
  return true;
}


bool isParton(int pdgId) {
  uint pdg = abs(pdgId);
  return ((pdg >= PDGID::DOWN_QUARK && pdg <= PDGID::TOP_QUARK) || pdg == PDGID::GLUON);
}


LorentzVector jet_constit_4vec(const Jet & jet, const std::vector<PFParticle> & pfparticles, bool doPuppi) {
  LorentzVector sumv4;
  for (const uint i : jet.pfcand_indexs()) {
    if (!doPuppi) {
      sumv4 += pfparticles.at(i).v4();
    } else {
      LorentzVectorXYZE v4XYZ = toXYZ(pfparticles.at(i).v4());
      v4XYZ *= pfparticles.at(i).puppiWeight();
      sumv4 += toPtEtaPhi(v4XYZ);
    }
  }
  return sumv4;
}


LorentzVector genjet_constit_4vec(const GenJet & jet, const std::vector<GenParticle> & genparticles) {
  LorentzVector sumv4;
  for (const uint i : jet.genparticles_indices()) {
    sumv4 += genparticles.at(i).v4();
  }
  return sumv4;
}


GenJetClusterer::GenJetClusterer(uhh2::Context & ctx,
                                 const std::string & genjet_coll_name,
                                 float radius,
                                 const std::string & genparticles_coll_name,
                                 const std::string & genparticles_exclude_coll_name):
  genjet_handle_(ctx.get_handle<std::vector<GenJet>>(genjet_coll_name)),
  genparticle_handle_(ctx.get_handle<std::vector<GenParticle>>(genparticles_coll_name)),
  jet_def_(antikt_algorithm, radius),
  ghostScaling_(1E-16)
{
  if (genparticles_exclude_coll_name != "") {
    genparticle_exclude_handle_ = ctx.get_handle<std::vector<GenParticle>>(genparticles_exclude_coll_name);
  }
  gpId_ = AndId<GenParticle>(NoNeutrinoCut(), FinalStateCut());
  partonId_ = IsPartonCut();
}

PseudoJet GenJetClusterer::convert_uhh_genparticle_to_pseudojet(const GenParticle & particle) {
  LorentzVectorXYZE lv = toXYZ(particle.v4());
  PseudoJet outputParticle(lv.px(), lv.py(), lv.pz(), lv.E());
  return outputParticle;
}

GenJet GenJetClusterer::convert_pseudojet_to_uhh_genjet(const PseudoJet & jet, const std::vector<GenParticle> & genparticles) {
  GenJet genjet;
  // use XYZE v4 as safer?
  genjet.set_v4(toPtEtaPhi(LorentzVectorXYZE(jet.px(), jet.py(), jet.pz(), jet.E())));
  int pdgid = 0;
  float maxGhostPt = 0;
  // cout << "GenJet " << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << endl;
  for (const auto & dItr : jet.constituents()) {
    // find ghost partons, use them to assign flavour
    // dont store in main jet
    if (dItr.pt() < (1E5 * ghostScaling_)) {
      if (dItr.pt() > maxGhostPt) {
        maxGhostPt = dItr.pt();
        // nb pseudojet doesn't have pdgid
        pdgid = genparticles.at(dItr.user_index()).pdgId();
        // const GenParticle & gp = genparticles.at(dItr.user_index());
        // cout << "Found parton " << gp.pt() << " : " << gp.pdgId() << " : " << gp.status() << endl;
      }
    } else {
      genjet.add_genparticles_index(dItr.user_index());
    }
  }
  genjet.set_pdgId(pdgid);
  genjet.set_partonFlavour(pdgid);
  return genjet;
}

bool GenJetClusterer::process(uhh2::Event & event) {
  // Get GenParticles, convert to Pseudojets, ready for clustering
  vector<GenParticle> gps = event.get(genparticle_handle_);
  vector<GenParticle> * gps_reject = nullptr;
  if (event.is_valid(genparticle_exclude_handle_)) {
    gps_reject = & event.get(genparticle_exclude_handle_);
  }

  vector<PseudoJet> pjs; // Get all relevant gen particles for genjet itself
  vector<PseudoJet> pjsGhosts; // Get Gen Partons for ghost clustering
  int gpCounter = -1;
  for (auto gp : gps) {
    gpCounter++;

    // select partons for ghosts
    if (partonId_(gp, event)) {
      PseudoJet pj = convert_uhh_genparticle_to_pseudojet(gp);
      pj.set_user_index(gpCounter); // keep link to original GP
      // // Check for duplicates
      // if (std::find_if(pjs.begin(), pjs.end(), [&pj] (const PseudoJet &arg) {
      //   return (fabs(arg.pt()-pj.pt()) < 1E-1 && arg.delta_R(pj)<1E-2);
      // } ) != pjs.end()) continue;
      pj *= ghostScaling_;
      // pjsGhosts.push_back(pj);
      pjs.push_back(pj);
    }

    // Keep status = 1, no neutrinos
    if (!gpId_(gp, event)) continue;
    // Reject any veto particles
    if (gps_reject && gps_reject->size() > 0) {
      if (std::find(gps_reject->begin(), gps_reject->end(), gp) != gps_reject->end()) {
        // cout << "Found my reject " << gp.pt() << " : " << gp.eta() << " : " << gp.status() << " : " << gp.pdgId() << endl;
        continue;
      }
    }
    PseudoJet pj = convert_uhh_genparticle_to_pseudojet(gp);
    pj.set_user_index(gpCounter); // keep link to original GP
    // // ignore duplicates - bug in how I stored genparticles originally
    // // needs custom predicate, since PseudoJet doesn't have ==
    // if (std::find_if(pjs.begin(), pjs.end(), [&pj] (const PseudoJet &arg) {
    //     return (fabs(arg.pt()-pj.pt()) < 1E-1 && arg.delta_R(pj)<1E-2);
    //   } ) != pjs.end()) continue;
    pjs.push_back(pj);
  }

  // run the jet clustering with the above jet definition
  fastjet::ClusterSequence clust_seq(pjs, jet_def_);

  // get the resulting jets ordered in pt
  double ptmin = 15.0;
  vector<fastjet::PseudoJet> akJets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  // TODO: flavour assignment?

  // Convert back to GenJet
  vector<GenJet> genJets;
  for (auto fjjet : akJets) {
    GenJet genjet = convert_pseudojet_to_uhh_genjet(fjjet, gps);
    genJets.push_back(genjet);

    // if (genjet.pdgId() == 0) {
    //   cout << " NO FLAVOUR" << endl;
    //   cout << genjet.pt() << " : " << genjet.eta() << " : " << genjet.phi() << endl;
    //   for (auto & gp : gps) {
    //     if (!partonId_(gp, event)) continue;
    //     cout << gp.pdgId() << " : " << gp.status() << " : " << gp.pt() << " : " << gp.eta() << " : " << gp.phi() << endl;
    //   }
    // }
  }
  event.set(genjet_handle_, genJets);

  return true;
}


GenJetSelector::GenJetSelector(uhh2::Context & ctx,
                               float pt_min,
                               float y_max,
                               const std::string & genjet_coll_name,
                               const std::string & output_genjet_coll_name):
  jet_pt_min_(pt_min),
  jet_y_max_(y_max),
  genjet_handle_(ctx.get_handle<std::vector<GenJet>>(genjet_coll_name)),
  out_genjet_handle_(ctx.get_handle<std::vector<GenJet>>(output_genjet_coll_name))
{}

bool GenJetSelector::process(uhh2::Event & event) {
  std::vector<GenJet> genjets_out;
  for (const auto jet : event.get(genjet_handle_)) {
    // check constituents
    // occasionally get one with 0 pt so skip those
    // Get constituents
    // std::vector<GenParticle> constits = get_jet_genparticles(jet, event);
    // clean_collection<GenParticle>(constits, event, PtEtaCut(1E-8, 10.)); // basic cut to remove weird 0 pt constits

    if ((jet.pt() > jet_pt_min_) && (fabs(jet.Rapidity()) < jet_y_max_)){ // && (constits.size()>1)) {
      genjets_out.push_back(jet);
    }
  }
  sort_by_pt(genjets_out);
  event.set(out_genjet_handle_, genjets_out);
  return true;
}

std::vector<GenParticle> GenJetSelector::get_jet_genparticles(const GenJet & genjet, uhh2::Event & event) {
  std::vector<GenParticle> * genparticles = event.genparticles;
  std::vector<GenParticle> gp;
  for (const uint i : genjet.genparticles_indices()) {
    gp.push_back(genparticles->at(i)); // TODO store copy incase we shift it?
  }
  return gp;
}



FinalStateGenObjSelector::FinalStateGenObjSelector(uhh2::Context & ctx,
                                                   const GenId & genId,
                                                   const std::string & output_coll_name):
  genId_(genId),
  out_handle_(ctx.get_handle<std::vector<GenParticle>>(output_coll_name))
{}

bool FinalStateGenObjSelector::process(uhh2::Event & event) {
  std::vector<GenParticle> objs;
  // Do in reverse order to pick up most evolved objs first
  for (auto itr = event.genparticles->rbegin(); itr != event.genparticles->rend(); ++itr){
    // We check to see if we already have a very similar, but not exact, object
    // since the MC "evolves" the particle and slightly changes pt/eta/phi
    // Status check may not be reliable for e.g. herwig
    bool alreadyFound = std::any_of(objs.begin(), objs.end(), [&itr] (const GenParticle & mtr) { return deltaR(*itr, mtr) < 0.05; });
    if (FinalStateCut()(*itr, event)
        && genId_(*itr)
        && !alreadyFound) {
      objs.push_back(*itr);
    }
  }
  sort_by_pt(objs);
  event.set(out_handle_, objs);
  return true;
}


GenMuonSelector::GenMuonSelector(uhh2::Context & ctx,
                                 float pt_min,
                                 float eta_max,
                                 const std::string & output_coll_name):
  pt_min_(pt_min),
  eta_max_(eta_max),
  out_handle_(ctx.get_handle<std::vector<GenParticle>>(output_coll_name))
{}

bool GenMuonSelector::process(uhh2::Event & event) {
  std::vector<GenParticle> muons;
  // Do in reverse order to pick up most evolved muons first
  for (auto itr = event.genparticles->rbegin(); itr != event.genparticles->rend(); ++itr){
    // We check to see if we already have a very similar, but not exact, muon
    // since the MC "evolves" the particle and slightly changes pt/eta/phi
    // Status check may not be reliable for e.g. herwig
    bool alreadyFound = std::any_of(muons.begin(), muons.end(), [&itr] (const GenParticle & mtr) { return deltaR(*itr, mtr) < 0.05 && itr->charge() == mtr.charge(); });
    if ((abs(itr->pdgId()) == PDGID::MUON)
        && (itr->status() == 1)
        && (itr->pt() > pt_min_)
        && (fabs(itr->eta()) < eta_max_)
        && !alreadyFound) {
      muons.push_back(*itr);
    }
  }
  sort_by_pt(muons);
  event.set(out_handle_, muons);
  return true;
}


VariableBinning::VariableBinning(bins genBinning):
genBinning_(genBinning)
{
  nbins_gen_ = genBinning_.size() - 1;
  recoBinning_ = calculate_fine_binning(genBinning);
  nbins_reco_ = recoBinning_.size() - 1;
}

bins VariableBinning::calculate_fine_binning(const bins & coarse_bin_edges)
{
  bins fine_bin_edges;
  for (uint i=0; i<coarse_bin_edges.size()-1; i++) {
    fine_bin_edges.push_back(coarse_bin_edges[i]);
    fine_bin_edges.push_back(0.5*(coarse_bin_edges[i+1] + coarse_bin_edges[i]));
  }
  fine_bin_edges.push_back(coarse_bin_edges[coarse_bin_edges.size()-1]);
  return fine_bin_edges;
}

// initialise all static members - can't do in ctor
// Initialise pt bins for dijet
// -----------------------------------------------------------------------------
bins Binning::pt_bin_edges_gen = {
  50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000 // maximum should be 13 TeV / 2, but get weird empty reco binning
};
int Binning::nbins_pt_gen = pt_bin_edges_gen.size() - 1;

// separate "underflow" bin edges for underflow binning region
// see comments in QGAnalysisHists.cxx
bins Binning::pt_bin_edges_gen_underflow = {
  15, 30, 38, 50
};
int Binning::nbins_pt_gen_underflow = pt_bin_edges_gen_underflow.size() - 1;

// calculate reco (fine binned) bins using gen bins and making half as wide
bins Binning::pt_bin_edges_reco = Binning::calculate_fine_binning(pt_bin_edges_gen);
int Binning::nbins_pt_reco = pt_bin_edges_reco.size() - 1;

bins Binning::pt_bin_edges_reco_underflow = Binning::calculate_fine_binning(pt_bin_edges_gen_underflow);
int Binning::nbins_pt_reco_underflow = pt_bin_edges_reco_underflow.size() - 1;

// for non-TUnfold things, need all bins together, plus extra underflow bins
// also need to remove duplicate bin at start of signal region,
// hence the inplace ctor with begin()+1
bins Binning::pt_bin_edges_gen_all = Binning::sum_vectors({0.},
                                                          Binning::sum_vectors(pt_bin_edges_gen_underflow,
                                                                               bins(pt_bin_edges_gen.begin()+1, pt_bin_edges_gen.end())
                                                                              )
);
int Binning::nbins_pt_gen_all = pt_bin_edges_gen_all.size() - 1;

bins Binning::pt_bin_edges_reco_all = Binning::sum_vectors({0.},
                                                            Binning::sum_vectors(pt_bin_edges_reco_underflow,
                                                                                 bins(pt_bin_edges_reco.begin()+1, pt_bin_edges_reco.end())
                                                                                )
);
int Binning::nbins_pt_reco_all = pt_bin_edges_reco_all.size() - 1;

// Separate pt binning for Z+jets
// ------------------------------
// lower last big bin for Z+jets - dont want many empty bins for tunfold
bins Binning::pt_bin_edges_zpj_gen = {
  50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 2000 // maximum should be 13 TeV / 2, but no events in the last reco bin once split into 2, so go for 2TeV
};
int Binning::nbins_pt_zpj_gen = pt_bin_edges_zpj_gen.size() - 1;

bins Binning::pt_bin_edges_zpj_gen_underflow = {
  15, 30, 38, 50
};
int Binning::nbins_pt_zpj_gen_underflow = pt_bin_edges_zpj_gen_underflow.size() - 1;

bins Binning::pt_bin_edges_zpj_reco = Binning::calculate_fine_binning(pt_bin_edges_zpj_gen);
int Binning::nbins_pt_zpj_reco = pt_bin_edges_zpj_reco.size() - 1;

bins Binning::pt_bin_edges_zpj_reco_underflow = Binning::calculate_fine_binning(pt_bin_edges_zpj_gen_underflow);
int Binning::nbins_pt_zpj_reco_underflow = pt_bin_edges_zpj_reco_underflow.size() - 1;

bins Binning::pt_bin_edges_zpj_gen_all = Binning::sum_vectors({0.},
                                                              Binning::sum_vectors(pt_bin_edges_zpj_gen_underflow,
                                                                                   bins(pt_bin_edges_zpj_gen.begin()+1, pt_bin_edges_zpj_gen.end())
                                                                                  )
);
int Binning::nbins_pt_zpj_gen_all = pt_bin_edges_zpj_gen_all.size() - 1;

bins Binning::pt_bin_edges_zpj_reco_all = Binning::sum_vectors({0.},
                                                               Binning::sum_vectors(pt_bin_edges_zpj_reco_underflow,
                                                                                    bins(pt_bin_edges_zpj_reco.begin()+1, pt_bin_edges_zpj_reco.end())
                                                                                   )
);
int Binning::nbins_pt_zpj_reco_all = pt_bin_edges_zpj_reco_all.size() - 1;


// Setup lambda var bins - only need to specify gen, class internals auto calculate reco
// ---------------------------------------------------------------------------
// let the compiler do the conversion/init - anything more complex,
// put into a function that returns a map, then assign that to the static var
VarBinningMap Binning::var_binning_map_ungroomed = {
  {"LHA",
    // target 0.5, cen+fwd ungroomed, WTA axis, pt > 0, fixLambda
    VariableBinning({0.0, 0.17, 0.25, 0.32, 0.38, 0.45, 0.52, 0.59, 0.66, 1.0})
  },
  {"LHA_charged",
    // target 0.6, spike smoothing, cen+fwd ungroomed, WTA axis, pt > 0, fixLambda, charged clustering
    VariableBinning({0, 0.06, 0.11, 0.15, 0.19, 0.23, 0.27, 0.31, 0.35, 0.39, 0.44, 0.49, 0.54, 0.6, 1})
  },

  {"puppiMultiplicity",
    // based on target 0.5, cen+fwd, ensuring even interval between bins.
    // values above 22 added by hand otherwise no granularity.
    // subtract 0.5 as mostly integers
    VariableBinning({-0.5, 9.5, 15.5, 21.5, 29.5, 39.5, 59.5, 99.5, 149.5})
  },
  {"puppiMultiplicity_charged",
    // target 0.5, cen+fwd ungroomed, highPt, AK axis, pt > 1, ensuring even interval between bins.
    // 59.5 added manually. upper limit set to 99.5 manually
    VariableBinning({-0.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 17.5, 21.5, 25.5, 29.5, 35.5, 41.5, 59.5, 99.5})
  },

  {"pTD",
    // target 0.5, cen+fwd groomed, highPt, AK axis, pt > 1, >=2 constits, ensuring even interval between bins.
    // 60 added manually. upper limit set to 100 manually
    VariableBinning({0.0, 0.06, 0.09, 0.13, 0.19, 0.3, 1.0})
  },
  {"pTD_charged",
    // target 0.6, cen+fwd ungroomed, pt > 0, charged clustering
    VariableBinning({0, 0.07, 0.09, 0.12, 0.15, 0.19, 0.24, 0.31, 0.4, 0.52, 0.69, 1})
  },

  {"thrust",
    // target 0.5, cen+fwd, antiKT axis, pt > 0, fixLambda, manually set end to 0.6
    VariableBinning({0.0, 0.05, 0.09, 0.15, 0.205, 0.26, 0.6})
  },
  {"thrust_charged",
    // target 0.6, spike smoothing, cen+fwd ungroomed, antiKT axis, pt > 0, fixLambda, charged clustering, manually set end to 0.6
    VariableBinning({0, 0.005, 0.0125, 0.0225, 0.035, 0.05, 0.07, 0.0925, 0.12, 0.152, 0.188, 0.228, 0.273, 0.6})
  },

  {"width",
    // target 0.5, cen+fwd, WTA axis, pt > 0, fixLambda
    VariableBinning({0.0, 0.105, 0.165, 0.23, 0.305, 0.38, 0.46, 0.55, 1.0})
  },
  {"width_charged",
    // target 0.6, spike smoothing, cen+fwd ungroomed, WTA axis, pt > 0, fixLambda, charged clustering
    VariableBinning({0, 0.0225, 0.04, 0.0575, 0.0775, 0.1, 0.125, 0.152, 0.185, 0.22, 0.26, 0.307, 0.362, 0.425, 0.497, 1})
  }
};


VarBinningMap Binning::var_binning_map_groomed = {
  {"LHA",
    // target 0.5, cen+fwd groomed midPt, WTA axis, pt > 0, fixLambda
    VariableBinning({0, 0.1, 0.18, 0.26, 0.34, 0.42, 0.5, 0.57, 0.64, 1})
  },
  {"LHA_charged",
    // target 0.6, cen+fwd groomed midPt, spike smoothing, WTA axis, pt > 0, fixLambda, charged clustering
    VariableBinning({0, 0.06, 0.09, 0.12, 0.15, 0.19, 0.23, 0.27, 0.32, 0.37, 0.42, 0.48, 0.54, 0.6, 1})
  },

  {"puppiMultiplicity",
    // target 0.5, cen+fwd groomed midPt, pt > 1
    // add 39.5, 49.5, etc manually otherwise no bins
    VariableBinning({-0.5, 7.5, 13.5, 19.5, 29.5, 39.5, 49.5, 75.5, 99.5, 149.5})
  },
  {"puppiMultiplicity_charged",
    // target 0.6, cen+fwd groomed highPt, pt > 1
    // add 49.5 manually otherwise no bins
    VariableBinning({-0.5, 3.5, 5.5, 9.5, 13.5, 17.5, 21.5, 27.5, 35.5, 49.5, 99.5})
  },

  {"pTD",
    // target 0.5, cen+fwd groomed midPt, pt > 0
    VariableBinning({0, 0.06, 0.09, 0.14, 0.22, 1})
  },
  {"pTD_charged",
    // target 0.6, cen+fwd groomed midPt, spike smoothing, pt > 0, fixLambda, charged clustering
    VariableBinning({0, 0.07, 0.09, 0.12, 0.16, 0.21, 0.27, 0.35, 0.44, 0.54, 0.66, 1})
  },

  {"thrust",
    // target 0.5, cen+fwd groomed midPt, AK axis, pt > 0, fixLambda, manually set end to 0.6
    VariableBinning({0, 0.0025, 0.01, 0.025, 0.06, 0.12, 0.177, 0.23, 0.285, 0.6})
  },
  {"thrust_charged",
    // target 0.65, cen+fwd groomed midPt, spike smoothing, WTA axis, pt > 0, fixLambda, charged clustering, manually set end to 0.6
    VariableBinning({0, 0.0025, 0.005, 0.0075, 0.0125, 0.02, 0.0325, 0.05, 0.0775, 0.115, 0.16, 0.21, 0.265, 0.6})
  },

  {"width",
    // target 0.5, cen+fwd groomed midPt, WTA axis, pt > 0, fixLambda
    VariableBinning({0, 0.02, 0.05, 0.095, 0.147, 0.225, 0.307, 0.388, 0.468, 0.56, 1})
  },
  {"width_charged",
    // target 0.65, cen+fwd groomed midPt, spike smoothing, WTA axis, pt > 0, fixLambda, charged clustering
    VariableBinning({0, 0.0125, 0.0225, 0.035, 0.05, 0.07, 0.095, 0.128, 0.17, 0.225, 0.29, 0.365, 0.45, 1})
  }
};


const bins & Binning::var_bin_edges(const std::string & variable, bool is_groomed, bool is_reco) {
  VarBinningMap & this_map = (is_groomed) ? Binning::var_binning_map_groomed : Binning::var_binning_map_ungroomed;
  return this_map[variable].getBinning(is_reco);
}

int Binning::nbins_var(const std::string & variable, bool is_groomed, bool is_reco) {
  VarBinningMap & this_map = (is_groomed) ? Binning::var_binning_map_groomed : Binning::var_binning_map_ungroomed;
  return this_map[variable].getNbins(is_reco);
}

bins Binning::Binning::calculate_fine_binning(const bins & coarse_bin_edges)
{
  bins fine_bin_edges;
  for (uint i=0; i<coarse_bin_edges.size()-1; i++) {
    fine_bin_edges.push_back(coarse_bin_edges[i]);
    fine_bin_edges.push_back(0.5*(coarse_bin_edges[i+1] + coarse_bin_edges[i]));
  }
  fine_bin_edges.push_back(coarse_bin_edges[coarse_bin_edges.size()-1]);
  return fine_bin_edges;
}

bins Binning::sum_vectors(const bins & vec1, const bins & vec2)
{
  bins result(vec1);
  result.insert(result.end(), vec2.begin(), vec2.end());
  return result;
}

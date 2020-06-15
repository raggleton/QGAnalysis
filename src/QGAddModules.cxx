#include "UHH2/QGAnalysis/include/QGAddModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"

#include "UHH2/core/include/Utils.h"


using namespace std;
using namespace uhh2;
using namespace fastjet;
// using namespace fastjet::contrib;


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
  std::string resolutionFilename = "2016/Summer16_25nsV1_MC_PtResolution_" + jet_cone + "PF" + puName + ".txt";
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

  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("globalTightHalo2016Filter", "Flag_globalTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
  PrimaryVertexId pvid = StandardPrimaryVertexId();
  metfilters_selection->add<NPVSelection>("1 good PV", 1, -1, pvid);

  pv_cleaner.reset(new PrimaryVertexCleaner(pvid));

  electron_cleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_Summer16_medium, PtEtaCut(20.0, 2.5))));

  muon_cleaner.reset(new MuonCleaner(AndId<Muon>(MuonID(Muon::CutBasedIdMedium), PtEtaCut(26.0, 2.4), MuonID(Muon::PFIsoMedium))));
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


RecoJetSetup::RecoJetSetup(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, float jet_radius, float jet_pt_min, float jet_y_max, bool doJetID) {
  bool is_mc = ctx.get("dataset_type") == "MC";

  if (doJetID) jet_pf_id.reset(new JetCleaner(ctx, JetPFID(Cuts::RECO_JET_ID)));

  if (is_mc) {
    jet_met_corrector.reset(new MCJetMetCorrector(ctx, pu_removal, jet_cone));
  } else {
    jet_met_corrector.reset(new DataJetMetCorrector(ctx, pu_removal, jet_cone));
  }

  // eta cut to allow for all of the jet to fall inside tracker (<2.5)
  // jet_cleaner.reset(new JetCleaner(ctx, PtEtaCut(jet_pt_min, 2.5-jet_radius)));

  // eta cut to allow for both AK4 & AK8 jets to fall inside tracker (<2.5)
  // also note we cut on y not eta, since our jets can be massive
  jet_cleaner.reset(new JetCleaner(ctx, PtYCut(jet_pt_min, jet_y_max)));

  jet_ele_cleaner.reset(new JetElectronOverlapRemoval(jet_radius));

  jet_mu_cleaner.reset(new JetMuonOverlapRemoval(jet_radius));
}

bool RecoJetSetup::process(uhh2::Event & event) {

  if (jet_pf_id) jet_pf_id->process(event);

  jet_met_corrector->process(event);

  jet_cleaner->process(event);

  jet_ele_cleaner->process(event);

  jet_mu_cleaner->process(event);

  sort_by_pt(*event.jets);

  return true;
}

GenZllFinder::GenZllFinder(uhh2::Context & ctx, bool onlyStatus23, const std::string & genZhandle, const std::string & genZLeptonhandle):
  onlyStatus23_(onlyStatus23),
  gen_z_handle(ctx.get_handle<GenParticle>(genZhandle)),
  gen_z_leptons_handle(ctx.get_handle<std::vector<GenParticle>>(genZLeptonhandle))
{
}

bool GenZllFinder::process(uhh2::Event & event) {
  bool foundFirstLepton = false;
  GenParticle firstLepton, secondLepton;
  for (const auto & itr : *event.genparticles) {
    bool goodPDGID = (abs(itr.pdgId()) == PDGID::ELECTRON || abs(itr.pdgId()) == PDGID::MUON || abs(itr.pdgId()) == PDGID::TAU);
    bool goodStatus = onlyStatus23_ ? itr.status() == 23 : true;
    if (goodPDGID && goodStatus) {
      if (!foundFirstLepton) {
        foundFirstLepton = true;
        firstLepton.set_pdgId(itr.pdgId());
        firstLepton.set_v4(itr.v4());
        firstLepton.set_status(itr.status());
      } else if (itr.pdgId() == -1*firstLepton.pdgId()) {
        secondLepton.set_pdgId(itr.pdgId());
        secondLepton.set_v4(itr.v4());
        secondLepton.set_status(itr.status());
        break;
      }
    }
  }

  if (firstLepton.pt() != 0 && secondLepton.pt() != 0) {
    // int counter = 0;
    // for (const auto & itr : *event.genparticles) {
    //   counter++;
    //   cout << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << itr.v4().px() << " : " << itr.v4().py() << " : " << itr.v4().pz() << endl;
    //   if (counter == 100) break;
    // }
    // cout << " 1st lepton: " << firstLepton.pdgId() << " : " << firstLepton.status() << " : " << firstLepton.pt() << " : " << firstLepton.v4().px() << " : " << firstLepton.v4().py() << " : " << firstLepton.v4().pz() << endl;
    // cout << " 2nd lepton: " << secondLepton.pdgId() << " : " << secondLepton.status() << " : " << secondLepton.pt() << " : " << secondLepton.v4().px() << " : " << secondLepton.v4().py() << " : " << secondLepton.v4().pz() << endl;
    GenParticle genZ;
    genZ.set_v4(firstLepton.v4() + secondLepton.v4());
    genZ.set_pdgId(PDGID::Z); // in reality could be a photon
    event.set(gen_z_handle, genZ);
    std::vector<GenParticle> leptons = {firstLepton, secondLepton};
    event.set(gen_z_leptons_handle, leptons);
    return true;
  }

  cout << "**** No Z event: ****" << endl;
  int counter = 0;
  for (const auto & itr : *event.genparticles) {
    counter++;
    cout << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << itr.v4().px() << " : " << itr.v4().py() << " : " << itr.v4().pz() << endl;
    if (counter == 100) break;
  }
  return false;
}

ZkFactorReweight::ZkFactorReweight(uhh2::Context & ctx, const std::string & weightFilename_, const std::string & genmuon_name):
  z_weight_handle(ctx.get_handle<double>("z_weight"))
{
  if (weightFilename_ != "") {
    zReweight.reset(new ZllKFactor(weightFilename_));
  }
  if (genmuon_name != "") {
    gen_muon_handle = ctx.get_handle<vector<GenParticle>>(genmuon_name);
  }
}

bool ZkFactorReweight::process(uhh2::Event & event) {
  if (zReweight) {
    double realZPt = 0;
    // the k factor uses the Z as reconstructed by muons
    // TODO: dress them with photons?
    vector<GenParticle> genMuons = event.get(gen_muon_handle);
    if (genMuons.size() >= 2) {
      const auto & muon1 = genMuons.at(0);
      const auto & muon2 = genMuons.at(1);
      // cout << "muon1 mum1: " << muon1.mother1() << endl;
      // cout << "muon1 mum2: " << muon1.mother2() << endl;
      // cout << "muon2 mum1: " << muon2.mother1() << endl;
      // cout << "muon2 mum2: " << muon2.mother2() << endl;
      // TODO: some way to check if these are actually from the Z?
      const auto z_cand = muon1.v4() + muon2.v4();
      float z_pt = z_cand.pt();
      float z_mass = z_cand.M();
      // cout << "Z pt: " << z_pt << " Z mass: " << z_mass << endl;
      // z_mass and z_pt cuts are as per in paper, to avoid phase space with
      // large Sudakov contributions
      if (muon1.charge() != muon2.charge() && z_mass > 30 && z_pt > 30) {
        realZPt = z_pt;
      }
    }

    double zWeight = 1;
    if (realZPt > 0) zWeight = zReweight->getKFactor(realZPt);
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
    } else {
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


MCReweighting::MCReweighting(uhh2::Context & ctx, const std::string & genjet_name, const std::string & genmuon_name) {
  gen_weight_handle = ctx.get_handle<double>("gen_weight");

  lumi_weighter.reset(new MCLumiWeight(ctx));

  pileup_reweighter.reset(new MCPileupReweight(ctx, ctx.get("pileup_direction", "central")));

  std::string datasetVersion = ctx.get("dataset_version");

  is_DY = isSubstring(datasetVersion, "DYJetsToLL", true);
  doMuons = string2bool(ctx.get("isZPlusJets")); // this is set in MCModule
  cout << "doMuons: " << doMuons << endl;

  if (doMuons) {
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

  if (is_DY) z_reweighter.reset(new ZkFactorReweight(ctx, ctx.get("z_reweight_file", ""), genmuon_name));

  std::string pt_filename = ctx.get("pt_reweight_file", "");
  if (doMuons && ! is_DY) {
    // don't apply to non-DY Z+Jets MC (e.g. ttbar)
    pt_filename = "";
  }
  std::string region = is_DY ? "zplusjets" : "dijet";
  pt_reweighter.reset(new PtReweight(ctx, genjet_name, pt_filename, region));

  mc_scalevar.reset(new MCScaleVariation(ctx));
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


ZFinder::ZFinder(uhh2::Context & ctx, const std::string & inputLabel_, const std::string & outputLabel_):
  // Would like more generic FlavourParticle handle, but may need to do additional declare_event_input?
  hndlInput(ctx.get_handle<vector<Muon>>(inputLabel_)),
  hndlZ(ctx.get_handle<vector<Muon>>(outputLabel_))
{}

bool ZFinder::process(uhh2::Event & event) {
  // Reweight to higher order cross-section using gen level Z pT
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


float calcGenHT(const std::vector<GenParticle> & genparticles) {
  // Find scalar sum of all matrix-element partons
  // Only works for Pythia8 samples due to status check
  float ht = 0.;
  for (const auto & itr: genparticles) {
    if (abs(itr.status()) != 23) continue;
    uint pdg = abs(itr.pdgId());
    if (( pdg <= PDGID::TOP_QUARK && pdg >= PDGID::DOWN_QUARK) || pdg == PDGID::GLUON) {
      ht += itr.pt();
    }
  }
  return ht;
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
  ptSum_(0),
  jetVector_(jet_vector),
  wtaVector_(wta_vector),
  usePuppiWeight_(usePuppiWeight)
{
  // cache pt_sum
  for (auto & dtr : constits_) {
    float weight = usePuppiWeight_ ? dtr.puppiWeight() : 1.;
    ptSum_ += weight*dtr.pt();
  }
}

template<>
LambdaCalculator<GenParticle>::LambdaCalculator(std::vector<GenParticle> & constits, float jet_radius, const LorentzVector & jet_vector, const LorentzVector & wta_vector, bool usePuppiWeight):
  constits_(constits),
  jetRadius_(jet_radius),
  ptSum_(0),
  jetVector_(jet_vector),
  wtaVector_(wta_vector),
  usePuppiWeight_(usePuppiWeight)
{
  // cache pt_sum
  for (auto & dtr : constits_) {
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
      for (auto dtr : constits_) {
        result += dtr.puppiWeight();
      }
    } else {
      result = constits_.size();
    }
    resultsCache_[thisArgs] = result;
    return result;
  }

  for (auto dtr : constits_) {
    float weight = usePuppiWeight_ ? dtr.puppiWeight() : 1.;
    float z = (kappa != 0) ? dtr.pt() / ptSum_ : 1.;
    z *= weight;
    // for the reference jet vector for deltaR, we use the WTA vector
    // if beta <=1, and the normal jet vector otherwise
    float theta = (beta != 0) ? deltaRUsingY(dtr.v4(), (beta <= 1) ? wtaVector_ : jetVector_) / jetRadius_ : 1.; // 1 as puppi doesn't change direction
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
  for (auto dtr : constits_) {
    float z = (kappa != 0) ? dtr.pt() / ptSum_ : 1.;
    // for the reference jet vector for deltaR, we use the WTA vector
    // if beta <=1, and the normal jet vector otherwise
    float theta = (beta != 0) ? deltaRUsingY(dtr.v4(), (beta <= 1) ? wtaVector_ : jetVector_) / jetRadius_ : 1.;
    // cout << "    adding z= " << z << " theta: " << theta << " for constit " << dtr.pdgId() << " : " << dtr.pt() << " : " << dtr.eta() << " : " << dtr.phi() << endl;
    // cout << "    => " << pow(z, kappa) << " * " << pow(theta, beta) << " = " << pow(z, kappa) * pow(theta, beta) << endl;
    result += (pow(z, kappa) * pow(theta, beta));
  }
  resultsCache_[thisArgs] = result;
  // cout << "  result = " << result << endl;
  return result;
}

template<class T>
void LambdaCalculator<T>::clearCache()
{
  resultsCache_.clear();
}



QGAnalysisJetLambda::QGAnalysisJetLambda(uhh2::Context & ctx,
                                         float jetRadius,
                                         int nJetsMax,
                                         bool doPuppi,
                                         const PFParticleId & pfId,
                                         const std::string & jet_coll_name,
                                         const std::string & output_coll_name):
  wta_cluster_(Recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme), false, Recluster::keep_only_hardest)),
  // ca_cluster_(Recluster(JetDefinition(cambridge_algorithm, JetDefinition::max_allowable_R), false, Recluster::keep_only_hardest)),
  mmdt_(contrib::ModifiedMassDropTagger(0.1)),
  jetRadius_(jetRadius),
  nJetsMax_(nJetsMax),
  doPuppi_(doPuppi),
  pfId_(pfId),
  chargedHadronShift_(0), // these are the fractional shift, e.g. 0.2 for 20%
  neutralHadronShift_(0), // these are the fractional shift, e.g. 0.2 for 20%
  photonShift_(0), // must init to 0 otherwise bad things will happen, will go haywire!
  jet_handle_(ctx.get_handle<std::vector<Jet>>(jet_coll_name)),
  output_handle_(ctx.get_handle<std::vector<JetLambdaBundle>>(output_coll_name))
{
  mmdt_.set_grooming_mode();
  mmdt_.set_reclustering(true);
}


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
    // - get constits, apply energy shifts if need be
    // - Calculate the WTA axis
    // - Apply grooming to create subset of groomed constits
    // - Calculate WTA axis for groomed jet
    // - Create subsets of charged constits (groomed & ungroomed)
    // - Apply any other cuts (e.g. pt)
    // - Create LambdaCalculator for each constit collection, save to struct
    // {ungroomed, groomed} x {neutral+charged, charged-only}

    // Get constituents
    // don't apply puppi weight here, since we already account for it in the LambdaCalculator
    std::vector<PFParticle> constits = get_jet_pfparticles(jet, event, false);
    clean_collection<PFParticle>(constits, event, PtEtaCut(1E-8, 5)); // basic cut to remove weird 0 pt constits

    if (constits.size() < 2) {
      throw runtime_error("QGAnalysisJetLambda: constits.size() < 2, jet filtering not done properly!");
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
      PseudoJet thisDau = convert_uhh_pfparticle_to_pseudojet(dau, doPuppi_); // apply puppi weight here to recluster
      thisDau.set_user_index(dauCounter);
      pjconstits.push_back(thisDau);
      dauCounter++;
    }

    // Cluster using AKx to make pseudojet with right history
    JetDefinition jet_def(antikt_algorithm, JetDefinition::max_allowable_R); // don't use jet radius, as it can lead to several jets instead of the original
    std::vector<PseudoJet> akJets = jet_def(pjconstits);
    if (akJets.size() > 1) {
      cout << " >1 ak jets" << endl;
      cout << "Original jet: " << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << endl;
      for (const auto & subjet : akJets) {
        cout << "ak jet: " << subjet.pt() << " : " << subjet.eta() << " : " << subjet.phi() << endl;
      }
    } else if (akJets.size() == 0) {
      cout << "uhh jet constits:" << endl;
      for (const auto & dau : constits) {
        cout << "    " << dau.pt() << " : " << dau.eta() << " :" << dau.phi() << endl;
      }
      cout << "pjconstits:" << endl;
      for (const auto & dau : pjconstits) {
        cout << "    " << dau.pt() << " : " << dau.eta() << " :" << dau.phi() << endl;
      }
      throw std::runtime_error("AKx reclustering failed QGAnalysisJetLambda - no jets");
    }
    // cout << "Original jet: " << jet.pt()*jet.JEC_factor_raw() << " : " << jet.eta() << " : " << jet.phi() << endl;
    // cout << "ak jet: " << akJets[0].pt() << " : " << akJets[0].eta() << " : " << akJets[0].phi() << endl;

    // Now recluster using WTA
    PseudoJet wtaJet = wta_cluster_(akJets[0]);

    // if (wtaJet.constituents().size() != constits.size()) {
    //   throw std::runtime_error("WTA constits != genjet constits");
    // }

    // Apply grooming to jet (does CA reclustering internally)
    // Get WTA axis, and constits left after grooming
    // PseudoJet caJet = ca_cluster_(akJets[0]);
    // cout << "CA jet: " << caJet.pt() << " : " << caJet.eta() << " : " << caJet.phi() << endl;
    PseudoJet mmdtJet = mmdt_(akJets[0]);
    // cout << "MMDT jet: " << mmdtJet.pt() << " : " << mmdtJet.eta() << " : " << mmdtJet.phi() << endl;
    PseudoJet mmdtJetWTA = wta_cluster_(mmdtJet); // want the AK WTA axis for lambda variables
    // create collection with only those in groomed jet
    std::vector<PFParticle> groomedConstits;
    for (const auto & dItr : mmdtJet.constituents()) {
      groomedConstits.push_back(constits.at(dItr.user_index()));
    }

    // create charged-only subset of constits
    std::vector<PFParticle> chargedConstits(constits);
    std::vector<PFParticle> groomedChargedConstits(groomedConstits);
    clean_collection<PFParticle>(chargedConstits, event, ChargedCut());
    clean_collection<PFParticle>(groomedChargedConstits, event, ChargedCut());

    // Finally apply any pt cuts etc
    if (pfId_) {
      clean_collection<PFParticle>(constits, event, pfId_);
      clean_collection<PFParticle>(chargedConstits, event, pfId_);
      clean_collection<PFParticle>(groomedConstits, event, pfId_);
      clean_collection<PFParticle>(groomedChargedConstits, event, pfId_);
    }

    // cout << "Original jet axis: " << jet.v4().px() << " : " << jet.v4().py() << " : " << jet.v4().pz() << endl;
    // cout << "Original jet axis: " << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << endl;
    // cout << "WTA jet axis: " << wtaJetAxis.px() << " : " << wtaJetAxis.py() << " : " << wtaJetAxis.pz() << endl;
    // cout << "WTA jet axis: " << wtaJetAxis.pt() << " : " << wtaJetAxis.eta() << " : " << wtaJetAxis.phi() << endl;
    // cout << "WTA groomed jet axis: " << wtaGroomedJetAxis.px() << " : " << wtaGroomedJetAxis.py() << " : " << wtaGroomedJetAxis.pz() << endl;
    // cout << "WTA groomed jet axis: " << wtaGroomedJetAxis.pt() << " : " << wtaGroomedJetAxis.eta() << " : " << wtaGroomedJetAxis.phi() << endl;

    auto wtaJetAxis = toPtEtaPhi(LorentzVectorXYZE(wtaJet.px(), wtaJet.py(), wtaJet.pz(), wtaJet.E()));
    LambdaCalculator<PFParticle> recoJetCalc(constits, jetRadius_, jet.v4(), wtaJetAxis, doPuppi_);
    LambdaCalculator<PFParticle> recoJetCalcCharged(chargedConstits, jetRadius_, jet.v4(), wtaJetAxis, doPuppi_);

    auto jetAxisGroomed = toPtEtaPhi(LorentzVectorXYZE(mmdtJet.px(), mmdtJet.py(), mmdtJet.pz(), mmdtJet.E()));
    auto wtaJetAxisGroomed = toPtEtaPhi(LorentzVectorXYZE(mmdtJetWTA.px(), mmdtJetWTA.py(), mmdtJetWTA.pz(), mmdtJetWTA.E()));
    LambdaCalculator<PFParticle> recoJetCalcGroomed(groomedConstits, jetRadius_, jetAxisGroomed, wtaJetAxisGroomed, doPuppi_);
    LambdaCalculator<PFParticle> recoJetCalcGroomedCharged(groomedChargedConstits, jetRadius_, jetAxisGroomed, wtaJetAxisGroomed, doPuppi_);

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


QGAnalysisGenJetLambda::QGAnalysisGenJetLambda(uhh2::Context & ctx,
                                               float jetRadius,
                                               int nJetsMax,
                                               const GenParticleId & genId,
                                               const std::string & jet_coll_name,
                                               const std::string & output_coll_name):
  wta_cluster_(Recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme), false, Recluster::keep_only_hardest)),
  // ca_cluster_(Recluster(JetDefinition(cambridge_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme), false, Recluster::keep_only_hardest)),
  mmdt_(contrib::ModifiedMassDropTagger(0.1)),
  jetRadius_(jetRadius),
  nJetsMax_(nJetsMax),
  genId_(genId),
  neutralHadronShift_(0), // these are the fractional shift, e.g. 0.2 for 20%
  photonShift_(0),
  genjet_handle_(ctx.get_handle<std::vector<GenJet>>(jet_coll_name)),
  output_handle_(ctx.get_handle<std::vector<GenJetLambdaBundle>>(output_coll_name))
{
  mmdt_.set_grooming_mode();
  mmdt_.set_reclustering(true);
}


PseudoJet QGAnalysisGenJetLambda::convert_uhh_genparticle_to_pseudojet(const GenParticle & particle) {
  LorentzVectorXYZE lv = toXYZ(particle.v4());
  PseudoJet outputParticle(lv.px(), lv.py(), lv.pz(), lv.E());
  return outputParticle;
}


bool QGAnalysisGenJetLambda::process(uhh2::Event & event) {
  std::vector<GenJet> jets = event.get(genjet_handle_);
  std::vector<GenJetLambdaBundle> outputs;
  int nJetCounter = 0;
  // cout << "---- Doing QGAnalysisGenJetLambda::process ----" << endl;
  for (auto & jet : jets) {
    if (nJetsMax_ > 0 && nJetCounter == nJetsMax_) break;
    // For each jet, we:
    // - get constits, apply energy shifts if need be
    // - Calculate the WTA axis
    // - Apply grooming to create subset of groomed constits
    // - Calculate WTA axis for groomed jet
    // - Create subsets of charged constits (groomed & ungroomed)
    // - Apply any other cuts (e.g. pt)
    // - Create LambdaCalculator for each constit collection, save to struct
    // {ungroomed, groomed} x {neutral+charged, charged-only}

    // Get constituents
    std::vector<GenParticle> constits = get_jet_genparticles(jet, event);
    clean_collection<GenParticle>(constits, event, PtEtaCut(1E-8, 5)); // basic cut to remove weird 0 pt constits

    if (constits.size() < 2) {
      throw runtime_error("QGAnalysisGenJetLambda: constits.size() < 2, jet filtering not done properly!");
    }

    // FIXME: handle when 0 leftover constituents

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
    std::vector<PseudoJet> akJets = jet_def(pjconstits);
    if (akJets.size() > 1) {
      cout << " >1 ak jets" << endl;
      cout << "Original jet: " << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << endl;
      for (const auto & subjet : akJets) {
        cout << "ak jet: " << subjet.pt() << " : " << subjet.eta() << " : " << subjet.phi() << endl;
      }
    } else if (akJets.size() == 0) {
      cout << "uhh jet constits:" << endl;
      for (const auto & dau : constits) {
        cout << "    " << dau.pt() << " : " << dau.eta() << " :" << dau.phi() << endl;
      }
      cout << "pjconstits:" << endl;
      for (const auto & dau : pjconstits) {
        cout << "    " << dau.pt() << " : " << dau.eta() << " :" << dau.phi() << endl;
      }
      throw std::runtime_error("AKx reclustering failed QGAnalysisGenJetLambda - no jets");
    }

    // Now recluster using WTA
    PseudoJet wtaJet = wta_cluster_(akJets[0]);

    // if (wtaJet.constituents().size() != constits.size()) {
    //   throw std::runtime_error("WTA constits != genjet constits");
    // }

    // Apply grooming to jet, internally does CA reclustering
    // Get WTA axis, and constits left after grooming
    // PseudoJet caJet = ca_cluster_(akJets[0]);
    PseudoJet mmdtJet = mmdt_(akJets[0]);
    PseudoJet mmdtJetWTA = wta_cluster_(mmdtJet); // want the AK WTA axis for lambda variables
    // create collection with only those in groomed jet
    std::vector<GenParticle> groomedConstits;
    for (const auto & dItr : mmdtJet.constituents()) {
      groomedConstits.push_back(constits.at(dItr.user_index()));
    }

    // create charged-only subset of constits
    std::vector<GenParticle> chargedConstits(constits);
    std::vector<GenParticle> groomedChargedConstits(groomedConstits);
    clean_collection<GenParticle>(chargedConstits, event, ChargedCut());
    clean_collection<GenParticle>(groomedChargedConstits, event, ChargedCut());

    // Finally apply any pt cuts etc
    if (genId_) {
      clean_collection<GenParticle>(constits, event, genId_);
      clean_collection<GenParticle>(chargedConstits, event, genId_);
      clean_collection<GenParticle>(groomedConstits, event, genId_);
      clean_collection<GenParticle>(groomedChargedConstits, event, genId_);
    }

    // cout << "Original jet axis: " << jet.v4().px() << " : " << jet.v4().py() << " : " << jet.v4().pz() << endl;
    // cout << "Original jet axis: " << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << endl;
    // cout << "WTA jet axis: " << wtaJetAxis.px() << " : " << wtaJetAxis.py() << " : " << wtaJetAxis.pz() << endl;
    // cout << "WTA jet axis: " << wtaJetAxis.pt() << " : " << wtaJetAxis.eta() << " : " << wtaJetAxis.phi() << endl;
    // cout << "WTA groomed jet axis: " << wtaGroomedJetAxis.px() << " : " << wtaGroomedJetAxis.py() << " : " << wtaGroomedJetAxis.pz() << endl;
    // cout << "WTA groomed jet axis: " << wtaGroomedJetAxis.pt() << " : " << wtaGroomedJetAxis.eta() << " : " << wtaGroomedJetAxis.phi() << endl;

    auto wtaJetAxis = toPtEtaPhi(LorentzVectorXYZE(wtaJet.px(), wtaJet.py(), wtaJet.pz(), wtaJet.E()));
    LambdaCalculator<GenParticle> genJetCalc(constits, jetRadius_, jet.v4(), wtaJetAxis, false);
    LambdaCalculator<GenParticle> genJetCalcCharged(chargedConstits, jetRadius_, jet.v4(), wtaJetAxis, false);

    auto jetAxisGroomed = toPtEtaPhi(LorentzVectorXYZE(mmdtJet.px(), mmdtJet.py(), mmdtJet.pz(), mmdtJet.E()));
    auto wtaJetAxisGroomed = toPtEtaPhi(LorentzVectorXYZE(mmdtJetWTA.px(), mmdtJetWTA.py(), mmdtJetWTA.pz(), mmdtJetWTA.E()));
    LambdaCalculator<GenParticle> genJetCalcGroomed(groomedConstits, jetRadius_, jetAxisGroomed, wtaJetAxisGroomed, false);
    LambdaCalculator<GenParticle> genJetCalcGroomedCharged(groomedChargedConstits, jetRadius_, jetAxisGroomed, wtaJetAxisGroomed, false);

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
  std::vector<GenParticle> gp;
  for (const uint i : genjet.genparticles_indices()) {
    gp.push_back(genparticles->at(i)); // TODO store copy incase we shift it?
  }
  return gp;
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


GenJetClusterer::GenJetClusterer(uhh2::Context & ctx, const std::string & genjet_coll_name, float radius, const std::string & genparticles_coll_name):
  genjet_handle_(ctx.get_handle<std::vector<GenJet>>(genjet_coll_name)),
  genparticle_handle_(ctx.get_handle<std::vector<GenParticle>>(genparticles_coll_name)),
  jet_def_(antikt_algorithm, radius)
{
  gpId_ = AndId<GenParticle>(NoNeutrinoCut(), FinalStateCut());
}

PseudoJet GenJetClusterer::convert_uhh_genparticle_to_pseudojet(const GenParticle & particle) {
  LorentzVectorXYZE lv = toXYZ(particle.v4());
  PseudoJet outputParticle(lv.px(), lv.py(), lv.pz(), lv.E());
  return outputParticle;
}

GenJet GenJetClusterer::convert_pseudojet_to_uhh_genjet(const PseudoJet & jet) {
  GenJet genjet;
  // use XYZE v4 as safer?
  genjet.set_v4(toPtEtaPhi(LorentzVectorXYZE(jet.px(), jet.py(), jet.pz(), jet.E())));
  for (const auto & dItr : jet.constituents()) {
    genjet.add_genparticles_index(dItr.user_index());
  }
  return genjet;
}

bool GenJetClusterer::process(uhh2::Event & event) {
  // Get GenParticles, convert to Pseudojets, ready for clustering
  vector<GenParticle> gps = event.get(genparticle_handle_);

  vector<PseudoJet> pjs;
  int gpCounter = -1;
  for (auto gp : gps) {
    gpCounter++;
    // Keep status = 1, no neutrinos
    if (!gpId_(gp, event)) continue;
    PseudoJet pj = convert_uhh_genparticle_to_pseudojet(gp);
    pj.set_user_index(gpCounter); // keep link to original GP
    // ignore duplicates - bug in how I stored genparticles originally
    // needs custom predicate, since PseudoJet doesn't have ==
    if (std::find_if(pjs.begin(), pjs.end(), [&pj] (const PseudoJet &arg) {
        return (fabs(arg.pt()-pj.pt()) < 1E-1 && arg.delta_R(pj)<1E-2);
      } ) != pjs.end()) continue;
    pjs.push_back(pj);
  }
  // cout << "Had " << gps.size() << " now have " << pjs.size() << endl;

  // run the jet clustering with the above jet definition
  fastjet::ClusterSequence clust_seq(pjs, jet_def_);

  // get the resulting jets ordered in pt
  double ptmin = 15.0;
  vector<fastjet::PseudoJet> akJets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  // TODO: flavour assignment?

  // Convert back to GenJet
  vector<GenJet> genJets;
  for (auto fjjet : akJets) {
    GenJet genjet = convert_pseudojet_to_uhh_genjet(fjjet);
    genJets.push_back(genjet);
  }
  event.set(genjet_handle_, genJets);

  return true;
}


GenJetSelector::GenJetSelector(uhh2::Context & ctx,
                               float pt_min,
                               float y_max,
                               float lepton_overlap_dr,
                               const std::string & genjet_coll_name,
                               const std::string & output_genjet_coll_name,
                               const std::string & genparticles_coll_name):
  jet_pt_min_(pt_min),
  jet_y_max_(y_max),
  lepton_overlap_dr_(lepton_overlap_dr),
  genjet_handle_(ctx.get_handle<std::vector<GenJet>>(genjet_coll_name)),
  out_genjet_handle_(ctx.get_handle<std::vector<GenJet>>(output_genjet_coll_name)),
  genparticle_handle_(ctx.get_handle<std::vector<GenParticle>>(genparticles_coll_name))
{}

bool GenJetSelector::process(uhh2::Event & event) {
  std::vector<GenJet> genjets_out;
  const std::vector<GenParticle> & genparticles = event.get(genparticle_handle_);
  for (const auto jet : event.get(genjet_handle_)) {
      bool found = (std::find(genjets_out.begin(), genjets_out.end(), jet) != genjets_out.end()); // avoid duplicate genjets
      // avoid jets that are just leptons + a few spurious gluons,
      // i.e. look for overlapping ones, and check their pT fraction
      bool leptonOverlap = false;
      for (const auto & gp : genparticles) {
        bool isLepton = ((abs(gp.pdgId()) == PDGID::MUON) || (abs(gp.pdgId()) == PDGID::ELECTRON) || (abs(gp.pdgId()) == PDGID::TAU));
        if (!isLepton) continue;
        bool thisLeptonOverlap = isLepton && (deltaRUsingY(gp.v4(), jet.v4()) < lepton_overlap_dr_) && ((gp.pt() / jet.pt()) > 0.5);
        leptonOverlap = leptonOverlap || thisLeptonOverlap;
      }
      // check constituents
      // occasionally get one with 0 pt so skip those
      // Get constituents
      std::vector<GenParticle> constits = get_jet_genparticles(jet, event);
      clean_collection<GenParticle>(constits, event, PtEtaCut(1E-8, 5)); // basic cut to remove weird 0 pt constits

      if ((jet.pt() > jet_pt_min_) && (fabs(jet.Rapidity()) < jet_y_max_) && !found && !leptonOverlap && (constits.size()>1)) {
        genjets_out.push_back(jet);
      }
    }
    sort_by_pt(genjets_out);
    event.set(out_genjet_handle_, genjets_out);
    return true;
};

std::vector<GenParticle> GenJetSelector::get_jet_genparticles(const GenJet & genjet, uhh2::Event & event) {
  std::vector<GenParticle> * genparticles = event.genparticles;
  std::vector<GenParticle> gp;
  for (const uint i : genjet.genparticles_indices()) {
    gp.push_back(genparticles->at(i)); // TODO store copy incase we shift it?
  }
  return gp;
}

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

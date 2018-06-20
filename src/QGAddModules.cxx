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
using namespace uhh2examples;


DataJetMetCorrector::DataJetMetCorrector(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone){
  std::vector<std::string> JEC_BCD, JEC_EFearly, JEC_FlateG, JEC_H;
  if (pu_removal == "CHS") {
    if (jet_cone == "AK4") {
      JEC_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA;
      JEC_EFearly = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA;
      JEC_FlateG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA;
      JEC_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA;
    } else if (jet_cone == "AK8") {
      JEC_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFchs_DATA;
      JEC_EFearly = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFchs_DATA;
      JEC_FlateG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFchs_DATA;
      JEC_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFchs_DATA;
    } else {
      throw runtime_error("CHS must have jet_cone of AK4 or AK8");
    }
  } else if (pu_removal == "PUPPI") {
    if (jet_cone == "AK4") {
      JEC_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFPuppi_DATA;
      JEC_EFearly = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFPuppi_DATA;
      JEC_FlateG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFPuppi_DATA;
      JEC_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFPuppi_DATA;
    } else if (jet_cone == "AK8") {
      JEC_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFPuppi_DATA;
      JEC_EFearly = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFPuppi_DATA;
      JEC_FlateG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFPuppi_DATA;
      JEC_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFPuppi_DATA;
    } else {
      throw runtime_error("PUPPI must have jet_cone of AK4 or AK8");
    }
  }

  jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_BCD));
  jet_corrector_EFearly.reset(new JetCorrector(ctx, JEC_EFearly));
  jet_corrector_FlateG.reset(new JetCorrector(ctx, JEC_FlateG));
  jet_corrector_H.reset(new JetCorrector(ctx, JEC_H));
}

bool DataJetMetCorrector::process(uhh2::Event & event) {
  if (event.run <= runnr_BCD) {
    jet_corrector_BCD->process(event);
    jet_corrector_BCD->correct_met(event);
  } else if (event.run < runnr_EFearly) { //< is correct, not <=
    jet_corrector_EFearly->process(event);
    jet_corrector_EFearly->correct_met(event);
  } else if (event.run <= runnr_FlateG) {
    jet_corrector_FlateG->process(event);
    jet_corrector_FlateG->correct_met(event);
  } else if (event.run > runnr_FlateG) {
    jet_corrector_H->process(event);
    jet_corrector_H->correct_met(event);
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
      JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
      resolutionFilename = "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt";
    } else if (jet_cone == "AK8") {
      JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK8PFchs_MC;
      resolutionFilename = "Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt";  // actually doesn't matter, they're all the same
    } else {
      throw runtime_error("CHS must have jet_cone of AK4 or AK8");
    }
  } else if (pu_removal == "PUPPI") {
    if (jet_cone == "AK4") {
      JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFPuppi_MC;
      resolutionFilename = "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt";
    } else if (jet_cone == "AK8") {
      JEC_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK8PFPuppi_MC;
      resolutionFilename = "Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt";
    } else {
      throw runtime_error("PUPPI must have jet_cone of AK4 or AK8");
    }
  }
  jet_corrector.reset(new JetCorrector(ctx, JEC_MC));
  jet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, jet_coll_name, genjet_coll_name, true, JERSmearing::SF_13TeV_2016_03Feb2017, resolutionFilename));
}

bool MCJetMetCorrector::process(uhh2::Event & event) {
  jet_corrector->process(event);
  jet_corrector->correct_met(event);
  jet_resolution_smearer->process(event);
  return true;
}


GeneralEventSetup::GeneralEventSetup(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, float jet_radius) {
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

  jet_pf_id.reset(new JetCleaner(ctx, JetPFID(JetPFID::wp::WP_LOOSE)));

  if (is_mc) {
    jet_met_corrector.reset(new MCJetMetCorrector(ctx, pu_removal, jet_cone));
  } else {
    jet_met_corrector.reset(new DataJetMetCorrector(ctx, pu_removal, jet_cone));
  }

  jet_cleaner.reset(new JetCleaner(ctx, PtEtaCut(30.0, 2.4)));

  jet_ele_cleaner.reset(new JetElectronOverlapRemoval(jet_radius));

  jet_mu_cleaner.reset(new JetMuonOverlapRemoval(jet_radius));
}

bool GeneralEventSetup::process(uhh2::Event & event) {

  if(event.isRealData && !lumi_selection->passes(event)) return false;

  if(!metfilters_selection->passes(event)) return false;

  pv_cleaner->process(event);

  electron_cleaner->process(event);

  muon_cleaner->process(event);

  jet_pf_id->process(event);

  jet_met_corrector->process(event);

  jet_cleaner->process(event);

  jet_ele_cleaner->process(event);

  jet_mu_cleaner->process(event);

  sort_by_pt(*event.jets);

  return true;
}


MCReweighting::MCReweighting(uhh2::Context & ctx) {
  lumi_weighter.reset(new MCLumiWeight(ctx));

  pileup_reweighter.reset(new MCPileupReweight(ctx, "central"));

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
  if (event.isRealData) return true;

  lumi_weighter->process(event);

  pileup_reweighter->process(event);

  muon_id_reweighter_pt_eta->process(event);

  // muon_id_reweighter_vtx->process(event);

  muon_trg_reweighter->process(event);

  muon_trk_reweighter->process(event);

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
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
    lumi_weighter.reset(new MCLumiWeight(ctx));
    pileup_reweighter.reset(new MCPileupReweight(ctx, "central"));
  } else {
    jet_met_corrector.reset(new DataJetMetCorrector(ctx, pu_removal, jet_cone));
  }

  jet_cleaner.reset(new JetCleaner(ctx, PtEtaCut(30.0, 4.7)));

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

  if (is_mc) {
    lumi_weighter->process(event);
    pileup_reweighter->process(event);
  }

  jet_cleaner->process(event);

  jet_ele_cleaner->process(event);

  jet_mu_cleaner->process(event);

  sort_by_pt(*event.jets);

  return true;
}
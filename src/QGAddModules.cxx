#include "UHH2/QGAnalysis/include/QGAddModules.h"

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

std::ostream& Color::operator<<(std::ostream& os, Code code) {
      return os << "\033[" << static_cast<int>(code) << "m";
}


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

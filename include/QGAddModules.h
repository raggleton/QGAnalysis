#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/common/include/JetCorrections.h"


namespace Color {
  enum Code {
      FG_RED      = 31,
      FG_GREEN    = 32,
      FG_YELLOW   = 33,
      FG_BLUE     = 34,
      FG_MAGENTA  = 35,
      FG_CYAN     = 36,
      FG_DEFAULT  = 39,
      BG_RED      = 41,
      BG_GREEN    = 42,
      BG_BLUE     = 44,
      BG_DEFAULT  = 49
  };
  std::ostream& operator<<(std::ostream& os, Code code);
  //  {
  //     return os << "\033[" << static_cast<int>(code) << "m";
  // };
};

namespace uhh2examples{

// Easy way to refer to PDGIDs
enum PDGID {
  UNKNOWN = 0,
  DOWN_QUARK = 1,
  UP_QUARK = 2,
  STRANGE_QUARK = 3,
  CHARM_QUARK = 4,
  BOTTOM_QUARK = 5,
  TOP_QUARK = 6,
  ELECTRON = 11,
  MUON = 13,
  TAU = 15,
  GLUON = 21
};


/**
 * Easy way to do JEC, JER, MET corrections in Data. Takes care of run period dependency.
 */
class DataJetMetCorrector : public uhh2::AnalysisModule {
public:
  explicit DataJetMetCorrector(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<JetCorrector> jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
  const int runnr_BCD = 276811;
  const int runnr_EFearly = 278802;
  const int runnr_FlateG = 280385;
};


/**
 * Easy way to do JEC, JER, MET corrections in MC
 */
class MCJetMetCorrector : public uhh2::AnalysisModule {
public:
  explicit MCJetMetCorrector(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, const std::string & jet_coll_name="jets", const std::string & genjet_coll_name="genjets");
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<JetCorrector> jet_corrector;
  std::unique_ptr<GenericJetResolutionSmearer> jet_resolution_smearer;
};


// Helper funcs
void printGenParticles(const std::vector<GenParticle> & gps, const std::string & info="", Color::Code color=Color::FG_DEFAULT);

std::vector<GenParticle*> print_genjet_genparticles(const GenJetWithParts & jet, std::vector<GenParticle>* genparticles);

void printGenJets(const std::vector<GenJetWithParts> & gps, std::vector<GenParticle>* genparticles, const std::string & info="", Color::Code color=Color::FG_BLUE);

void printJets(const std::vector<Jet> & jets, const std::string & info="", Color::Code color=Color::FG_GREEN);

void printMuons(const std::vector<Muon> & muons, const std::string & info="", Color::Code color=Color::FG_RED);

void printElectrons(const std::vector<Electron> & electrons, const std::string & info="", Color::Code color=Color::FG_YELLOW);

};
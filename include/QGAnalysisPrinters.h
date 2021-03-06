#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"

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
    BG_YELLOW   = 43,
    BG_BLUE     = 44,
    BG_MAGENTA  = 45,
    BG_CYAN     = 46,
    BG_WHITE    = 47,
    BG_DEFAULT  = 49
  };
  std::ostream& operator<<(std::ostream& os, Code code);
};

namespace uhh2examples {

// Helper funcs
void printGenParticles(const std::vector<GenParticle> & gps, const std::string & label="", Color::Code color=Color::FG_DEFAULT);

void print_genjet_genparticles(const GenJetWithParts & jet, const std::vector<GenParticle>* genparticles);

void printPFParticles(const std::vector<PFParticle> & pfparticles, const std::string & label="", Color::Code color=Color::FG_DEFAULT);

void print_jet_pfparticles(const Jet & jet, const std::vector<PFParticle>* pfparticles);

void printGenJets(const std::vector<GenJetWithParts> & gps, const std::string & label="", Color::Code color=Color::FG_BLUE);

void printGenJetsWithParts(const std::vector<GenJetWithParts> & gps, const std::vector<GenParticle>* genparticles, const std::string & label="", Color::Code color=Color::FG_BLUE);

void printJets(const std::vector<Jet> & jets, const std::string & label="", Color::Code color=Color::FG_GREEN);

void printJetsWithParts(const std::vector<Jet> & jets, const std::vector<PFParticle> * pfparticles, const std::string & label="", Color::Code color=Color::FG_GREEN);

void printMuons(const std::vector<Muon> & muons, const std::string & label="", Color::Code color=Color::FG_RED);

void printElectrons(const std::vector<Electron> & electrons, const std::string & label="", Color::Code color=Color::FG_YELLOW);

};
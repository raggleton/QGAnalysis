#include "UHH2/QGAnalysis/include/QGAnalysisPrinters.h"

using namespace std;
using namespace uhh2;
using namespace uhh2examples;


std::ostream& Color::operator<<(std::ostream& os, Code code) {
  return os << "\033[" << static_cast<int>(code) << "m";
}



void uhh2examples::printGenParticles(const std::vector<GenParticle> & gps, const std::string & label, Color::Code color) {
  for (auto & itr: gps) {
    // if (itr.status() != 1) continue;
    cout << color << "GP";
    if (label != "") {
      cout << " [" << label << "]";
    }
    cout << ": " << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << std::endl;
  }
}

std::vector<GenParticle*> uhh2examples::print_genjet_genparticles(const GenJetWithParts & jet, std::vector<GenParticle>* genparticles) {
  std::vector<GenParticle*> gp;
  for (const uint i : jet.genparticles_indices()) {
    gp.push_back(&(genparticles->at(i)));
  }
  for (const auto *itr : gp) {
    std::cout << itr->pdgId() << " : " << itr->pt() << " : " << deltaR(jet.v4(), itr->v4()) << std::endl;
  }
  return gp;
}

void uhh2examples::printGenJets(const std::vector<GenJetWithParts> & gps, const std::string & label, Color::Code color) {
  for (auto & itr: gps) {
    std::cout << color << "GenJet";
    if (label != "") {
        std::cout << " [" << label << "]";
    }
    std::cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << std::endl;
  }
}

void uhh2examples::printGenJetsWithParts(const std::vector<GenJetWithParts> & gps, std::vector<GenParticle>* genparticles, const std::string & label, Color::Code color) {
  for (auto & itr: gps) {
    std::cout << color << "GenJet";
    if (label != "") {
        std::cout << " [" << label << "]";
    }
    std::cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << std::endl;
    print_genjet_genparticles(itr, genparticles);  // requires the print* funcs to have namespace explicitly written
  }
}

void uhh2examples::printJets(const std::vector<Jet> & jets, const std::string & label, Color::Code color) {
  for (auto & itr: jets) {
    std::cout << color << "jet";
    if (label != "") {
        std::cout << " [" << label << "]";
    }
    std::cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << " : GenJet: " << itr.genjet_index() << Color::FG_DEFAULT << std::endl;
  }
}

void uhh2examples::printMuons(const std::vector<Muon> & muons, const std::string & label, Color::Code color) {
  for (auto & itr: muons) {
    std::cout << color << "muon";
    if (label != "") {
        std::cout << " [" << label << "]";
    }
    std::cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << " : " << itr.relIso() << Color::FG_DEFAULT << std::endl;
  }
}

void uhh2examples::printElectrons(const std::vector<Electron> & electrons, const std::string & label, Color::Code color) {
  for (auto & itr: electrons) {
    std::cout << color << "electron";
    if (label != "") {
        std::cout << " [" << label << "]";
    }
    std::cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << std::endl;
  }
}
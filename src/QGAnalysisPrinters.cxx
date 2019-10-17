#include "UHH2/QGAnalysis/include/QGAnalysisPrinters.h"
#include "UHH2/core/include/Utils.h"

using namespace std;
using namespace uhh2;
using namespace uhh2examples;


std::ostream& Color::operator<<(std::ostream& os, Code code) {
  return os << "\033[" << static_cast<int>(code) << "m";
}



void uhh2examples::printGenParticles(const std::vector<GenParticle> & gps, const std::string & label, Color::Code color) {
  int counter = 0;
  for (auto & itr: gps) {
    // if (itr.status() != 1) continue;
    cout << color << "GP" << std::to_string(counter);
    if (label != "") {
      cout << " [" << label << "]";
    }
    cout << ": " << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << std::endl;
    counter++;
  }
}

void uhh2examples::print_genjet_genparticles(const GenJetWithParts & jet, const std::vector<GenParticle>* genparticles) {
  for (const uint i : jet.genparticles_indices()) {
    std::cout << genparticles->at(i).pdgId() << " : " << genparticles->at(i).pt() << " : " << deltaR(jet.v4(), genparticles->at(i).v4()) << std::endl;
  }
}

void uhh2examples::print_jet_pfparticles(const Jet & jet, const std::vector<PFParticle>* pfparticles) {
  for (const uint i : jet.daughterIndices()) {
    std::cout << pfparticles->at(i).pt() << " : " << pfparticles->at(i).charge() << " : " << deltaR(jet.v4(), pfparticles->at(i).v4()) << std::endl;
  }
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

void uhh2examples::printGenJetsWithParts(const std::vector<GenJetWithParts> & gps, const std::vector<GenParticle>* genparticles, const std::string & label, Color::Code color) {
  for (auto & itr: gps) {
    std::cout << color << "GenJet";
    if (label != "") {
        std::cout << " [" << label << "]";
    }
    std::cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << Color::FG_DEFAULT << std::endl;
    uhh2examples::print_genjet_genparticles(itr, genparticles);  // requires the print* funcs to have namespace explicitly written
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


void uhh2examples::printJetsWithParts(const std::vector<Jet> & jets, const std::vector<PFParticle> * pfparticles, const std::string & label, Color::Code color) {
  for (auto & itr: jets) {
    std::cout << color << "jet";
    if (label != "") {
        std::cout << " [" << label << "]";
    }
    std::cout << ": " << itr.pt() << " : " << itr.eta() << " : " << itr.phi() << " : GenJet: " << itr.genjet_index() << Color::FG_DEFAULT << std::endl;
    uhh2examples::print_jet_pfparticles(itr, pfparticles);
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
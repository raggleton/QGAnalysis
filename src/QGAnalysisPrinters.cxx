#include "UHH2/QGAnalysis/include/QGAnalysisPrinters.h"
#include "UHH2/core/include/Utils.h"

using namespace std;
using namespace uhh2;
using namespace uhh2examples;


std::ostream& Color::operator<<(std::ostream& os, Code code) {
  return os << "\033[" << static_cast<int>(code) << "m";
}



void uhh2examples::printGenParticles(const vector<GenParticle> & gps, const string & label, Color::Code color) {
  if (gps.size() == 0) return;

  int counter = 0;
  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
  cout << color << "GenParticle : pdgId() : status() : pt() : y() : phi() " << Color::FG_DEFAULT << endl;
  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;

  for (auto & itr: gps) {
    // if (itr.status() != 1) continue;
    cout << color << "GP" << to_string(counter);
    if (label != "") {
      cout << "[" << label << "]";
    }
    cout << ": " << itr.pdgId() << " : " << itr.status() << " : " << itr.pt() << " : " << itr.Rapidity() << " : " << itr.phi() << Color::FG_DEFAULT << endl;
    counter++;
  }

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
}

void uhh2examples::print_genjet_genparticles(const GenJet & jet, const vector<GenParticle>* genparticles) {
  if (jet.genparticles_indices().size() == 0) return;

  cout << "  ................................................................................" << endl;
  cout << "  GenJet's GenParticle : pdgId() : pt() : deltaR(jet, genparticle)" << endl;
  cout << "  ................................................................................" << endl;

  LorentzVector sumv4;
  for (const uint i : jet.genparticles_indices()) {
    cout << "  " << genparticles->at(i).pdgId() << " : " << genparticles->at(i).pt() << " : " << deltaRUsingY(jet.v4(), genparticles->at(i).v4()) << endl;
    sumv4 += genparticles->at(i).v4();
  }
  cout << "SUM: " << sumv4.pt() << " : " << sumv4.Rapidity() << " : " << sumv4.phi() << endl;
  cout << "  ................................................................................" << endl;
}

void uhh2examples::printPFParticles(const vector<PFParticle> & pfparticles, const string & label, Color::Code color) {
  if (pfparticles.size() == 0) return;

  cout << "--------------------------------------------------------------------------------" << endl;
  cout << "PFParticle : charge() : particleID() : pt() : y() : phi() : puppiWeight()" << endl;
  cout << "--------------------------------------------------------------------------------" << endl;

  int counter = 0;
  for (auto & itr: pfparticles) {
    cout << color << "PF" << to_string(counter);
    if (label != "") {
      cout << "[" << label << "]";
    }
    cout << ": " << itr.charge() << " : " << itr.particleID() << " : " << itr.pt() << " : " << itr.Rapidity() << " : " << itr.phi() << " : " << itr.puppiWeight() << Color::FG_DEFAULT << endl;
    counter++;
  }

  cout << "--------------------------------------------------------------------------------" << endl;
}

void uhh2examples::print_jet_pfparticles(const Jet & jet, const vector<PFParticle>* pfparticles) {
  if (jet.pfcand_indexs().size() == 0) return;

  cout << "  ................................................................................" << endl;
  cout << "  Jet's PFParticle : pt() : charge() : puppiWeight() : deltaR(jet, PFparticle)" << endl;
  cout << "  ................................................................................" << endl;

  LorentzVector sumv4;
  LorentzVector sumv4Puppi;
  for (const uint i : jet.pfcand_indexs()) {
    cout << "  " << pfparticles->at(i).pt() << " : " << pfparticles->at(i).charge() << " : " << pfparticles->at(i).puppiWeight() << " : " << deltaRUsingY(jet.v4(), pfparticles->at(i).v4()) << endl;
    sumv4 += pfparticles->at(i).v4();
    LorentzVectorXYZE v4XYZ = toXYZ(pfparticles->at(i).v4());
    v4XYZ *= pfparticles->at(i).puppiWeight();
    sumv4Puppi += toPtEtaPhi(v4XYZ);
  }
  cout << "SUM (no Puppi Weight): " << sumv4.pt() << " : " << sumv4.Rapidity() << " : " << sumv4.phi() << endl;
  cout << "SUM (Puppi Weight): " << sumv4Puppi.pt() << " : " << sumv4Puppi.Rapidity() << " : " << sumv4Puppi.phi() << endl;
  cout << "  ................................................................................" << endl;
}

void uhh2examples::printGenJets(const vector<GenJet> & gps, const string & label, Color::Code color) {
  if (gps.size() == 0) return;

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
  cout << color << "GenJet : pt() : eta() : y() : phi()" << Color::FG_DEFAULT << endl;
  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;

  for (auto & itr: gps) {
    cout << color;
    if (label != "") {
        cout << "[" << label << "]: ";
    }
    cout << itr.pt() << " : " << itr.eta() << " : " << itr.Rapidity() << " : " << itr.phi() << Color::FG_DEFAULT << endl;
  }

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
}

void uhh2examples::printGenJetsWithParts(const vector<GenJet> & gps, const vector<GenParticle>* genparticles, const string & label, Color::Code color) {
  if (gps.size() == 0) return;

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
  cout << color << "GenJet : pt() : eta() : y() : phi() : genparticles_indices().size()" << Color::FG_DEFAULT << endl;
  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;

  for (auto & itr: gps) {
    cout << color;
    if (label != "") {
        cout << "[" << label << "]: ";
    }
    cout << itr.pt() << " : " << itr.eta() << " : " << itr.Rapidity() << " : " << itr.phi() << " : " << itr.genparticles_indices().size() << Color::FG_DEFAULT << endl;
    uhh2examples::print_genjet_genparticles(itr, genparticles);  // requires the print* funcs to have namespace explicitly written
  }

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
}

void uhh2examples::printJets(const vector<Jet> & jets, const string & label, Color::Code color) {
  if (jets.size() == 0) return;

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
  cout << color << "Jet : pt() : pt() (no JEC) : eta() : y() : phi()" << Color::FG_DEFAULT << endl;
  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;

  for (auto & itr: jets) {
    cout << color;
    if (label != "") {
        cout << "[" << label << "]: ";
    }
    cout << itr.pt() << " : " << itr.pt() * itr.JEC_factor_raw() << " : " << itr.eta() << " : " << itr.Rapidity() << " : " << itr.phi() << " : GenJet: " << itr.genjet_index() << Color::FG_DEFAULT << endl;
  }

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
}


void uhh2examples::printJetsWithParts(const vector<Jet> & jets, const vector<PFParticle> * pfparticles, const string & label, Color::Code color) {
  if (jets.size() == 0) return;

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
  cout << color << "Jet : pt() : pt() (no JEC) : eta() : y() : phi() : numberOfDaughters()" << Color::FG_DEFAULT << endl;
  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;

  for (auto & itr: jets) {
    cout << color;
    if (label != "") {
        cout << "[" << label << "]: ";
    }
    cout << itr.pt() << " : " << itr.pt() * itr.JEC_factor_raw() << " : " << itr.eta() << " : " << itr.Rapidity() << " : " << itr.phi() << " : " << itr.numberOfDaughters() << " : GenJet: " << itr.genjet_index() << Color::FG_DEFAULT << endl;
    uhh2examples::print_jet_pfparticles(itr, pfparticles);
  }

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
}

void uhh2examples::printMuons(const vector<Muon> & muons, const string & label, Color::Code color) {
  if (muons.size() == 0) return;

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
  cout << color << "Muon : pt() : y() : phi() : relIso() : selection bools" << Color::FG_DEFAULT << endl;
  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;

  for (auto & itr: muons) {
    cout << color;
    if (label != "") {
        cout << "[" << label << "]: ";
    }
    // create string with selection bools
    string selection = "";
    for (int i=0; i<=Muon::TkIsoTight; i++) {
      selection += itr.get_selector(static_cast<Muon::Selector>(i)) ? "1" : "0";
    }
    cout << itr.pt() << " : " << itr.Rapidity() << " : " << itr.phi() << " : " << itr.relIso() << " : " << selection << Color::FG_DEFAULT << endl;
  }

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
}

void uhh2examples::printElectrons(const vector<Electron> & electrons, const string & label, Color::Code color) {
  if (electrons.size() == 0) return;

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
  cout << color << "Electron : pt() : y() : phi() : selection bools" << Color::FG_DEFAULT << endl;
  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;

  for (auto & itr: electrons) {
    cout << color;
    if (label != "") {
        cout << "[" << label << "]: ";
    }
    // create string with selection bools
    string selection = "";
    for (int i=0; i<=Electron::mvaEleID_Fall17_iso_V2_wpLoose; i++) {
      auto tag = static_cast<Electron::tag>(i);
      if (!itr.has_tag(tag)) {
        selection += "0";
        continue;
      }
      selection += (itr.get_tag(tag) == 1) ? "1" : "0";
    }
    cout << itr.pt() << " : " << itr.Rapidity() << " : " << itr.phi() << " : " << selection << Color::FG_DEFAULT << endl;
  }

  cout << color << "--------------------------------------------------------------------------------" << Color::FG_DEFAULT << endl;
}
#pragma once

#include <limits>       // std::numeric_limits

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/NtupleObjects.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MCWeight.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#pragma GCC diagnostic pop

#include "TFile.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TH3.h"


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
  ELECTRON_NEUTRINO = 12,
  MUON = 13,
  MUON_NEUTRINO = 14,
  TAU = 15,
  TAU_NEUTRINO = 16,
  GLUON = 21,
  PHOTON = 22,
  Z = 23,
  W = 24
};


/**
 * Easy way to do JEC, JER, MET corrections in Data. Takes care of run period dependency.
 */
class DataJetMetCorrector : public uhh2::AnalysisModule {
public:
  explicit DataJetMetCorrector(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<JetCorrector> jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_GH;
  const int runnr_BCD = 276811;
  const int runnr_EFearly = 278801;
  const int runnr_FlateG = 280385;
  // G&H last set
  // const int runnr_FlateG = 280385;
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


/**
 * Standardised way of cleaning objects, applying IDs, corrections, etc
 * But don't do jets, they are done separately
 */
class GeneralEventSetup : public uhh2::AnalysisModule {
public:
  explicit GeneralEventSetup(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<Selection> lumi_selection;
  std::unique_ptr<AndSelection> metfilters_selection;
  std::unique_ptr<PrimaryVertexCleaner> pv_cleaner;
  std::unique_ptr<ElectronCleaner> electron_cleaner;
  std::unique_ptr<MuonCleaner> muon_cleaner;
};


/**
 * RecoJet-specific bits of setup e.g. JEC, JER
 */
class RecoJetSetup : public uhh2::AnalysisModule {
public:
  RecoJetSetup(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, float jet_radius, float jet_pt_min);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<JetCleaner> jet_pf_id;
  std::unique_ptr<AnalysisModule> jet_met_corrector;
  std::unique_ptr<JetCleaner> jet_cleaner;
  std::unique_ptr<JetElectronOverlapRemoval> jet_ele_cleaner;
  std::unique_ptr<JetMuonOverlapRemoval> jet_mu_cleaner;
};


/**
 * Find generator dilepton invariant object - could be Z, or gamma
 */
class GenZllFinder : public uhh2::AnalysisModule {
public:
  GenZllFinder(uhh2::Context & ctx, bool onlyStatus23, const std::string & genZhandle="genZ", const std::string & genZLeptonhandle="genZleptons");
  virtual bool process(uhh2::Event & event) override;
private:
  bool onlyStatus23_;
  uhh2::Event::Handle<GenParticle> gen_z_handle;
  uhh2::Event::Handle<std::vector<GenParticle>> gen_z_leptons_handle;
};


/**
 * Get k factor for Z(ll)+jets
 */
class ZllKFactor {
public:
  ZllKFactor(const std::string & weightFilename_);
  virtual float getKFactor(float zPt);
private:
  std::unique_ptr<TFile> file;
  std::unique_ptr<TGraph> grNNLO;

};

/**
 * Analysis module to actually apply k factor weight
 */
class ZkFactorReweight : public uhh2::AnalysisModule {
public:
  ZkFactorReweight(uhh2::Context & ctx, const std::string & weightFilename_="", const std::string & genMuonName="");
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<double> z_weight_handle;
  uhh2::Event::Handle<std::vector<GenParticle>> gen_muon_handle;
  std::unique_ptr<ZllKFactor> zReweight;
};


/**
 * Analysis module to reweight using pt spectrum, e.g. to make Herwig++ match Pythia8
 */
class PtReweight : public uhh2::AnalysisModule {
public:
  PtReweight(uhh2::Context & ctx, const std::string & genjets_name, const std::string & weightFilename_="", const std::string & region_="");
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<TH1F> reweightHist;
  std::string method;
  Event::Handle<std::vector<GenJetWithParts>> genjets_handle;
};


/**
 * Apply track SFs to MC, with/without uncertianties. Adds collections of promoted & dropped PF particles to the event
 */
class MCTrackScaleFactor : public uhh2::AnalysisModule {
public:
  MCTrackScaleFactor(uhh2::Context & ctx, const std::string & direction="central");
  virtual bool process(uhh2::Event & event) override;
  PFParticle::EParticleID pdgid_to_pfid(int pdgid, int charge);
private:
  const std::vector<float> eta_regions;
  const std::vector<std::vector<float>> SF;
  const std::vector<std::vector<float>> SF_uncert;
  float run_BtoF_lumi, run_GtoH_lumi, total_lumi;
  const std::string direction_;
  float drMax_;
  std::unique_ptr<TRandom3> random_;
  Event::Handle<std::vector<PFParticle>> dropped_pf_handle;
  Event::Handle<std::vector<PFParticle>> promoted_pf_handle;
  std::map<int, TH3D*> matching_pf_hists;

};


/**
 * Update event pfparticles and jets with dropped/promoted particles from MCTrackScaleFactor
 */
class JetPFUpdater : public uhh2::AnalysisModule {
public:
  JetPFUpdater(uhh2::Context & ctx, const std::string & jet_coll_name="jets");
  virtual bool process(uhh2::Event & event) override;
private:
  Event::Handle<std::vector<Jet>> jet_handle;
  Event::Handle<std::vector<PFParticle>> dropped_pf_handle;
  Event::Handle<std::vector<PFParticle>> promoted_pf_handle;
};


/**
 * Correct for tracking efficiency, adding/removing PF particles & updating reco jets
 */
class TrackingEfficiency : public uhh2::AnalysisModule {
public:
  explicit TrackingEfficiency(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<AnalysisModule> track_sf, jet_updater;
};


/**
 * Apply weights/SF specifically for MC
 */
class MCReweighting : public uhh2::AnalysisModule {
public:
  MCReweighting(uhh2::Context & ctx, const std::string & genjet_name="", const std::string & genmuon_name="");
  virtual bool process(uhh2::Event & event) override;
private:
  Event::Handle<double> gen_weight_handle;
  std::unique_ptr<MCLumiWeight> lumi_weighter;
  std::unique_ptr<MCPileupReweight> pileup_reweighter;
  std::unique_ptr<MCMuonScaleFactor> muon_id_reweighter_pt_eta, muon_id_reweighter_vtx, muon_trg_reweighter;
  std::unique_ptr<MCMuonTrkScaleFactor> muon_trk_reweighter;
  std::unique_ptr<ZkFactorReweight> z_reweighter;
  std::unique_ptr<PtReweight> pt_reweighter;
  std::unique_ptr<MCScaleVariation> mc_scalevar;
  bool doMuons, is_DY;
};


/**
 * find the Z->mumu
 */
class ZFinder : public uhh2::AnalysisModule {
public:
  ZFinder(uhh2::Context & ctx, const std::string & inputLabel_, const std::string & outputLabel_);
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<std::vector<Muon>> hndlInput, hndlZ;
};

/**
 * reweight event based on jet pt
 * useful for e.g. reweighting herwig to match pythia
 */
// class PtReweight : public uhh2::AnalysisModule {
// public:
//   PtReweight(uhh2::Context & ctx, const std::string & selection, const std::string & weightFilename);
//   bool process(uhh2::Event & event, float value);
// private:
//   std::unique_ptr<TFile> f_weight;
//   std::unique_ptr<TH1F> reweightHist;
//   uhh2::Event::Handle<double> gen_weight_handle;
// };

/**
 * Calculate any LesHouches variable from jet constituents
 *
 * Note that we explicity ask for the constituents, and don't just take a jet
 * as an argument. This is so we can customise the constituents separately
 * (e.g. shift constituents' energy, rejecting some).
 *
 * We also have 2 jet axes: the standard axis (jet_vector),
 * and the WTA axis (wta_vector).
 */
template <class T> class LambdaCalculator {
public:
  LambdaCalculator(std::vector<T> & constits, float jet_radius, const LorentzVector & jet_vector, const LorentzVector & wta_vector, bool usePuppiWeight);
  float getLambda(float kappa, float beta);
  void clearCache();
  std::vector<T> constits() { return constits_; }
  float getPtSum() { return ptSum_; }
private:
  std::vector<T> constits_;
  float jetRadius_, ptSum_;
  LorentzVector jetVector_;
  LorentzVector wtaVector_;
  bool usePuppiWeight_;
  std::map<std::pair<float, float>, float> resultsCache_;
};

template class LambdaCalculator<PFParticle>;
template class LambdaCalculator<GenParticle>;

/**
 *
 * Structure to hold Jet + lambda variables
 * Holds:
 * - ungroomed, charged+neutral constits
 * - ungroomed, charged-only constits
 * - groomed, charged+neutral constits
 * - groomed, charged-only constits
 */
struct JetLambdaBundle {
  Jet jet;
  LambdaCalculator<PFParticle> lambda;
  LambdaCalculator<PFParticle> chargedLambda;
  LambdaCalculator<PFParticle> groomedLambda;
  LambdaCalculator<PFParticle> groomedChargedLambda;

  LambdaCalculator<PFParticle> getLambdaCalculator(const bool isCharged, const bool isGroomed) const {
    if (isCharged) {
      if (isGroomed) {
        return groomedChargedLambda;
      } else {
        return chargedLambda;
      }
    } else {
      if (isGroomed) {
        return groomedLambda;
      } else {
        return lambda;
      }
    }
  }
};

/**
 * Structure to hold Jet + lambda variables
 * Holds:
 * - ungroomed, charged+neutral constits
 * - ungroomed, charged-only constits
 * - groomed, charged+neutral constits
 * - groomed, charged-only constits
 *
 * It's easier to create & store all these LambdaCalculators in once go,
 * because then we only need to run fastjet once per jet for e.g. WTA axis
 */
struct GenJetLambdaBundle {
  GenJetWithParts jet; // original AKx jet
  LambdaCalculator<GenParticle> lambda;
  LambdaCalculator<GenParticle> chargedLambda;
  LambdaCalculator<GenParticle> groomedLambda;
  LambdaCalculator<GenParticle> groomedChargedLambda;

  LambdaCalculator<GenParticle> getLambdaCalculator(const bool isCharged, const bool isGroomed) const {
    if (isCharged) {
      if (isGroomed) {
        return groomedChargedLambda;
      } else {
        return chargedLambda;
      }
    } else {
      if (isGroomed) {
        return groomedLambda;
      } else {
        return lambda;
      }
    }
  }

};


float get_jet_radius(const std::string & jet_cone);

float calcGenHT(const std::vector<GenParticle> & gps);

float calcJetKt(const std::vector<GenParticle> & genparticles);

/**
* Class to create Jet, LambdaCalculator structs
*/
class QGAnalysisJetLambda : public uhh2::AnalysisModule {
public:
  QGAnalysisJetLambda(uhh2::Context & ctx,
                      float jetRadius,
                      int nJetsMax,
                      bool doPuppi,
                      const PFParticleId & pfId, // applied to all lambdas
                      const std::string & jet_coll_name,
                      const std::string & output_coll_name);
  bool process(uhh2::Event & event);
  fastjet::PseudoJet convert_uhh_pfparticle_to_pseudojet(const PFParticle & particle, bool applyPuppiWeight);
  std::vector<PFParticle> get_jet_pfparticles(const Jet & jet, uhh2::Event & event, bool applyPuppiWeight);

  void set_charged_hadron_shift(int direction, float rel_shift);
  void shift_charged_hadron_pfparticles(std::vector<PFParticle> & pfparticles, float shift);

  void set_neutral_hadron_shift(int direction, float rel_shift);
  void shift_neutral_hadron_pfparticles(std::vector<PFParticle> & pfparticles, float shift);

  void set_photon_shift(int direction, float rel_shift);
  void shift_photon_pfparticles(std::vector<PFParticle> & pfparticles, float shift);

private:
  fastjet::Recluster wta_cluster_;
  fastjet::Recluster ca_cluster_;
  fastjet::contrib::ModifiedMassDropTagger mmdt_;
  float jetRadius_;
  int nJetsMax_;
  bool doPuppi_;
  PFParticleId pfId_;
  float chargedHadronShift_, neutralHadronShift_, photonShift_;
  uhh2::Event::Handle<std::vector<Jet>> jet_handle_;
  uhh2::Event::Handle<std::vector<JetLambdaBundle>> output_handle_;
};

/**
 * Class to create GenJet, LambdaCalculator structs
 */
class QGAnalysisGenJetLambda : public uhh2::AnalysisModule {
public:
  QGAnalysisGenJetLambda(uhh2::Context & ctx,
                         float jetRadius,
                         int nJetsMax,
                         const GenParticleId & genId, // applied to all lambdas
                         const std::string & jet_coll_name,
                         const std::string & output_coll_name);
  bool process(uhh2::Event & event);
  fastjet::PseudoJet convert_uhh_genparticle_to_pseudojet(const GenParticle & particle);
  std::vector<GenParticle> get_jet_genparticles(const GenJetWithParts & genjet, uhh2::Event & event);

private:
  fastjet::Recluster wta_cluster_;
  fastjet::Recluster ca_cluster_;
  fastjet::contrib::ModifiedMassDropTagger mmdt_;
  float jetRadius_;
  int nJetsMax_;
  GenParticleId genId_;
  float neutralHadronShift_, photonShift_;
  uhh2::Event::Handle<std::vector<GenJetWithParts>> genjet_handle_;
  uhh2::Event::Handle<std::vector<GenJetLambdaBundle>> output_handle_;
};


/**
 * Check if object is charged
 */
class ChargedCut {
public:
  ChargedCut() {}
  bool operator()(const Particle & p, const uhh2::Event & ) const {
    return p.charge() != 0;
  }

};


/**
 * Do reco-genjet matching
 */
class JetMatcher : public uhh2::AnalysisModule {
public:
  JetMatcher(uhh2::Context & ctx, const std::string & jet_coll_name, const std::string & genjet_coll_name, float matchRadius, bool uniqueMatch);
  bool process(uhh2::Event & event);
private:
  uhh2::Event::Handle<std::vector<Jet>> recojet_handle_;
  uhh2::Event::Handle<std::vector<GenJetWithParts>> genjet_handle_;
  float matchRadius_;
  bool uniqueMatch_;
};



/**
 * Check object isn't a neutrino
 */
class NoNeutrinoCut {
public:
  NoNeutrinoCut() {}
  bool operator()(const GenParticle & p, const uhh2::Event & ) const {
    int pid = abs(p.pdgId());
    return ((pid != 12) && (pid != 14) && (pid != 16));
  }
};


/**
 * Check object is status = 1
 */
class Status1Cut {
public:
  Status1Cut() {}
  bool operator()(const GenParticle & p, const uhh2::Event & ) const {
    return (p.status() == 1);
  }
};

/**
 * Check object is final state: has status = 1 and has no daughters
 */
class FinalStateCut {
public:
  FinalStateCut() {}
  bool operator()(const GenParticle & p, const uhh2::Event & event) const {
    return (Status1Cut()(p, event) &&
            // ensure no daughters ie final state
            // note that daughter1/2 are unsigned shorts, so may be 65535!
            (p.daughter1() == -1 || p.daughter1() == std::numeric_limits<unsigned short>::max()) &&
            (p.daughter2() == -1 || p.daughter2() == std::numeric_limits<unsigned short>::max()));
  }
};


/**
 * Cluster GenJets with anti-kT
 */
class GenJetClusterer : public uhh2::AnalysisModule {
public:
  GenJetClusterer(uhh2::Context & ctx, const std::string & genjet_coll_name, float radius, const std::string & genparticles_coll_name="genparticles");
  virtual bool process(uhh2::Event & event) override;
  fastjet::PseudoJet convert_uhh_genparticle_to_pseudojet(const GenParticle & particle);
  GenJetWithParts convert_pseudojet_to_uhh_genjet(const fastjet::PseudoJet & jet);
  inline bool floatMatch(float a, float b) { return (fabs(a-b) < (1E-6 * std::max(a, b))); }
private:
  uhh2::Event::Handle<std::vector<GenJetWithParts>> genjet_handle_;
  uhh2::Event::Handle<std::vector<GenParticle>> genparticle_handle_;
  fastjet::JetDefinition jet_def_;
  GenParticleId gpId_;
};


namespace Binning {
  // For unfolding, the "coarse" binning is the binning used for the final unfolded plots
  // The "fine" binning is each of the coarse bins divided by 2, which is better for TUnfold

  std::vector<double> calculate_fine_binning(const std::vector<double> & coarse_bin_edges);
  std::vector<double> sum_vectors(const std::vector<double> & vec1, const std::vector<double> & vec2);

  // pt bins
  // ---------
  const std::vector<double> pt_bin_edges_gen = {
    50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000 // maximum should be 13 TeV / 2, but get weird empty reco binning
  };
  const int nbins_pt_gen(pt_bin_edges_gen.size() - 1);

  // separate "underflow" bin edges for underflow binning region
  // see comments in QGAnalysisHists.cxx
  const std::vector<double> pt_bin_edges_gen_underflow = {
    15, 30, 38, 50
  };
  const int nbins_pt_gen_underflow(pt_bin_edges_gen_underflow.size() - 1);

  // calculate reco (fine binned) bins using gen bins and making half as wide
  const std::vector<double> pt_bin_edges_reco = calculate_fine_binning(pt_bin_edges_gen); // MUST be const otherwise get mulitple defintion issues
  const int nbins_pt_reco(pt_bin_edges_reco.size() - 1);

  const std::vector<double> pt_bin_edges_reco_underflow = calculate_fine_binning(pt_bin_edges_gen_underflow);
  const int nbins_pt_reco_underflow(pt_bin_edges_reco_underflow.size() - 1);

  // for non-TUnfold things, need all bins together, plus extra underflow bins
  // also need to remove duplicate bin at start of signal region,
  // hence the inplace ctor with begin()+1
  const std::vector<double> pt_bin_edges_gen_all = sum_vectors({0., 15.},
                                                               sum_vectors(pt_bin_edges_gen_underflow,
                                                                           std::vector<double>(pt_bin_edges_gen.begin()+1, pt_bin_edges_gen.end()) )
                                                              );
  const int nbins_pt_gen_all(pt_bin_edges_gen_all.size() - 1);

  const std::vector<double> pt_bin_edges_reco_all = sum_vectors({0., 15.},
                                                                sum_vectors(pt_bin_edges_reco_underflow,
                                                                            std::vector<double>(pt_bin_edges_reco.begin()+1, pt_bin_edges_reco.end()) )
                                                               );
  const int nbins_pt_reco_all(pt_bin_edges_reco_all.size() - 1);

  // Separate pt binning for Z+jets
  // ------------------------------
  // lower last big bin for Z+jets - dont want many empty bins for tunfold
  const std::vector<double> pt_bin_edges_zpj_gen = {
    50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 2000 // maximum should be 13 TeV / 2, but no events in the last reco bin once split into 2, so go for 2TeV
  };
  const int nbins_pt_zpj_gen(pt_bin_edges_zpj_gen.size() - 1);

  const std::vector<double> pt_bin_edges_zpj_gen_underflow = {
    15, 30, 38, 50
  };
  const int nbins_pt_zpj_gen_underflow(pt_bin_edges_zpj_gen_underflow.size() - 1);

  const std::vector<double> pt_bin_edges_zpj_reco = calculate_fine_binning(pt_bin_edges_zpj_gen);
  const int nbins_pt_zpj_reco(pt_bin_edges_zpj_reco.size() - 1);

  const std::vector<double> pt_bin_edges_zpj_reco_underflow = calculate_fine_binning(pt_bin_edges_zpj_gen_underflow);
  const int nbins_pt_zpj_reco_underflow(pt_bin_edges_zpj_reco_underflow.size() - 1);

  const std::vector<double> pt_bin_edges_zpj_gen_all = sum_vectors({0., 15.},
                                                                   sum_vectors(pt_bin_edges_zpj_gen_underflow,
                                                                               std::vector<double>(pt_bin_edges_zpj_gen.begin()+1, pt_bin_edges_zpj_gen.end()) )
                                                                  );
  const int nbins_pt_zpj_gen_all(pt_bin_edges_zpj_gen_all.size() - 1);

  const std::vector<double> pt_bin_edges_zpj_reco_all = sum_vectors({0., 15.},
                                                                    sum_vectors(pt_bin_edges_zpj_reco_underflow,
                                                                                std::vector<double>(pt_bin_edges_zpj_reco.begin()+1, pt_bin_edges_zpj_reco.end()) )
                                                                   );
  const int nbins_pt_zpj_reco_all(pt_bin_edges_zpj_reco_all.size() - 1);


  // LHA binning
  // -----------
  const std::vector<double> lha_bin_edges_gen = {
    // 0.0, 0.16, 0.23, 0.3, 0.36, 0.42, 0.49, 0.56, 0.63, 0.72, 1.0 // target 0.5, across both dijets
    // 0.0, 0.19, 0.29, 0.38, 0.49, 0.6, 0.75, 1.0 // target 0.6, forward
    // 0.0, 0.14, 0.22, 0.29, 0.35, 0.42, 0.49, 0.56, 0.64, 1.0 //target 0.5, cen+fwd
    0.0, 0.14, 0.21, 0.28, 0.34, 0.4, 0.47, 0.54, 0.61, 1.0 //target 0.5, cen+fwd, AK axis
  };
  const int nbins_lha_gen(lha_bin_edges_gen.size() - 1);

  const std::vector<double> lha_bin_edges_reco= calculate_fine_binning(lha_bin_edges_gen);
  const int nbins_lha_reco(lha_bin_edges_reco.size() - 1);

  // Charged LHA binning
  // -----------
  const std::vector<double> lha_charged_bin_edges_gen = {
    // 0.0, 0.09, 0.14, 0.19, 0.24, 0.3, 0.36, 0.43, 0.51, 0.6, 0.71, 1.0 // target 0.5, across both dijets
    // 0.0, 0.11, 0.17, 0.24, 0.32, 0.42, 0.54, 0.69, 1.0 // target 0.6, forward
    // 0.0, 0.05, 0.09, 0.12, 0.15, 0.18, 0.22, 0.26, 0.31, 0.36, 0.42, 0.49, 0.56, 0.64, 1.0 // target 0.5, cen+fwd
    0.0, 0.05, 0.09, 0.14, 0.19, 0.25, 0.32, 0.39, 0.47, 0.55, 0.63, 0.73, 1.0 // target 0.5, cen+fwd groomed (as coarser), AK axis
  };
  const int nbins_lha_charged_gen(lha_charged_bin_edges_gen.size() - 1);

  const std::vector<double> lha_charged_bin_edges_reco= calculate_fine_binning(lha_charged_bin_edges_gen);
  const int nbins_lha_charged_reco(lha_charged_bin_edges_reco.size() - 1);

  // Multiplicity binning
  // --------------------
  const std::vector<double> puppiMultiplicity_bin_edges_gen = {
    // 0, 9, 15, 22, 35, 50, 75, 100, 150 // target 0.5, across both dijets
    // 0.0, 11.0, 18.0, 25.0, 150.0 // target 0.6, central
    // 0.0, 10, 15, 20, 27, 50, 75, 100, 150 // target 0.5, cen+fwd
    0.0, 10, 15, 20, 30, 50, 75, 100, 150 // target 0.5, cen+fwd, AK axis. 30, 50, 75, 100 added by hand otherwise no granulairty
  };
  const int nbins_puppiMultiplicity_gen(puppiMultiplicity_bin_edges_gen.size() - 1);

  const std::vector<double> puppiMultiplicity_bin_edges_reco = calculate_fine_binning(puppiMultiplicity_bin_edges_gen);
  const int nbins_puppiMultiplicity_reco(puppiMultiplicity_bin_edges_reco.size() - 1);

  // Charged Multiplicity binning
  // --------------------
  const std::vector<double> puppiMultiplicity_charged_bin_edges_gen = {
    // 0, 9, 15, 22, 35, 50, 75, 100, 150 // target 0.5, across both dijets
    // 0.0, 4.0, 6.0, 9.0, 12.0, 16.0, 23.0, 150.0 // target 0.6, forward
    // 0.0, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 18.0, 21.0, 25.0, 32.0, 91, 150.0 // target 0.5, cen+fwd, 91 added as halfway between 32 and 150
   0.0, 3.0, 5.0, 7.0, 10.0, 13.0, 16.0, 20.0, 30, 50, 75, 100, 150.0 // target 0.5, cen+fwd groomed, AK axis, added 30, 50, 75, 100 manually
  };
  const int nbins_puppiMultiplicity_charged_gen(puppiMultiplicity_charged_bin_edges_gen.size() - 1);

  const std::vector<double> puppiMultiplicity_charged_bin_edges_reco = calculate_fine_binning(puppiMultiplicity_charged_bin_edges_gen);
  const int nbins_puppiMultiplicity_charged_reco(puppiMultiplicity_charged_bin_edges_reco.size() - 1);

  // pTD binning
  // --------------------
  const std::vector<double> pTD_bin_edges_gen = {
    // 0.0, 0.09, 0.14, 0.25, 1.0 // target 0.5, across both dijets
    // 0.0, 0.1, 0.17, 0.4, 1.0 // target 0.6, forward
    // 0.0, 0.07, 0.1, 0.15, 0.24, 0.45, 1.0 // target 0.5, cen+fwd
    0.0, 0.07, 0.1, 0.15, 0.24, 1.0 // target 0.5, cen+fwd, AK axis
  };
  const int nbins_pTD_gen(pTD_bin_edges_gen.size() - 1);

  const std::vector<double> pTD_bin_edges_reco = calculate_fine_binning(pTD_bin_edges_gen);
  const int nbins_pTD_reco(pTD_bin_edges_reco.size() - 1);

  // Charged pTD binning
  // --------------------
  const std::vector<double> pTD_charged_bin_edges_gen = {
    // 0.0, 0.09, 0.12, 0.15, 0.19, 0.24, 0.31, 0.4, 0.53, 0.73, 1.0 // target 0.5, across both dijets
    // 0.0, 0.09, 0.12, 0.16, 0.21, 0.29, 0.41, 0.6, 1.0 // target 0.6, forward
    // 0.0, 0.08, 0.1, 0.12, 0.14, 0.17, 0.2, 0.24, 0.29, 0.34, 0.4, 0.47, 0.55, 0.65, 0.76, 0.89, 1.0 // target 0.5, cen+fwd
    0.0, 0.07, 0.09, 0.11, 0.14, 0.18, 0.23, 0.3, 0.39, 0.51, 0.64, 1.0 // target 0.5, cen+fwd groomed (as ungroomed had too large drop in purity/stab), AK axis
  };
  const int nbins_pTD_charged_gen(pTD_charged_bin_edges_gen.size() - 1);

  const std::vector<double> pTD_charged_bin_edges_reco = calculate_fine_binning(pTD_charged_bin_edges_gen);
  const int nbins_pTD_charged_reco(pTD_charged_bin_edges_reco.size() - 1);

  // thrust binning
  // --------------------
  const std::vector<double> thrust_bin_edges_gen = {
    // 0.0, 0.04, 0.07, 0.12, 0.19, 0.27, 0.37, 0.52, 1.0 // target 0.5, across both dijets
    // 0.0, 0.06, 0.16, 0.3, 0.51, 1.0 // target 0.6, forward
    // 0.0, 0.04, 0.08, 0.145, 0.225, 0.32, 0.445, 1.0 // target 0.5, cen+fwd
    0.0, 0.04, 0.08, 0.14, 0.195, 0.25, 1.0 // target 0.5, cen+fwd, AK axis
  };
  const int nbins_thrust_gen(thrust_bin_edges_gen.size() - 1);

  const std::vector<double> thrust_bin_edges_reco = calculate_fine_binning(thrust_bin_edges_gen);
  const int nbins_thrust_reco(thrust_bin_edges_reco.size() - 1);

  // Charged thrust binning
  // --------------------
  const std::vector<double> thrust_charged_bin_edges_gen = {
    // 0.0, 0.02, 0.04, 0.07, 0.11, 0.17, 0.25, 0.36, 0.52, 1.0 // target 0.5, across both dijets
    // 0.0, 0.01, 0.02, 0.04, 0.07, 0.11, 0.18, 0.29, 0.48, 1.0 // target 0.6, forward
    // 0.0, 0.005, 0.01, 0.015, 0.025, 0.035, 0.05, 0.065, 0.085, 0.115, 0.15, 0.2, 0.265, 0.355, 0.48, 0.725, 1.0 // target 0.5, cen+fwd
    0.0, 0.005, 0.015, 0.03, 0.06, 0.105, 0.16, 0.23, 0.315, 0.425, 1.0 // target 0.5, cen+fwd groomed, AK axis
  };
  const int nbins_thrust_charged_gen(thrust_charged_bin_edges_gen.size() - 1);

  const std::vector<double> thrust_charged_bin_edges_reco = calculate_fine_binning(thrust_charged_bin_edges_gen);
  const int nbins_thrust_charged_reco(thrust_charged_bin_edges_reco.size() - 1);

  // width binning
  // --------------------
  const std::vector<double> width_bin_edges_gen = {
    // 0.0, 0.08, 0.13, 0.18, 0.24, 0.31, 0.39, 0.47, 0.58, 1.0 // target 0.5, across both dijets
    // 0.0, 0.12, 0.2, 0.32, 0.45, 0.64, 1.0 // target 0.6, forward
    // 0.0, 0.09, 0.145, 0.205, 0.28, 0.36, 0.445, 0.545, 1.0 // target 0.5, cen+fwd
    0.0, 0.09, 0.145, 0.205, 0.28, 0.355, 0.43, 0.515, 1.0 // target 0.5, cen+fwd, AK axis
  };
  const int nbins_width_gen(width_bin_edges_gen.size() - 1);

  const std::vector<double> width_bin_edges_reco = calculate_fine_binning(width_bin_edges_gen);
  const int nbins_width_reco(width_bin_edges_reco.size() - 1);

  // Charged width binning
  // --------------------
  const std::vector<double> width_charged_bin_edges_gen = {
    // 0.0, 0.04, 0.07, 0.1, 0.14, 0.18, 0.23, 0.29, 0.37, 0.46, 0.58, 0.99, 1.0 // target 0.5, across both dijets
    // 0.0, 0.04, 0.07, 0.11, 0.16, 0.22, 0.3, 0.41, 0.57, 1.0 // target 0.6, forward
    // 0.0, 0.015, 0.025, 0.035, 0.05, 0.065, 0.085, 0.105, 0.13, 0.16, 0.195, 0.235, 0.28, 0.335, 0.4, 0.48, 0.585, 1.0 // target 0.5, cen+fwd
    0.0, 0.01, 0.025, 0.045, 0.075, 0.11, 0.155, 0.21, 0.275, 0.34, 0.415, 0.495, 0.595, 1.0 // target 0.5, cen+fwd groomed, AK axis
  };
  const int nbins_width_charged_gen(width_charged_bin_edges_gen.size() - 1);

  const std::vector<double> width_charged_bin_edges_reco = calculate_fine_binning(width_charged_bin_edges_gen);
  const int nbins_width_charged_reco(width_charged_bin_edges_reco.size() - 1);
}

#pragma once

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

#include "TFile.h"
#include "TGraph.h"


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
 */
class GeneralEventSetup : public uhh2::AnalysisModule {
public:
  GeneralEventSetup(uhh2::Context & ctx, const std::string & pu_removal, const std::string & jet_cone, float jet_radius, float jet_pt_min=30);
  virtual bool process(uhh2::Event & event) override;
private:
  bool is_mc;
  std::unique_ptr<Selection> lumi_selection;
  std::unique_ptr<AndSelection> metfilters_selection;
  std::unique_ptr<PrimaryVertexCleaner> pv_cleaner;
  std::unique_ptr<ElectronCleaner> electron_cleaner;
  std::unique_ptr<MuonCleaner> muon_cleaner;
  std::unique_ptr<JetCleaner> jet_pf_id;
  std::unique_ptr<AnalysisModule> jet_met_corrector;
  std::unique_ptr<JetCleaner> jet_cleaner;
  std::unique_ptr<JetElectronOverlapRemoval> jet_ele_cleaner;
  std::unique_ptr<JetMuonOverlapRemoval> jet_mu_cleaner;
};


/**
 * Apply weights/SF specifically for MC
 */
class MCReweighting : public uhh2::AnalysisModule {
public:
  MCReweighting(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  Event::Handle<double> gen_weight_handle;
  std::unique_ptr<MCLumiWeight> lumi_weighter;
  std::unique_ptr<MCPileupReweight> pileup_reweighter;
  std::unique_ptr<MCMuonScaleFactor> muon_id_reweighter_pt_eta, muon_id_reweighter_vtx, muon_trg_reweighter;
  std::unique_ptr<MCMuonTrkScaleFactor> muon_trk_reweighter;
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
 * find the Z->mumu
 */
class ZFinder : public uhh2::AnalysisModule {
public:
  ZFinder(uhh2::Context & ctx, const std::string & inputLabel_, const std::string & outputLabel_, const std::string & weightFilename_="");
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<std::vector<Muon>> hndlInput, hndlZ;
  uhh2::Event::Handle<double> z_weight_handle;
  std::unique_ptr<ZllKFactor> zReweight;
};

/**
 * reweight event based on jet pt
 * useful for e.g. reweighting herwig to match pythia
 */
class PtReweight : public uhh2::AnalysisModule {
public:
  PtReweight(uhh2::Context & ctx, const std::string & selection, const std::string & weightFilename);
  bool process(uhh2::Event & event, float value);
private:
  std::unique_ptr<TFile> f_weight;
  std::unique_ptr<TH1F> reweightHist;
  uhh2::Event::Handle<double> gen_weight_handle;
};

/**
 * Calculate any LesHouches variable from jet constituents
 */
template <class T> class LambdaCalculator {
public:
  LambdaCalculator(std::vector<T> & daughters, float jet_radius, const LorentzVector & jet_vector, bool usePuppiWeight);
  float getLambda(float kappa, float beta);
  void clearCache();
  std::vector<T> daughters() { return daughters_; }
private:
  std::vector<T> daughters_;
  float jetRadius_, ptSum_;
  LorentzVector jetVector_;
  bool usePuppiWeight_;
  std::map<std::pair<float, float>, float> resultsCache_;
};

template class LambdaCalculator<PFParticle>;
template class LambdaCalculator<GenParticle>;

/**
 *
 * Structure to hold Jet + lambda variable
 */
struct JetLambdaBundle {
  Jet jet;
  LambdaCalculator<PFParticle> lambda;
};

/**
 * Structure to hold Jet + lambda variable
 */
struct GenJetLambdaBundle {
  GenJetWithParts jet;
  LambdaCalculator<GenParticle> lambda;
};


float get_jet_radius(const std::string & jet_cone);

float calcGenHT(const std::vector<GenParticle> & gps);


/**
* Class to create Jet, LambdaCalculator structs
*/
class QGAnalysisJetLambda : public uhh2::AnalysisModule {
public:
  QGAnalysisJetLambda(uhh2::Context & ctx,
                      float jetRadius,
                      int nJetsMax,
                      bool doPuppi,
                      const PFParticleId & pfId,
                      const std::string & jet_coll_name="jets",
                      const std::string & output_coll_name="jetlambdas");
  bool process(uhh2::Event & event);
  std::vector<PFParticle> get_jet_pfparticles(const Jet & jet, uhh2::Event & event);

  void set_neutral_hadron_shift(int direction, float rel_shift);
  void shift_neutral_hadron_pfparticles(std::vector<PFParticle> & pfparticles, float shift);

  void set_photon_shift(int direction, float rel_shift);
  void shift_photon_pfparticles(std::vector<PFParticle> & pfparticles, float shift);

private:
  float jetRadius_;
  int nJetsMax_;
  bool doPuppi_;
  PFParticleId pfId_;
  float neutralHadronShift_, photonShift_;
  uhh2::Event::Handle<std::vector<Jet>> jet_handle;
  uhh2::Event::Handle<std::vector<JetLambdaBundle>> output_handle;
};

/**
 * Class to create GenJet, LambdaCalculator structs
 */
class QGAnalysisGenJetLambda : public uhh2::AnalysisModule {
public:
  QGAnalysisGenJetLambda(uhh2::Context & ctx,
                         float jetRadius,
                         int nJetsMax,
                         const GenParticleId & genId,
                         const std::string & jet_coll_name="genjets",
                         const std::string & output_coll_name="genjetlambdas");
  bool process(uhh2::Event & event);
  std::vector<GenParticle> get_jet_genparticles(const GenJetWithParts & genjet, uhh2::Event & event);

  // Is shifting genjet constit energies sensible?
  // void set_neutral_hadron_shift(float direction, float rel_shift);
  // void shift_neutral_hadron_genparticles(std::vector<GenParticle*> genparticles, float shift);

  // void set_photon_shift(float direction, float rel_shift);
  // void shift_photon_genparticles(std::vector<GenParticle*> genparticles, float shift);

private:
  float jetRadius_;
  int nJetsMax_;
  GenParticleId genId_;
  float neutralHadronShift_, photonShift_;
  uhh2::Event::Handle<std::vector<GenJetWithParts>> genjet_handle;
  uhh2::Event::Handle<std::vector<GenJetLambdaBundle>> output_handle;
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


namespace Binning {
  // For unfolding, the "coarse" binning is the binning used for the final unfolded plots
  // The "fine" binning is each of the coarse bins divided by 2, which is better for TUnfold

  std::vector<double> calculate_fine_binning(const std::vector<double> & coarse_bin_edges);
  std::vector<double> sum_vectors(const std::vector<double> & vec1, const std::vector<double> & vec2);

  // pt bins
  // ---------
  const std::vector<double> pt_bin_edges_gen = {
    50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 1500, 2000, 10000
  };
  const int nbins_pt_gen(pt_bin_edges_gen.size() - 1);

  // separate "underflow" bin edges for underflow binning region
  // see comments in QGAnalysisHists.cxx
  const std::vector<double> pt_bin_edges_gen_underflow = {
    30, 38, 50
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
    50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 10000
  };
  const int nbins_pt_zpj_gen(pt_bin_edges_zpj_gen.size() - 1);

  const std::vector<double> pt_bin_edges_zpj_gen_underflow = {
    30, 38, 50
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
    0.0, 0.12, 0.25, 0.33, 0.39, 0.45, 0.51, 0.57, 0.62, 0.66, 0.7, 0.76, 1.0
    // 0.0, 0.35, 0.5, 0.64, 1.0
  };
  const int nbins_lha_gen(lha_bin_edges_gen.size() - 1);

  const std::vector<double> lha_bin_edges_reco= calculate_fine_binning(lha_bin_edges_gen);
  const int nbins_lha_reco(lha_bin_edges_reco.size() - 1);

  // Charged LHA binning
  // -----------
  const std::vector<double> lha_charged_bin_edges_gen = {
    0.0, 0.1, 0.2, 0.26, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.66, 0.7, 0.75, 1.0
  // 0.0, 0.29, 0.41, 0.53, 0.66, 1.0
  };
  const int nbins_lha_charged_gen(lha_charged_bin_edges_gen.size() - 1);

  const std::vector<double> lha_charged_bin_edges_reco= calculate_fine_binning(lha_charged_bin_edges_gen);
  const int nbins_lha_charged_reco(lha_charged_bin_edges_reco.size() - 1);

  // Multiplicity binning
  // --------------------
  const std::vector<double> puppiMultiplicity_bin_edges_gen = {
    0, 9, 15, 22, 35, 50, 75, 100, 150
    // 0, 22, 44, 150
  };
  const int nbins_puppiMultiplicity_gen(puppiMultiplicity_bin_edges_gen.size() - 1);

  const std::vector<double> puppiMultiplicity_bin_edges_reco = calculate_fine_binning(puppiMultiplicity_bin_edges_gen);
  const int nbins_puppiMultiplicity_reco(puppiMultiplicity_bin_edges_reco.size() - 1);

  // Charged Multiplicity binning
  // --------------------
  const std::vector<double> puppiMultiplicity_charged_bin_edges_gen = {
    0, 9, 15, 22, 35, 50, 75, 100, 150
    // 0, 16, 32, 150
  };
  const int nbins_puppiMultiplicity_charged_gen(puppiMultiplicity_charged_bin_edges_gen.size() - 1);

  const std::vector<double> puppiMultiplicity_charged_bin_edges_reco = calculate_fine_binning(puppiMultiplicity_charged_bin_edges_gen);
  const int nbins_puppiMultiplicity_charged_reco(puppiMultiplicity_charged_bin_edges_reco.size() - 1);

  // pTD binning
  // --------------------
  const std::vector<double> pTD_bin_edges_gen = {
    0.0, 0.09, 0.14, 0.25, 1.0
    // 0.0, 0.15, 1.0
  };
  const int nbins_pTD_gen(pTD_bin_edges_gen.size() - 1);

  const std::vector<double> pTD_bin_edges_reco = calculate_fine_binning(pTD_bin_edges_gen);
  const int nbins_pTD_reco(pTD_bin_edges_reco.size() - 1);

  // Charged pTD binning
  // --------------------
  const std::vector<double> pTD_charged_bin_edges_gen = {
    0.0, 0.09, 0.12, 0.15, 0.19, 0.24, 0.31, 0.4, 0.53, 0.73, 1.0
    // 0.0, 0.15, 0.26, 0.6, 1.0
  };
  const int nbins_pTD_charged_gen(pTD_charged_bin_edges_gen.size() - 1);

  const std::vector<double> pTD_charged_bin_edges_reco = calculate_fine_binning(pTD_charged_bin_edges_gen);
  const int nbins_pTD_charged_reco(pTD_charged_bin_edges_reco.size() - 1);

  // thrust binning
  // --------------------
  const std::vector<double> thrust_bin_edges_gen = {
    0.0, 0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.31, 0.47, 1.0
    // 0.0, 0.06, 0.17, 0.31, 1.0
  };
  const int nbins_thrust_gen(thrust_bin_edges_gen.size() - 1);

  const std::vector<double> thrust_bin_edges_reco = calculate_fine_binning(thrust_bin_edges_gen);
  const int nbins_thrust_reco(thrust_bin_edges_reco.size() - 1);

  // Charged thrust binning
  // --------------------
  const std::vector<double> thrust_charged_bin_edges_gen = {
    0.0, 0.02, 0.04, 0.07, 0.11, 0.16, 0.22, 0.28, 0.37, 1.0
    // 0.0, 0.05, 0.15, 0.45, 1.0
  };
  const int nbins_thrust_charged_gen(thrust_charged_bin_edges_gen.size() - 1);

  const std::vector<double> thrust_charged_bin_edges_reco = calculate_fine_binning(thrust_charged_bin_edges_gen);
  const int nbins_thrust_charged_reco(thrust_charged_bin_edges_reco.size() - 1);

  // width binning
  // --------------------
  const std::vector<double> width_bin_edges_gen = {
    0.0, 0.09, 0.14, 0.19, 0.25, 0.31, 0.37, 0.43, 0.48, 0.54, 1.0
  // 0.0, 0.16, 0.3, 0.45, 1.0
  };
  const int nbins_width_gen(width_bin_edges_gen.size() - 1);

  const std::vector<double> width_bin_edges_reco = calculate_fine_binning(width_bin_edges_gen);
  const int nbins_width_reco(width_bin_edges_reco.size() - 1);

  // Charged width binning
  // --------------------
  const std::vector<double> width_charged_bin_edges_gen = {
    0.0, 0.06, 0.09, 0.13, 0.17, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.52, 0.6, 1.0
    // 0.0, 0.11, 0.2, 0.33, 0.49, 1.0
  };
  const int nbins_width_charged_gen(width_charged_bin_edges_gen.size() - 1);

  const std::vector<double> width_charged_bin_edges_reco = calculate_fine_binning(width_charged_bin_edges_gen);
  const int nbins_width_charged_reco(width_charged_bin_edges_reco.size() - 1);
}

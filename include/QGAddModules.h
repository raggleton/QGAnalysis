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
  uhh2::Event::Handle<double> gen_weight_handle;
  std::unique_ptr<uhh2examples::ZllKFactor> zReweight;
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
  LambdaCalculator(std::vector<T*> daughters, float jet_radius, const LorentzVector & jet_vector, bool usePuppiWeight);
  float getLambda(float kappa, float beta);
private:
  std::vector<T*> daughters_;
  float jetRadius_, ptSum_;
  LorentzVector jetVector_;
  bool usePuppiWeight_;
};

template class LambdaCalculator<PFParticle>;
template class LambdaCalculator<GenParticle>;
};

float get_jet_radius(const std::string & jet_cone);

float calcGenHT(const std::vector<GenParticle> & gps);

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
    0.0, 0.14, 0.29, 0.37, 0.44, 0.5, 0.56, 0.62, 0.68, 0.75, 1.0
  };
  const int nbins_lha_gen(lha_bin_edges_gen.size() - 1);

  const std::vector<double> lha_bin_edges_reco= calculate_fine_binning(lha_bin_edges_gen);
  const int nbins_lha_reco(lha_bin_edges_reco.size() - 1);

  // Charged LHA binning
  // -----------
  const std::vector<double> lha_charged_bin_edges_gen = {
    0.0, 0.14, 0.29, 0.37, 0.44, 0.5, 0.56, 0.62, 0.68, 0.75, 1.0
  };
  const int nbins_lha_charged_gen(lha_charged_bin_edges_gen.size() - 1);

  const std::vector<double> lha_charged_bin_edges_reco= calculate_fine_binning(lha_charged_bin_edges_gen);
  const int nbins_lha_charged_reco(lha_charged_bin_edges_reco.size() - 1);

  // Multiplicity binning
  // --------------------
  const std::vector<double> puppiMultiplicity_bin_edges_gen = {
    1, 5, 10, 13, 19, 25, 35, 50, 75, 100, 150
  };
  const int nbins_puppiMultiplicity_gen(puppiMultiplicity_bin_edges_gen.size() - 1);

  const std::vector<double> puppiMultiplicity_bin_edges_reco = calculate_fine_binning(puppiMultiplicity_bin_edges_gen);
  const int nbins_puppiMultiplicity_reco(puppiMultiplicity_bin_edges_reco.size() - 1);

  // Charged Multiplicity binning
  // --------------------
  const std::vector<double> puppiMultiplicity_charged_bin_edges_gen = {
    1, 5, 10, 13, 19, 25, 35, 50, 75, 100, 150
  };
  const int nbins_puppiMultiplicity_charged_gen(puppiMultiplicity_charged_bin_edges_gen.size() - 1);

  const std::vector<double> puppiMultiplicity_charged_bin_edges_reco = calculate_fine_binning(puppiMultiplicity_charged_bin_edges_gen);
  const int nbins_puppiMultiplicity_charged_reco(puppiMultiplicity_charged_bin_edges_reco.size() - 1);

  // pTD binning
  // --------------------
  const std::vector<double> pTD_bin_edges_gen = {
    0.0, 0.09, 0.12, 0.16, 0.21, 0.29, 0.43, 0.7, 1.0
  };
  const int nbins_pTD_gen(pTD_bin_edges_gen.size() - 1);

  const std::vector<double> pTD_bin_edges_reco = calculate_fine_binning(pTD_bin_edges_gen);
  const int nbins_pTD_reco(pTD_bin_edges_reco.size() - 1);

  // Charged pTD binning
  // --------------------
  const std::vector<double> pTD_charged_bin_edges_gen = {
    0.0, 0.09, 0.12, 0.16, 0.21, 0.29, 0.43, 0.7, 1.0
  };
  const int nbins_pTD_charged_gen(pTD_charged_bin_edges_gen.size() - 1);

  const std::vector<double> pTD_charged_bin_edges_reco = calculate_fine_binning(pTD_charged_bin_edges_gen);
  const int nbins_pTD_charged_reco(pTD_charged_bin_edges_reco.size() - 1);

  // thrust binning
  // --------------------
  const std::vector<double> thrust_bin_edges_gen = {
    0.0, 0.04, 0.08, 0.12, 0.17, 0.24, 0.33, 0.66, 1.0
  };
  const int nbins_thrust_gen(thrust_bin_edges_gen.size() - 1);

  const std::vector<double> thrust_bin_edges_reco = calculate_fine_binning(thrust_bin_edges_gen);
  const int nbins_thrust_reco(thrust_bin_edges_reco.size() - 1);

  // Charged thrust binning
  // --------------------
  const std::vector<double> thrust_charged_bin_edges_gen = {
    0.0, 0.04, 0.08, 0.12, 0.17, 0.24, 0.33, 0.66, 1.0
  };
  const int nbins_thrust_charged_gen(thrust_charged_bin_edges_gen.size() - 1);

  const std::vector<double> thrust_charged_bin_edges_reco = calculate_fine_binning(thrust_charged_bin_edges_gen);
  const int nbins_thrust_charged_reco(thrust_charged_bin_edges_reco.size() - 1);

  // width binning
  // --------------------
  const std::vector<double> width_bin_edges_gen = {
    0.0, 0.11, 0.17, 0.23, 0.29, 0.35, 0.42, 0.6, 1.0
  };
  const int nbins_width_gen(width_bin_edges_gen.size() - 1);

  const std::vector<double> width_bin_edges_reco = calculate_fine_binning(width_bin_edges_gen);
  const int nbins_width_reco(width_bin_edges_reco.size() - 1);

  // Charged width binning
  // --------------------
  const std::vector<double> width_charged_bin_edges_gen = {
    0.0, 0.11, 0.17, 0.23, 0.29, 0.35, 0.42, 0.6, 1.0
  };
  const int nbins_width_charged_gen(width_charged_bin_edges_gen.size() - 1);

  const std::vector<double> width_charged_bin_edges_reco = calculate_fine_binning(width_charged_bin_edges_gen);
  const int nbins_width_charged_reco(width_charged_bin_edges_reco.size() - 1);
}

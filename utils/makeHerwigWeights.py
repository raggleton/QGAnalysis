#!/usr/bin/env python

"""Construct weights for Herwig to ensure jet spectra matches Pythia

Need for both Dijet & Z+jets regions, for both Reco & Gen jets
"""

from itertools import product
import ROOT


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


def make_weights_file(pythia_dir, herwig_dir, output_file):
    """Create file with histogram of weights for each selection region.
    The weights reweight the herwig jet pt spectrum to that of pythia.
    """
    fout = ROOT.TFile(output_file, "RECREATE")
    
    # Do QCD sample
    f_pythia = ROOT.TFile("../Selection/"+pythia_dir+"/uhh2.AnalysisModuleRunner.MC.MC_QCD.root")
    f_herwig = ROOT.TFile("../Selection/"+herwig_dir+"/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root")

    if (f_pythia.IsZombie()):
        print "Cannot open Pythia QCD file"
    if (f_herwig.IsZombie()):
        print "Cannot open Herwig QCD file"
    if not (f_pythia.IsZombie() or f_herwig.IsZombie()):
        h_dijet_p = f_pythia.Get("Dijet_gen/pt_dijet_ave").Clone("dijet")
        h_dijet_h = f_herwig.Get("Dijet_gen/pt_dijet_ave")
        h_dijet_p.SetDirectory(0)
        h_dijet_h.SetDirectory(0)
        h_dijet_p.Scale(1./h_dijet_p.Integral())
        h_dijet_h.Scale(1./h_dijet_h.Integral())
        h_dijet_p.Divide(h_dijet_h)

        fout.cd()
        h_dijet_p.Write()

    exit()

    # Do Z+Jets
    f_pythia = ROOT.TFile("../Selection/"+pythia_dir+"/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
    f_herwig = ROOT.TFile("../Selection/"+herwig_dir+"/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
    if (f_pythia.IsZombie()):
        print "Cannot open Pythia DY file"
    if (f_herwig.IsZombie()):
        print "Cannot open Herwig DY file"
    if not (f_pythia.IsZombie() or f_herwig.IsZombie()):
        h_zpj_reco_p = f_pythia.Get("ZPlusJets_QG/jet_pt").Clone("zpj_reco")
        h_zpj_reco_h = f_herwig.Get("ZPlusJets_QG/jet_pt")
        h_zpj_reco_p.SetDirectory(0)
        h_zpj_reco_h.SetDirectory(0)
        # h_zpj_reco_p.Rebin(2)
        # h_zpj_reco_h.Rebin(2)
        h_zpj_reco_p.Scale(1./h_zpj_reco_p.Integral())
        h_zpj_reco_h.Scale(1./h_zpj_reco_h.Integral())
        h_zpj_reco_p.Divide(h_zpj_reco_h)

        # h_zpj_gen_p = f_pythia.Get("ZPlusJets_genjet/genjet_pt").Clone("zpj_gen")
        # h_zpj_gen_h = f_herwig.Get("ZPlusJets_genjet/genjet_pt")
        # h_zpj_gen_p.SetDirectory(0)
        # h_zpj_gen_h.SetDirectory(0)
        # # h_zpj_gen_p.Rebin(2)
        # # h_zpj_gen_h.Rebin(2)
        # h_zpj_gen_p.Scale(1./h_zpj_gen_p.Integral())
        # h_zpj_gen_h.Scale(1./h_zpj_gen_h.Integral())
        # h_zpj_gen_p.Divide(h_zpj_gen_h)
        
        fout.cd()
        h_zpj_reco_p.Write()
        # h_zpj_gen_p.Write()


if __name__ == "__main__":

    ALGOS = ["ak4", "ak8"][0:1]
    PUS = ["chs", "puppi"][1:2]
    for algo, pu, in product(ALGOS, PUS):
        SETUP = algo + pu
        # Check your directories are in the right places
        PYTHIA_DIR = "MGPythia/workdir_%s_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto" % SETUP
        HERWIG_DIR = "Herwig/workdir_%s_herwig_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut" % SETUP
        # POWHEG_DIR = "Powheg/workdir_%s_powheg_newFlav_withPtRatio_0p94_withLeptonOverlapVeto_eta2p4_ssEta_dEta1p2" % SETUP
        # PYTHIAONLY_DIR = "PythiaOnlyFlat/workdir_%s_pythiaOnlyFlat_newFlav_withPtRatio_0p94_withLeptonOverlapVeto_eta2p4_ssEta_dEta1p2" % SETUP
        make_weights_file(PYTHIA_DIR, HERWIG_DIR, "herwig_dijet_only_weights_%s.root" % SETUP)
        # make_weights_file(PYTHIA_DIR, PYTHIAONLY_DIR, "pythiaOnlyFlat_newflav_weights_%s.root" % SETUP)
        # make_weights_file(PYTHIA_DIR, POWHEG_DIR, "powheg_newflav_weights_%s.root" % SETUP)

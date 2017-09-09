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
    The wieghts reweight the herwig jet pt spectrum to that of pythia.
    """
    # Do QCD sample
    f_pythia = ROOT.TFile("../Selection/"+pythia_dir+"/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root")
    f_herwig = ROOT.TFile("../Selection/"+herwig_dir+"/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root")

    if (f_pythia.IsZombie()):
        raise RuntimeError("Cannot open Pythia file")
    if (f_herwig.IsZombie()):
        raise RuntimeError("Cannot open Herwig file")

    h_dijet_reco_p = f_pythia.Get("Dijet/pt_jet1").Clone("dijet_reco")
    h_dijet_reco_h = f_herwig.Get("Dijet/pt_jet1")
    h_dijet_reco_p.SetDirectory(0)
    h_dijet_reco_h.SetDirectory(0)
    h_dijet_reco_p.Scale(1./h_dijet_reco_p.Integral())
    h_dijet_reco_h.Scale(1./h_dijet_reco_h.Integral())
    h_dijet_reco_p.Divide(h_dijet_reco_h)

    h_dijet_gen_p = f_pythia.Get("Dijet_genjet/genjet_pt").Clone("dijet_gen")
    h_dijet_gen_h = f_herwig.Get("Dijet_genjet/genjet_pt")
    h_dijet_gen_p.SetDirectory(0)
    h_dijet_gen_h.SetDirectory(0)
    h_dijet_gen_p.Scale(1./h_dijet_gen_p.Integral())
    h_dijet_gen_h.Scale(1./h_dijet_gen_h.Integral())
    h_dijet_gen_p.Divide(h_dijet_gen_h)

    # Do Z+Jets
    f_pythia = ROOT.TFile("../Selection/"+pythia_dir+"/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root")
    f_herwig = ROOT.TFile("../Selection/"+herwig_dir+"/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root")
    h_zpj_reco_p = f_pythia.Get("ZPlusJets_QG/jet_pt").Clone("zpj_reco")
    h_zpj_reco_h = f_herwig.Get("ZPlusJets_QG/jet_pt")
    h_zpj_reco_p.SetDirectory(0)
    h_zpj_reco_h.SetDirectory(0)
    # h_zpj_reco_p.Rebin(2)
    # h_zpj_reco_h.Rebin(2)
    h_zpj_reco_p.Scale(1./h_zpj_reco_p.Integral())
    h_zpj_reco_h.Scale(1./h_zpj_reco_h.Integral())
    h_zpj_reco_p.Divide(h_zpj_reco_h)

    h_zpj_gen_p = f_pythia.Get("ZPlusJets_genjet/genjet_pt").Clone("zpj_gen")
    h_zpj_gen_h = f_herwig.Get("ZPlusJets_genjet/genjet_pt")
    h_zpj_gen_p.SetDirectory(0)
    h_zpj_gen_h.SetDirectory(0)
    # h_zpj_gen_p.Rebin(2)
    # h_zpj_gen_h.Rebin(2)
    h_zpj_gen_p.Scale(1./h_zpj_gen_p.Integral())
    h_zpj_gen_h.Scale(1./h_zpj_gen_h.Integral())
    h_zpj_gen_p.Divide(h_zpj_gen_h)

    fout = ROOT.TFile(output_file, "RECREATE")
    h_dijet_reco_p.Write()
    h_dijet_gen_p.Write()
    h_zpj_reco_p.Write()
    h_zpj_gen_p.Write()


if __name__ == "__main__":

    ALGOS = ["ak4", "ak8"]
    PUS = ["chs", "puppi"]
    for algo, pu, in product(ALGOS, PUS):
        SETUP = algo + pu
        # Check your directories are in the right places
        PYTHIA_DIR = "MGPythia/workdir_%s_mgpythia" % SETUP
        HERWIG_DIR = "Herwig/workdir_%s_herwig" % SETUP
        make_weights_file(PYTHIA_DIR, HERWIG_DIR, "herwig_weights_%s.root" % SETUP)

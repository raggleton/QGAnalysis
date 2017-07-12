#!/usr/bin/env python

"""Construct weights for Herwig to ensure jet spectra matches Pythia

Need for both Dijet & Z+jets regions, for both Reco & Gen jets
"""

import ROOT


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


if __name__ == "__main__":
    f_pythia = ROOT.TFile("../Selection/100K/workdir_ak4chs/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root")
    f_herwig = ROOT.TFile("../Selection/Herwig/workdir_ak4chs_herwig/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root")

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

    f_pythia = ROOT.TFile("../Selection/100K/workdir_ak4chs/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root")
    f_herwig = ROOT.TFile("../Selection/Herwig/workdir_ak4chs_herwig/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL.root")
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

    fout = ROOT.TFile("herwig_weights.root", "RECREATE")
    h_dijet_reco_p.Write()
    h_dijet_gen_p.Write()
    h_zpj_reco_p.Write()
    h_zpj_gen_p.Write()
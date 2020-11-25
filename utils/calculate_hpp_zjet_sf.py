#!/usr/bin/env python

"""
Calculate the scale factor needed for high-pT Heriwg++ Z+jet sample.

This is because the cross section doesn't match the inclusive sample,
so instead we just normalise in a given pT range to the inclusive sample.
"""

from __future__ import print_function, division

from math import hypot
import numpy as np
from array import array
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


def calc_sf(inclusive_hist, high_pt_hist, pt_min, pt_max):
    """Calculate the scale factor to be applied to the high_pt_hist,
    based on the ratio of integrals between pt_min and pt_max"""
    xax = inclusive_hist.GetXaxis()
    first_bin = xax.FindBin(pt_min)
    last_bin = xax.FindBin(pt_max)
    inclusive_err = array('d', [0])  # do this way as expecting reference in C++
    high_pt_err = array('d', [0])
    inclusive_integral = inclusive_hist.IntegralAndError(first_bin, last_bin, inclusive_err)
    high_pt_integral = high_pt_hist.IntegralAndError(first_bin, last_bin, high_pt_err)
    ratio = (inclusive_integral / high_pt_integral)
    err = ratio * hypot((inclusive_err[0] / inclusive_integral), (high_pt_err[0] / high_pt_integral))
    return ratio, err


def plot_hists(inclusive_hist, high_pt_hist, output_filename, title=""):
    canv = ROOT.TCanvas("c", "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLogy()
    inclusive_hist.SetLineColor(ROOT.kBlue)
    inclusive_hist.SetMarkerColor(ROOT.kBlue)
    high_pt_hist.SetLineColor(ROOT.kRed)
    high_pt_hist.SetMarkerColor(ROOT.kRed)
    hst = ROOT.THStack("hst", "%s;Parton k_{T} [GeV];N" % title)
    hst.Add(inclusive_hist)
    hst.Add(high_pt_hist)
    hst.Draw("NOSTACK HIST E0")
    hst.SetMinimum(1E-1)
    # hst.SetMaximum(1E6)
    leg = ROOT.TLegend(0.56, 0.7, 0.88, 0.88)
    leg.AddEntry(inclusive_hist, "Inclusive sample", "L")
    leg.AddEntry(high_pt_hist, "High-k_{T} sample", "L")
    leg.Draw()
    canv.SaveAs(output_filename)


def plot_pt_Z_hists(inclusive_hist, high_pt_hist, output_filename, title=""):
    canv = ROOT.TCanvas("c", "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLogy()
    inclusive_hist.SetLineColor(ROOT.kBlue)
    inclusive_hist.SetMarkerColor(ROOT.kBlue)
    high_pt_hist.SetLineColor(ROOT.kRed)
    high_pt_hist.SetMarkerColor(ROOT.kRed)
    hst = ROOT.THStack("hst", "%s;Matrix-Element p_{T}^{ll} [GeV];N" % title)
    hst.Add(inclusive_hist)
    hst.Add(high_pt_hist)
    hst.Draw("NOSTACK HIST E0")
    hst.SetMinimum(1E-1)
    # hst.SetMaximum(1E6)
    leg = ROOT.TLegend(0.56, 0.7, 0.88, 0.88)
    leg.AddEntry(inclusive_hist, "Inclusive sample", "L")
    leg.AddEntry(high_pt_hist, "High-k_{T} sample", "L")
    leg.Draw()
    canv.SaveAs(output_filename)


def plot_sf_vs_pt(pts, sfs, output_filename):
    canv = ROOT.TCanvas("c", "", 800, 600)
    canv.SetTicks(1, 1)
    print(len(pts))

    graph = ROOT.TGraphErrors(len(pts), array('d', [p for p in pts]), array('d', [s[0] for s in sfs]), array('d', [0 for s in sfs]), array('d', [s[1] for s in sfs]))
    graph.Draw("ALP")
    graph.SetTitle("Scale factor vs minimum parton k_{T};Minimum parton k_{T} [GeV];Scale factor")
    canv.SaveAs(output_filename)


if __name__ == "__main__":
    inclusive_filename = "/nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Selection/Herwig/workdir_ak4puppi_herwig_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_passGenRsp_sameGenCuts_mySamples/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_Incl.root"
    high_pt_filename = "/nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Selection/Herwig/workdir_ak4puppi_herwig_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_passGenRsp_sameGenCuts_mySamples/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_JetKtMin170.root"

    inclusive_filename = "/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_Incl_all.root"
    high_pt_filename = "/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_JetKtMin170_PartonKtMin300.root"

    inclusive_filename = "/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/Herwig/workdir_102X_v2_ak4puppi_herwig_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_Incl.root"
    high_pt_filename = "/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/Herwig/workdir_102X_v2_ak4puppi_herwig_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_JetKtMin170_PartonKtMin300.root"

    hist_name = "ZPlusJets_gen/jet_kt"

    inclusive_tfile = ROOT.TFile(inclusive_filename)
    inclusive_hist = inclusive_tfile.Get(hist_name)

    high_pt_tfile = ROOT.TFile(high_pt_filename)
    high_pt_hist = high_pt_tfile.Get(hist_name)

    rebin = 1
    append = "_rebin" if rebin > 1 else ""
    inclusive_hist.Rebin(rebin)
    high_pt_hist.Rebin(rebin)

    output_filename = "herwigpp_unscaled%s.pdf" % (append)
    plot_hists(inclusive_hist=inclusive_hist,
               high_pt_hist=high_pt_hist,
               output_filename=output_filename,
               title="Unscaled histograms")

    pt_min, pt_max = 300, 700
    sf, sf_err = calc_sf(inclusive_hist, high_pt_hist, pt_min, pt_max)
    print("Integral:", pt_min, "to", pt_max, "gives scale factor for high pT sample:", sf, "+/-", sf_err)
    print("You should scale your cross section by this value")
    print("i.e. you should **divide** your Lumi argument by this value")

    # Test convergence of the scale factor,
    # by testing different lower bounds on the integral range
    # The idea being that at a certain point it should ~converge
    pt_mins = np.arange(200, 520, 20)
    sfs = [calc_sf(inclusive_hist, high_pt_hist, pt, pt_max) for pt in pt_mins]
    print(pt_mins)
    print(sfs)
    plot_sf_vs_pt(pt_mins, sfs, output_filename="herwigpp_sf_vs_pt%s.pdf" % (append))

    high_pt_hist_scaled = high_pt_hist.Clone()
    high_pt_hist_scaled.Scale(sf)
    output_filename = "herwigpp_scaled_highpt%s.pdf" % (append)
    plot_hists(inclusive_hist=inclusive_hist,
               high_pt_hist=high_pt_hist_scaled,
               output_filename=output_filename,
               title="Normalised to %g < parton k_{T} < %g GeV" % (pt_min, pt_max))

    inclusive_hist_ptZ = inclusive_tfile.Get("ZPlusJets_gen/pt_z")
    high_pt_hist_ptZ = high_pt_tfile.Get("ZPlusJets_gen/pt_z")
    inclusive_hist_ptZ.Rebin(rebin)
    high_pt_hist_ptZ.Rebin(rebin)
    high_pt_hist_ptZ.Scale(sf)
    output_filename = "herwigpp_scaled_highpt_ptZ%s.pdf" % (append)
    plot_pt_Z_hists(inclusive_hist=inclusive_hist_ptZ,
                    high_pt_hist=high_pt_hist_ptZ,
                    output_filename=output_filename,
                    title="Normalised to %g < parton k_{T} < %g GeV" % (pt_min, pt_max))

    inclusive_tfile.Close()
    high_pt_tfile.Close()

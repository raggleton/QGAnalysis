#!/usr/bin/env python

"""
Calculate correct weights, create XML fragments
"""

import math
import sys
import os
import pandas as pd
import xml.etree.ElementTree as ET
import ROOT


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

# Load CSV files as DataFrames
df_bkg = pd.read_csv("/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Samples_CSV/Ntuple_production_(Run_II,_80X,_Moriond17)_-_Background_samples.csv", header=1)


# DO ALL THE FOLLOWING, the short name must correspond to that in the CSV
DO_THESE = [
('MC_DYJetsToLL_M-50_HT-70to100', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/DYJetsToLL_HT70to100.xml'),
('MC_DYJetsToLL_M-50_HT-100to200', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/DYJetsToLL_HT100to200.xml'),
('MC_DYJetsToLL_M-50_HT-200to400', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/DYJetsToLL_HT200to400.xml'),
('MC_DYJetsToLL_M-50_HT-400to600', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/DYJetsToLL_HT400to600.xml'),
('MC_DYJetsToLL_M-50_HT-600to800', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/DYJetsToLL_HT600to800.xml'),
('MC_DYJetsToLL_M-50_HT-800to1200', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/DYJetsToLL_HT800to1200.xml'),
('MC_DYJetsToLL_M-50_HT-1200to2500', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/DYJetsToLL_HT1200to2500.xml'),
('MC_DYJetsToLL_M-50_HT-2500toInf', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/DYJetsToLL_HT2500toInf.xml'),
('MC_QCD_HT50to100', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT50to100.xml'),
('MC_QCD_HT100to200', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT100to200.xml'),
('MC_QCD_HT200to300', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT200to300.xml'),
('MC_QCD_HT300to500', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT300to500.xml'),
('MC_QCD_HT500to700', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT500to700.xml'),
('MC_QCD_HT700to1000', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT700to1000.xml'),
('MC_QCD_HT1000to1500', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT1000to1500.xml'),
('MC_QCD_HT1500to2000', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT1500to2000.xml'),
('MC_QCD_HT2000toInf', '/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/datasets/MG-Pythia/QCD_HT2000toInf.xml'),
]


TEMPLATE = """        <InputData Lumi="{LUMI}" NEventsMax="&NEVT;" Type="MC" Version="{VERSION}" Cacheable="False">
        {XMLPATH}
        <InputTree Name="AnalysisTree"/>
        </InputData>

"""



def get_xsec(df, short_name):
    this_df = df[df['Short name'] == short_name]
    if len(this_df.index) == 0:
        raise RuntimeError("No matches for dataset: %s / short_name: %s" % (dataset, short_name))
    elif len(this_df.index) > 1:
        print this_df
        raise RuntimeError("More than 1 match for dataset: %s / short_name: %s" % (dataset, short_name))
    return float(this_df['Cross section [pb]'])


def get_num_events(list_of_files):
    """Get the total number of events in the files

    NB tree name is hard wired!
    """
    chain = ROOT.TChain("AnalysisTree")
    for fname in list_of_files:
        chain.Add(fname)
    return chain.GetEntries()


def get_lumiweight(xml, shortname):
    # get list of files from XML
    with open(xml) as f:
        contents = f.read()

    # hack to parse as proper XML
    tree = ET.fromstring("<root>"+contents+"</root>")
    filenames = [child.attrib['FileName'] for child in tree]
    num_events = get_num_events(filenames)
    xsec = get_xsec(df_bkg, shortname)
    # print "LumiWeight: ", num_events/xsec
    return num_events/xsec


if __name__ == "__main__":
    for shortname, xml in DO_THESE:
        weight = get_lumiweight(xml, shortname)
        if not math.isnan(weight):
            # print shortname, "=", weight
            print TEMPLATE.format(VERSION=shortname, XMLPATH=os.path.abspath(xml), LUMI=weight)

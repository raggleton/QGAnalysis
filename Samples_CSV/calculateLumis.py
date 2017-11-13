#!/usr/bin/env python


"""Script to calculate the necessary lumis for each sample when you don't use every event"""


import os
import numpy as np
import pandas as pd
import ROOT


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)


# Load CSV files as DataFrames
df_bkg = pd.read_csv("/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Samples_CSV/Ntuple_production_(Run_II,_80X,_Moriond17)_-_Background_samples.csv", header=1)


def parse_sample(df, num_events=0, dataset=None, xml_filename=None):
    """Calculate correct lumi weight for one sample, identified either by XML filename or full dataset."""
    if dataset != None:
        this_df = df[df['Sample name'] == dataset]
    elif xml_filename != None:
        this_df = df[df['Short name'] == xml_filename]
    else:
        raise RuntimeError("Need dataset or xml_filename to parse sample")

    if len(this_df.index) == 0:
        raise RuntimeError("No matches for dataset: %s / xml_filename: %s" % (dataset, xml_filename))
    elif len(this_df.index) > 1:
        print this_df
        raise RuntimeError("More than 1 match for dataset: %s / xml_filename: %s" % (dataset, xml_filename))

    return num_events / float(this_df['Cross section [pb]'])


def get_num_events(list_of_files):
    """Get the total number of events in the files

    NB tree name is hard wired!
    """
    total_num_events = 0
    for fname in list_of_files:
        # print fname
        if not os.path.isfile(fname):
            raise RuntimeError("No such file %s" % fname)
        f = ROOT.TFile(fname)
        tree = f.Get("AnalysisTree")
        total_num_events += tree.GetEntries()
        f.Close()
    return total_num_events


def parse_samples(samples_dicts):
    """Go through all user samples, calc new lumi weight, printout to screen"""
    for sd in samples_dicts:
        num_events = get_num_events(sd['local_filenames'])
        sd['num_events'] = num_events
        new_lumi = parse_sample(df=sd['df'],
                                num_events=num_events,
                                dataset=sd.get("dataset", None),
                                xml_filename=sd.get("xml_filename", None))
        sample_label = sd.get("xml_filename", None) or sd.get("dataset", None)
        print "%s (%d events): %10f" % (sample_label, num_events, new_lumi)


if __name__ == "__main__":
    # List all your samples here
    ntuple_dir = "../Ntuples/100K/"
    my_setup = [
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_DYJetsToLL_M-50_HT70to100.root"],
            "xml_filename": "MC_DYJetsToLL_M-50_HT-70to100"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_DYJetsToLL_M-50_HT100to200.root"],
            "xml_filename": "MC_DYJetsToLL_M-50_HT-100to200"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_DYJetsToLL_M-50_HT200to400.root"],
            "xml_filename": "MC_DYJetsToLL_M-50_HT-200to400"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_DYJetsToLL_M-50_HT400to600.root"],
            "xml_filename": "MC_DYJetsToLL_M-50_HT-400to600"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_DYJetsToLL_M-50_HT600to800.root"],
            "xml_filename": "MC_DYJetsToLL_M-50_HT-600to800"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_QCD_HT100to200.root"],
            "xml_filename": "MC_QCD_HT100to200"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_QCD_HT200to300.root"],
            "xml_filename": "MC_QCD_HT200to300"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_QCD_HT300to500.root"],
            "xml_filename": "MC_QCD_HT300to500"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_QCD_HT500to700.root"],
            "xml_filename": "MC_QCD_HT500to700"
        },
        {
            "df": df_bkg,
            "local_filenames" : [ntuple_dir + "Ntuple_MC_QCD_HT700to1000.root"],
            "xml_filename": "MC_QCD_HT700to1000"
        }

    ]

    parse_samples(my_setup)

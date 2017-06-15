#!/usr/bin/env python

"""
Auto create the necessary XML snippets from my custom Ntuples, with correct weights.
"""

import sys
sys.path.append("../../core/python")
from wrapper import samples_dict
import os
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import HTMLParser
from collections import OrderedDict
import ROOT


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

# CHANGE ME - dir of NTuples
ntuple_dir = "/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Ntuples/100K"

# CHANGE ME - output dir for XML file
xml_dir = "/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/config/100K"

# For all XML input snippets:
common_attr = dict(NEventsMax="&NEVT;", Type="MC", Cacheable="False")


# Load CSV files as DataFrames
df_bkg = pd.read_csv("/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Samples_CSV/Ntuple_production_(Run_II,_80X,_Moriond17)_-_Background_samples.csv", header=1)
df_bkg_names = df_bkg['Short name'].tolist()


def get_lumi_weight(df, num_events, short_name):
    """Calculate correct lumi weight for one sample, identified either by XML filename (short_name) or full dataset."""
    this_df = df[df['Short name'] == short_name]

    if len(this_df.index) == 0:
        raise RuntimeError("No matches for dataset: %s / short_name: %s" % (dataset, short_name))
    elif len(this_df.index) > 1:
        print this_df
        raise RuntimeError("More than 1 match for dataset: %s / short_name: %s" % (dataset, short_name))

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


def create_xml_snippet(sdict, common_attr):
    """Create an XML snippet for this sample dictionary"""
    new_attrib = {"Lumi" : "%.8f" % sdict['lumi_weight'],
                  "Type": sdict['type'],
                  "Version": sdict['short_name']}
    new_attrib.update(common_attr)
    root = ET.Element("InputData", attrib=new_attrib)
    for f in sdict['local_filenames']:
        fattrib = {"FileName": f, "Lumi": "0.0"}
        in_el = ET.SubElement(root, "In", attrib=fattrib)
    in_tree = ET.SubElement(root, "InputTree", attrib={"Name": "AnalysisTree"})
    # need to convert &amp; back to & etc
    contents = HTMLParser.HTMLParser().unescape(ET.tostring(root))
    return contents


if __name__ == "__main__":
    xml_snippets = []
    # sort by sample name
    ordered_samples_dict = OrderedDict(sorted(samples_dict.items(), key=lambda t: t[0]))
    for sname, sdict in ordered_samples_dict.iteritems():
        sdict['df'] = df_bkg  # TODO: select
        sdict['type'] = "MC"  # TODO: select
        sdict['local_filenames'] = [os.path.join(ntuple_dir, sdict['output'])]  # TODO: why is this a list?
        sdict['short_name'] = sname
        num_events = get_num_events(sdict['local_filenames'])
        sdict['num_events'] = num_events
        lumi_weight = get_lumi_weight(df=sdict['df'],
                                      num_events=num_events,
                                      short_name=sname)
        sdict['lumi_weight'] = lumi_weight
        print "%s (%d events): %10f" % (sname, num_events, lumi_weight)
        snip = create_xml_snippet(sdict, common_attr)
        xml_snippets.append(snip)

    if not os.path.isdir(xml_dir):
        os.makedirs(xml_dir)
    xml_filename = os.path.join(xml_dir, "input.xml")
    with open(xml_filename, "w") as f:
        f.write("\n".join(xml_snippets))

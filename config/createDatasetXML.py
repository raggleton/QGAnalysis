#!/usr/bin/env python

"""Create XML files with <In> elements of input files. Loops over many crab output dirs

Usage: update DO_THESE with the correct crab output dirs and XML filenames
"""

import sys
import glob
import os


def do_xml_file(directories, xml_filename):
    xml_filename = os.path.abspath(xml_filename)
    this_dir = os.path.dirname(xml_filename)
    if not os.path.isdir(this_dir):
        os.makedirs(this_dir)
    with open(xml_filename, "w") as outfile:
        for directory in directories:
            glob_query = os.path.join(directory, "*.root")
            file_list = glob.glob(glob_query)
            if len(file_list) == 0:
                raise RuntimeError("O files found for %s" % glob_query)
            else:
                print "%d files for %s" % (len(file_list), glob_query)
            for file in file_list:
                outfile.write('<In FileName="%s" Lumi="0.0"/>\n' % os.path.abspath(file))


if __name__ == "__main__":
    # For each output XML, we provide a list of dirs to glob over
    DO_THESE = [
        # MG-Pythia DY
        ('../datasets/MG-Pythia/DYJetsToLL_HT70to100.xml', 
            [
                'DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-70to100_mg-pythia_18_Nov_17_newRecoJetFlav_v4'
            ]
        ),
        ('../datasets/MG-Pythia/DYJetsToLL_HT100to200.xml', 
            [
                'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-100to200_mg-pythia_16_Nov_17_newRecoJetFlav_v2',
                'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-100to200_ext_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/DYJetsToLL_HT200to400.xml', 
            [
                'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-200to400_mg-pythia_16_Nov_17_newRecoJetFlav_v2',
                'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-200to400_ext_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/DYJetsToLL_HT400to600.xml', 
            [
                'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-400to600_mg-pythia_16_Nov_17_newRecoJetFlav_v2',
                'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-400to600_ext_mg-pythia_18_Nov_17_newRecoJetFlav_v4'
            ]
        ),
        ('../datasets/MG-Pythia/DYJetsToLL_HT600to800.xml', 
            [
                'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-600to800_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/DYJetsToLL_HT800to1200.xml', 
            [
                'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-800to1200_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/DYJetsToLL_HT1200to2500.xml', 
            [
                'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-1200to2500_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/DYJetsToLL_HT2500toInf.xml', 
            [
                'DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-2500toInf_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        # MG-Pythia QCD
        ('../datasets/MG-Pythia/QCD_HT50to100.xml', 
            [
                'QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT50to100_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/QCD_HT100to200.xml', 
            [
                'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT100to200_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/QCD_HT200to300.xml', 
            [
                'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT200to300_mg-pythia_18_Nov_17_newRecoJetFlav_v4'
            ]
        ),
        ('../datasets/MG-Pythia/QCD_HT300to500.xml', 
            [
                'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT300to500_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/QCD_HT500to700.xml', 
            [
                'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT500to700_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/QCD_HT700to1000.xml', 
            [
                'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT700to1000_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/QCD_HT1000to1500.xml', 
            [
                'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1000to1500_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/QCD_HT1500to2000.xml', 
            [
                'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1500to2000_mg-pythia_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/MG-Pythia/QCD_HT2000toInf.xml', 
            [
                'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT2000toInf_mg-pythia_18_Nov_17_newRecoJetFlav_v4'
            ]
        ),
        # Herwig
        ('../datasets/MG-Herwig/DYJetsToLL.xml', 
            [
                'DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/crab_DYJetsToLL_M-50_mg-herwig_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        ('../datasets/HerwigOnly/QCD_Pt15to7000_Flat.xml', 
            [
                'QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/crab_QCD_Pt-15to7000_herwig_16_Nov_17_newRecoJetFlav_v2'
            ]
        ),
        # Pythia only
        ('../datasets/PythiaOnly/QCD_Pt15to7000_Flat.xml', 
            [
                'QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/crab_QCD_Pt-15to7000_pythia_16_Nov_17_newRecoJetFlav_v2/'
            ]
        ),
        # Powheg
        ('../datasets/Powheg-Pythia/Dijet.xml', 
            [
                'Dijet_NNPDF30_powheg_pythia8_TuneCUETP8M1_13TeV_bornktmin150/crab_Dijet_NNPDF30_powheg_pythia_18_Nov_17_newRecoJetFlav_v4/'
            ]
        ),
    ]

    BASE = "/pnfs/desy.de/cms/tier2/store/user/raggleto/"
    for xml_filename, directories in DO_THESE:
        do_xml_file([os.path.join(BASE, d, "*/*/") for d in directories], xml_filename)

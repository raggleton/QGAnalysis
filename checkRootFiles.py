#!/usr/bin/env python

"""Check all ROOT files in directory to see if valid

Can optionally delete them
"""

import os
import sys
import argparse
from glob import glob

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gErrorIgnoreLevel = ROOT.kWarning


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory with ROOT files") 
    parser.add_argument("--rm",
                        action='store_true',
                        help='Remove bad files')
    args = parser.parse_args()

    if not os.path.isdir(args.source):
        raise IOError("%s not valid directory" % args.source)

    bad_filenames = []
    for filename in glob(args.source+"/*.root"):
        f = ROOT.TFile(filename)
        if f.IsZombie() or not f or f.TestBit(ROOT.TFile.kRecovered):
            print 'Bad file:', filename
            bad_filenames.append(filename)
            if args.rm:
                os.remove(filename)
        else:
            f.Close()
    print bad_filenames
    print 'Found', len(bad_filenames), 'bad files'

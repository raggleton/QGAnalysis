#!/usr/bin/env python

"""Remove zombie ROOT file - the ones that won't open properly

Usage:

    ./removeZombies.py <directory-with-root-files>

"""


from __future__ import print_function

import os
import sys
import ROOT
from glob import glob

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)


if __name__ == "__main__":
    num_del = 0
    num_updated = 0
    for f in glob(os.path.join(sys.argv[1], "*.root")):
        # print(f)
        rootf = ROOT.TFile(f)
        if not rootf or rootf is None or rootf.IsZombie():
            print("Bad file", f)
            rootf.Close()
            num_del += 1
            os.remove(f)
        elif rootf.TestBit(ROOT.TFile.kRecovered):
            # to solve this in future, write file again with recovered keys
            print("Fixing keys", f)
            rootf.Close()
            rootf = ROOT.TFile(f, "UPDATE")
            rootf.Write()
            rootf.Close()
            num_updated += 1
        else:
            rootf.Close()
    print("Finished")
    print("num deleted:", num_del)
    print("num updated:", num_updated)


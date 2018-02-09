#!/usr/bin/env python

"""
Print prescales for triggers from csv file

e.g. 

brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_DATACERT.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt  -u /pb --hltpath 'HLT_PFJet*_v*' -o pfjet_trig.csv
"""


import csv
import sys
from collections import defaultdict
import re


def tryint(s):
    try:
        return int(s)
    except:
        return s
     
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]


results = defaultdict(float)

with open(sys.argv[1]) as f:
	reader = csv.reader(f)
	for row in reader:
		# print row
		if row[0].startswith("#"):
			continue
		trigname = re.sub('_v.*', '', row[3])
		lumi = float(row[-1])
		results[trigname] += lumi

total_lumi = 35863818448.025

for k in sorted(results.keys(), key=alphanum_key):
	print k, ' = ', total_lumi / results[k]

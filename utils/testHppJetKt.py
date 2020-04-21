#!/usr/bin/env python

"""
Test to see if jet KT is same at GEN and ntuple
"""

from __future__ import print_function, division

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


f_gen = ROOT.TFile("/pnfs/desy.de/cms/tier2/store/user/raggleto/PrivateHerwigppDYJetsToLL/DYJetsToLL_M-50_JetKtMin175_TuneCUETHS1_13TeV-herwigpp/herwigpp_zjet_jetktmin175_flat_pow2p5_500K_2020_03_12/200312_151951/0000/JME-RunIISummer15GS-00016_982.root")
tree_gen = f_gen.Get("Events")

f_ntuple = ROOT.TFile("/nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src/UHH2/core/python/Ntuple_JetKtMin170_982.root")
tree_ntuple = f_ntuple.Get("AnalysisTree")

print("******Gen tree:")
gen_dict = {}
for ievt, evt in enumerate(tree_gen):
    evt_num = evt.generator.GetEvent().event_number()
    largest_pt = -99
    # print("Event:", evt_num)
    for igp, gp in enumerate(evt.genParticles):
        # if gp.status() != 11:
        #     continue
        abs_id = abs(gp.pdgId())
        if abs_id < 6 or abs_id == 21:
            largest_pt = max(largest_pt, gp.pt())
            # print("using", largest_pt, gp.pdgId(), gp.status())
    gen_dict[evt_num] = largest_pt

print("******Ntuple tree:")
ntuple_dict = {}
for ievt, evt in enumerate(tree_ntuple):
    evt_num = evt.event
    # print("Event:", evt_num)
    largest_pt = -99
    for igp, gp in enumerate(evt.GenParticles):
        if gp.status() != 11:
            continue
        abs_id = abs(gp.pdgId())
        if abs_id < 6 or abs_id == 21:
            largest_pt = max(largest_pt, gp.pt())
            # print("using", largest_pt, gp.pdgId(), gp.status())
    ntuple_dict[evt_num] = largest_pt


gen_evts = set(gen_dict.keys())
ntuple_evts = set(ntuple_dict.keys())

missing = gen_evts ^ ntuple_evts
print("Missing", missing)
common = gen_evts & ntuple_evts
for k in common:
    alert = ''
    gen_val = gen_dict[k]
    ntuple_val = ntuple_dict[k]
    if gen_val != ntuple_val:
        alert = ' ********** '
    print(k, 'Gen:', gen_val, "ntuple:", ntuple_val, alert)

print("min:", min(ntuple_dict.values()))

f_gen.Close()
f_ntuple.Close()

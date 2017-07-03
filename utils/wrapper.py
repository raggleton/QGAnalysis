import FWCore.ParameterSet.Config as cms
from ntuplewriter import *
import sys
import os

NUM_EVT = 100000

# samples in /pnfs/desy.de/cms/tier2/

# output = Ntuple_<key>.root
samples_dict = {
    # DY MC
    # "MY_DYJetsToLL_M-50_1": {
    #     "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00099D43-77ED-E611-8889-5065F381E1A1.root', ],
    #     "output": "Ntuple_MC_DYJetsToLL_M-50_1.root",
    #     "nevents": -1
    # },
    # "MY_DYJetsToLL_M-50_2": {
    #     "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/008B4E7F-8FED-E611-897F-5065F381F291.root'],
    #     "output": "Ntuple_MC_DYJetsToLL_M-50_2.root",
    #     "nevents": -1
    # },
    # "MC_DYJetsToLL_M-50_HT-70to100": {
    #     "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/52D899A9-92D0-E611-857A-001E674FCAE9.root'],
    #     "output": "Ntuple_MC_DYJetsToLL_M-50_HT-70to100.root",
    #     "nevents": NUM_EVT
    # },
    "MC_DYJetsToLL_M-50_HT-100to200": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/0EEDEFFF-58CE-E611-8ED3-0025905A612E.root'],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT-100to200.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_M-50_HT-200to400": {
        "input": [
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/76738514-DAD2-E611-8F95-D067E5F9156C.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/5A1CAB38-FFD1-E611-B1FA-A0000420FE80.root'
            ],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT-200to400.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_M-50_HT-400to600": {
        "input": [
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/52D899A9-92D0-E611-857A-001E674FCAE9.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/48A0E465-46D0-E611-94D2-001E674DA1AD.root'
        ],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT-400to600.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_M-50_HT-600to800": {
        "input": [
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/0A0083AF-EFBD-E611-ABF7-F04DA275BF11.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/98CE9891-55BD-E611-89FD-02163E017703.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/1206115B-44BD-E611-A3FF-02163E00B1A5.root'
        ],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT-600to800.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_M-50_HT-800to1200": {
        "input": [
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/0456978B-A2C0-E611-9368-02163E0119D7.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/C0344D02-B8C0-E611-B8BC-0CC47A4D7630.root'
        ],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT-800to1200.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_M-50_HT-1200to2500": {
        "input": [
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/1CE9292C-7CBD-E611-AE79-1866DAEA79A4.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/8E3E0DF0-94BD-E611-B665-001EC94BA119.root'
        ],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT-1200to2500.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_M-50_HT-2500toInf": {
        "input": [
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/947C3809-B4C0-E611-8344-0025905A609E.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/CAD47782-AEC0-E611-A8D3-02163E00B829.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/5C8C9696-AFC0-E611-BF54-0025905AA9CC.root',
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/0C99DA93-79C1-E611-8D0B-FA163EEA85D4.root'
        ],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT-2500toInf.root",
        "nevents": NUM_EVT
    },
    # QCD MC
    # "MC_QCD_HT50to100": {
    #     "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/128D78DF-AFB6-E611-BB25-14187764197C.root'],
    #     "output": "Ntuple_MC_QCD_HT50to100.root",
    #     "nevents": NUM_EVT
    # },
    "MC_QCD_HT100to200": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/AA7E2F95-FDBD-E611-8BEB-70106F4A929C.root'],
        "output": "Ntuple_MC_QCD_HT100to200.root",
        "nevents": NUM_EVT
    },
    "MC_QCD_HT200to300": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/02FDEBBA-9ABD-E611-888C-001E67E6F7F6.root'],
        "output": "Ntuple_MC_QCD_HT200to300.root",
        "nevents": NUM_EVT
    },
    "MC_QCD_HT300to500": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/04B528C6-7FBD-E611-AD61-3417EBE6446E.root'],
        "output": "Ntuple_MC_QCD_HT300to500.root",
        "nevents": NUM_EVT
    },
    "MC_QCD_HT500to700": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/000316AF-9FBE-E611-9761-0CC47A7C35F8.root'],
        "output": "Ntuple_MC_QCD_HT500to700.root",
        "nevents": NUM_EVT
    },
    "MC_QCD_HT700to1000": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00D17FD4-8EBD-E611-B17D-002590D0AFC2.root'],
        "output": "Ntuple_MC_QCD_HT700to1000.root",
        "nevents": NUM_EVT
    },
    "MC_QCD_HT1000to1500": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0E5C4A43-04B7-E611-B3D9-A0000420FE80.root'],
        "output": "Ntuple_MC_QCD_HT1000to1500.root",
        "nevents": NUM_EVT
    },
    "MC_QCD_HT1500to2000": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/1C1FD02A-21BB-E611-8071-6C3BE5B50180.root'],
        "output": "Ntuple_MC_QCD_HT1500to2000.root",
        "nevents": NUM_EVT
    },
    "MC_QCD_HT2000toInf": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/02FBE5E0-1AB6-E611-85B1-0CC47A7DFF82.root'],
        "output": "Ntuple_MC_QCD_HT2000toInf.root",
        "nevents": NUM_EVT
    },

}

# auto setup output names
for k, v in samples_dict.iteritems():
    v['output'] = 'Ntuple_'+k+".root"


# For Herwig samples
Herwig_DYJetsToLL_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/3ED57B25-F7DE-E611-810A-24BE05C46B01.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/F4C2FE80-F7DE-E611-8C93-A0369F7FC934.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/029CE79C-05E0-E611-9ABB-1866DAEA79D0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/600CE904-F3E7-E611-B106-0025905A60D0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/F85A36FD-19E0-E611-8512-001E67A42026.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/18DC1706-3EE5-E611-BA51-782BCB67A0F8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/D2CCF66A-DEE2-E611-BDAF-0025902008EC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/24CDCD06-6DDF-E611-ADEF-FA163E7F05CC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/60A7FC29-32DF-E611-AB1E-0090FAA57430.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/2E1E7311-44DF-E611-945F-0090FAA58124.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/110000/1A1BA35C-21E4-E611-956A-FA163EAC6AD1.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/5A8365B3-2DE4-E611-A6C9-002590A370FE.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/383D3C45-3CDF-E611-B532-ECF4BBE16230.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/C07D78AE-FADF-E611-B407-1418774105B6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/B0D5737B-68DF-E611-A9CF-001E67792888.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/BA04A44C-61E3-E611-B7D8-001E677925F6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/110000/F66E11D7-22E3-E611-8A6E-1866DAEECFDC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/84156A86-F9DE-E611-8299-5065F381A2F1.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/3686E1D4-E4E2-E611-9C4D-008CFA0516BC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/FC83ED14-2DE0-E611-A357-001E67F11FC7.root"
]
Herwig_QCD_files = [
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/0401C3EE-60D2-E611-BBB0-1866DAEA7E28.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/06912832-62D2-E611-927A-0025901D0C52.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/087192DE-EFD3-E611-A677-1C6A7A26C8DF.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/0A6F941E-67CF-E611-9644-0025905B85B2.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1A95AB56-75D2-E611-98A9-002590FD5A48.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1E51A681-5DD2-E611-B946-02163E00B926.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/22D4BBF0-19D3-E611-8DC8-008CFA197CA0.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/22FC9BB3-4FD2-E611-AACE-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/246E47CF-69D2-E611-800C-0242AC130003.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/24C6B2CF-73D2-E611-9EA7-02163E013C82.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/2CE66D77-38D3-E611-93B0-0025905B8576.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/3E161A19-6DD2-E611-A3E6-0CC47A78A440.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/40260F97-4CD2-E611-93FA-0CC47A4D75EE.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/4C130C59-3FD2-E611-896B-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/4C49E7ED-43D2-E611-9614-4C79BA320D7D.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/4E289EF0-44D2-E611-B2BC-0CC47A4C8EEA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/4E8449DA-5BCF-E611-8D61-0025905A60B6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/50DF3854-63CF-E611-9B16-0025905A612C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/54529306-58CF-E611-A063-0CC47A7C3638.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/58A51ADC-69D2-E611-8D05-0025905A48C0.root"
]

samples_dict = {}
NUM_EVT = 100000

for ind, f in enumerate(Herwig_DYJetsToLL_files):
    name = "MC_HERWIG_DYJetsToLL_M-50_%d" % ind
    samples_dict[name] = {
        "input": [f],
        "output": "Ntuple_"+name+".root",
        "nevents": NUM_EVT
    }

for ind, f in enumerate(Herwig_QCD_files):
    name = "MC_HERWIG_QCD_%d" % ind
    samples_dict[name] = {
        "input": [f],
        "output": "Ntuple_"+name+".root",
        "nevents": NUM_EVT
    }


def setup_sample(sample_info):
    process.source.fileNames = cms.untracked.vstring(sample_info['input'])
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(sample_info['nevents']))
    process.MyNtuple.fileName = cms.string(sample_info['output'])


if __name__ == "__main__":
    arr_id = int(os.environ['SGE_TASK_ID']) - 1  # -1 as array IDs start at 1
    key = sorted(samples_dict.keys())[arr_id]
    print 'Doing', key
    setup_sample(samples_dict[key])
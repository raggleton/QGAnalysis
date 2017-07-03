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


def setup_sample(sample_info):
    process.source.fileNames = cms.untracked.vstring(sample_info['input'])
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(sample_info['nevents']))
    process.MyNtuple.fileName = cms.string(sample_info['output'])


if __name__ == "__main__":
    arr_id = int(os.environ['SGE_TASK_ID']) - 1  # -1 as array IDs start at 1
    key = sorted(samples_dict.keys())[arr_id]
    print 'Doing', key
    setup_sample(samples_dict[key])
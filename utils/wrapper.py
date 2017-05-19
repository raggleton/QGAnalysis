import FWCore.ParameterSet.Config as cms
from ntuplewriter import *
import sys
import os

NUM_EVT = 20000

samples_dict = {
    # DY MC
    "MC_DYJetsToLL_HT70to100": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/52D899A9-92D0-E611-857A-001E674FCAE9.root'],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT70to100.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_HT100to200": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/362BD21E-47CE-E611-9D61-0025905A612E.root'],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT100to200.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_HT200to400": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/76738514-DAD2-E611-8F95-D067E5F9156C.root'],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT200to400.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_HT400to600": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/52D899A9-92D0-E611-857A-001E674FCAE9.root'],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT400to600.root",
        "nevents": NUM_EVT
    },
    "MC_DYJetsToLL_HT600to800": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/0A0083AF-EFBD-E611-ABF7-F04DA275BF11.root'],
        "output": "Ntuple_MC_DYJetsToLL_M-50_HT600to800.root",
        "nevents": NUM_EVT
    },

    # QCD MC
    "MC_QCD_HT100to200": {
        "input": ['/store/mc/RunIISummer16MiniAODv2/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00DED9E0-F7BD-E611-94B8-02163E015EF5.root'],
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
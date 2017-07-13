import FWCore.ParameterSet.Config as cms
from ntuplewriter import *
import sys
import os
import re


NUM_EVT = 100000


def produce_sample_dict(file_list, name_stem):
    sd = {}
    for ind, f in enumerate(file_list):
        name = "%s_%d" % (name_stem, ind)
        sd[name] = {
            "input": [f],
            "output": "Ntuple_"+name+".root",
            "nevents": NUM_EVT
        }
    return sd


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    """For "natural" sorting using numerical values. Use as key arg in sorted() etc

    Stolen from: https://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    """
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]


# QCD samples
samples_dict = {
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
    }
}

# samples in /pnfs/desy.de/cms/tier2/
samples_dict = {}

# HT binned pythia samples
Pythia_DYJetsToLL_100to200_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/4848FEE1-C8CE-E611-9E62-0CC47A4C8E86.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/601C84AE-60D0-E611-8A73-0025905A48D6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/94E1541E-A4D0-E611-BA5C-0CC47A745294.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E269B108-66CE-E611-8FC5-0025905A60D6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/362BD21E-47CE-E611-9D61-0025905A612E.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/8294F98A-4ECE-E611-8BE8-0CC47A78A458.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4CC85607-9ACE-E611-807D-0025905A48EC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A6A204BB-5ECE-E611-8B18-0CC47A4D769A.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/50A8E072-57CE-E611-9177-0CC47A4C8F10.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/041BA6C0-8ED0-E611-AD0B-0025905A608E.root"
]
samples_dict.update(produce_sample_dict(Pythia_DYJetsToLL_100to200_files, "MC_DYJetsToLL_M-50_HT-100to200"))

Pythia_DYJetsToLL_200to400_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/3E3C8E53-1FD3-E611-AFA2-001E67457A5D.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/6E3A407E-BAD2-E611-B989-B083FECF837B.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/5E97F1F8-04D3-E611-9E11-549F3525DB98.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6493464C-0AD3-E611-BFB4-001E674FAF23.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/76738514-DAD2-E611-8F95-D067E5F9156C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/C0A4F58C-C6D2-E611-8433-5065F3820351.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/7058283A-01D2-E611-B723-E0071B7A48A0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/0472BF59-3DD2-E611-883B-D067E5FA0623.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/50E3F38E-02D2-E611-89C7-0CC47A13D176.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/F6FEFB6F-11D2-E611-93D0-0CC47A0AD69E.root",
]
samples_dict.update(produce_sample_dict(Pythia_DYJetsToLL_200to400_files, "MC_DYJetsToLL_M-50_HT-200to400"))

# --limit=10 --idx=2
Pythia_DYJetsToLL_400to600_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/1E0BBFFA-AED0-E611-B452-001E674FCA99.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/3494FC7C-64D0-E611-8B6F-001E674FB063.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/F00B3200-9BD0-E611-83B9-001E673476AA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1C6C4217-42D0-E611-984C-0CC47A537688.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/3CF1FDC0-FAD0-E611-A8E5-001E67397CB5.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/DA120E92-06CF-E611-8395-001E674DA347.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/46EF88CB-49D0-E611-A4AD-001E67348055.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/BC37B8A5-95D0-E611-B739-001E674FBFC2.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/48DB01D6-15CF-E611-91DC-0CC47A537688.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/ECF12F2F-1ACF-E611-AE84-001E674DA1AD.root",
]
samples_dict.update(produce_sample_dict(Pythia_DYJetsToLL_400to600_files, "MC_DYJetsToLL_M-50_HT-400to600"))

# --idx=250 --limit=10
Pythia_DYJetsToLL_600to800_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/2E0728C0-23BD-E611-B283-001E6779238C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/0A0083AF-EFBD-E611-ABF7-F04DA275BF11.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/6A63D7AE-6ABD-E611-91C0-FA163EA451AB.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/6A4FB05C-0BBE-E611-83EB-001E67C7AF3F.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/265552A7-E9BD-E611-9488-0CC47A7C340C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/78AE4D07-D3BE-E611-A22F-FA163E074AF8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/FC63A375-C2BD-E611-BBD5-02163E017665.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/10FF29D2-A8BD-E611-9719-B499BAAC0068.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/92AB1C8E-99BD-E611-9C9D-7845C4FC3B8D.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/987CD64E-B6BD-E611-8375-00259022277E.root",
]
samples_dict.update(produce_sample_dict(Pythia_DYJetsToLL_600to800_files, "MC_DYJetsToLL_M-50_HT-600to800"))

#  --limit=10 --idx=60
Pythia_DYJetsToLL_800to1200_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/96E2ECFF-5FC0-E611-BADC-782BCB20D86D.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/7A17E1A0-ACC0-E611-8E8F-02163E00BCE4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/D85FB7C3-ADC0-E611-A207-0CC47A1E0C76.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/5CE931E7-6BC0-E611-82F1-008CFA582C28.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/56762ABE-92C0-E611-8044-FA163E616349.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/16AF32AC-A8C0-E611-B282-0025905A60C6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/B68B5E99-9AC0-E611-A0DB-002618FDA208.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/CA520016-88C0-E611-8BBD-0090FAA57260.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/ACF987CC-98C0-E611-9FAF-02163E00B2D9.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/4EDD4524-99C0-E611-9200-02163E0146AE.root",
]
samples_dict.update(produce_sample_dict(Pythia_DYJetsToLL_800to1200_files, "MC_DYJetsToLL_M-50_HT-800to1200"))

# --idx=20
Pythia_DYJetsToLL_1200to2500_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8AA311D6-EEBD-E611-8455-0CC47AD98F70.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/7E08BD8B-93BE-E611-890D-0CC47A7C34E6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/F877A452-5BBE-E611-B76D-B083FED42488.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/E6F44563-CFBE-E611-853F-0CC47A78A4A0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/248B8E65-7FBD-E611-926B-00266CF9BDFC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/72066A7C-C8BD-E611-A1B9-0025905A612C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/2AE9984E-D2BD-E611-8213-0025905A6066.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/EAA7B2BB-95BD-E611-AAFA-0CC47AD98D08.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/3ADC4A69-8FBD-E611-8739-B083FED42FC4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/E2A0DD9A-77BD-E611-92B4-7845C4F92E7F.root",
]
samples_dict.update(produce_sample_dict(Pythia_DYJetsToLL_1200to2500_files, "MC_DYJetsToLL_M-50_HT-1200to2500"))

# --idx=10
Pythia_DYJetsToLL_2500toInf_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/F2CBBFBF-AAC0-E611-8F2A-FA163E4BD49B.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/844FF1A4-B5C0-E611-90C3-0CC47A4D7630.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/BE1438F4-B1C0-E611-A7B5-FA163E6EC48D.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/308C289B-B8C0-E611-A7C6-0025905A606A.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/5C8C9696-AFC0-E611-BF54-0025905AA9CC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/CAD47782-AEC0-E611-A8D3-02163E00B829.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/947C3809-B4C0-E611-8344-0025905A609E.root",
]
samples_dict.update(produce_sample_dict(Pythia_DYJetsToLL_2500toInf_files, "MC_DYJetsToLL_M-50_HT-2500toInf"))


# For Herwig samples
# --limit=80 --idx=175
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

"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/168F3697-6FDF-E611-9AA3-FA163EFDA7C4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/74B272A7-B1E3-E611-9AEF-0CC47A7FC7C8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/08831EF2-9BDF-E611-94C2-008CFAFBFEBE.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/80069042-FCDF-E611-859A-001E67E71DFD.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/0273B7FF-83DF-E611-A3B3-FA163E1DA94B.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/2AFFBF07-5AE0-E611-A221-FA163E979FFE.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/A29DA66C-44DF-E611-BCA1-485B3919F0B5.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/6440F8AF-CCDF-E611-B8A1-001E67E5A39A.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/D6375BC0-FFDE-E611-A303-0CC47AD99238.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/564C0808-A3DF-E611-BC8C-008CFA002830.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/72150FF0-E3DE-E611-9968-B083FED42A6E.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/02456AB6-BBDF-E611-97A4-0CC47AD98F70.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/9051B6E4-6FDF-E611-9BD6-001E67792608.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/72EF507E-FDDE-E611-9612-0025901F8740.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/110000/6ACC9402-2AE3-E611-91B0-0CC47ADAF3DA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/C4B3AADA-26DF-E611-B3AC-A4BF01011BF7.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/A82F9B52-21E0-E611-928C-4C79BA320DB9.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/04A78E35-F6DE-E611-B822-001E67A402B2.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/9EB9E28D-1ADF-E611-A26A-3417EBE706ED.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/C0DADECA-93DF-E611-80EF-0CC47AD98C88.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/C6DDD9CF-B6E5-E611-8B53-FA163E0DF1DA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/0C8F8079-FDDE-E611-9806-001E67DFF7CB.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/302C3A69-E9E7-E611-A376-24BE05BD4FC2.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/5E869B34-F1DF-E611-A556-02163E013B28.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/E48D2AAD-FDDE-E611-B8DF-0090FAA57FC4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/C6388305-06DF-E611-9008-A0369F7F92FC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/6A500C3A-19DF-E611-946D-047D7BD6DF5E.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/58B40AFA-53DF-E611-B55A-20CF3027A5B5.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/4459F107-C4E9-E611-AD87-90B11C0BD467.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/7641DF43-EFDE-E611-B27F-0090FAA582E4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/B6E5AB2F-33E0-E611-A222-FA163E0E1A14.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/D6739E25-ACDF-E611-9E34-A4BADB1E6B36.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/FA720CC6-DCE1-E611-92B9-FA163E4D51D9.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/28AD2DAB-DCE1-E611-9ADD-FA163E8928E8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/64731DB7-DEE1-E611-808E-02163E013CD5.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/D6C5793A-3CE0-E611-9469-001E67E95BD8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/9AC4F292-DDE1-E611-9395-02163E012EBF.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/3E7C3735-0CDF-E611-A415-D8D385FF35E1.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/029AF739-B5DF-E611-8A77-141877411C7F.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/A08014EB-26DF-E611-8C8A-0090FAA57600.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/0CD44D44-0CDF-E611-9E89-0CC47A7E6A6C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/46BF4209-E4E1-E611-BA94-0025904B1F90.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/8CCA03EF-3EDF-E611-B167-E41D2D08E110.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/AE5ABD7C-DAE5-E611-9B96-0025905A6110.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/808F3C9A-75E3-E611-9D27-001E67E34165.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/CC8C528A-43DF-E611-96C1-70106F49CDF4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/586E5B8F-90DF-E611-B1CD-848F69FD2955.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/110000/58487268-FFE2-E611-98E7-0CC47A4D765A.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/32E6A2AB-26DF-E611-AE51-002590FCF738.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/42E6987F-F3E7-E611-895A-141877411D83.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/5840E9DB-01E0-E611-B088-E0071B7A58B0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/5EB2DF7A-EBDE-E611-BF6A-002590747E28.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/3EE6B37E-E0DF-E611-A3CD-B083FED04CAB.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/AC7B9565-B0DF-E611-963F-848F69FD45A7.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/AA34FB03-F6DE-E611-A10A-00266CF3E0FC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/C2EF26A4-EDDE-E611-95ED-D4AE526A0B47.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/1E694E8C-EFDE-E611-8957-002590DD7E50.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/D2022D1D-45E0-E611-88EB-FA163E9B5AE6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/1AA29BC0-19E0-E611-802F-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/1CB9D4A4-75E9-E611-BB01-A4BF01027688.root",

"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/110000/A0498243-CAE3-E611-9455-001E67A3FD26.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/3AEC7843-82E9-E611-A1C8-008CFA111330.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/5A24060D-73EA-E611-967A-001E673CC77F.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/C02C8177-E2E5-E611-8AE7-008CFA197E58.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/BAD258B5-E6E5-E611-9CB5-0242AC130004.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/BE8CE8D9-F1E0-E611-8816-008CFAFBFB94.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/02003B4D-E7E1-E611-BF78-FA163EC0471C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/D4CE9C1A-E1E5-E611-B20D-0025905A611C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/A2C0EF2E-D4DF-E611-87D9-001E67DDD0AA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/967B7C9F-59E0-E611-9D2E-FA163EC78C09.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/96FD0413-DFE1-E611-9339-0025904B2C2C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/C4089704-01E0-E611-A603-001E67444EAC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/C634E712-D8E0-E611-AA3A-0CC47A745298.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/7278F5CA-BAE8-E611-9B7D-001E674FAF23.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/2804C089-F7E0-E611-9DAE-C4346BC84780.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/70C069A3-12E0-E611-88AD-00266CFFC7B0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/7C911CF5-F8E0-E611-875F-00259073E466.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/2AC5FBF8-E9DF-E611-8C01-008CFA11114C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/200640B7-79E4-E611-A519-008CFAF750B6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/DC6E9371-CAE8-E611-8D0F-0CC47A4D76B2.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/B0814511-82E9-E611-B920-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/BC9A5509-60E3-E611-A45D-FA163E4A6ED3.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/54338C90-38E5-E611-AA8D-141877411FED.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/1A2F20C9-BAE8-E611-9133-0CC47A57D086.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/F6CAA3A4-E1E5-E611-9015-44A84224BE51.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/72DBB0A5-06E0-E611-B3BF-70106F4A91E8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/ACA22895-85DF-E611-8569-E0071B7A8590.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/7CD80360-CADF-E611-AAE5-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/110000/720F854C-C7E3-E611-B091-00259075D6B4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/4808294B-B5E8-E611-A2CB-008CFA111330.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/80916F1F-E1E5-E611-B518-0CC47A745298.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/0A841C92-BAE3-E611-9229-008CFA152104.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/7059879D-99DF-E611-92D8-0CC47AD98D00.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/F46C0133-CEE0-E611-9A3E-24BE05C68671.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/7A88E614-BBE0-E611-B2AA-001E67C9B465.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/E02EF1A0-77E0-E611-825B-A0369FC5ED30.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/5E8442CD-67E0-E611-9BBA-008CFA197E58.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/CAF01998-3DE0-E611-AF62-001E67E6F8AA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/1C6CF104-86E0-E611-9C49-0CC47AD98D08.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/7A2D14E4-F6E9-E611-B951-0025905A6094.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/340D565D-58E3-E611-8BA2-6CC2173BC0B0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/A2208EEF-57E3-E611-9159-FA163E4BA6AB.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/26E404AE-E3DF-E611-963C-002590DD61E0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/2CB65936-C5DF-E611-A84E-002590D9D8B6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/7E84996E-B9DF-E611-AE9E-00259074AE82.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/5C41EC5F-E4DF-E611-8B61-0025901D0C52.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/54AB9758-F7E9-E611-B370-001E674FB149.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/0A278ACF-F6E0-E611-B6FB-0025905B85DE.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/08A4D558-D4E0-E611-8EE3-002590AC4C72.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/9C27FC69-86E4-E611-A101-002590AC4FC8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/DA340939-EAE0-E611-8C46-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/5ACDBB39-CEE0-E611-BE49-14187741208F.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/023538A4-B6E8-E611-9A1B-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/4826CE69-CDDF-E611-B761-00259074AE54.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/90000/884FF8A7-E1E5-E611-B013-002590AB3A82.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/2820B309-EBE0-E611-AF2B-0025905B857C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/46DE450E-C0E8-E611-9C7D-001E67D80528.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/F0B2F66E-C5DF-E611-9DF4-0019B9CB0300.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/9EE96CDB-B1E0-E611-B075-047D7B881D1E.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/FEC894E6-D0E0-E611-830B-0CC47A78A496.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/8E1B214A-BDE0-E611-81D8-008CFA110C78.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/925A2E1A-57E0-E611-81FB-0025904B1F92.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/789046E4-E0DF-E611-9DA7-003048322C46.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/38929ED7-A1DF-E611-98CB-0CC47A7E6A3E.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/E4D91DE8-22E0-E611-A4BC-0090FAA57AE0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/110000/6422F165-DDE3-E611-BA5E-0025905B8572.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/9266D308-A9DF-E611-A5D0-008CFA111244.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/020D5832-7CE3-E611-9B31-1CC1DE053050.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/EE37ECD8-ABE9-E611-8629-ECF4BBE1BE70.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/50000/448E9E80-F6DE-E611-B6B1-0CC47AA989C6.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/38B95FA9-05E0-E611-B205-002590E3A224.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/E23DF795-2AE0-E611-8853-D067E5F914D3.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/147EC226-44DF-E611-9E91-008CFA0518D4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/4A9BA054-2AE0-E611-A933-001E67792872.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/FCD3D146-22DF-E611-9583-A0369FC5E22C.root",
]

# das_client.py --query="file dataset=/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM | sort file.size | grep file.name" --idx=29 --limit=80
Herwig_QCD_files = [
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/8EAD8B76-38D3-E611-A34E-0025905A60E4.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/CAB7C0F5-77D5-E611-A2B7-02163E011DBE.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/68C3AAD0-48D2-E611-B366-002590D9D8D4.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/A42256D0-48D2-E611-B4D0-001E67A3FDF8.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/DA783774-53CF-E611-AC42-0CC47A4D7616.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/B8680DFD-A3CF-E611-A20B-0CC47A7C353E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/4E8449DA-5BCF-E611-8D61-0025905A60B6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/BAF29B73-BDCF-E611-9762-0CC47A4D768E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1A339F33-96CF-E611-A8C0-0025905A609E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/14CA1357-FECF-E611-B09B-0025905A613C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/9C73E145-5BD0-E611-A734-0CC47A78A4A0.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/24C6B2CF-73D2-E611-9EA7-02163E013C82.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/4C130C59-3FD2-E611-896B-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/E626D86C-6BD2-E611-8359-B083FED406AD.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/90B0CE5B-57D2-E611-A8F8-001E67DFFB86.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/AABF2807-AECF-E611-83FF-0025905A6070.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/5AEF0FE9-98D2-E611-9FA3-02163E013268.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/348BFB0B-B8CF-E611-B371-0CC47A78A3B4.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/FEE57F90-A8CF-E611-964F-0025905B8592.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/866D2C05-7DCF-E611-928B-0025905A60D6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/CE812454-FCCF-E611-A52A-00248C55CC40.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/DC1B20A3-13D0-E611-A7FC-0025905A6110.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/6025067D-5AD2-E611-9197-02163E013A09.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/F6DE5C53-56D2-E611-8FB8-0CC47AA989BE.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/EEFE70B5-C0CF-E611-8133-0CC47A7C3410.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/AC0EEC50-78CF-E611-95DA-0CC47A7C345E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/4E289EF0-44D2-E611-B2BC-0CC47A4C8EEA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/FA6E7F70-A8CF-E611-803A-0CC47A745294.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/8CE25629-51D2-E611-BE30-0025905B856C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/04C750B7-A6CF-E611-BA92-0CC47A4C8E20.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/5ACE6E78-4ED2-E611-BB88-002590E2F9D4.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/B2A55617-A5CF-E611-A9F8-00248C55CC3C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CA4A3955-81CF-E611-ABDB-0CC47A4D769A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1E51A681-5DD2-E611-B946-02163E00B926.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/84724259-58D2-E611-A281-842B2B42BC3A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/0A6F941E-67CF-E611-9644-0025905B85B2.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/06912832-62D2-E611-927A-0025901D0C52.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/C4A5867C-4ECF-E611-B54B-0CC47A7C3420.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/54529306-58CF-E611-A063-0CC47A7C3638.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/A89B3957-57D2-E611-97DF-02163E0176D4.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1A95AB56-75D2-E611-98A9-002590FD5A48.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/844DA74B-67D2-E611-9E24-0025905B85B6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/10CCEDB3-A1CF-E611-B567-003048FFD798.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4AF53AD0-A2CF-E611-A73E-0025905A612C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/68B9B6C6-55D2-E611-850A-C81F66B7EAE5.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/0401C3EE-60D2-E611-BBB0-1866DAEA7E28.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2C3023E1-6ECF-E611-A56A-0CC47A4D7650.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/D0A33AB7-5CD2-E611-8586-141877344D39.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/B815ED5B-5BD2-E611-89AF-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/F6D39503-58CF-E611-8605-0CC47A7C340E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/DE4A6EED-51D2-E611-A160-0CC47A4D767C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/50DF3854-63CF-E611-9B16-0025905A612C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/B4E7B27B-A0CF-E611-B994-0025905A60BE.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6A66EB4E-9CCF-E611-B94C-0CC47A78A4BA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4CB2A17F-A0CF-E611-98AD-0025905B8600.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E6FE12A2-9ECF-E611-8076-0025905A60DA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/EA460D32-78CF-E611-95DE-0CC47A4C8E28.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/ACCFF00E-54CF-E611-9EDD-0025905B85BA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/54DCB20C-76CF-E611-808C-0CC47A78A3EE.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/C2470C37-48D2-E611-8EE6-782BCB20DC03.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/945B1955-7BCF-E611-BB99-0CC47A74527A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/D0523749-9FCF-E611-8A14-0CC47A78A3EE.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/BA6C3CBF-A2CF-E611-AD93-0025905B8566.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/0A5A443F-A6CF-E611-837E-0025905A60DA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5A7B666D-85CF-E611-B1B6-0CC47A4C8E98.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/9448A80C-93CF-E611-BAB5-0CC47A4D7602.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/869DB2C0-C1D2-E611-A3D9-02163E019D42.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/243F94A6-6DCF-E611-9696-0025905A60B0.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/72A7C164-85CF-E611-A85C-0CC47A4C8F0A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/54E1A7FF-D5D2-E611-87F2-02163E014563.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/D484DC72-B3CF-E611-88E8-0CC47A4D767E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/642F0552-ADCF-E611-8DB5-0CC47A7C3604.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/7266BC0D-AFCF-E611-9BBB-0025905B8564.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/22FC9BB3-4FD2-E611-AACE-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/9090C38A-83CF-E611-B09C-0CC47A4D7640.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/A41881D1-B5CF-E611-B1E0-0CC47A4C8E46.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/B260102A-52D2-E611-BDD4-02163E00ADDD.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/24B20AF0-94CF-E611-BBEA-0CC47A4C8F10.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/58A51ADC-69D2-E611-8D05-0025905A48C0.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/06B3FF01-B5CF-E611-BCD9-0CC47A4D768E.root",
]

samples_dict = {}
# NUM_EVT = 100000

samples_dict.update(produce_sample_dict(Herwig_DYJetsToLL_files, 'MC_HERWIG_DYJetsToLL_M-50'))
# samples_dict.update(produce_sample_dict(Herwig_QCD_files, "MC_HERWIG_QCD"))

# inclusive pythia samples
# --limit=100 --idx=25
Pythia_DYJetsToLL_files = [
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/70E1BAEB-CDC4-E611-AC82-008CFA19741C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/26461FA6-3AC4-E611-9AA8-008CFA11136C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/4A7546B7-77C4-E611-8F5D-008CFA197904.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/E05826A9-86C4-E611-A6AF-C4346BBCB6A8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/B206DF7C-34C4-E611-8298-C4346BC7EDD8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/A8FD91A4-76C4-E611-82E7-001E67DBE700.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/4EFFE6C7-DEC3-E611-823D-001E67F337BC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/B27F2ADE-2DC3-E611-917C-00215AA65030.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/80DC6B18-BBC6-E611-B6D4-001E67DBE3EF.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/9AD9DB3A-AAC4-E611-BF13-008CFA1112B0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/9457B4BD-50C3-E611-9165-001A648F1C36.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/DA8D0892-44C6-E611-B425-6CC2173BBF00.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/6645E291-A1C4-E611-9905-008CFA197B74.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/F83B71BE-C1C4-E611-BA9D-00237DF244A8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/7AE0D8B0-91C4-E611-BEF8-008CFA197C34.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/6AD3B402-A9C4-E611-B3EF-6CC2173BC0A0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/DC0E703D-B1C3-E611-85B1-001A648F1DC2.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/E6B38172-2AC3-E611-9FD6-002481ACEA28.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/862DCFB8-6DC3-E611-9EFD-002481AD3168.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/3ECFE8DA-FEC4-E611-99A1-C4346BC75558.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/E659E662-F7C3-E611-8127-008CFA110C78.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/40B69D14-60C3-E611-9048-001A648F153E.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/5A133CC2-91C4-E611-A5B5-6CC2173BBA40.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/BA414361-2CC3-E611-BF87-00237DF28460.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/0CDDA2E8-4CC4-E611-816D-001E0BC1ECC8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/9AFA8857-B0C4-E611-999A-008CFA197AF4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/10ADD2A8-9CC4-E611-9A1D-00266CFF0034.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/9CFAB56E-3CC5-E611-B47C-001E67F3332A.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/0CF33E90-B4C6-E611-9B5A-008CFA197BD4.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/D228D640-DBC5-E611-87A7-008CFA197B74.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/30F48B9A-68C3-E611-A945-00266CFFBE5C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/38D061B1-C8C4-E611-AECB-6CC2173BB980.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/DCD89426-4CC4-E611-A4E5-001E0BC198D8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/9A1E501B-C0C4-E611-8050-00266CFFC948.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/4AB71EB5-41C4-E611-8E49-001E67E5E889.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/6CD7353D-0BC4-E611-85F9-001A648F1ACA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/5AAE228B-A1C4-E611-8F83-008CFA1979FC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/065F4D50-B4C4-E611-98E7-008CFA110DD8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/2CF3885D-6CC3-E611-B9D4-001F29EACF54.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/AEF1D97D-7AC3-E611-8C6D-001A648F7EF2.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/BCB3923C-33C3-E611-AE00-001E0BC1E34C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/F0643443-D6C4-E611-9A6B-008CFA05206C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/382EF287-41C3-E611-9432-00237DF29490.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/40C8A837-4CC4-E611-8EB4-00215E2EB6A0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/C2D0ADAC-42C4-E611-AB04-008CFA110A98.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/302A9628-9BC7-E611-AF16-001E67F33348.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/2840DEBD-17C4-E611-9481-6CC2173BBA60.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/9CC2E5FA-73C3-E611-968F-001E0BC1ECC8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/0646F754-D8C4-E611-8D2A-00215E2EB72A.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/C8DBD5EF-F2C4-E611-B41D-008CFA197B74.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/FE203AE8-08C3-E611-9A26-6CC2173BC0A0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/6C70B38B-5AC4-E611-8DBA-001E0BC1DE40.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/36B4EFDB-79C4-E611-ACE2-001E67DBE700.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/3877AFE4-28C4-E611-8D0A-001E0BC1DDEA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/4C1E45A4-72C4-E611-AABF-C4346BC808B8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/AE57FFCB-54C4-E611-96C7-6CC2173BBD70.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/94DC48A5-6DC3-E611-8E5E-C4346BC808B8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/0AE16771-E3C4-E611-842D-6CC2173BBED0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/188BC1B5-BBC3-E611-9BF5-001E67F120A8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/C87FFBB2-BBC4-E611-9FF1-00266CFFBCD0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/1E854471-31C4-E611-AA4B-00266CFFC948.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/BA1B6F53-E3C3-E611-B114-001E67F3332A.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/04D1AC84-DEC3-E611-8301-001E67F33348.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/78CD6016-DCC4-E611-A33A-C4346BC84780.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/A0841E0B-96C4-E611-8BCE-008CFA1112B0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/3A010999-DEC4-E611-9C3A-008CFA11131C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/B2275F85-DBC3-E611-93B8-001A648F1A56.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/7AEA0309-5AC4-E611-9FC0-00266CFF090C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/ACD5108C-DEC4-E611-97C1-008CFA110DD8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/444E4ED0-B0C4-E611-884E-001E67E5E889.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/14C4A28E-4DC5-E611-9AD2-C4346BC75558.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/8CA67EE3-4AC4-E611-9208-00266CFFCA1C.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/7248D349-B8C3-E611-88B6-001E67E5EAF0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/024D0A89-34C4-E611-AEF9-008CFA111348.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/48ABA96A-3CC4-E611-839D-002481ACEC90.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/EE0E5A01-E1C4-E611-8694-00266CFFBC38.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/B0C44FA4-BFC6-E611-A699-001E67F120A8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/B4A6EC5D-80C4-E611-ABFC-008CFA111314.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/683C8BB2-B2C3-E611-B62A-001A648F1DEA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/405D889E-39C4-E611-82D8-6CC2173BBA60.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/720EF489-F0C3-E611-9783-00266CFFCCBC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/26649648-8CC7-E611-93B2-00266CFFBC60.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/6877DC40-51C4-E611-A897-C4346BC87798.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/2C95CA8D-DEC4-E611-8C68-008CFA197904.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/F283191C-11C4-E611-973D-00215E2EB74E.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/D2287568-4CC5-E611-AA9D-001E67DBDF5D.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/0840E216-1DC3-E611-8753-00266CF89498.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/B6DFD257-E3C3-E611-AA33-00266CFEFDEC.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/46124E04-21C4-E611-A31E-001E0BC1DE40.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/E6D778F4-8BC4-E611-88CE-C4346BBCB6A8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/6E5577BC-1AC4-E611-9BB1-00215E2EB6A0.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/360532E5-2DC4-E611-840C-002481ACEC90.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/D2A4CB2C-2FC4-E611-BCB7-6CC2173BBA60.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/CEAB3688-1CC7-E611-8BC3-C4346BBCB6A8.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/C4FDDA7B-A9C3-E611-A23B-001A648F1CAA.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/70000/F623FA51-8DC4-E611-A7E7-008CFA197C34.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/A08B1790-FDC3-E611-828B-6CC2173BBD70.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/1A8BEA4B-C5C6-E611-A82A-001E67F333BB.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/585B2BBC-A9C3-E611-8DF9-002481FF6531.root",
"/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/8801FDF4-68C4-E611-A0B9-00266CFEFDE0.root",
]

# das_client.py --query="file dataset=/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM | sort file.size | grep file.name" --limit=75 --idx=20
Pythia_only_flat_QCD_files = [
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/C6FC2F30-63BE-E611-A730-001E674FBFC2.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/F0FE6C72-67BE-E611-9D19-90B11C2C93B1.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/0C2141B5-62BE-E611-82EA-1866DAEA6D08.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/F0168BF2-6FBE-E611-8BA4-002590D8C68A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/0C61D691-6ABE-E611-8E67-002590D8C68A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/B8D10797-6CBE-E611-A835-002590791D3C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/8291DA18-75BE-E611-9DEF-0242AC130002.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/1EDFE14B-6ABE-E611-9873-782BCB47D0B6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/06E658EF-65BE-E611-BA87-0025907D2212.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/F82578FF-72BE-E611-8235-20CF3027A582.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A8238E31-64BE-E611-8F3F-D067E5F914D3.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/AAADCB4E-73BE-E611-A675-001EC94B4EF5.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/64B984B1-69BE-E611-927E-0025905C2CE6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/CE26F592-68BE-E611-A0A8-FA163E8B40DC.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/088CE532-70BE-E611-8264-002590791D3C.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/8C8B17EA-67BE-E611-BC81-001E67DFFB86.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/5E639D3F-6ABE-E611-AC14-00259048A860.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E828EB9B-70BE-E611-A05F-0CC47A78A45A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/CCE2D52B-62BE-E611-85D1-0CC47A78A3B4.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A23B7C9D-70BE-E611-9AB1-0CC47A4C8E0E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/947F1B79-65BE-E611-9F2E-0CC47A4C8E5E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/B65F9C2F-63BE-E611-8390-0025907D2212.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/54B35CD4-66BE-E611-B8E9-02163E017630.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/46B8FBBE-63BE-E611-A577-002590E7DDE6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/14D3659F-62BE-E611-BDDB-0CC47A4C8E5E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E2D8B0AD-62BE-E611-BDEF-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/2CA2F01C-6EBE-E611-A431-D4AE527F2883.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/6ED4FC75-6CBE-E611-B58C-0CC47A4D7674.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A0F754FE-63BE-E611-A93B-02163E013035.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/12E82276-61BE-E611-B400-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/F6AFD94A-6BBE-E611-BF10-002590E7E004.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/58C8F64D-62BE-E611-B7C5-0025905A6066.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DABED09D-6ABE-E611-AAF1-90B11C2CB121.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/5E128058-6FBE-E611-9163-FA163E8B7DCF.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/90B5973F-64BE-E611-B8AE-0025904C66F6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/722DC8C3-5FBE-E611-A027-0CC47A13D16E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/AC1906A9-61BE-E611-8018-20CF3027A59F.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A0F450C5-62BE-E611-AB04-0025905C4300.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/6EE69CA5-5DBE-E611-A456-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/02E84E91-61BE-E611-AB67-B083FED03632.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/B056021F-6DBE-E611-9FBB-002590D8C68A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/60D73874-60BE-E611-82C5-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/72A5934D-65BE-E611-88D9-0025905A60D0.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/B872C727-67BE-E611-99AB-0025905C543A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/9CF7A17C-6ABE-E611-8615-0CC47A4D76D6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/F678DACC-65BE-E611-BF1D-00259073E504.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/4E3C12DD-6BBE-E611-A74E-0025905A60CA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/5E11DA30-61BE-E611-B18B-02163E013035.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/CE772BB1-63BE-E611-9256-FA163E479835.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/58D38CE8-67BE-E611-8FFD-001E674FBA1D.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/CEEF6B2D-61BE-E611-8071-FA163E747C1B.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/467C1CF2-5FBE-E611-88B2-02163E0130D0.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/9471ABF0-61BE-E611-A98B-90B11C08CDC7.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/36D1F109-65BE-E611-BBA2-0025905C53DA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E8C1FB0B-61BE-E611-B64F-0025905A60CA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/949AE742-5FBE-E611-ACBF-0242AC130002.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E441E702-71BE-E611-B3CC-D067E5F914D3.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/20083B7E-68BE-E611-BFE2-0CC47A13D3A8.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A27FD12D-65BE-E611-A139-002590FD5122.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/D6A68248-6BBE-E611-95BF-00259073E3E0.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/840CBC45-62BE-E611-B95D-02163E0131D5.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/42C501A7-62BE-E611-A8E3-001E67A40451.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/30F9ECA6-69BE-E611-98A4-0CC47A74524E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/2456B812-64BE-E611-9318-002590D8C68A.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/62A6FB9C-67BE-E611-A506-0090FAA57310.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/981F1294-64BE-E611-9619-FA163E4BEB03.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/F48ABE8D-6FBE-E611-96B2-0025905C2CBA.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/9CA860D0-7EBE-E611-8BB3-001E675A681F.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A637E08A-6ABE-E611-8D2F-0025905B860E.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DE6FF30E-63BE-E611-8ED0-A0000420FE80.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/30E3A7B2-63BE-E611-BF81-0025905B85B6.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/342AD212-66BE-E611-86A0-0242AC130002.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/E62BABD0-66BE-E611-AB35-FA163E90EFC1.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/2C6C48A6-68BE-E611-85B5-90B11C0BD63B.root",
"/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/A0D6B55A-63BE-E611-9EFE-0090FAA57AF0.root",
]

samples_dict = {}
# samples_dict.update(produce_sample_dict(Pythia_DYJetsToLL_files, "MC_PYTHIA_DYJetsToLL_M-50"))
samples_dict.update(produce_sample_dict(Pythia_only_flat_QCD_files, "MC_PYTHIA_ONLY_FLAT_QCD"))

def setup_sample(sample_info):
    process.source.fileNames = cms.untracked.vstring(sample_info['input'])
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(sample_info['nevents']))
    process.MyNtuple.fileName = cms.string(sample_info['output'])


if __name__ == "__main__":
    arr_id = int(os.environ['SGE_TASK_ID']) - 1  # -1 as array IDs start at 1
    # arr_id = 0
    key = sorted(samples_dict.keys(), key=natural_sort_key)[arr_id]
    print 'Doing', key
    setup_sample(samples_dict[key])

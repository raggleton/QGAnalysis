#!/usr/bin/env python

"""
This script is designed to submit lots of CRAB jobs.
"""


import subprocess
import os
from time import strftime
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from httplib import HTTPException
from DasQuery import autocomplete_Datasets
from multiprocessing import Process, Pool


inputDatasets = [
# MG+PYTHIA
'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v*/MINIAODSIM',
'/DYJetsToLL_M-50_HT-*to*_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_*/MINIAODSIM',
'/QCD_HT*to*_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v*/MINIAODSIM',
# MG+HERWIG
'/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v*/MINIAODSIM',
'/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v*/MINIAODSIM',
# PYTHIA ONLY
'/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v*/MINIAODSIM',
#
][-1:]

filter_keywords = [
'GenJets5',
'BGenFilter',
]

requestNameCustom = "_test8"


def filter_datasets(input_datasets):
    return [x for x in input_datasets if not any(y in x for y in filter_keywords)]


def create_request_name(input_dataset):
    name = input_dataset.split('/')[1]
    modified_name = (name.replace('_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','')
                         .replace('_TuneCUETHS1_13TeV-madgraphMLM-herwigpp', '')
                         .replace('_TuneCUETP8M1_FlatP6_13TeV_pythia8', '')
                         .replace('_TuneCUETHS1_Flat_13TeV_herwigpp', ''))
    if 'ext' in input_dataset:
        modified_name += '_ext'
    if 'madgraphMLM-pythia8' in input_dataset:
        modified_name += '_mg-pythia'
    elif 'madgraphMLM-herwigpp' in input_dataset:
        modified_name += '_mg-herwig'
    elif '_pythia8' in input_dataset:
        modified_name += '_pythia'
    elif '_herwigpp' in input_dataset:
        modified_name += '_herwig'
    modified_name += "_" + strftime('%d_%b_%y') + requestNameCustom
    return modified_name


def create_crab_config(request_name, input_dataset):
    conf = config()
    conf.General.workArea = 'crab_test'
    conf.General.transferOutputs = True
    conf.General.transferLogs = True
    conf.General.requestName = request_name

    conf.JobType.pluginName = 'Analysis'
    conf.JobType.psetName = '../../core/python/ntuplewriter.py'
    conf.JobType.outputFiles = ["Ntuple.root"]
    conf.JobType.maxMemoryMB = 3000
    #conf.JobType.inputFiles = ['/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_4_15_patch1/src/UHH2/core/python/Summer15_25nsV2_MC.db']

    conf.Data.inputDataset = input_dataset
    conf.Data.inputDBS = 'global'
    conf.Data.splitting = 'EventAwareLumiBased'
    conf.Data.unitsPerJob = 50000
    # conf.Data.totalUnits = 
    # conf.Data.outLFNDirBase = '/store/user/%s/QGNtuples/' % (getUsernameFromSiteDB())
    conf.Data.publication = False
    conf.JobType.sendExternalFolder = True
    #conf.Data.allowNonValidInputDataset = True

    conf.Site.storageSite = 'T2_DE_DESY'

    return conf


def write_crab_config_file(conf, filename):
    contents = conf.pythonise_()
    if not os.path.isdir(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    with open(filename, "w") as f:
        f.write(contents)


def submit_config(config_filename):
    # Why not use crabCommand? Well it does weird thing importing the analysis config
    # So you end up with duplciated modules sometimes
    subprocess.call('crab submit -c %s' % os.path.abspath(config_filename), shell=True)


def main():
    print "Getting datasets..."
    input_datasets = filter_datasets(autocomplete_Datasets(inputDatasets))
    request_names = [create_request_name(in_data) for in_data in input_datasets]
    print "Creating Configuration objects..."
    configs = [create_crab_config(rn, in_data) for rn, in_data in zip(request_names, input_datasets)]
    config_filenames = [os.path.join('crab_configs', c.General.requestName, 'crab.py') for c in configs]
    for c, fn in zip(configs, config_filenames):
        write_crab_config_file(c, fn)
    print "Submitting..."
    try:
        p = Pool(5)
        # result = p.map_async(submit_config, configs)
        result = p.map_async(submit_config, config_filenames)
        print result.get()
        p.close()
    except KeyboardInterrupt:
        print 'Got ^C while pool mapping, terminating the pool'
        p.terminate()
        print 'Pool is terminated'
    except Exception, e:
        print 'Got exception: %r, terminating the pool' % (e,)
        p.terminate()
        print 'Pool is terminated'
    finally:
        print 'Joining pool processes'
        p.join()
        print 'Join complete'
    print 'The end'


if __name__ == '__main__':
    main()

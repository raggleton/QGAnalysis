#!/usr/bin/env python

"""Create copies of XML configs with different jet radii and systematics varied, (re)submit, check, hadd"""

from __future__ import print_function

import os
import sys
import argparse
import re
from collections import namedtuple
import subprocess
from itertools import product
from time import sleep, strftime
import shutil
import stat


# structure to hold info about a given systematic
Systematic = namedtuple("Systematic", ['names', 'values', 'onData', 'onMC'])


def get_workdir(line):
    match = re.search(r'Workdir="(.*)"/>', line)
    if not match:
        raise RuntimeError("Couldn't find Workdir")
    return match.group(1)


def write_updated_file(contents, new_xml_filename, radius, systematic_names=None, systematic_values=None, workdir_append=""):
    # Should really generalise this
    new_workdir = None
    output_dir = None
    with open(new_xml_filename, 'w') as f:
        for line in contents:
            if "Workdir=" in line:
                # update workdir
                orig_workdir = get_workdir(line)
                new_workdir = orig_workdir.replace("_ak4puppi_", "_%s_" % (radius.lower()))
                if systematic_names:
                    if workdir_append is "":
                        # auto generate workdir append
                        workdir_append = "_".join(["%s%s" % (n, v.capitalize()) for n, v in zip(systematic_names, systematic_values)])
                    new_workdir += "_" + workdir_append
                else:
                    if workdir_append:
                        new_workdir += "_" + workdir_append

                new_line = line.replace(orig_workdir, new_workdir)

                # Update RAM if necessary
                if systematic_names and "PDFvariations" in systematic_names:
                    match = re.search(r'RAM *= *"[0-9]+"', new_line)
                    if not match:
                        raise RuntimeError("Can't find RAM requirement")
                    current_RAM = match.group(0)
                    new_RAM = 'RAM="8"'
                    new_line = new_line.replace(current_RAM, new_RAM)

                f.write(new_line)

            elif "&AK4PUPPI;" in line:  # cos the default is AK4PUPPI
                f.write(line.replace("AK4PUPPI", radius))

            # elif systematic_names and "<!--@SYST-->" in line:
            elif "<!--@SYST-->" in line:
                # Turn off the standard hists, we only care about lambda & unfolding ones (currently)
                new_line = ""
                if systematic_names:
                    for key in ["DO_PU_BINNED_HISTS", "DO_FLAVOUR_HISTS", "DO_KINEMATIC_HISTS", "DO_LAMBDA_HISTS"][:-1]:
                        new_line = '            <Item Name="%s" Value="False"/>\n' % (key)
                        f.write(new_line)
                else:
                    pass
                    # new_line = '            <Item Name="DO_KINEMATIC_HISTS" Value="False"/>\n'
                    # new_line = '            <Item Name="DO_LAMBDA_HISTS" Value="False"/>\n'
                    # new_line += '            <Item Name="DO_UNFOLD_HISTS" Value="False"/>\n'
                    # new_line += '            <Item Name="DO_WEIGHT_HISTS" Value="False"/>\n'
                    # f.write(new_line)

                # Write out systematics
                if systematic_values:
                    for sn, val in zip(systematic_names, systematic_values):
                        new_line = '            <Item Name="%s" Value="%s"/>\n' % (sn, val)
                        f.write(new_line)

            elif '<!--' not in line and "<Item Name=" in line:
                if systematic_names:
                    # remove any existing entries for this systematic
                    if any(['<Item Name="%s"' % (sn) in line for sn in systematic_names]):
                        continue
                    elif "JackknifeVariations" in line:
                        # turn off Jackknife for systematics
                        f.write('            <Item Name="JackknifeVariations" Value="False"/>\n')
                    else:
                        f.write(line)
                else:
                    if "JackknifeVariations" in line:
                        # turn off Jackknife for all
                        f.write('            <Item Name="JackknifeVariations" Value="False"/>\n')
                    else:
                        f.write(line)

            # elif "OutputDirectory=" in line:
            #     match = re.search(r'OutputDirectory ?= ?"(.*?)"', line)
            #     if not match:
            #         raise RuntimeError("Can't find OutputDirectory")
            #     output_dir = match.group(0)
            else:
                f.write(line)
    return new_workdir


def submit_xml(xml_filename, local=False, grid_control=False):
    local_opt = "--local" if local else ""
    gc_opt = "--gc" if grid_control else ""
    cmd = 'sframe_batch.py -s %s %s %s' % (local_opt, gc_opt, xml_filename)
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        print(result)
    except Exception as err:
        print(err)
        print(err.output)
    sleep(5)


def resubmit_xml(xml_filename, local=False, grid_control=False):
    local_opt = "--local" if local else ""
    gc_opt = "--gc" if grid_control else ""
    cmd = 'sframe_batch.py -r %s %s %s' % (local_opt, gc_opt, xml_filename)
    try:
        if not local:
            result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
            print(result)
        else:
            subprocess.Popen(cmd, shell=True)
    except Exception as err:
        print(err)
        print(err.output)
    sleep(5)


def check_xml(xml_filename):
    cmd = 'sframe_batch.py %s' % (xml_filename)
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        print(result)
    except Exception as err:
        print(err)
        print(err.output)
    sleep(5)


def hadd_job_gc(workdir, sample, xml_name, cmd, local=False, dry_run=False):
    """Create & submit hadd job using grid-control

    Returns grid-control config filename
    """
    stem = os.path.splitext(os.path.basename(xml_name))[0]
    conf_basename = "hadd_%s_%s.conf" % (sample, stem)
    conf_name = os.path.join(workdir, conf_basename)
    gc_workdir = os.path.join(workdir, 'work.%s' % os.path.splitext(conf_basename)[0])
    if os.path.isdir(gc_workdir):
        # print("Deleting old workdir", gc_workdir)
        shutil.rmtree(gc_workdir)
    print("Writing new conf:", conf_name)
    with open(conf_name, 'w') as f:
        f.write("[global]\n")
        f.write('workdir create = True\n')
        if local:
            f.write('backend = host\n')  # run locally
        else:
            f.write('backend = local\n')  # Send to local batch system

        f.write('[workflow]\n')
        f.write('task = UserTask\n')

        f.write('[jobs]\n')
        f.write('wall time = 3:00\n')  # in hours
        f.write('max retry = 2\n')
        f.write('memory =1000\n')  # in MB
        f.write('jobs = 1\n')  # need this for exactly 1 job
        if local:
            f.write('in flight = 6\n')

        f.write("[UserTask]\n")
        f.write('executable = %s\n' % (os.path.abspath("hadd_gc.sh")))
        # need to transfer env vars ourselves
        f.write("constants = LD_LIBRARY_PATH_STORED PATH_STORED CMSSW_BASE\n")
        f.write('LD_LIBRARY_PATH_STORED = %s\n' % os.getenv('LD_LIBRARY_PATH'))
        f.write('PATH_STORED = %s\n' % os.getenv('PATH'))
        f.write('CMSSW_BASE = %s\n' % os.getenv('CMSSW_BASE'))
        # add XML filenames as arguments
        f.write("parameters = CMD STARTDIR\n")
        f.write('STARTDIR = "%s"\n' % workdir)
        f.write('CMD = "%s"\n' % cmd)

    if dry_run:
        print("Please submit with:")
        print("    grid-control", conf_name)
    else:
        subprocess.check_output(["grid-control", conf_name])
        print("Check with:")
        print("grid-control -c", conf_name)

    return conf_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-s", "--submit", action='store_true', help='submit on sframe_batch')
    group.add_argument("-r", "--resubmit", action='store_true', help='resubmit on sframe_batch')
    group.add_argument("-c", "--check", action='store_true', help='check with sframe_batch')
    group.add_argument("--hadd", action='store_true', help='hadd output with grid-control jobs')
    parser.add_argument("--onlySubmitSystematics", action='store_true', help='Only submit systematic variation jobs, not nominal')
    parser.add_argument("--local", action='store_true', help='Run locally')
    parser.add_argument("--gc", action='store_true', help='Use grid-control')
    args = parser.parse_args()

    base_files = [
        # 'QGAnalysisHerwig.xml',
        # 'QGAnalysisHerwigDYIncl.xml',
        # 'QGAnalysisHerwigDYJetKt170.xml',

        # 'QGAnalysisPythia.xml',

        'QGAnalysisMGPythiaQCD.xml',
        'QGAnalysisMGPythiaDY.xml',
        'QGAnalysisMGPythiaDYInclHTMax70.xml',

        # 'QGAnalysisMGPythiaDYIncl.xml',

        # 'QGAnalysisDataJetHT.xml',
        # 'QGAnalysisDataZeroBias.xml',
        # 'QGAnalysisDataSingleMu.xml',

        # 'QGAnalysisDataZeroBiasTest.xml',
    ]

    # These are entities that are used
    radii = [
        # "AK4PUPPI",
        "AK8PUPPI",
    ]

    # map of hadd command for each XML, and its general output dir (can't infer from XML yet)
    hadd_cmds = {
        'QGAnalysisHerwig.xml': {
            "cmd": "hadd -f uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_*",
            "sample_dir": "Herwig",
        },
        'QGAnalysisHerwigDYIncl.xml': {
            "cmd": "hadd -f uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_merged_PartonKtMin300.root uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_*PartonKt*_*[0-9].root",
            "sample_dir": "Herwig",
        },
        'QGAnalysisHerwigDYJetKt170.xml': {
            "cmd": None,  # dont hadd this, as done as part of QGAnalysisHerwigDYIncl.xml
            "sample_dir": None,
        },

        'QGAnalysisPythia.xml': {
            "cmd": "hadd -f uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD-Pt*",
            "sample_dir": "Pythia",
        },


        'QGAnalysisMGPythiaQCD.xml': {
            "cmd": "hadd -f uhh2.AnalysisModuleRunner.MC.MC_QCD.root uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT*",
            "sample_dir": "MGPythia",
        },
        'QGAnalysisMGPythiaDY.xml': {
            "cmd": "hadd -f uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-[0-9]*.root",
            "sample_dir": "MGPythia",
        },
        'QGAnalysisMGPythiaDYInclHTMax70.xml': {
            "cmd": None, # done as part of QGAnalysisMGPythiaDY.xml
            "sample_dir": None,
        },

        'QGAnalysisDataJetHT.xml': {
            "cmd": "hadd -f uhh2.AnalysisModuleRunner.DATA.Data_JetHT.root uhh2.AnalysisModuleRunner.DATA.Data_JetHT_*.root",
            "sample_dir": "JetHT",
        },
        'QGAnalysisDataZeroBias.xml': {
            "cmd": "hadd -f uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias.root uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias_*.root",
            "sample_dir": "ZeroBias",
        },
        'QGAnalysisDataSingleMu.xml': {
            "cmd": "hadd -f uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_*.root",
            "sample_dir": "SingleMu",
        },
    }

    # Systematics to be changed
    # Note that the original XML should have the "nominal" settings if needed
    # - here are just the variations
    # The values are the possible settings for the names,
    # and each entry in values should match the length of the names
    # - they get done zip-wise, not product-wise
    systematics = [
        # Systematic(names=['chargedHadronShift'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        # Systematic(names=['neutralHadronShift'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        # Systematic(names=['photonShift'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        # Systematic(names=['jecsmear_direction'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        # Systematic(names=['jersmear_direction'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        # Systematic(names=['track_direction'], onMC=True, onData=False,
        #           # values=[('up'), ('down')]),
        Systematic(names=['pileup_direction'], onMC=True, onData=False,
                   # values=[('up'), ('down')]),
                   values=[('down')]),
        # Systematic(names=['PDFvariations'], onMC=True, onData=False,
        #            values=[ ('true') ]),
        # Systematic(names=['ScaleVariationMuR', 'ScaleVariationMuF'], onMC=True, onData=False,
        #            values=[
        #                    ('up', 'up'),
        #                    ('nominal', 'up'),
        #                    ('up', 'nominal'),
        #                    ('down', 'down'),
        #                    ('down', 'nominal'),
        #                    ('nominal', 'down'),
        #                    ]),
    ]

    # JEC individual source systematic
    directions = ['up', 'down']
    jec_sources = [
        "AbsoluteStat",
        "AbsoluteScale",
        "AbsoluteFlavMap",
        "AbsoluteMPFBias",
        "Fragmentation",
        "SinglePionECAL",
        "SinglePionHCAL",
        "FlavorQCD",
        "TimePtEta",
        "RelativeJEREC1",
        "RelativeJEREC2",
        "RelativeJERHF",
        "RelativePtBB",
        "RelativePtEC1",
        "RelativePtEC2",
        "RelativePtHF",
        "RelativeBal",
        "RelativeSample",
        "RelativeFSR",
        "RelativeStatFSR",
        "RelativeStatEC",
        "RelativeStatHF",
        "PileUpDataMC",
        "PileUpPtRef",
        "PileUpPtBB",
        "PileUpPtEC1",
        "PileUpPtEC2",
        "PileUpPtHF",
        "PileUpMuZero",
        "PileUpEnvelope",
        "FlavorZJet",
        "FlavorPhotonJet",
        "FlavorPureGluon",
        "FlavorPureQuark",
        "FlavorPureCharm",
        "FlavorPureBottom",
        "TimeRunBCD",
        "TimeRunEF",
        "TimeRunGH",
        "CorrelationGroupMPFInSitu",
        "CorrelationGroupIntercalibration",
        "CorrelationGroupbJES",
        "CorrelationGroupFlavor",
        "CorrelationGroupUncorrelated",
    ]
    jec_systematics = [
        Systematic(names=['jecsmear_direction', 'jecsmear_source'], onMC=True, onData=True,
                   values=[direction, src_name])
        for direction, src_name in product(directions, jec_sources)
    ]

    hadd_filenames = []

    for xml_filename in base_files:
        with open(xml_filename) as f:
            contents = f.readlines()

        xml_stem = os.path.splitext(xml_filename)[0]

        for radius in radii:
            # no systematics
            new_xml_filename = xml_stem + '_' + radius + '.xml'

            # if radius is not "AK4PUPPI":
            if not args.onlySubmitSystematics:
                new_workdir = write_updated_file(contents, new_xml_filename, radius)
                if args.submit:
                    print("Submitting", new_xml_filename)
                    submit_xml(new_xml_filename, local=args.local, grid_control=args.gc)
                elif args.resubmit:
                    print("Re-Submitting", new_xml_filename)
                    resubmit_xml(new_xml_filename, local=args.local, grid_control=args.gc)
                elif args.check:
                    print("Checking", new_xml_filename)
                    check_xml(new_xml_filename)
                elif args.hadd:
                        cmd = hadd_cmds[xml_filename]['cmd']
                        sample_dir = hadd_cmds[xml_filename]['sample_dir']
                        if cmd is None:
                            continue
                        workdir = os.path.abspath("/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/%s/%s" % (sample_dir, new_workdir))
                        conf_filename = hadd_job_gc(workdir=workdir, sample=sample_dir, xml_name=new_xml_filename, cmd=cmd)
                        hadd_filenames.append(conf_filename)

            # For now, do one systematic shift at a time
            for syst in systematics:
                is_data = "data" in xml_filename.lower()
                if (is_data and not syst.onData) or (not is_data and not syst.onMC):
                    continue

                # generate all combos of systematics, and do each one
                for syst_values in syst.values:
                    if isinstance(syst_values, str):
                        # important - handle case where single entry in list gets treated as str
                        # instead of list - turn into list
                        syst_values = [syst_values]
                    name = "_".join(["%s_%s" % (n, v) for n, v in zip(syst.names, syst_values)])
                    new_xml_filename = xml_stem + '_' + radius + '_' + name + '.xml'

                    append = "_".join(["%s%s" % (n, v.capitalize()) for n, v in zip(syst.names, syst_values)])
                    new_workdir = write_updated_file(contents, new_xml_filename, radius, systematic_names=syst.names, systematic_values=syst_values, workdir_append=append)
                    if args.submit:
                        print("Submitting", new_xml_filename)
                        submit_xml(new_xml_filename, local=args.local, grid_control=args.gc)
                    elif args.resubmit:
                        print("Re-Submitting", new_xml_filename)
                        resubmit_xml(new_xml_filename, local=args.local, grid_control=args.gc)
                    elif args.check:
                        print("Checking", new_xml_filename)
                        check_xml(new_xml_filename)
                    elif args.hadd:
                        cmd = hadd_cmds[xml_filename]['cmd']
                        sample_dir = hadd_cmds[xml_filename]['sample_dir']
                        if cmd is None:
                            continue
                        workdir = os.path.abspath("/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/%s/%s" % (sample_dir, new_workdir))
                        conf_filename = hadd_job_gc(workdir=workdir, sample=sample_dir, xml_name=new_xml_filename, cmd=cmd)
                        hadd_filenames.append(conf_filename)

    if args.hadd:
        # file to easily check all hadd jobs
        hadd_script = "hadd-check-%s.sh" % strftime("%Y%m%d-%H%M%S")
        with open(hadd_script, "w") as f:
            for conf_name in hadd_filenames:
                f.write("grid-control -c %s\n" % conf_name)
        os.chmod(hadd_script,
            stat.S_IRUSR |
            stat.S_IWUSR |
            stat.S_IXUSR |
            stat.S_IRGRP |
            stat.S_IROTH )
        print("Check all hadd jobs:")
        print("./%s" % hadd_script)

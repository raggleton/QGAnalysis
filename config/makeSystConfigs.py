#!/usr/bin/env python

"""Create copies of XML configs with different jet radii and systematics varied"""

import os
import sys
import argparse
import re
from collections import namedtuple
import subprocess
from itertools import product
from time import sleep


# structure to hold info about a given systematic
Systematic = namedtuple("Systematic", ['names', 'values', 'onData', 'onMC'])


def get_workdir(line):
    match = re.search(r'Workdir="(.*)"/>', line)
    if not match:
        raise RuntimeError("Couldn't find Workdir")
    return match.group(1)


def write_updated_file(contents, new_xml_filename, radius, systematic_names=None, systematic_values=None, workdir_append=""):
    # Should really generalise this
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
                if systematic_names:
                    for key in ["DO_PU_BINNED_HISTS", "DO_FLAVOUR_HISTS", "DO_KINEMATIC_HISTS"]:
                        new_line = '            <Item Name="%s" Value="False"/>\n' % (key)
                        f.write(new_line)

                new_line = '            <Item Name="DO_LAMBDA_HISTS" Value="True"/>\n'
                new_line += '            <Item Name="DO_UNFOLD_HISTS" Value="True"/>\n'
                f.write(new_line)

                # Write out systematics
                if systematic_values:
                    for sn, val in zip(systematic_names, systematic_values):
                        new_line = '            <Item Name="%s" Value="%s"/>\n' % (sn, val)
                        f.write(new_line)

            elif systematic_names and '<!--' not in line and "<Item Name=" in line:
                # remove any existing entries for this systematic
                if any(['<Item Name="%s"' % (sn) in line for sn in systematic_names]):
                    continue
                elif "JackknifeVariations" in line:
                    # turn off Jackknife for systematics
                    f.write('            <Item Name="JackknifeVariations" Value="False"/>\n')
                else:
                    f.write(line)

            else:
                f.write(line)


def submit_xml(xml_filename, local=False, el7_worker=False):
    local_opt = "--local" if local else ""
    worker_opt = "--el7worker" if el7_worker else ""
    cmd = 'sframe_batch.py -s %s %s %s' % (local_opt, worker_opt, xml_filename)
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        print(result)
    except Exception as err:
        print(err)
        print(err.output)
    sleep(5)


def resubmit_xml(xml_filename, local=False, el7_worker=False):
    local_opt = "--local" if local else ""
    worker_opt = "--el7worker" if el7_worker else ""
    cmd = 'sframe_batch.py -r %s %s %s' % (local_opt, worker_opt, xml_filename)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument( "-s", "--submit", action='store_true', help='submit on sframe_batch')
    group.add_argument("-r", "--resubmit", action='store_true', help='resubmit on sframe_batch')
    group.add_argument("-c", "--check", action='store_true', help='check with sframe_batch')
    parser.add_argument("--onlySubmitSystematics", action='store_true', help='Only submit systematic variation jobs, not nominal')
    parser.add_argument("--local", action='store_true', help='Run locally')
    parser.add_argument("--el7worker", action='store_true', help='Run on EL7 worker nodes')
    args = parser.parse_args()

    base_files = [
        # 'QGAnalysisHerwig.xml',
        # 'QGAnalysisPythia.xml',
        # 'QGAnalysisMGPythiaQCD.xml',
        'QGAnalysisMGPythiaDY.xml',
        'QGAnalysisMGPythiaDYInclHTMax70.xml',
        'QGAnalysisMGPythiaDYIncl.xml',

        # 'QGAnalysisDataJetHT.xml',
        # 'QGAnalysisDataZeroBias.xml',
        # 'QGAnalysisDataSingleMu.xml',
    ]

    # These are entities that are used
    radii = [
        "AK4PUPPI",
        # "AK8PUPPI",
    ]

    # Systematics to be changed
    # Note that the original XML should have the "nominal" settings if needed
    # - here are just the variations
    # The values are the possible settings for the names,
    # and each entry in values should match the length of the names
    # - they get done zip-qise, not product-wise
    systematics = [
        Systematic(names=['chargedHadronShift'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        Systematic(names=['neutralHadronShift'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        Systematic(names=['photonShift'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        Systematic(names=['jecsmear_direction'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        Systematic(names=['jersmear_direction'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        Systematic(names=['track_direction'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        Systematic(names=['pileup_direction'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        Systematic(names=['PDFvariations'], onMC=True, onData=False,
                   values=[ ('true') ]),
        Systematic(names=['ScaleVariationMuR', 'ScaleVariationMuF'], onMC=True, onData=False,
                   values=[('up', 'up'),
                           ('nominal', 'up'),
                           ('up', 'nominal'),
                           ('down', 'down'),
                           ('down', 'nominal'),
                           ('nominal', 'down'),
                           ]),
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

    for xml_filename in base_files:
        with open(xml_filename) as f:
            contents = f.readlines()

        xml_stem = os.path.splitext(xml_filename)[0]

        for radius in radii:
            # no systematics
            new_xml_filename = xml_stem + '_' + radius + '.xml'

            # if radius is not "AK4PUPPI":
            write_updated_file(contents, new_xml_filename, radius)
            if not args.onlySubmitSystematics:
                if args.submit:
                    print "Submitting", new_xml_filename
                    submit_xml(new_xml_filename, local=args.local, el7_worker=args.el7worker)
                elif args.resubmit:
                    print "Re-Submitting", new_xml_filename
                    resubmit_xml(new_xml_filename, local=args.local, el7_worker=args.el7worker)
                elif args.check:
                    print "Checking", new_xml_filename
                    check_xml(new_xml_filename)

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
                    write_updated_file(contents, new_xml_filename, radius, systematic_names=syst.names, systematic_values=syst_values, workdir_append=append)
                    if args.submit:
                        print "Submitting", new_xml_filename
                        submit_xml(new_xml_filename, local=args.local, el7_worker=args.el7worker)
                    elif args.resubmit:
                        print "Re-Submitting", new_xml_filename
                        resubmit_xml(new_xml_filename, local=args.local, el7_worker=args.el7worker)
                    elif args.check:
                        print "Checking", new_xml_filename
                        check_xml(new_xml_filename)

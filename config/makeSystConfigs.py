#!/usr/bin/env python

"""Create copies of XML configs with different jet radii and systematics varied"""

import os
import sys
import argparse
import re
from collections import namedtuple
import subprocess


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
                f.write(new_line)

            elif "&AK4PUPPI;" in line:  # cos the default is AK4PUPPI
                f.write(line.replace("AK4PUPPI", radius))

            elif systematic_names and "<!--@SYST-->" in line:
                for sn, val in zip(systematic_names, systematic_values):
                    new_line = '            <Item Name="%s" Value="%s"/>\n' % (sn, val)
                    f.write(new_line)

            elif systematic_names and '<!--' not in line and "<Item Name=" in line:
                # remove any existing entries for this systematic
                for sn in systematic_names:
                    if '<Item Name="%s"' % (sn) in line:
                        continue
            else:
                f.write(line)


def submit_xml(xml_filename, local=False):
    local_opt = "--local" if local else ""
    cmd = 'sframe_batch.py -s %s %s' % (local_opt, xml_filename)
    subprocess.check_output(cmd, shell=True)

def resubmit_xml(xml_filename, local=False):
    local_opt = "--local" if local else ""
    cmd = 'sframe_batch.py -r %s %s' % (local_opt, xml_filename)
    subprocess.Popen(cmd, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--submit", action='store_true', help='submit on sframe_batch')
    group.add_argument("--resubmit", action='store_true', help='resubmit on sframe_batch')
    parser.add_argument("--onlySubmitSystematics", action='store_true', help='Only submit systematic variation jobs, not nominal')
    parser.add_argument("--local", action='store_true', help='Run locally')
    args = parser.parse_args()

    base_files = [
        'QGAnalysisHerwig.xml',
        # 'QGAnalysisPythia.xml',
        'QGAnalysisMGPythia.xml',
        'QGAnalysisMGPythiaDYIncl.xml',

        'QGAnalysisDataJetHT.xml',
        'QGAnalysisDataZeroBias.xml',
        'QGAnalysisDataSingleMu.xml',
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
    systematics = [
        Systematic(names=['chargedHadronShift'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        # Systematic(names=['neutralHadronShift'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        # Systematic(names=['photonShift'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        Systematic(names=['jecsmear_direction'], onMC=True, onData=False,
                   values=[('up'), ('down')]),
        # Systematic(names=['jersmear_direction'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        # Systematic(names=['pileup_direction'], onMC=True, onData=False,
        #            values=[('up'), ('down')]),
        # Systematic(names=['PDFvariations'], onMC=True, onData=False,
        #            values=[ ('true') ]),
        # Systematic(names=['ScaleVariationMuR', 'ScaleVariationMuF'], onMC=True, onData=False,
        #            values=[('up', 'up'),
        #                    ('nominal', 'up'),
        #                    ('up', 'nominal'),
        #                    ('down', 'down'),
        #                    ('down', 'nominal'),
        #                    ('nominal', 'down'),
                           # ]),
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
                    submit_xml(new_xml_filename, local=args.local)
                elif args.resubmit:
                    print "Re-Submitting", new_xml_filename
                    resubmit_xml(new_xml_filename, local=args.local)

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
                    write_updated_file(contents, new_xml_filename, radius, syst.names, syst_values, append)
                    if args.submit:
                        print "Submitting", new_xml_filename
                        submit_xml(new_xml_filename, local=args.local)
                    elif args.resubmit:
                        print "Re-Submitting", new_xml_filename
                        resubmit_xml(new_xml_filename, local=args.local)

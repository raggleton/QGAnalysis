#!/usr/bin/env python

"""Create copies of XML configs with different jet radii and systematics varied"""

import os
import sys
import argparse
import re
from collections import namedtuple
import subprocess


# structure to hold info about a given systematic
Systematic = namedtuple("Systematic", ['name', 'values', 'onData', 'onMC'])


def get_workdir(line):
    match = re.search(r'Workdir="(.*)"/>', line)
    if not match:
        raise RuntimeError("Couldn't find Workdir")
    return match.group(1)


def write_updated_file(contents, new_xml_filename, jet_config, systematic_name=None, systematic_value=None):
    # Should really generalise this
    with open(new_xml_filename, 'w') as f:
        for line in contents:
            if "Workdir=" in line:
                # update workdir
                orig_workdir = get_workdir(line)
                new_workdir = orig_workdir.replace("_ak4puppi_", "_%s_" % (radius.lower()))
                if systematic_name:
                    new_workdir += "_" + systematic_name + systematic_value.capitalize()
                new_line = line.replace(orig_workdir, new_workdir)
                f.write(new_line)
            elif "&AK4PUPPI;" in line:  # cos the default is AK4PUPPI
                f.write(line.replace("AK4PUPPI", radius))
            elif systematic_name and "<!--@SYST-->" in line:
                new_line = '            <Item Name="%s" Value="%s"/>\n' % (systematic_name, systematic_value)
                f.write(new_line)
            elif systematic_name and '<Item Name="%s"' % (systematic_name) in line and '<!--' not in line:
                # remove any existing entries for this systematic
                continue
            else:
                f.write(line)


def submit_xml(xml_filename):
    cmd = 'sframe_batch.py -s %s' % (xml_filename)
    subprocess.check_call(cmd.split())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submit", action='store_true', help='submit on sframe_batch')
    args = parser.parse_args()

    base_files = [
        'QGAnalysisHerwig.xml',
        'QGAnalysisPythia.xml',
        'QGAnalysisMGPythia.xml',
        'QGAnalysisMGPythiaDYIncl.xml',

        'QGAnalysisDataJetHT.xml',
        'QGAnalysisDataZeroBias.xml',
        'QGAnalysisDataSingleMu.xml',
    ]

    # These are entities that are used
    radii = [
        "AK4PUPPI",
        "AK8PUPPI",
    ]

    systematics = [
        # Systematic(name='neutralHadronShift', onMC=True, onData=False,
        #            values=['up', 'down']),
        Systematic(name='photonShift', onMC=True, onData=False,
                   values=['up', 'down']),
        # Systematic(name='jecsmear_direction', onMC=True, onData=False,
        #            values=['up', 'down']),
        # Systematic(name='jersmear_direction', onMC=True, onData=False,
        #            values=['up', 'down']),
        # Systematic(name='pileup_direction', onMC=True, onData=False,
        #            values=['up', 'down']),
        # Systematic(name='ScaleVariationMuR', onMC=True, onData=False,
        #            values=['up', 'down']),
        # Systematic(name='ScaleVariationMuF', onMC=True, onData=False,
        #            values=['up', 'down']),
    ]

    for xml_filename in base_files:
        with open(xml_filename) as f:
            contents = f.readlines()

        xml_stem = os.path.splitext(xml_filename)[0]

        for radius in radii:
            # no systematics
            new_xml_filename = xml_stem + '_' + radius + '.xml'

            write_updated_file(contents, new_xml_filename, radius, None, None)
            if args.submit:
                submit_xml(new_xml_filename)

            # For now, do one systematic shift at a time
            for syst in systematics:
                for variation in syst.values:

                    is_data = "data" in xml_filename.lower()
                    if (is_data and not syst.onData) or (not is_data and not syst.onMC):
                        continue

                    syst_name = syst.name
                    new_xml_filename = xml_stem + '_' + radius + '_'+syst_name+'_'+variation+'.xml'

                    write_updated_file(contents, new_xml_filename, radius, syst.name, variation)
                    if args.submit:
                        submit_xml(new_xml_filename)

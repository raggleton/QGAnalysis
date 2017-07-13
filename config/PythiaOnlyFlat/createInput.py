#!/usr/bin/env python


"""Create input XML file from all ROOT ntuples in a dir"""


import xml.etree.ElementTree as ET
import os
import HTMLParser


ntuple_dir = "/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Ntuples/PythiaOnly/FlatQCD"

# For all XML input snippets:
common_attr = dict(NEventsMax="&NEVT;", Type="MC", Cacheable="False", Lumi= '1.0')


def create_xml_snippet(data_attrib, in_attribs):
    """Create an XML snippet for this sample dictionary"""
    root = ET.Element("InputData", attrib=data_attrib)
    for ia in in_attribs:
        in_el = ET.SubElement(root, "In", attrib=ia)
    in_tree = ET.SubElement(root, "InputTree", attrib={"Name": "AnalysisTree"})
    # need to convert &amp; back to & etc
    contents = HTMLParser.HTMLParser().unescape(ET.tostring(root))
    return contents


if __name__ == "__main__":
    xml_snippets = []

    for f in os.listdir(ntuple_dir):
        # if not os.path.isfile(f):
        #     continue

        data_attrib = {
            'Version': f.replace("Ntuple_", "").replace(".root", ""),
        }
        data_attrib.update(common_attr)
        in_attrib = {"FileName": os.path.join(ntuple_dir, f), "Lumi": "0.0"}
        snip = create_xml_snippet(data_attrib, [in_attrib])
        xml_snippets.append(snip)

    xml_filename = "input.xml"

    with open(xml_filename, "w") as f:
        f.write("\n".join(xml_snippets))

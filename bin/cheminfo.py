#!/usr/bin/env python

import sys
import urllib.request, json


compound = sys.argv[1]

query = "https://data.rcsb.org/rest/v1/core/chemcomp/%s" % compound
imgurl_base = "http://hulab.rxnfinder.org/smi2img/"

with urllib.request.urlopen(query) as url:
    data = json.load(url)
    json_formatted_str = json.dumps(data, indent=2)
    # print(json_formatted_str)
    print("Name: %s" % data["chem_comp"]["name"])
    print("Formula: %s" % data["chem_comp"]["formula"])
    print("Molecular Weight: %s" % data["chem_comp"]["formula_weight"])
    for entry in data["pdbx_chem_comp_descriptor"]:
        if entry["program"] == "OpenEye OEToolkits" and entry["type"] == "SMILES_CANONICAL":
            smiles = entry["descriptor"]
            print("OpenEyeo SMILES: %s" % smiles)
            print("Molecule image: %s%s/" % (imgurl_base, smiles))
            break





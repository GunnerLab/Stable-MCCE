#!/usr/bin/env python

import sys

CONVERT = {"1H   HOH": " H1  HOH",
           "2H   HOH": " H2  HOH",
           "1HG1 VAL": "HG11 VAL",
           "2HG1 VAL": "HG12 VAL",
           "3HG1 VAL": "HG13 VAL",
           "1HG2 VAL": "HG21 VAL",
           "2HG2 VAL": "HG22 VAL",
           "3HG2 VAL": "HG23 VAL",
           "1HA  GLY": " HA2 GLY",
           "2HA  GLY": " HA3 GLY",
           "1HB  ALA": " HB1 ALA",
           "2HB  ALA": " HB2 ALA",
           "3HB  ALA": " HB3 ALA",
           "1HB  LEU": " HB2 LEU",
           "2HB  LEU": " HB3 LEU",
           "1HD1 LEU": "HD11 LEU",
           "2HD1 LEU": "HD12 LEU",
           "3HD1 LEU": "HD13 LEU",
           "1HD2 LEU": "HD21 LEU",
           "2HD2 LEU": "HD22 LEU",
           "3HD2 LEU": "HD23 LEU",
           "1HB  TRP": " HB2 TRP",
           "2HB  TRP": " HB3 TRP",
           }

if __name__ == "__main__":

    lines = open(sys.argv[1]).readlines()
    for line in lines:
        if line[:6] == "ATOM  " or line[:6] == "HETATM":
            if line[12:20] in CONVERT:
                print(line[:12]+CONVERT[line[12:20]]+line[20:].rstrip())
            else:
                print(line, end="")
        else:
            print(line, end="")


#!/usr/bin/env python

from mfe import *
import pandas as pd
import numpy as np

def get_resids(fname):
    resid_pkas = []
    lines = open(fname).readlines()
    lines.pop(0)
    for line in lines:
        fields = line.split()
        if len(fields) >= 2:
            resid = fields[0]
            try:
                pka = float(fields[1])
            except ValueError:
                pka = None

            resid_pkas.append((resid, pka))

    return resid_pkas

if __name__ == "__main__":
    # Get the command arguments
    helpmsg = "Calculate mean field energy ionization energy on all residues at a specific pH/eH."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-p", metavar="pH/Eh", default="m",
                        help="pH or Eh mfe analysis is carried out, or \'m\' for midpoint")
    parser.add_argument("-x", metavar="TS_correction", default="r",
                        help="f: False, t: True, or r: determined by run.prm (default)")
    parser.add_argument("-f", metavar="format", default="xls", help="Output format: xls (default), csv, or txt")
    parser.add_argument("-o", metavar="filename", default="mfe_all", help="Output file name")

    args = parser.parse_args()

    args.c = -0.001
    resid_pkas = get_resids("pK.out")


    # Index column
    index = {"Terms": ["vdw0", "vdw1", "tors", "ebkb", "dsol", "offset", "pH&pK0", "Eh&Em0", "-TS", "residues", "TOTAL"]}
    df = pd.DataFrame(index)
    for resid_pka in resid_pkas:
        resid, pka = resid_pka
        args.residue = [resid]
        report = mfe(args)
        names = {"Terms": []}
        analysis = {resid: []}
        if len(report) <= 1:  # out of range
            df[resid] = analysis
        else:
            for line in report:
                fields = line.split()
                if len(fields) > 2:
                    names["Terms"].append(fields[0])
                    analysis[resid].append(float(fields[1]))
            if index["Terms"][:11] == names["Terms"][:11]:
                print("OK")
            else:
                print("Mismatch")

        #break
        #add_output(report)

    print(df)

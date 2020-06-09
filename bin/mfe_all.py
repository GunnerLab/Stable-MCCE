#!/usr/bin/env python

from mfe import *
import pandas as pd
import numpy as np
import os

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
    parser.add_argument("-f", metavar="format", default="xls", help="Output format: xls (default), csv, or html")
    parser.add_argument("-o", metavar="filename", default="mfe_all", help="Output file name")

    args = parser.parse_args()

    args.c = -0.001
    resid_pkas = get_resids("pK.out")


    # Index column
    index = ["Terms", "pKa/Eh", "vdw0", "vdw1", "tors", "ebkb", "dsol", "offset", "pH&pK0", "Eh&Em0", "-TS", "residues", "TOTAL"]
    current_len_of_index = len(index)
    res_matrix = []
    for resid_pka in resid_pkas:
        resid, pka = resid_pka
        args.residue = [resid]
        report = mfe(args)
        res_column_index = ["Terms", "pKa/Eh"]
        if pka:
            res_column_value = [resid, float(pka)]
        else:
            res_column_value = [resid, ""]

        if len(report) > 1:  # not "titration out of range titration"
            report = report[4:]
            for line in report:
                fields = line.split()
                if len(fields) > 2:
                    #print(fields)
                    res_column_index.append(fields[0])
                    res_column_value.append(float(fields[1]))
        if len(res_column_index) > current_len_of_index:
            current_len_of_index = len(res_column_index)
            index = res_column_index

        res_matrix.append(res_column_value)

    # Now we have a table, with some cells missing, due to midpoint not specified, we will complete the table with the longest column
    c1 = {index[0]: index[1:]}
    df = pd.DataFrame(c1)
    for column in res_matrix:
        if len(column) < current_len_of_index:  # complete the empty cells with empty strings
            tail = [""] * (current_len_of_index - 2)
            column.extend(tail)
        df[column[0]] = column[1:]

    # print(df)

    if args.f == "csv":
        fname = args.o
        filename, extension = os.path.splitext(fname)
        if extension == "csv":
            fname = "%s.csv" % filename
        else:
            fname = "%s.csv" % fname

        df.to_csv(fname, float_format="%.3f")
        print("MFE output saved in file %s" % fname)
    elif args.f == "html":
        fname = args.o
        filename, extension = os.path.splitext(fname)
        if extension == "html":
            fname = "%s.html" % filename
        else:
            fname = "%s.html" % fname

        df.to_html(fname, float_format="%.3f")
        print("MFE output saved in file %s" % fname)
    elif args.f == "xls" or args.f == "xlsx":
        fname = args.o
        filename, extension = os.path.splitext(fname)
        if extension == "xls" or extension == "xlsx":
            fname = "%s.xlsx" % filename
        else:
            fname = "%s.xlsx" % fname

        df.to_excel(fname, float_format="%.3f")
        print("MFE output saved in file %s" % fname)
    else:
        print(df)

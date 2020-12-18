#!/usr/bin/env python
# Subtract two Monte Carlo sampling result reports. The reports can be fort.38, sum_crg.out, or anything with matching columns.

import os, argparse

def read_matrix(fname):
    "Read a file into a matrix, str as str, numbers to float"
    result = []
    lines = open(fname).readlines()
    for line in lines:
        fields = line.strip().split()
        ldata = []
        for f in fields:
            try:
                ldata.append(float(f))
            except ValueError:
                ldata.append(f)
        result.append(ldata)

    return result


def msubtract(m1, m2):
    nrow_m1 = len(m1)
    nrow_m2 = len(m2)
    result = []
    if nrow_m1 != nrow_m2:
        print("Rows in two data files don't match")
        return None

    result = [m1[0]]
    for irow in range(1, nrow_m1):
        if len(m1[irow]) != len(m2[irow]):
            print("Columns %s in two data files don't match" % m1[irow][0])
            return None
        row_result = []
        for icol in range(len(m1[irow])):
            if isinstance(m1[irow][icol], str):  # keep
                row_result.append(m1[irow][icol])
            else:
                row_result.append(m1[irow][icol] - m2[irow][icol])
        result.append(row_result)

    return result


def diff(files):
    data1 = read_matrix(files[0])
    data2 = read_matrix(files[1])
    data_diff = msubtract(data1, data2)

    for row_data in data_diff:
        print("%14s %s" % (row_data[0], " ".join(["%6.3f" % x for x in row_data[1:]])))
    return

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Subtract two sum_crg.out or fort.38 files."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("files", metavar="file1 file2", nargs=2)
    args = parser.parse_args()

    diff(args.files)

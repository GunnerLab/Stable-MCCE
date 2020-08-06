#!/usr/bin/env python
import argparse

def cutnmr(fname):
    lines = [x for x in open(fname).readlines() if x[:6] == "ATOM  " or x[:6] == "HETATM" or x[:6] == "MODEL " or x[:6] == "ENDMDL"]

    # collect models
    models = []

    for i_line in range(len(lines)):
        line = lines[i_line]
        if line[:6] == "MODEL ":
            _, name = line.split()[:2]
            i_start = i_line
        if line[:6] == "ENDMDL":
            i_end = i_line
            models.append((int(name), i_start, i_end))

    if models:
        for model in models:
            outname = "model%02d.pdb" % model[0]
            outlines = lines[model[1]: model[2]+1]
            open(outname, "w").writelines(outlines)
    else:
        print("No NMR models detected.")


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Cut pdb file NMR models into model01.pdb, model02.pdb ..."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("prot", metavar="pdb", nargs=1)
    args = parser.parse_args()

    inputfile = args.prot[0]
    cutnmr(inputfile)

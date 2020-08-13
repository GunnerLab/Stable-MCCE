#!/usr/bin/env python
import argparse

def missing_ftpl(line):
    fields = line.split()
    missing = fields[-1]
    print("Warning: The residue \"%s\" in pdb file doesn't have parameter file. If this resuisue has charge, you should prepare a ftpl file." % (missing))
    return True


def asNTR(line):
    w = False
    # these residues can be termini
    aminoacids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    # Labeling "LYS A   1" as NTR
    fields = line.split('"')
    res, _, _ = fields[1].split()
    if res not in aminoacids:
        print("Error: Residue \"%s\" was mistakenly as NTR. Report this bug to developers." % res)
        w = True

    return w




def analyze(log):
    warning = False
    lines = open(log).readlines()
    for line in lines:
        if line.find("Error! premcce_confname(): Can't get conformer list of this residue") >= 0:
            missing_ftpl(line)
            warning = True
        elif line.find("as NTR") >= 0:
            w = asNTR(line)
            if not warning:
                warning = w
        elif line.find("has non integer charge"):
            w = non_int_crg(line)
            if not warning:
                warning = w

    return warning


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Analyze log file and report potential problems and suggestions."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("log", metavar="log", nargs=1)
    args = parser.parse_args()

    log = args.log[0]


    if not analyze(log):
        print("No warning and errors detected.")
#!/usr/bin/env python
"""
This program set up and run MCCE step 1.

Input:
 * PDB file

Output:
 * run.prm
 * step1_out.pdb
 * list_rot.gold
 * head1.lst

Usage examples:

1. Run step 1
    step1.py prot.pdb

2. Write run.prm for step 1, do not actually run step 1. (No run)
    step1.py prot.pdb --norun

3. Run step 1 with specified mcce
    step1.py prot.pdb -e /path/to/mcce

4. Run step 1 with other customized parameters
    step1.py prot.pdb -u HOME_MCCE=/path/to/mcce_home,H2O_SASCUTOFF=0.05
"""

import os, argparse, shutil
import subprocess
from mccesteps import *

def fix_format(fname):
    # make sure NT atoms appear before CTR atoms.
    # truncate ending special characters from Windows.
    NTR_atoms = [" CA ", " N  ", " HA ", " H  ", " H2 ", " H3 "]
    CTR_atoms = [" C  ", " O  ", " HO ", " OXT", " HXT"]

    pdblines = open(fname).readlines()

    new_pdblines = []
    cur_resid = ""
    ntr_atoms = []
    mid_atoms = []
    ctr_atoms = []
    for line in pdblines:
        line = line.rstrip() + "\n"   # remove any special char and replace it by Linux end of line
        if line[:6] == "ATOM  " or line[:6] == "HETATM":
            # Remove alternative location
            if line[16] == "A":
                line = line[:16]+" "+line[17:]
            elif line[16] != " ":
                continue

            resid = line[17:27]
            if resid != cur_resid:
                new_pdblines += ntr_atoms
                new_pdblines += mid_atoms
                new_pdblines += ctr_atoms
                ntr_atoms = []
                mid_atoms = []
                ctr_atoms = []
                cur_resid = resid
            if line[12:16] in NTR_atoms:
                ntr_atoms.append(line)
            elif line[12:16] in CTR_atoms:
                ctr_atoms.append(line)
            else:
                mid_atoms.append(line)
        else:
            new_pdblines.append(line)

    # last one
    new_pdblines += ntr_atoms
    new_pdblines += mid_atoms
    new_pdblines += ctr_atoms

    return new_pdblines


def write_runprm(args):
    runprm = {}

    path = str(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.dirname(path)
    #print(base_path)
    runprm["INPDB"] = "step0_out.pdb"
    runprm["DO_PREMCCE"] = "t"
    runprm["MCCE_HOME"] = base_path
    runprm["MINIMIZE_SIZE"] = "t"
    if args.noter:
        runprm["TERMINALS"] = "f"
    else:
        runprm["TERMINALS"] = "t"
    runprm["CLASH_DISTANCE"] = "2.0"
    runprm["H2O_SASCUTOFF"] = "0.05"
    runprm["IGNORE_INPUT_H"] = "t"
    runprm["RENAME_RULES"] = "%s/name.txt" % base_path
    runprm["EPSILON_PROT"] = args.d

    if args.dry:
        runprm["H2O_SASCUTOFF"] = "-0.01"

    if args.u:
        fields = args.u.split(",")
        for field in fields:
            try:
                key, value = field.split("=")
                runprm[key] = value
            except ValueError:
                print("Argument must be \"KEY=VALUE\" format, but got \"%s\" instead" % field)


    # write run.prm
    export_runprm(runprm)
    record_runprm(runprm, "#STEP1")
    return


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Run mcce step 1, premcce to format PDB file to MCCE PDB format."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("--norun", default=False, help="Create run.prm but do not run step 1", action="store_true")
    parser.add_argument("--noter", default=False, help="Do not label terminal residues (for making ftpl).", action="store_true")
    parser.add_argument("-e", metavar="/path/to/mcce", default="mcce", help="mcce executable location, default to \"mcce\"")
    parser.add_argument("-u", metavar="Key=Value", default="", help="User customized variables")
    parser.add_argument("-d", metavar="epsilon", default="4.0", help="protein dielectric constant for delphi, default to 4.0")
    parser.add_argument("--dry", default=False, help="Delete all water molecules.", action="store_true")
    parser.add_argument("prot", metavar="pdb", nargs=1)
    args = parser.parse_args()

    #print(args)

    write_runprm(args)
    if not args.norun:
        new_pdblines = fix_format(args.prot[0])
        open("step0_out.pdb", "w").writelines(new_pdblines)

        mcce = args.e
        process = subprocess.Popen(["mcce"], stdout=subprocess.PIPE)
        for line in process.stdout:
            print("%s" % line.decode(), end="")

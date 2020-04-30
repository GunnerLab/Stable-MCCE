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

2. Write run.prm for step 1, do not actually run step 1. (Dry run)
    step1.py prot.pdb --dry

3. Run step 1 with specified mcce
    step1.py prot.pdb -e /path/to/mcce

4. Run step 1 with other customized parameters
    step1.py prot.pdb -u HOME_MCCE=/path/to/mcce_home,H2O_SASCUTOFF=0.05
"""

import os, argparse, shutil
import subprocess

def write_runprm(args):
    runprm = {}

    path = str(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.dirname(path)
    #print(base_path)
    runprm["INPDB"] = args.prot[0]
    runprm["DO_PREMCCE"] = "t"
    runprm["MCCE_HOME"] = base_path
    runprm["MINIMIZE_SIZE"] = "t"
    runprm["TERMINALS"] = "t"
    runprm["CLASH_DISTANCE"] = "2.0"
    runprm["H2O_SASCUTOFF"] = "0.05"
    runprm["IGNORE_INPUT_H"] = "t"
    runprm["RENAME_RULES"] = "%s/name.txt" % base_path

    if args.u:
        fields = args.u.split(",")
        for field in fields:
            try:
                key, value = field.split("=")
                runprm[key] = value
            except ValueError:
                print("Argument must be \"KEY=VALUE\" format, but got \"%s\" instead" % field)


    # write
    lines = []
    for key in runprm:
        line = "%-20s    (%s)\n" % (runprm[key], key)
        lines.append(line)
    open("run.prm", "w").writelines(lines)

    return

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Run mcce step 1, premcce to format PDB file to MCCE PDB format."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("--dry", default=False, help="Dry run, create run.prm but do not run step 1", action="store_true")
    parser.add_argument("-e", metavar="/path/to/mcce", default="mcce", help="mcce executable location, default to \"mcce\"")
    parser.add_argument("-u", metavar="Key=Value", default="", help="User customized variables")
    parser.add_argument("prot", metavar="pdb", nargs=1)
    args = parser.parse_args()

    #print(args)

    write_runprm(args)
    if not args.dry:
        mcce = args.e
        process = subprocess.Popen(["mcce"], stdout=subprocess.PIPE)
        for line in process.stdout:
            print("%s" % line.decode(), end="")

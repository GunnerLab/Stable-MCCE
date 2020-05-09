#!/usr/bin/env python
"""
This program set up and run MCCE step 2.

Input:
 * step1_out.pdb
 * list_rot.gold

Output:
 * step2_out.pdb

Usage examples:

1. Run step 1 with default (quick conformers)
    step2.py

2. Write run.prm for step 2, do not actually run step 2.
    step2.py --norun

3. Run step 2 with quick (default), medium, and comprehensive conformer making
    step2.py -l 1   # quick
    step2.py -l 2   # medium
    step2.py -l 3   # comprehensive

4. Run step 2 using specific mcce executable
    step2.py -e /path/to/mcce

5. Run step 2 using specific mcce executable
    step2.py -d 4.0

6. Run step 2 with other customized parameters
    step1.py -u HOME_MCCE=/path/to/mcce_home,H2O_SASCUTOFF=0.05
"""

import os, argparse
import subprocess
import sys
from mccesteps import *

def write_runprm(args):
    runprm = {}

    path = str(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.dirname(path)
    #print(base_path)
    runprm["DO_ROTAMERS"] = "t"
    runprm["MCCE_HOME"] = base_path
    runprm["EXTRA"] = "%s/extra.tpl" % base_path
    runprm["EPSILON_PROT"] = args.d
    runprm["ROT_SPECIF"] = "f"
    runprm["ROT_SWAP"] = "t"
    runprm["REPACKS"] = "5000"
    runprm["REPACK_CUTOFF"] = "0.01"
    runprm["VDW_CUTOFF"] = "10.00"
    runprm["HDIRECTED"] = "t"
    runprm["HDIRDIFF"] = "1.0"
    runprm["HDIRLIMT"] = "36"
    runprm["RELAX_H"] = "t"
    runprm["RELAX_E_THR"] = "-1.0"
    runprm["RELAX_NSTATES"] = "100"
    runprm["RELAX_CLASH_THR"] = "5.0"
    runprm["RELAX_PHI"] = "1.0"
    runprm["RELAX_NITER"] = "300"
    runprm["RELAX_TORQ_THR"] = "0.5"

    runprm["NCONF_LIMIT"] = "999"

    runprm["MINIMIZE_SIZE"] = "t"


    if args.l == 1:
        runprm["PACK"] = "f"
        runprm["ROTATIONS"] = "0"
        runprm["SAS_CUTOFF"] = "1.00"

        runprm["HV_RELAX_NCYCLE"] = "0"
        runprm["HV_RELAX_DT"] = "4"
        runprm["HV_RELAX_NITER"] = "100"
        runprm["HV_RELAX_VDW_THR"] = "2.0"
        runprm["HV_RELAX_HV_VDW_THR"] = "10.0"
        runprm["HV_TORS_SCALE"] = "20.0"
        runprm["HV_RELAX_N_SHAKE"] = "10000"
        runprm["HV_RELAX_CONSTRAINT"] = "1.0"
        runprm["HV_RELAX_CONSTRAINT_FRC"] = "20.0"
        runprm["HV_RELAX_ELEC_THR"] = "-2.0"
        runprm["HV_RELAX_ELEC_CRG_THR"] = "0.3"
        runprm["HV_RELAX_ELEC_DIST_THR"] = "2.0"

        runprm["RELAX_N_HYD"] = "36"

        runprm["PRUNE_THR"] = "0.01"
        runprm["PRUNE_RMSD"] = "2.0"
        runprm["PRUNE_ELE"] = "2.0"
        runprm["PRUNE_VDW"] = "2.0"
    elif args.l == 2:
        runprm["PACK"] = "t"
        runprm["ROTATIONS"] = "6"
        runprm["SAS_CUTOFF"] = "0.2"

        runprm["HV_RELAX_NCYCLE"] = "3"
        runprm["HV_RELAX_DT"] = "4"
        runprm["HV_RELAX_NITER"] = "100"
        runprm["HV_RELAX_VDW_THR"] = "2.0"
        runprm["HV_RELAX_HV_VDW_THR"] = "10.0"
        runprm["HV_TORS_SCALE"] = "20.0"
        runprm["HV_RELAX_N_SHAKE"] = "10000"
        runprm["HV_RELAX_CONSTRAINT"] = "1.0"
        runprm["HV_RELAX_CONSTRAINT_FRC"] = "20.0"
        runprm["HV_RELAX_ELEC_THR"] = "-2.0"
        runprm["HV_RELAX_ELEC_CRG_THR"] = "0.3"
        runprm["HV_RELAX_ELEC_DIST_THR"] = "2.0"

        runprm["RELAX_N_HYD"] = "72"

        runprm["PRUNE_THR"] = "0.02"
        runprm["PRUNE_RMSD"] = "2.0"
        runprm["PRUNE_ELE"] = "2.0"
        runprm["PRUNE_VDW"] = "8.0"
    elif args.l == 3:
        runprm["PACK"] = "t"
        runprm["ROTATIONS"] = "6"
        runprm["SAS_CUTOFF"] = "4.0"

        runprm["SWING"] = "f"
        runprm["PHI_SWING"] = "10.0"

        runprm["RELAX_WAT"] = "t"
        runprm["WATER_RELAX_THR"] = "3.2"

        runprm["HV_RELAX_NCYCLE"] = "10"
        runprm["HV_RELAX_DT"] = "10"
        runprm["HV_RELAX_NITER"] = "100"
        runprm["HV_RELAX_VDW_THR"] = "2.0"
        runprm["HV_RELAX_HV_VDW_THR"] = "5.0"
        runprm["HV_TORS_SCALE"] = "20.0"
        runprm["HV_RELAX_N_SHAKE"] = "10000"
        runprm["HV_RELAX_CONSTRAINT"] = "1.0"
        runprm["HV_RELAX_CONSTRAINT_FRC"] = "20.0"
        runprm["HV_RELAX_ELEC_THR"] = "-1.0"
        runprm["HV_RELAX_ELEC_CRG_THR"] = "0.1"
        runprm["HV_RELAX_ELEC_DIST_THR"] = "3.0"

        runprm["RELAX_N_HYD"] = "72"

        runprm["PRUNE_THR"] = "0.02"
        runprm["PRUNE_RMSD"] = "1.0"
        runprm["PRUNE_ELE"] = "1.0"
        runprm["PRUNE_VDW"] = "8.0"
    else:
        print("This level %d of rotamer making is not supported." % args.l)
        sys.exit()


    if args.u:
        fields = args.u.split(",")
        for field in fields:
            try:
                key, value = field.split("=")
                runprm[key] = value
            except ValueError:
                print("Argument must be \"KEY=VALUE\" format, but got \"%s\" instead" % field)

    # write
    # write run.prm
    export_runprm(runprm)
    record_runprm(runprm, "#STEP2")

    return

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Run mcce step 2, make side chain conformers from step1_out.pdb."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("--norun", default=False, help="Create run.prm but do not run step 2", action="store_true")
    parser.add_argument("-d", metavar="epsilon", default="4.0", help="dielectric constant for optimizing conformers")
    parser.add_argument("-e", metavar="/path/to/mcce", default="mcce", help="mcce executable location, default to \"mcce\"")
    parser.add_argument("-u", metavar="Key=Value", default="", help="User customized variables")
    parser.add_argument("-l", metavar="level", default=1, help="conformer level 1: quick, 2: medium, 3: comnprehensive", type=int)
    args = parser.parse_args()

    #print(args)

    write_runprm(args)
    if not args.norun:
        mcce = args.e
        process = subprocess.Popen([mcce], stdout=subprocess.PIPE)
        for line in process.stdout:
            print("%s" % line.decode(), end="")

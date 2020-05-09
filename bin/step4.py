#!/usr/bin/env python
"""
This program set up and run MCCE step 4.

Input:
 * energies/
 * head3.lst

Output:
 * fort.38

Usage examples:

1. Run step 4 with default (pH titration from 0.0 to 14.0 at 15 points, without entropy correction)
    step4.py

2. Write run.prm for step 4, do not actually run step 4. (Dry run)
    step4.py --dry

3. Run step 4 at defined points
    step4.py -i 0.0 -d 1 -n 15

4. Run step 4 with entropy correction
    step2.py --xts

5. Run step 4 using specific mcce executable
    step2.py -e /path/to/mcce

6. Run step 4 using Eh titration
    step2.py -t eh

7. Run step 4 with other customized parameters
    step1.py -u EXTRA=./extra.tpl
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
    runprm["DO_MONTE"] = "t"
    runprm["EXTRA"] = "%s/extra.tpl" % base_path
    runprm["TITR_TYPE"] = args.t.lower()
    runprm["TITR_PH0"] = args.i
    runprm["TITR_EH0"] = args.i
    runprm["TITR_PHD"] = args.d
    runprm["TITR.EHD"] = args.d
    runprm["TITR_STEPS"] = args.n
    runprm["BIG_PAIRWISE"] = "5.0"
    runprm["MONTE_SEED"] = "-1"
    runprm["MONTE_T"] = "298.15"
    runprm["MONTE_FLIPS"] = "3"
    runprm["MONTE_NSTART"] = "100"
    runprm["MONTE_NEQ"] = "300"
    runprm["MONTE_REDUCE"] = "0.001"
    runprm["MONTE_RUNS"] = "6"
    runprm["MONTE_NITER"] = "2000"
    runprm["MONTE_TRACE"] = "50000"
    runprm["NSTATE_MAX"] = "1000000"
    if args.xts:
        runprm["MONTE_TSX"] = "t"
    else:
        runprm["MONTE_TSX"] = "f"
    runprm["MFE_POINT"] = "f"
    runprm["MFE_CUTOFF"] = "-1.0"
    if args.ms:
        runprm["MS_OUT"] = "t"
    else:
        runprm["MS_OUT"] = "f"


    if args.u:
        fields = args.u.split(",")
        for field in fields:
            try:
                key, value = field.split("=")
                runprm[key] = value
            except ValueError:
                print("Argument must be \"KEY=VALUE\" format, but got \"%s\" instead" % field)

    # write
    export_runprm(runprm)
    record_runprm(runprm, "#STEP4")

    return

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Run mcce step 4, Monte Carlo sampling to simulate a titration."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("--norun", default=False, help="Create run.prm but do not run step 4", action="store_true")
    parser.add_argument("-i", metavar="initial ph/eh", default="0.0", help="Initial pH/Eh of titration")
    parser.add_argument("-d", metavar="interval", default="1.0", help="titration interval in pJ or mV")
    parser.add_argument("-n", metavar="steps", default="15", help="number of steps of titration")
    parser.add_argument("--xts", default=False, help="Enable entropy correction, default is false", action="store_true")
    parser.add_argument("--ms", default=False, help="Enable microstate output", action="store_true")
    parser.add_argument("-e", metavar="/path/to/mcce", default="mcce", help="mcce executable location, default to \"mcce\"")
    parser.add_argument("-t", metavar="ph or eh", default="ph", help="titration type, pH or Eh")
    parser.add_argument("-u", metavar="Key=Value", default="", help="User customized variables")
    args = parser.parse_args()

    #print(args)

    write_runprm(args)
    if not args.norun:
        mcce = args.e
        process = subprocess.Popen([mcce], stdout=subprocess.PIPE)
        for line in process.stdout:
            print("%s" % line.decode(), end="")

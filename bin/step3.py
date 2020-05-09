#!/usr/bin/env python
"""
This program set up and run MCCE step 3 in multiple threads.

Usage examples:

1. Run step 3 with 6 threads
    step3.py -p 6 [-e mcce] [-d delphi]

2. Run step 3 to recreate head3.lst
    step3.py -r [-e mcce] [-d delphi]

3. Run partial conformers. conformer 1 to 100
    step3.py -c 1, 100 [-e mcce] [-d delphi]

-p number: run with number of processes
-r: refresh opp files vdw and ele, and head3.lst
-c start, end: run conformer from start to end
-d: the optional path to delphi program

Assuming the underlying c code will:
1. Write only opp files that are calculated (so that multi-threads won't overwrite each other). Dummies and
un-calculated conformers don't have an opp file.
2. When writing opp file, no previous run is loaded.
3. head3.lst is generated at the end of step 3.


Mechanism:
 * mcce is modified to run part of the conformers, as controlled by run.prm.
 * if an opp file exists, a conformer is considered to be already computed and therefore loaded.
 * no matter how many conformers are calculated, the averaging and vdw terms are always recalculated.
 * by modifying run.prm, we can run partial step 3.
 * a dummy run after all partial runs will update the final head3.lst and vdw terms.
"""

import sys, argparse, shutil
import logging
import threading
import time
import os
import subprocess
from mccesteps import *

def mcce_thread(mcce, name):
    logging.info("Thread %s: starting", name)
    subprocess.run(mcce)
    #time.sleep(name+2)
    logging.info("Thread %s: finishing", name)

    return


def write_runprm(runprm):
    lines = []
    for key in runprm:
        line = "%-20s    (%s)\n" % (runprm[key], key)
        lines.append(line)
    open("run.prm", "w").writelines(lines)
    return


def process(args):
    runprm = {}

    path = str(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.dirname(path)

    runprm["DO_ENERGY"] = "t"
    runprm["PBE_SOLVER"] = "delphi"
    runprm["EPSILON_PROT"] = args.d
    runprm["EXTRA"] = "%s/extra.tpl" % base_path
    runprm["EPSILON_SOLV"] = "80.0"
    runprm["GRIDS_DELPHI"] = "65"
    runprm["GRIDS_PER_ANG"] = "2.0"
    runprm["RADIUS_PROBE"] = "1.4"
    runprm["IONRAD"] = "2.0"
    runprm["SALT"] = "0.15"
    runprm["DELPHI_EXE"] = "%s/bin/delphi" % base_path
    runprm["DELPHI_CLEAN"] = "t"
    runprm["PBE_FOLDER"] = args.f

    if args.u:
        fields = args.u.split(",")
        for field in fields:
            try:
                key, value = field.split("=")
                runprm[key] = value
            except ValueError:
                print("Argument must be \"KEY=VALUE\" format, but got \"%s\" instead" % field)


    fname = "step2_out.pdb"
    lines = open(fname).readlines()
    conformers = []
    for line in lines:
        # residue name, type and chain/sequence/icode/confNumber
        if line[80:82] != "BK":
            uniqID = line[17:20] + line[80:82] + line[21:30]
            if uniqID not in conformers:
                conformers.append(uniqID)
    nconf = len(conformers)

    # write run.prm.record as a single thread run.prm
    runprm["PBE_START"] = "1"
    runprm["PBE_END"] = "99999"
    record_runprm(runprm, "#STEP3")

    if args.norun:
        export_runprm(runprm)
        record_runprm(runprm, "#STEP3")
    else:
        if args.r:
            runprm["PBE_START"] = "1"
            runprm["PBE_END"] = "-11"
            write_runprm(runprm)

            logging.info("Create a single thread to run this job.")
            x = threading.Thread(target=mcce_thread, args=(args.e, 1))
            x.start()
        else:
            nrange = args.c[1] - args.c[0] + 1
            if nconf >= nrange:
                nconf = nrange

            # loop over number of threads
            threads = []
            for i in range(args.p):
                start = int(args.c[0] + i * nconf/args.p)
                end = int(args.c[0] -1 + (i+1) * nconf/args.p)
                runprm["PBE_START"] = str(start)
                runprm["PBE_END"] = str(end)
                write_runprm(runprm)

                x = threading.Thread(target=mcce_thread, args=(args.e, i))
                threads.append(x)
                x.start()
                time.sleep(5)

            for i, x in enumerate(threads):
                x.join()
                logging.info("Main: join thread %d" % i)

            logging.info("Main: Done")

    return


if __name__ == "__main__":
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")
    helpmsg = "Run mcce step 3, energy calculations, with multiple threads."

    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-c", metavar=('start', 'end'), default=[1, 99999], nargs=2, help="starting and ending "
                                                                                         "conformer, default to 1 and 9999", type=int)
    parser.add_argument("-d", metavar="epsilon", default="4.0", help="protein dielectric constant for delphi, default to 4.0")
    parser.add_argument("-e", metavar="/path/to/mcce", default="mcce", help="mcce executable location, default to \"mcce\"")
    parser.add_argument("-f", metavar="tmp folder", default="/tmp", help="delphi temporary folder, default to /tmp")
    parser.add_argument("-p", metavar="processes", default=1, help="run mcce with number of processes, default to 1", type=int)
    parser.add_argument("-r", default=False, help="refresh opp files and head3.lst without running delphi", action="store_true")
    parser.add_argument("-u", metavar="Key=Value", default="", help="User customized variables")
    parser.add_argument("-x", metavar="/path/to/delphi", default="delphi", help="delphi executable location, default to \"delphi\"")
    parser.add_argument("--norun", default=False, help="Create run.prm but do not run step 3", action="store_true")

    args = parser.parse_args()

    process(args)



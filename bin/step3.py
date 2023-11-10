#!/usr/bin/env python
"""
Step 3 computes energy look up table. 
The input is step2_out.pdb.
The output is energies/*.opp, energies/*.oppl and head3.lst

Command line options:
-d number: dielectric constant (default 4)
-s PBSolver: PB solver name (default delphi)
-p number: run with number of threads (default 1)
-c start end: run conformer from start to end (default 0 to 99999)
-t path: temporary folder path (default /tmp) 
--vdw: run vdw calculation only
--fly: on-the-fly rxn0 calculation
--head3: recreate *.opp and head3.lst from step2_out.pdb and *.oppl files without doing calculation
-l file: load above options from a file

Usage examples:

1. Run step 3 with 6 threads
    step3.py -p 6

"""

import sys, argparse, shutil, logging, time, os
from vdw_pw import *



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



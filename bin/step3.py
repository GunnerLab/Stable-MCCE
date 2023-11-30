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
--refresh: recreate *.opp and head3.lst from step2_out.pdb and *.oppl files without doing calculation
-l file: load above options from a file

Usage examples:

1. Run step 3 with 6 threads
    step3.py -p 6

"""

import sys, argparse, shutil, logging, time, os, json
from pdbio import *


class RunOptions:
    def __init__(self, args):
        self.inputpdb = "step2_out.pdb"  # implicit input pdb file
        self.start = args.c[0]
        self.end = args.c[1]
        self.d = float(args.d)
        self.s = args.s
        self.p = args.p
        self.t = args.t
        self.ftpl = args.ftpl
        self.salt = args.salt
        self.vdw = args.vdw
        self.fly = args.fly
        self.refresh = args.refresh
        if args.l:  # load options from the specified file 
            lines=open(args.l).readlines()
            for line in lines:
                line = line.split("#")[0].strip()
                fields = [x.strip() for x in line.split()]  # Extract an option line and strip off the spaces
                if len(fields) > 0:
                    key = fields[0]
                    if key == "-c":
                        self.start = int(fields[1])
                        self.end = int(fields[2])
                    elif key == "-d":                        
                        self.d = float(fields[1])
                    elif key == "-s":
                        self.s = fields[1]
                    elif key == "-p":
                        self.p = int(fields[1])
                    elif key == "-t":
                        self.t = fields[1]
                    elif key == "--vdw":
                        self.vdw = True
                    elif key == "--fly":
                        self.fly = True
                    elif key == "--refresh":
                        self.refresh = True
                    elif key == "-ftpl":
                        self.ftpl = fields[1]

#    def toJSON(self):
#        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
                

class ExchangeAtom:
    def __init__(self):
        self.x = 0.0  
        self.y = 0.0
        self.z = 0.0
        self.r = 0.0  # radius
        self.c = 0.0  # charge
        return

class Exchange:
    # raw data and functions to process exchange data with PB wrapper
    def __init__(self, run_options, protein):
        #
        # Boundary set is composed by 3 atom collections
        # 1. Backbone atoms: these never change in position and charge, and match one to one from boundary to step2_out.pdb
        # 2. Side chain atoms: these may change depending on side chain selection.
        #
        # How to index atoms in the boundary list to conformer atom record and vice versa?
        # In conformer's atom record, i_compressed_bnd is the atom index to the compressed boundary set.
        # 
        #
        # Methods to operate on boundary sets
        # * delete a side chain
        # * add a side chain
        # * compress a boundary set by combining shared atoms (same coordinates) and removing atoms r=0 in the set and update ibound in atom record

        self.backbone = []
        self.all = []
        self.single = []
        self.compressed_bnd = []
        return

    def add_conformer(self, bnd, protein, ir, ic):
        for atom in protein.residue[ir].conf[ic].atom:




    # sites to receive potential and index to conformer atoms
        

        return
    
    def find_common_boundary(self, protein):
        boundary = []

        return boundary
    
    


if __name__ == "__main__":
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO, datefmt="%D %H:%M:%S")
    logging.info("Step 3 starts")

    helpmsg = "Run mcce step 3, energy lookup table calculations."

    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-c", metavar=('start', 'end'), default=[0, 99999], nargs=2, help="starting and ending "
                                                                                         "conformer, default to 0 and 9999", type=int)
    parser.add_argument("-d", metavar="epsilon", default="4.0", help="protein dielectric constant, default to 4.0")
    parser.add_argument("-s", metavar="pbs_name", default="delphi", help="PSE solver. Choices are delphi. default to \"delphi\"")
    parser.add_argument("-t", metavar="tmp folder", default="/tmp", help="PB solver temporary folder, default to /tmp")
    parser.add_argument("-p", metavar="processes", default=1, help="run step 3 with number of processes, default to 1", type=int)
    parser.add_argument("-ftpl", metavar="ftpl folder", default="", help="ftpl folder, default to \"param/\" of mcce exeuctable location")
    parser.add_argument("-salt", metavar="salt concentration", default=0.15, help="Salt concentration in moles/L. default to 0.15", type=float)
    parser.add_argument("--vdw", default=False, help="run vdw calculation only", action="store_true")
    parser.add_argument("--fly", default=False, help="don-the-fly rxn0 calculation", action="store_true")
    parser.add_argument("--refresh", default=False, help="recreate *.opp and head3.lst from step2_out.pdb and *.oppl files", action="store_true")
    parser.add_argument("-l", metavar="file", default="", help="load above options from a file")

    args = parser.parse_args()


    # Process run time options
    run_options = RunOptions(args)
    print(vars(run_options))

    # environment and ftpl
    if run_options.ftpl:
        env.runprm["FTPLDIR"] = run_options.ftpl     
    else:
        path = str(os.path.dirname(os.path.abspath(__file__)))
        base_path = os.path.dirname(path)
        env.runprm["FTPLDIR"] = base_path + "/param"
    env.load_ftpl()


    # read step2_out.pdb and convert to mcce structure
    protein = Protein()
    protein.loadpdb(run_options.inputpdb)


    # Prepare input for PB solver: common_boundary, sites to receive potential, and PB conditions


    boundary_conditions = BoundaryConditions(run_options, protein)

    # Set up parallel envrionment and run PB solver

    # Post-process electrostatic potential

    # Compute vdw

    # Assemble output files

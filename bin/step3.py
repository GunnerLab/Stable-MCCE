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
    def __init__(self, atom):
        self.x = atom.xyz[0]
        self.y = atom.xyz[1]
        self.z = atom.xyz[2]
        self.r = atom.r_bound  # radius
        self.c = atom.charge   # charge
        self.p = 0.0
        return

class Exchange:
    # We have to abadon the mop on and off mechanism when modifying the dielectric boundary.
    # In the parallel for loop, we don't want the boundary revising to be dependent on previous step.
    # Each step in the loop should be an addition of deletion from a static starting point.
    # Therefore, we will compose 
    # * backbone atoms 
    # * method to compose single side chain condition
    # * method to compose multi side chain condition
    # * method to compress boundary condition 

    def __init__(self, protein):
        self.backbone = self.add_backbone(protein)
        return

    def add_backbone(self, protein):
        backbone = []
        for res in protein.residue:
            if res.conformer:
                for atom in res.conformer[0]:
                    xyzrcp = ExchangeAtom(atom)
                    backbone.append(xyzrcp)
                    atom.ibound = len(backbone) - 1
        return backbone


    def compose_single(self, protein, ir, ic):
        """ Compose a single side chain boundary condition.
            The atoms are added in addition to backbone.
            Atoms other than in residue[ir], conformer[ic] are then appended.
            Atoms in residue[ir], conformer[ic] are appended last
            When appending an atom, this subroutine will check if the same atom (xyzrc identical) already exists. If yes, just update icound
        """

        return
    
    def compose_multi(self, protein, ir, ic):
        """ Compose a multi side chain boundary condition.
            The atoms are added in addition to backbone.
            Atoms other than in residue[ir], conformer[ic] are appended next.
            Atoms in residue[ir], conformer[ic] are appended last.
            When appending an atom, this subroutine will check if the same atom (xyzrc identical) already exists. If yes, just update icound
        """

        return
    
    


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


    boundary = Exchange()
    for res in protein.residue:
        for conf in res.conf:
            boundary.add_conformer(conf)


    # Set up parallel envrionment and run PB solver

    # Post-process electrostatic potential

    # Compute vdw

    # Assemble output files

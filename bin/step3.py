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
from multiprocess import Pool
from pdbio import *
from pbs_interfaces import *
import uuid


global protein, run_options

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

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
                

class ExchangeAtom:
    def __init__(self, atom):
        self.x = atom.xyz[0]
        self.y = atom.xyz[1]
        self.z = atom.xyz[2]
        self.r = atom.r_bound  # radius
        self.c = atom.charge   # charge
        self.p = 0.0
        return

class Exchange:  # This is the data passed to the PB wrapper, together with runoptions
    # We have to abadon the mop on and off mechanism when modifying the dielectric boundary.
    # In the parallel for loop, we don't want the boundary revising to be dependent on previous step.
    # Each step in the loop should be an addition of deletion from a static starting point.
    # Therefore, we will compose 
    # * backbone atoms
    # * index to match multiple atoms to boundary line number 
    # * method to compose single side chain condition
    # * method to compose multi side chain condition
    

    def __init__(self, protein):
        # initilaize backbone
        self.backbone_xyzrcp = []
        self.backbone_atom = []
        for res in protein.residue:
            if res.conf:
                for atom in res.conf[0].atom:
                    xyzrcp = ExchangeAtom(atom)
                    self.backbone_xyzrcp.append(xyzrcp)
                    self.backbone_atom.append([atom])   # the atom is in an array because it is allowed to have multiple atoms to match the same line in xyzrcp line


        return

    def compose_single(self, protein, ir, ic):
        """ Compose a single side chain boundary condition.
            The atoms are added in addition to backbone.
            Atoms other than in residue[ir], conformer[ic] are then appended.
            Atoms in residue[ir], conformer[ic] are appended last
        """
        self.single_bnd_xyzrcp = self.backbone_xyzrcp.copy()
        self.single_bnd_atom = self.backbone_atom.copy()


        for ires in range(len(protein.residue)):
            #print(protein.residue[ires].resID)
            if ires == ir:  # this is the residue we want to put desired side chain conf
                for atom in protein.residue[ires].conf[ic].atom:
                    xyzrcp = ExchangeAtom(atom)
                    self.single_bnd_xyzrcp.append(xyzrcp)
                    self.single_bnd_atom.append([atom])
                    
            else:  # find the first charged conformer if any, otherwise use the first
                if len(protein.residue[ires].conf) > 1:  # skip dummy or backbone only residue
                    i_useconf = 1  # defaul the first conformer
                    for iconf in range(1, len(protein.residue[ires].conf)):
                        # print(protein.residue[ires].conf[iconf].confID, protein.residue[ires].conf[iconf].crg)
                        if abs(protein.residue[ires].conf[iconf].crg) > 0.0001:
                            i_useconf = iconf
                            break

                    #print(protein.residue[ires].conf[i_useconf].confID, protein.residue[ires].conf[i_useconf].crg)
                    for atom in protein.residue[ires].conf[i_useconf].atom:
                        xyzrcp = ExchangeAtom(atom)
                        xyzrcp.c = 0.0   # boundary defining atom charge should be set as 0
                        self.single_bnd_xyzrcp.append(xyzrcp)
                        self.single_bnd_atom.append([atom])

        # Error checking, basic
        # for atom_list in self.ibound2atoms:
        #     print(len(atom_list), atom_list[0].atomID)
        #     if len(atom_list) != 1:
        #         print("ERROR")
        #         break

        logging.debug("%s Single-sidechain boundary record length should be equal: %d, %d" % (protein.residue[ir].conf[ic].confID, len(self.single_bnd_xyzrcp), len(self.single_bnd_atom)))

        return


    
    def compose_multi(self, protein, ir, ic):
        """ Compose a multi side chain boundary condition.
            The atoms are added in addition to backbone.
            Atoms other than in residue[ir], conformer[ic] are appended next.
            Atoms in residue[ir], conformer[ic] are appended last.
            When appending an atom, this subroutine will check if the same atom (xyzrc identical) already exists. If yes, just update icound
        """
        self.multi_bnd_xyzrcp = self.backbone_xyzrcp.copy()
        self.multi_bnd_atom = self.backbone_atom.copy()
        for ires in range(len(protein.residue)):
            #print(protein.residue[ires].resID)
            if ires == ir:  # this is the residue we want to put desired side chain conf
                for atom in protein.residue[ires].conf[ic].atom:
                    xyzrcp = ExchangeAtom(atom)
                    self.multi_bnd_xyzrcp.append(xyzrcp)
                    self.multi_bnd_atom.append([atom])
                    
            else:  # other residues will have all conformers with 0 charge
                if len(protein.residue[ires].conf) > 1:  # skip dummy or backbone only residue
                    residue_bnd_xyzrcp = []  # this is the xyzrcp record for all side chain atoms of this residue
                    residue_bnd_atom = []   # this points to the atom records of each line in residue_bnd_xyzrpc
                    for iconf in range(1, len(protein.residue[ires].conf)):
                        for atom in protein.residue[ires].conf[iconf].atom:
                            xyzrcp = ExchangeAtom(atom)
                            xyzrcp.c = 0.0   # boundary defining atom charge should be set as 0
                            # test if this atom existed within this residue already
                            found = False
                            for ib in range(len(residue_bnd_xyzrcp)):
                                if  abs(residue_bnd_xyzrcp[ib].x - xyzrcp.x) < 0.001 and \
                                    abs(residue_bnd_xyzrcp[ib].y - xyzrcp.y) < 0.001 and \
                                    abs(residue_bnd_xyzrcp[ib].z - xyzrcp.z) < 0.001 and \
                                    abs(residue_bnd_xyzrcp[ib].r - xyzrcp.r) < 0.001:  # identical atom 
                                    residue_bnd_atom[ib].append(atom)
                                    found = True
                                    break
                            
                            if not found:
                                residue_bnd_xyzrcp.append(xyzrcp)
                                residue_bnd_atom.append([atom])

                    # merge this residue to the multi-bnd
                    self.multi_bnd_xyzrcp += residue_bnd_xyzrcp
                    self.multi_bnd_atom += residue_bnd_atom


        # Basic error checking
        logging.debug("%s Multi-sidechain boundary record length should be equal: %d, %d" % (protein.residue[ir].conf[ic].confID, len(self.multi_bnd_xyzrcp), len(self.multi_bnd_atom)))

        return
    
    def write_single_bnd(self, fname):
        "This writes out both the xyzrcp and atom index files, for error checking and potentially as data exchange with PB wrapper."
        # write xyzrpc
        lines = []
        for xyzrpc in self.single_bnd_xyzrcp:
            line = "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n" % (xyzrpc.x,
                                                              xyzrpc.y,
                                                              xyzrpc.z,
                                                              xyzrpc.r,
                                                              xyzrpc.c,
                                                              xyzrpc.p)
            lines.append(line)
        open(fname+".xyzrcp", "w").writelines(lines)

        # write index file
        lines = []
        for matched_atoms in self.single_bnd_atom:
            atoms = " ".join([atom.atomID for atom in matched_atoms])
            xyzrc = "%8.3f %8.3f %8.3f %8.3f %8.3f" % (matched_atoms[0].xyz[0], 
                                           matched_atoms[0].xyz[1], 
                                           matched_atoms[0].xyz[2],
                                           matched_atoms[0].r_bound,
                                           matched_atoms[0].charge)

            line = "%s %s\n" % (xyzrc, atoms)
            lines.append(line)
        open(fname+".atoms", "w").writelines(lines)
    
    def write_multi_bnd(self, fname):
        "This writes out both the xyzrcp and atom index files, for error checking and potentially as data exchange with PB wrapper."
        # write xyzrpc
        lines = []
        for xyzrpc in self.multi_bnd_xyzrcp:
            line = "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n" % (xyzrpc.x,
                                                              xyzrpc.y,
                                                              xyzrpc.z,
                                                              xyzrpc.r,
                                                              xyzrpc.c,
                                                              xyzrpc.p)
            lines.append(line)
        open(fname+".xyzrcp", "w").writelines(lines)

        # write index file
        lines = []
        for matched_atoms in self.multi_bnd_atom:
            atoms = " ".join([atom.atomID for atom in matched_atoms])
            xyzrc = "%8.3f %8.3f %8.3f %8.3f %8.3f" % (matched_atoms[0].xyz[0], 
                                           matched_atoms[0].xyz[1], 
                                           matched_atoms[0].xyz[2],
                                           matched_atoms[0].r_bound,
                                           matched_atoms[0].charge)

            line = "%s %s\n" % (xyzrc, atoms)
            lines.append(line)
        open(fname+".atoms", "w").writelines(lines)


def def_boundary(ir, ic):
    boundary = Exchange(protein)

    boundary.compose_single(protein, ir, ic)
    boundary.compose_multi(protein, ir, ic)
    
    # Do not write out boundary condition except for debug purpose
    boundary.write_single_bnd("single_bnd")
    boundary.write_multi_bnd("multi_bnd")

    return boundary


def pbe(iric):
    ir = iric[0]
    ic = iric[1]
    bound = def_boundary(ir, ic)

    # switch to temporary unique directory
    cwd = os.getcwd()
    tmp_pbe = run_options.t + "/pbe_" + cwd.strip("/").replace("/", ".") + "/" + uuid.uuid4().hex[:6]
    os.makedirs(tmp_pbe)
    os.chdir(tmp_pbe)

    
    print(os.getcwd())

    # print(run_options.toJSON())
    # decide which pb solver, delphi = delphi legacy
    if run_options.s.upper() == "DELPHI":
        logging.info("Calling delphi to calulate conformer %s" % protein.residue[ir].conf[ic].confID)
        rxn = pbs_delphi(bound)
        # write raw files

    else:
        print("No compatible PBE solver detected, given pb solver is %s" % run_options.s)




    # switch back to the current directory
    os.rmdir(tmp_pbe)
    os.chdir(cwd)

    return(ir, ic)



if __name__ == "__main__":
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
    parser.add_argument("--debug", default=False, help="print debug information", action="store_true")
    parser.add_argument("-l", metavar="file", default="", help="load above options from a file")
    args = parser.parse_args()

    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=log_level, datefmt="%D %H:%M:%S")
    logging.info("Step 3 starts")


    # Process run time options
    run_options = RunOptions(args)
    # print(vars(run_options))

    # environment and ftpl
    if run_options.ftpl:
        env.runprm["FTPLDIR"] = run_options.ftpl     
    else:
        path = str(os.path.dirname(os.path.abspath(__file__)))
        base_path = os.path.dirname(path)
        env.runprm["FTPLDIR"] = base_path + "/param"
    env.load_ftpl()


    # read step2_out.pdb and convert to mcce structure
    logging.info("Read from step2_out.pdb")
    protein = Protein()
    protein.loadpdb(run_options.inputpdb)
    protein.update_confcrg()


    # Prepare input for PB solver: common_boundary, sites to receive potential, and PB conditions
    # make conformer list with their corresponding ir and ic. This list or (ir, ic) will be passed as an array 
    # that multiprocess module will take in as work load.
    work_load = []
    for ir in range(len(protein.residue)):
        if len(protein.residue[ir].conf) > 1:  # skip dummy or backbone only residue
            for ic in range(1, len(protein.residue[ir].conf)):
                work_load.append((ir,ic))
    logging.debug("work_load as (ir ic) list", str(work_load))


    # Set up parallel envrionment and run PB solver
    max_pool = run_options.p
    logging.info("Running PBE solver in %d threads" % max_pool)
    with Pool(max_pool) as process:
        work_out = process.imap(pbe, work_load)
        logging.debug("Done PDE solving on %s" % str(list(work_out)))
        

    # Post-process electrostatic potential

    # Compute vdw

    # Assemble output files

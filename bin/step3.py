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
from multiprocess import Pool, current_process
from pdbio import *
from pbs_interfaces import *

energy_folder = "energies"
PW_CUTOFF = 0.001   # cut off value for pairwise interaction to report



global protein, run_options


class ElePW:
    def __init__(self):
        self.mark = ""
        self.multi = 0.0
        self.single = 0.0
        self.scaled = 0.0
        self.averaged = 0.0
        return


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
        self.debug = args.debug
        self.refresh = args.refresh
        if args.l:  # load options from the specified file 
            lines = open(args.l).readlines()
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
                    elif key == "--debug":
                        self.debug = True
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
        self.c = 0.0  # default to boundary defining atom, charge should be set as 0
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
                    self.backbone_atom.append([
                                                  atom])  # the atom is in an array because it is allowed to have multiple atoms to match the same line in xyzrcp line

        self.float_bnd_xyzrcp = []
        self.float_bnd_atom = []
        self.single_bnd_xyzrcp = []
        self.single_bnd_atom = []
        self.multi_bnd_xyzrcp = []
        self.multi_bnd_atom = []
        return

    def compose_float(self, protein, ir, ic):
        """ Compose a floating side chain boundary condition for rxn0 calculation.
            Only atoms in residue[ir], conformer[ic] are in this list
        """
        self.float_bnd_xyzrcp = []
        self.float_bnd_atom = []

        for atom in protein.residue[ir].conf[ic].atom:
            xyzrcp = ExchangeAtom(atom)
            xyzrcp.c = atom.charge
            self.float_bnd_xyzrcp.append(xyzrcp)
            self.float_bnd_atom.append([atom])

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
            # print(protein.residue[ires].resID)
            if ires == ir:  # this is the residue we want to put desired side chain conf
                for atom in protein.residue[ires].conf[ic].atom:
                    xyzrcp = ExchangeAtom(atom)
                    xyzrcp.c = atom.charge
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

                    # just for boundary
                    for atom in protein.residue[ires].conf[i_useconf].atom:
                        xyzrcp = ExchangeAtom(atom)
                        self.single_bnd_xyzrcp.append(xyzrcp)
                        self.single_bnd_atom.append([atom])

        # Error checking, basic
        # for atom_list in self.ibound2atoms:
        #     print(len(atom_list), atom_list[0].atomID)
        #     if len(atom_list) != 1:
        #         print("ERROR")
        #         break

        logging.debug("%s Single-sidechain boundary record length should be equal: %d, %d" % (
        protein.residue[ir].conf[ic].confID, len(self.single_bnd_xyzrcp), len(self.single_bnd_atom)))

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
            # print(protein.residue[ires].resID)
            if ires == ir:  # this is the residue we want to put desired side chain conf
                for atom in protein.residue[ires].conf[ic].atom:
                    xyzrcp = ExchangeAtom(atom)
                    xyzrcp.c = atom.charge
                    self.multi_bnd_xyzrcp.append(xyzrcp)
                    self.multi_bnd_atom.append([atom])

            else:  # other residues will have all conformers with 0 charge
                if len(protein.residue[ires].conf) > 1:  # skip dummy or backbone only residue
                    residue_bnd_xyzrcp = []  # this is the xyzrcp record for all side chain atoms of this residue
                    residue_bnd_atom = []  # this points to the atom records of each line in residue_bnd_xyzrpc
                    for iconf in range(1, len(protein.residue[ires].conf)):
                        for atom in protein.residue[ires].conf[iconf].atom:
                            xyzrcp = ExchangeAtom(atom)
                            # test if this atom existed within this residue already
                            found = False
                            for ib in range(len(residue_bnd_xyzrcp)):
                                if abs(residue_bnd_xyzrcp[ib].x - xyzrcp.x) < 0.001 and \
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
        logging.debug("%s Multi-sidechain boundary record length should be equal: %d, %d" % (
        protein.residue[ir].conf[ic].confID, len(self.multi_bnd_xyzrcp), len(self.multi_bnd_atom)))

        return

    def write_float_bnd(self, fname):
        "This writes out both the xyzrcp and atom index files, for error checking and potentially as data exchange with PB wrapper."
        # write xyzrpc
        lines = []
        for xyzrcp in self.float_bnd_xyzrcp:
            line = "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n" % (xyzrcp.x,
                                                              xyzrcp.y,
                                                              xyzrcp.z,
                                                              xyzrcp.r,
                                                              xyzrcp.c,
                                                              xyzrcp.p)
            lines.append(line)
        open(fname + ".xyzrcp", "w").writelines(lines)

        # write index file
        lines = []
        for matched_atoms in self.float_bnd_atom:
            atoms = " ".join([atom.atomID for atom in matched_atoms])
            xyzrc = "%8.3f %8.3f %8.3f %8.3f %8.3f" % (matched_atoms[0].xyz[0],
                                                       matched_atoms[0].xyz[1],
                                                       matched_atoms[0].xyz[2],
                                                       matched_atoms[0].r_bound,
                                                       matched_atoms[0].charge)

            line = "%s %s\n" % (xyzrc, atoms)
            lines.append(line)
        open(fname + ".atoms", "w").writelines(lines)


    def write_single_bnd(self, fname):
        "This writes out both the xyzrcp and atom index files, for error checking and potentially as data exchange with PB wrapper."
        # write xyzrpc
        lines = []
        for xyzrcp in self.single_bnd_xyzrcp:
            line = "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n" % (xyzrcp.x,
                                                              xyzrcp.y,
                                                              xyzrcp.z,
                                                              xyzrcp.r,
                                                              xyzrcp.c,
                                                              xyzrcp.p)
            lines.append(line)
        open(fname + ".xyzrcp", "w").writelines(lines)

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
        open(fname + ".atoms", "w").writelines(lines)

    def write_multi_bnd(self, fname):
        "This writes out both the xyzrcp and atom index files, for error checking and potentially as data exchange with PB wrapper."
        # write xyzrpc
        lines = []
        for xyzrcp in self.multi_bnd_xyzrcp:
            line = "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n" % (xyzrcp.x,
                                                              xyzrcp.y,
                                                              xyzrcp.z,
                                                              xyzrcp.r,
                                                              xyzrcp.c,
                                                              xyzrcp.p)
            lines.append(line)
        open(fname + ".xyzrcp", "w").writelines(lines)

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
        open(fname + ".atoms", "w").writelines(lines)


def def_boundary(ir, ic):
    boundary = Exchange(protein)

    boundary.compose_float(protein, ir, ic)
    boundary.compose_single(protein, ir, ic)
    boundary.compose_multi(protein, ir, ic)

    # Do not write out boundary condition except for debug purpose
    if run_options.debug:
        boundary.write_float_bnd("float_bnd")
        boundary.write_single_bnd("single_bnd")
        boundary.write_multi_bnd("multi_bnd")

    return boundary


def get_conflist(protein):
    conf_list = []
    bkb_list = []
    for res in protein.residue:
        for conf in res.conf:
            if conf.confID[3:5] == "BK":
                bkb_list.append(conf.confID)
            else:
                conf_list.append(conf.confID)
    return conf_list, bkb_list


def pbe(iric):
    "Calculate electrostatic terms: pairwise and reaction filed energy, store results in energies/*.raw file"

    ir = iric[0]
    ic = iric[1]
    pid = current_process()  # Identify this worker
    confid = protein.residue[ir].conf[ic].confID
    resid = confid[:3] + confid[5:11]
    bound = def_boundary(ir, ic)
    rxn = 0.0
    

    # skip pbe if atoms in this conformer are all 0 charged
    all_0 = True
    for atom in protein.residue[ir].conf[ic].atom:
        if abs(atom.charge) > 0.001:
            all_0 = False
            break
    if all_0:  # skip
        logging.info("Skipping PBE solver for non-charge confortmer %s..." % confid)
    else:
        # switch to temporary unique directory
        cwd = os.getcwd()
        # tmp_pbe = run_options.t + "/pbe_" + cwd.strip("/").replace("/", ".") + "/" + uuid.uuid4().hex[:6]
        tmp_pbe = run_options.t + "/pbe_" + cwd.strip("/").replace("/", ".") + "/" + confid
        if not os.path.exists(tmp_pbe):
            os.makedirs(tmp_pbe)
        os.chdir(tmp_pbe)

        # decide which pb solver, delphi = delphi legacy
        if run_options.s.upper() == "DELPHI":
            logging.info("%s: Calling delphi to calulate conformer %s" % (pid.name, confid))
            pbs_delphi = PBS_DELPHI()
            rxn0, rxn = pbs_delphi.run(bound, run_options)

        else:
            print("No compatible PBE solver detected, given pb solver is %s" % run_options.s)

        # switch back to the current directory
        if not run_options.debug:
            shutil.rmtree(tmp_pbe)
        os.chdir(cwd)

    # write raw opp file
    fname = "%s/%s.raw" % (energy_folder, confid)

    if not os.path.exists(energy_folder):
        os.mkdir(energy_folder)

    if all_0:
        # create an empty file as a marker to indicate this conformer has been calculated
        with open(fname, 'w') as fp:
            pass
    else:
        # generate electrostatic interaction raw file
        raw_lines = []

        # Part 1: Method = run_options.s
        line = "[Method] %s\n" % run_options.s
        raw_lines.append(line)

        # Part 2: pw to other conformers, single and multi
        pw_single = {}
        pw_multi = {}
        for ia in range(len(bound.single_bnd_xyzrcp)):
            # print(confid, bound.single_bnd_atom[ia][0].confID, bound.single_bnd_atom[ia][0].name)
            pw_confname = bound.single_bnd_atom[ia][0].confID  # single conformer
            p = bound.single_bnd_xyzrcp[ia].p * bound.single_bnd_atom[ia][0].charge
            if abs(p) > 0.0001:
                if pw_confname in pw_single:
                    pw_single[pw_confname] += p
                else:
                    pw_single[pw_confname] = p

        # convert to Kcal and remove interaction within residue including backbone piece
        pw_single.update((key, value / KCAL2KT) for key, value in pw_single.items())

        # print(pw_single)
        # print(len(pw_single))

        pw_multi = {}
        for ia in range(len(bound.multi_bnd_xyzrcp)):
            for atom in bound.multi_bnd_atom[ia]:
                pw_confname = atom.confID
                if pw_confname[3:5] == "BK":
                    continue
                p = bound.multi_bnd_xyzrcp[ia].p * atom.charge
                if abs(p) > 0.0001:
                    if pw_confname in pw_multi:
                        pw_multi[pw_confname] += p
                    else:
                        pw_multi[pw_confname] = p

        # convert to Kcal
        pw_multi.update((key, value / KCAL2KT) for key, value in pw_multi.items())
        # for key, value in pw_multi.items():
        #     if abs(value) >= 0.001:
        #         print("%s %8.3f" % (key, value))

        # get conformer list and backbone segment list
        conf_list, bkb_list = get_conflist(protein)

        # write pw section
        line = "\n[PAIRWISE confID single multi flag, kcal/mol]\n"
        raw_lines.append(line)

        for pw_conf in conf_list:
            if resid == pw_conf[:3] + pw_conf[5:11]:
                continue  # skip conformer within residue
            single = 0.0
            multi = 0.0
            non0 = False
            reference = ""
            if pw_conf in pw_single:
                non0 = True
                single = pw_single[pw_conf]
                reference = "*"  # this is marked as a reference conformer for boundary correction
            if pw_conf in pw_multi:
                non0 = True
                multi = pw_multi[pw_conf]
            if non0 and (abs(single) >= PW_CUTOFF or abs(multi) >= PW_CUTOFF):
                line = "%s %8.3f %8.3f %s\n" % (pw_conf, single, multi, reference)
                raw_lines.append(line)

        # Part 3: backbone interaction total
        bkb_total = 0.0
        bkb_breakdown_lines = ["\n[BACKBONE breakdown, kcal/mol]\n"]
        for pw_conf in bkb_list:
            non0 = False
            bkb_pw = 0.0
            if pw_conf in pw_single:
                non0 = True
                bkb_pw = pw_single[pw_conf]
                # if resid != (pw_conf[:3] + pw_conf[5:11]):    # exclude the backbone piece when calculating the total
                #     bkb_total += bkb_pw
                bkb_total += bkb_pw  # do inclusive calculation for now

            if non0 and abs(bkb_pw) >= PW_CUTOFF:
                line = "%s %8.3f\n" % (pw_conf, bkb_pw)
                bkb_breakdown_lines.append(line)

        line = "\n[BACKBONE total including self, kcal/mol] %8.3f\n" % bkb_total
        raw_lines.append(line)

        # Part 4: backbone interaction breakdown
        line = "\n[BACKBONE breakdown]\n"
        raw_lines += bkb_breakdown_lines

        # Part 5: rxn
        if run_options.fly:
            line = "\n[RXN0, kcal/mol] %8.3f" % rxn0
            raw_lines.append(line)

        line = "\n[RXN, kcal/mol] %8.3f" % rxn
        raw_lines.append(line)

        open(fname, "w").writelines(raw_lines)

    return (ir, ic)


def is_conf_clash(conf1, conf2, use_r_bound=True):
    """Quick detection of the conformer to conformer clash without considering connectivity.
        If use_r_bound is True, then use the radius of dielectric boundary,
        otherwise use r_vdw, Van der Waals radius
    """
    clash = False
    for atom1 in conf1.atom:
        if clash:
            break
        if use_r_bound:
            r1 = atom1.r_bound * 0.5  # artificially allow less strict clash detection to match old mcce
        else:
            r1 = atom1.r_vdw
        for atom2 in conf2.atom:
            if use_r_bound:
                r2 = atom2.r_bound * 0.5  # artificially allow less strict clash detection to match old mcce
            else:
                r2 = atom2.r_vdw
            d2 = ddvv(atom1.xyz, atom2.xyz)
            if d2 < (r1 + r2) ** 2:
                clash = True
                # print(atom1.atomID, atom2.atomID, math.sqrt(d2))
                break

    return clash


def postprocess_ele():
    conf_list, _ = get_conflist(protein)
    ele_matrix = {}
    ele_k_s2m_byresid = {}

    # load raw
    for conf in conf_list:
        fname = "%s/%s.raw" % (energy_folder, conf)
        if os.path.isfile(fname):
            lines = open(fname).readlines()
            pw_start = False
            resid1 = conf[:3] + conf[5:11]

            for line in lines:
                if line[:9] == "[PAIRWISE":
                    pw_start = True
                    continue
                elif line[:9] == "[BACKBONE":
                    break
                if pw_start:
                    fields = line.split()
                    if len(fields) >= 3:
                        conf2 = fields[0]
                        single = float(fields[1])
                        multi = float(fields[2])
                        if len(fields) > 3:
                            mark = fields[3]
                        else:
                            mark = ""
                        ele_pw = ElePW()
                        ele_pw.multi = multi
                        ele_pw.single = single
                        ele_pw.mark = mark
                        ele_matrix[(conf, conf2)] = ele_pw

                        if "*" in mark:
                            resid2 = conf2[:3] + conf2[5:11]
                            if abs(multi) > 0.1:
                                k_single_multi = single / multi
                            else:
                                k_single_multi = 1.0
                            if k_single_multi > 1.0:
                                k_single_multi = 1.0
                            ele_k_s2m_byresid[(resid1, resid2)] = k_single_multi

    # scale multi by the k factor, or use single for reference point
    for conf_pair, ele_pw in ele_matrix.items():
        conf1, conf2 = conf_pair
        resid1 = conf1[:3] + conf1[5:11]
        resid2 = conf2[:3] + conf2[5:11]
        key = (resid1, resid2)
        if key in ele_k_s2m_byresid:
            k = ele_k_s2m_byresid[(resid1, resid2)]
        else:
            k = 1.0
        if "*" in ele_pw.mark:
            ele_pw.scaled = ele_pw.single  # use the reference directly
        else:
            ele_pw.scaled = k * ele_pw.multi

    # for conf_pair, ele_pw in ele_matrix.items():
    #     print("%s: %8.3f %8.3f %8.3f %s" % (conf_pair, ele_pw.scaled, ele_pw.single, ele_pw.multi, ele_pw.mark))

    # Find clashing dielectric boundary. This is because the single conformation boundary is composed by the selected
    # conformer + the native conformer of other residues. This selected conformer may have geometry conflict with the
    # native conformers. If a clash is identified, a question mark will be added to the ele_pw mark.
    logging.debug("Detecting conformer to conformer clashes ...")
    for res1 in protein.residue:
        if len(res1.conf) > 1:  # skip dummy or backbone only residue
            for conf1 in res1.conf[1:]:
                conf1_id = conf1.confID
                for res2 in protein.residue:
                    if res2 == res1:
                        continue
                    # find the reference conformer to res2
                    conf2_ref = None
                    if len(res2.conf) > 1:  # skip dummy or backbone only residue
                        for conf2 in res2.conf[1:]:
                            conf2_id = conf2.confID
                            key = (conf1_id, conf2_id)
                            if key in ele_matrix:
                                ele_pw = ele_matrix[key]
                                if "*" in ele_pw.mark:
                                    conf2_ref = conf2
                                    break

                    # detect the clash between conf1 and conf2_ref
                    clash = False
                    if conf2_ref:
                        clash = is_conf_clash(conf1, conf2_ref)

                    # mark res2 conformers with "?"
                    if clash:
                        # print(conf1.confID, conf2_ref.confID)
                        for conf2 in res2.conf[1:]:
                            conf2_id = conf2.confID
                            key = (conf1_id, conf2_id)
                            if key in ele_matrix:
                                ele_pw = ele_matrix[key]
                                ele_pw.mark = "?" + ele_pw.mark

    # ele correction and average
    # The correction starts with ele_pw.scaled, and follow the following rules to make correction
    # C-C: Charged to charged abnormal case: reduce by a factor of 1.5 of the smaller pw at multi condition
    # C-C: Charged to charged normal case: average scaled pws
    # D-D: dipole to dipole, reduce by a factor of 2.0
    # C-D: charged to dipole,
    for conf_pair, ele_pw in ele_matrix.items():
        if conf_pair[0][3] == "+" or conf_pair[0][3] == "-":
            conf1_type = "C"
        else:
            conf1_type = "D"
        if conf_pair[1][3] == "+" or conf_pair[1][3] == "-":
            conf2_type = "C"
        else:
            conf2_type = "D"
        reversed_key = (conf_pair[1], conf_pair[0])
        if reversed_key in ele_matrix:
            reversed_ele_pw = ele_matrix[reversed_key]
        else:
            reversed_ele_pw = ElePW()
        if conf1_type == "C" and conf2_type == "C":
            # abnormal case, opposite sign when scaled, or both directions have ? mark, scaled down by factor 1.5 from
            # the smaller of the two multi
            if ele_pw.scaled * ele_pw.multi < 0 or \
                    reversed_ele_pw.scaled * reversed_ele_pw.multi < 0 or \
                    ("?" in ele_pw.mark and "?" in reversed_ele_pw.mark):
                if abs(ele_pw.multi) < abs(reversed_ele_pw.multi):
                    ele_pw.averaged = ele_pw.multi / 1.5
                    reversed_ele_pw.averaged = ele_pw.multi / 1.5
                else:
                    ele_pw.averaged = reversed_ele_pw.multi / 1.5
                    reversed_ele_pw.averaged = reversed_ele_pw.multi / 1.5
                # debug
                # print("%s -> %s    |    %s -> %s" % (conf_pair[0], conf_pair[1], reversed_key[0], reversed_key[1]))
                # print("Multi:        %8.3f              |          %8.3f" % (ele_matrix[conf_pair].multi, ele_matrix[reversed_key].multi))
                # print("Single:       %8.3f              |          %8.3f" % (ele_matrix[conf_pair].single, ele_matrix[reversed_key].single))
                # print("Scaled:       %8.3f              |          %8.3f" % (ele_matrix[conf_pair].scaled, ele_matrix[reversed_key].scaled))
                # print("Averaged:     %8.3f              |          %8.3f" % (ele_matrix[conf_pair].averaged, ele_matrix[reversed_key].averaged))
                # print(vars(ele_pw))
                # print(vars(reversed_ele_pw))
            elif "?" in ele_pw.mark and "?" not in reversed_ele_pw.mark:
                ele_pw.averaged = reversed_ele_pw.scaled
                reversed_ele_pw.averaged = reversed_ele_pw.scaled
            elif "?" not in ele_pw.mark and "?" in reversed_ele_pw.mark:
                ele_pw.averaged = ele_pw.scaled
                reversed_ele_pw.averaged = ele_pw.scaled
            else:
                if abs(ele_pw.scaled) < abs(reversed_ele_pw.scaled):
                    ele_pw.averaged = reversed_ele_pw.averaged = ele_pw.scaled
                else:
                    ele_pw.averaged = reversed_ele_pw.averaged = reversed_ele_pw.scaled
        elif conf1_type == "D" and conf2_type == "D":
            ele_pw.averaged = reversed_ele_pw.averaged = (ele_pw.scaled + reversed_ele_pw.scaled)/2.0  # Originally averaging multi
        else:  # must be "D" to "C"
            if ele_pw.scaled * ele_pw.multi < 0 or \
                    reversed_ele_pw.scaled * reversed_ele_pw.multi < 0 or \
                    ("?" in ele_pw.mark and "?" in reversed_ele_pw.mark):
                ele_pw.averaged = reversed_ele_pw.averaged = (ele_pw.multi+reversed_ele_pw.multi)/2.0/1.5
            elif "?" in ele_pw.mark and "?" not in reversed_ele_pw.mark:
                ele_pw.averaged = reversed_ele_pw.averaged = reversed_ele_pw.scaled
            elif "?" not in ele_pw.mark and "?" in reversed_ele_pw.mark:
                ele_pw.averaged = reversed_ele_pw.averaged = ele_pw.scaled
            else:
                ele_pw.averaged = reversed_ele_pw.averaged = (ele_pw.scaled + reversed_ele_pw.scaled) / 2.0

    return ele_matrix


def compose_opp(protein, ele_matrix):

    epath = "energies"
    for res1 in protein.residue:
        for conf1 in res1.conf[1:]:
            fname = "%s/%s.raw" % (epath, conf1.confID)
            if os.path.isfile(fname):  # only create opp files when a raw file exists
                fname = "%s/%s.opp" % (epath, conf1.confID)
                lines = []
                for res2 in protein.residue:
                    if res2 != res1:
                        for conf2 in res2.conf[1:]:
                            conf_pair = (conf1.confID, conf2.confID)
                            i1 = conf1.i
                            i2 = conf2.i
                            if conf_pair in ele_matrix:
                                average = ele_matrix[conf_pair].averaged
                                scaled = ele_matrix[conf_pair].scaled
                                multi = ele_matrix[conf_pair].multi
                                mark = ele_matrix[conf_pair].mark
                            else:
                                average = scaled = multi = 0.0
                                mark = ""
                            if abs(protein.vdw_pw[i1, i2]) > PW_CUTOFF or abs(average) > PW_CUTOFF:
                                lines.append("%05d %s %8.3f %7.3f %7.3f %7.3f %s\n" % (conf2.serial, conf2.confID, average, protein.vdw_pw[i1, i2], scaled, multi, mark))
                open(fname, "w").writelines(lines)

    return


def compose_head3(protein):
    epath = "energies"
    head3lines = ["iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    history\n"]
    # read backbone ele interaction epol
    epol_all = {}
    rxn_all = {}
    rxn0_all = {}
    for res in protein.residue:
        for conf in res.conf[1:]:
            fname = "%s/%s.raw" % (epath, conf.confID)
            if os.path.isfile(fname):
                lines = open(fname).readlines()
                for line in lines:
                    if line.startswith("[BACKBONE total"):
                        fields = line.split("]")
                        epol_all[conf.confID] = float(fields[-1])
                    elif line.startswith("[RXN0,"):
                        fields = line.split("]")
                        rxn0_all[conf.confID] = float(fields[-1])
                    elif line.startswith("[RXN,"):
                        fields = line.split("]")
                        rxn_all[conf.confID] = float(fields[-1])

    # natom dictonary to determine dummy conformers
    natom_byconftype = {}
    for key, value in env.param.items():
        if key[0] == "CONFLIST":
            for conftype in env.param[("CONFLIST", key[1])]:
                natom_byconftype[conftype] = 0


    for key, value in env.param.items():
        if key[0] == "CONNECT":
            conftype = key[2]
            if conftype in natom_byconftype:  # there are cases CONNECT exists but conftype was commented out
                natom_byconftype[conftype] += 1

    serial = 1
    for res in protein.residue:

        # add dummy conformers
        for conftype in env.param[("CONFLIST", res.resID[:3])]:
            natom = natom_byconftype[conftype]
            if natom == 0 and conftype[-2:] != "BK":
                newconf = Conformer()
                n = len(res.conf)
                newconf.confID = "%s%s%s%s_%03d" % (res.resID[:3], conftype[-2:], res.resID[7], res.resID[3:7], n)
                newconf.history = "DM"
                newconf.mark = "d"
                res.conf.append(newconf)

        tors_confs = []
        for conf in res.conf[1:]:
            tors_confs.append(torsion(conf))
        if tors_confs:
            min_tors = min(tors_confs)
            tors_confs = [x - min_tors for x in tors_confs]
        count = 0
        for conf in res.conf[1:]:
            if conf.mark == "d":
                iconf = 0
            else:
                iconf = conf.i + 1
            confID = conf.confID
            flag = "f"
            occ = 0.0
            crg = conf.crg
            conftype = conf.confID[:5]
            em0 = env.param["CONFORMER", conftype].param["em0"]
            pka0 = env.param["CONFORMER", conftype].param["pka0"]
            ne = env.param["CONFORMER", conftype].param["ne"]
            nh = env.param["CONFORMER", conftype].param["nh"]
            vdw0 = conf.vdw0
            vdw1 = conf.vdw1

            tors = tors_confs[count]
            count += 1

            if confID in epol_all:
                epol = epol_all[conf.confID]
            else:
                epol = 0.0

            if confID in rxn_all:
                rxn = rxn_all[conf.confID]
            else:
                rxn = 0.0
            
            if run_options.fly:
                if conf.confID in rxn0_all:
                    rxn0 = rxn0_all[conf.confID]
                else:
                    rxn0 = 0.0
            else:
                epsilon = run_options.d
                key3 = "rxn%02d" % int(epsilon)
                rxn0 = env.param["CONFORMER", conftype].param[key3]
            dsolv = rxn - rxn0
            history = conf.history

            key = ("EXTRA", conftype)
            if key in env.param:
                extra = env.param["EXTRA", conftype]
            else:
                extra = 0.0

            fname = "%s/%s.raw" % (epath, conf.confID)
            if os.path.isfile(fname):  # only create opp files when a raw file exists
                mark = "t"
            else:
                mark = "f"
            if conf.mark == "d":
                mark = conf.mark  # inherit dummy mark from conformer

            head3lines.append("%05d %s %s %4.2f %6.3f %5d %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10s %s\n" % \
              (serial, confID, flag, occ, crg, em0, pka0, ne, nh, vdw0, vdw1, tors, epol, dsolv, extra, history, mark))
            conf.serial = serial            
            serial += 1
    open("head3.lst", "w").writelines(head3lines)
    return


if __name__ == "__main__":
    helpmsg = "Run mcce step 3, energy lookup table calculations."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-c", metavar=('start', 'end'), default=[1, 99999], nargs=2, help="starting and ending "
                                                                                          "conformer, default to 1 and 9999",
                        type=int)
    parser.add_argument("-d", metavar="epsilon", default="4.0", help="protein dielectric constant, default to 4.0")
    parser.add_argument("-s", metavar="pbs_name", default="delphi",
                        help="PSE solver. Choices are delphi. default to \"delphi\"")
    parser.add_argument("-t", metavar="tmp folder", default="/tmp", help="PB solver temporary folder, default to /tmp")
    parser.add_argument("-p", metavar="processes", default=1, help="run step 3 with number of processes, default to 1",
                        type=int)
    parser.add_argument("-ftpl", metavar="ftpl folder", default="",
                        help="ftpl folder, default to \"param/\" of mcce exeuctable location")
    parser.add_argument("-salt", metavar="salt concentration", default=0.15,
                        help="Salt concentration in moles/L. default to 0.15", type=float)
    parser.add_argument("--vdw", default=False, help="run vdw calculation only", action="store_true")
    parser.add_argument("--fly", default=False, help="don-the-fly rxn0 calculation", action="store_true")
    parser.add_argument("--refresh", default=False,
                        help="recreate *.opp and head3.lst from step2_out.pdb and *.oppl files", action="store_true")
    parser.add_argument("--debug", default=False, help="print debug information and keep pb solver tmp",
                        action="store_true")
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
    env.load_runprm()
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
    counter = 1
    for ir in range(len(protein.residue)):
        if len(protein.residue[ir].conf) > 1:  # skip dummy or backbone only residue
            for ic in range(1, len(protein.residue[ir].conf)):
                if run_options.end >= counter >= run_options.start:
                    work_load.append((ir, ic))
                counter += 1
    logging.debug("work_load as (ir ic) list [%s]" % ','.join(map(str, work_load)))

    # Set up parallel envrionment and run PB solver
    max_pool = run_options.p
    logging.info("Running PBE solver in %d threads" % max_pool)
    with Pool(max_pool) as process:
        work_out = process.imap(pbe, work_load)
        logging.debug("Done PDE solving on %s" % str(list(work_out)))

    cwd = os.getcwd()
    pbe_folder = run_options.t + "/pbe_" + cwd.strip("/").replace("/", ".")
    if not run_options.debug and os.path.exists(pbe_folder):
        shutil.rmtree(pbe_folder)

    # Post-process electrostatic potential
    logging.info("Processing ele pairwise interaction...")
    ele_matrix = postprocess_ele()
    logging.info("Processing ele pairwise interaction...Done")

    # Debug
    # for conf_pair, ele_pw in ele_matrix.items():
    #     reversed_key = (conf_pair[1], conf_pair[0])
    #     ele_pw = ele_matrix[conf_pair]
    #     if reversed_key in ele_matrix:
    #         reversed_ele_pw = ele_matrix[reversed_key]
    #     else:
    #         reversed_ele_pw = ElePW()
    #
    #     if abs(ele_pw.multi) < 0.1:
    #         continue
    #     print("%s -> %s    |    %s -> %s" % (conf_pair[0], conf_pair[1], reversed_key[0], reversed_key[1]))
    #     print("Multi:        %8.3f              |          %8.3f" % (ele_pw.multi, reversed_ele_pw.multi))
    #     print("Single:       %8.3f              |          %8.3f" % (ele_pw.single, reversed_ele_pw.single))
    #     print("Scaled:       %8.3f              |          %8.3f" % (ele_pw.scaled, reversed_ele_pw.scaled))
    #     print("Averaged:     %8.3f              |          %8.3f" % (ele_pw.averaged, reversed_ele_pw.averaged))
    #     print("Mark:         %8s              |          %8s" % (ele_pw.mark, reversed_ele_pw.mark))

    # Compute vdw, not doing parallelization at this moment
    logging.info("Making atom connectivity ...")
    protein.make_connect12()
    protein.make_connect13()
    protein.make_connect14()

    logging.info("Calculating vdw ...")
    protein.calc_vdw()
    # For efficiency reason, the vdw pairwise table is a matrix protein.vdw_pw[conf1.i, conf2.i]


    # Assemble output files, order sensitive as head3.lst subroutine will make serial for conformers later used by opp files
    logging.info("Composing head3.lst ...")
    compose_head3(protein)

    logging.info("Composing opp files ...")
    compose_opp(protein, ele_matrix)

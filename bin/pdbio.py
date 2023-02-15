#!/usr/bin/env python
"""Sub routines to input and write out MCCE PDB files."""

import math
import os
import glob

# bond distance scaling factor: cutoff = k*(r_vdw1 + r_vdw2)
BONDDISTANCE_scaling = 0.54  # calibrated by 1akk
#CUTOFF2 = 1.65*1.65

# scaling factor for vdw
VDW_SCALE14 = 0.5

# set any big conf vdw to 999
VDW_CUTOFF_FAR2 = 100     # set atom vdw to 0 if atoms are further than this value
VDW_CUTOFF_NEAR2 = 1      # set atom vdw to 999 if atoms are closer than this value
VDW_UPLIMIT = 320.0     # set conf vdw to 999 if bigger than this number


def ddvv(xyz1, xyz2):
    """Distance squared between two vectors."""
    dx=xyz1[0]-xyz2[0]
    dy=xyz1[1]-xyz2[1]
    dz=xyz1[2]-xyz2[2]
    return dx*dx+dy*dy+dz*dz

class Atom:
    def __init__(self):
        self.serial = 0
        self.name = "  X "
        self.altLoc = " "
        self.resName = "UNK"
        self.chainID = "A"
        self.resSeq = 0
        self.iCode = " "
        self.confNum = 0
        self.atomID = ""
        self.confID = ""
        self.confType = ""
        self.resID = ""
        self.xyz = (0.0, 0.0, 0.0)
        self.connectivity_param = ""
        self.r_bound = 0.0
        self.r_vdw = 0.0
        self.e_vdw= 0.0
        self.connect12 = []
        self.connect13 = []
        self.connect14 = []
        return

    def loadline(self, line):
        self.serial = int(line[6:11])
        self.name = line[12:16]
        self.altLoc = line[16]
        self.resName = line[17:20]
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.confNum = int(line[27:30])
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.confType = "%3s%2s" % (self.resName, line[80:82])
        self.atomID = "%4s%3s%04d%c%03d" % (self.name, self.resName, self.resSeq, self.chainID, self.confNum)
        self.confID = "%5s%c%04d%c%03d" % (self.confType, self.chainID, self.resSeq, self.iCode, self.confNum)
        self.resID = "%3s%04d%c" % (self.resName, self.resSeq, self.chainID)

        # extended records
        connect_key = ("CONNECT", self.name, self.confType)
        self.connectivity_param = env.param[connect_key]
        radius_key = ("RADIUS", self.confType, self.name)
        radius_values = env.param[radius_key]
        self.r_bound = radius_values.r_bound
        self.r_vdw = radius_values.r_vdw
        self.e_vdw = radius_values.e_vdw
        return

    def writeline(self):
        line = "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f\n" % \
               (self.serial, self.name, self.altLoc, self.resName, self.chainID, \
                self.resSeq, self.iCode, self.xyz[0], self.xyz[1], \
                self.xyz[2])
        return line


class Conformer:
    def __init__(self):
        self.confID = ""
        self.resID = ""
        self.atom = []
        self.vdw1 = 0.0
        return

class Residue:
    def __init__(self):
        self.resID = ""
        self.conf = []
        return

class Protein:
    def __init__(self):
        self.residue = []
        self.vdw_pw = {}
        return

    def loadpdb(self, fname):
        rawlines = open(fname).readlines()
        lines = [x.strip("\n") for x in rawlines if x[:6] == "ATOM  " or x[:6] == "HETATM"]
        for line in lines:
            atom = Atom()
            atom.loadline(line)
            # reverse search if this atom belongs to an existing conformer and residue
            found_conf = False
            found_res = False
            for i_res in range(len(self.residue) - 1, -1, -1):
                insert_res = i_res
                for i_conf in range(len(self.residue[i_res].conf) - 1, -1, -1):
                    if self.residue[i_res].conf[i_conf].confID == atom.confID:
                        found_conf = True
                        self.residue[i_res].conf[i_conf].atom.append(atom)
                        break
                if found_conf:
                    break
                elif self.residue[i_res].resID == atom.resID:  # not in conf but residue ID matches
                    conf = Conformer()
                    conf.confID = atom.confID
                    conf.resID = atom.resID
                    conf.atom.append(atom)
                    self.residue[i_res].conf.append(conf)
                    found_res = True

            if not found_conf and not found_res:  # new residue
                conf = Conformer()
                conf.confID = atom.confID
                conf.resID = atom.resID
                conf.atom.append(atom)
                res = Residue()
                res.resID = conf.resID
                res.conf.append(conf)
                self.residue.append(res)

        # Insert an empty conformer for cofactors that do not have backbone
        for res in self.residue:
            first_confID = res.conf[0].confID
            if first_confID[3:5] != "BK":
                conf = Conformer()
                conf.confID = "%sBK%s000" % (first_confID[:3], first_confID[5:11])
                conf.resID = res.resID
                res.conf.insert(0, conf)
        return

    def make_connect12(self):
        # make connect table, need to start from connect records!
        # an optimized version should search connect12 for backbone and side chain separately.
        # backbone only connects the residue before and after, while side chain may connect with all other side chains
        # Special rules for backbone atoms:
        #    if ligated, search only residue before and after
        #    if connected atom is missing (CA), search residue before and after in case it is terminal residue
        #    all side chain conformers

        for i_res in range(len(self.residue)):
            res = self.residue[i_res]
            for conf in res.conf:  # search all conformers including backbone
                for atom in conf.atom:
                    connect_key = ("CONNECT", atom.name, atom.confType)
                    connected_atoms = env.param[connect_key].connected
                    for c_atom in connected_atoms:
                        found = False
                        if "?" in c_atom:  # ligated
                            for res2 in self.residue:
                                if res != res2:
                                    for conf2 in res2.conf:
                                        for atom2 in conf2.atom:
                                            connect_key2 = ("CONNECT", atom2.name, atom2.confType)
                                            connected_atoms2 = env.param[connect_key2].connected
                                            ligated2 = False
                                            for ligated_atom2 in connected_atoms2:
                                                if "?" in ligated_atom2:
                                                    ligated2 = True
                                                    break
                                            if ligated2:
                                                r = (atom.r_vdw + atom2.r_vdw) * BONDDISTANCE_scaling
                                                CUTOFF2 = r * r
                                                # if " C  " in atom.atomID:
                                                #     print("%s <-> %s: d=%.3f cut=%.3f" % (atom.atomID, atom2.atomID, math.sqrt(ddvv(atom.xyz, atom2.xyz)), math.sqrt(CUTOFF2)))
                                                #     print(atom2.atomID, [x.atomID for x in atom.connect12])
                                                #     print("Here", found)
                                                if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                                    if atom2 not in atom.connect12:
                                                        atom.connect12.append(atom2)
                                                        found = True  # after ligand found, do not break, continue to search other conformers within residue

                                    if found:   # one "?" for one ligand
                                        # print(conf.confID, atom.atomID, c_atom, res2.resID)
                                        break

                            if not found and res.resID[:3] == "NTR" and atom.name == " CA ":   # NTR CA connects to CB of next residue
                                res2 = self.residue[i_res+1]
                                # find " CB "
                                found = False
                                for conf2 in res2.conf[1:]:
                                    for atom2 in conf2.atom:
                                        if atom2.name == " CB ":
                                            r = (atom.r_vdw + atom2.r_vdw) * BONDDISTANCE_scaling
                                            CUTOFF2 = r * r
                                            if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                                if atom2 not in atom.connect12:
                                                    atom.connect12.append(atom2)
                                                    found = True

                            if not found and res.resID[:3] == "CTR" and atom.name == " C  ":   # CTR C connects to CA of previous atom
                                res2 = self.residue[i_res-1]
                                # find " CA "
                                found = False
                                for atom2 in res2.conf[0].atom:
                                    if atom2.name == " CA ":
                                        r = (atom.r_vdw + atom2.r_vdw) * BONDDISTANCE_scaling
                                        CUTOFF2 = r * r
                                        if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                            if atom2 not in atom.connect12:
                                                atom.connect12.append(atom2)
                                                found = True

                            if not found:  # no actual CTR case
                                if atom.name == " C  ":
                                    #print("No CTR after atom \"%s\"?" % atom.atomID)
                                    found = True

                            if not found:
                                if not "CTR" in atom.atomID:  # ignore CTR due to CA not specified as ligand
                                    print("Warning: Ligand atom bond to \"%s\" was not found" % atom.atomID)

                        else:  # a named atom
                            # 1) backbone
                            for atom2 in res.conf[0].atom:
                                if atom2.name == c_atom:
                                    atom.connect12.append(atom2)
                                    found = True
                                    break
                            # 2) own conformer
                            if not found:
                                for atom2 in conf.atom:
                                    if atom2.name == c_atom:
                                        atom.connect12.append(atom2)
                                        found = True
                                        break
                            # 3) a backbone atom could connect to side chain conformer atoms
                            if not found:
                                if conf.confID[3:5] == "BK":
                                    for conf2 in res.conf[1:]:
                                        for atom2 in conf2.atom:
                                            if atom2.name == c_atom:
                                                atom.connect12.append(atom2)
                                                found = True
                            # 4) NTR case, residue before
                            if not found:
                                if atom.name == " C  " and c_atom == " CA ":  # NTR separated from this res
                                    for conf2 in self.residue[i_res-1].conf[1:]:
                                        for atom2 in conf2.atom:
                                            if atom2.name == c_atom:
                                                atom.connect12.append(atom2)
                                                found = True
                                                break

                            # 5) CTR case, " C  " bond to " CA " not found
                            if not found:
                                if atom.name == " CA " and c_atom == " C  ":  # CTR is separated from this residue
                                    if i_res + 1 >= len(self.residue): # last residue
                                        found = True
                                    else:
                                        for conf2 in self.residue[i_res+1].conf[1:]:
                                            for atom2 in conf2.atom:
                                                if atom2.name == c_atom:
                                                    atom.connect12.append(atom2)
                                                    found = True
                                                    break

                            # 6) " CB " connects to " CA " of NTR case
                            if not found:
                                if atom.name == " CB " and c_atom == " CA ":
                                    res2 = self.residue[i_res - 1]
                                    for conf2 in res2.conf:
                                        for atom2 in conf2.atom:
                                            if atom2.name == " CA ":
                                                r = (atom.r_vdw + atom2.r_vdw) * BONDDISTANCE_scaling
                                                CUTOFF2 = r * r
                                                if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                                    if atom2 not in atom.connect12:
                                                        atom.connect12.append(atom2)
                                                        found = True
                                                        break

                            if not found:
                                print("Warning: Atom \"%s\" bond to \"%s\" was not found" % (c_atom, atom.atomID))

        return

    def print_connect12(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print(atom.atomID)
                    for atom2 in atom.connect12:
                        print("   -> %s" % atom2.atomID)
        return

    def make_connect13(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    for atom2 in atom.connect12:
                        if atom2 != atom:
                            for atom3 in atom2.connect12:
                                if (atom3 != atom) and (atom3 not in atom.connect12) and (atom3 not in atom.connect13):
                                    atom.connect13.append(atom3)
                        else:
                            print("Warning: Atom \"%s\" has itself in connect12." % atom.atomID)
        return

    def make_connect14(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    for atom3 in atom.connect13:
                        if atom3 == atom: continue
                        for atom4 in atom3.connect12:
                            if (atom4 != atom) \
                                    and (atom4 not in atom.connect12) \
                                    and (atom4 not in atom.connect13) \
                                    and (atom4 not in atom.connect14):
                                atom.connect14.append(atom4)
        return

    def print_connect13(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print(atom.atomID)
                    for atom2 in atom.connect13:
                        print("   -X-> %s" % atom2.atomID)
        return

    def print_connect14(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print(atom.atomID)
                    for atom2 in atom.connect14:
                        print("   -X-X-> %s" % atom2.atomID)
        return

    def exportpdb(self, fname):
        lines = []
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    lines.append(atom.writeline())
        open(fname, "w").writelines(lines)
        return

    def print_atom_structure(self):
        for res in self.residue:
            print("Residue %s" % res.resID)
            for conf in res.conf:
                print("-->Conformer %s" % conf.confID)
                for atom in conf.atom:
                    print("---->Atom %s" % atom.atomID)
        return

    def calc_vdw(self):
        # do it on two sides so the two-way interaction numbers can be checked
        for res1 in self.residue:
            for conf1 in res1.conf:
                for res2 in self.residue:
                    if res1 == res2: # we need to do self to self vdw - vdw0
                        vdw = vdw_conf(conf1, conf1)
                        if abs(vdw) > 0.001:
                            self.vdw_pw[(conf1.confID, conf1.confID)] = vdw
                        continue
                    for conf2 in res2.conf:
                        vdw = vdw_conf(conf1, conf2)
                        if abs(vdw) > 0.001:
                            #print("%s - %s: %.3f" % (conf1.confID, conf2.confID, vdw))
                            self.vdw_pw[(conf1.confID, conf2.confID)] = vdw

                # compute vdw1, vdw to all backbone
                vdw1 = 0.0
                for res2 in self.residue:
                    conf2 = res2.conf[0]
                    vdw = vdw_conf(conf1, conf2)
                    vdw1 += vdw
                conf1.vdw1 = vdw1

    def connect_reciprocity_check(self):
        # connectivity should be reciprocal except backbone atoms
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    for atom2 in atom.connect12:
                        if atom not in atom2.connect12:
                            print("Atom %s in connect12 of atom %s but the other way is not true" % (atom2.atomID, atom.atomID))
                    for atom2 in atom.connect13:
                        if atom2.confType[-2:] == "BK":
                            continue
                        if atom not in atom2.connect13:
                            print("Atom %s in connect13 of atom %s but the other way is not true" % (atom2.atomID, atom.atomID))
                    for atom2 in atom.connect14:
                        if atom2.confType[-2:] == "BK":
                            continue
                        if atom not in atom2.connect14:
                            print("Atom %s in connect14 of atom %s but the other way is not true" % (atom2.atomID, atom.atomID))

    def vdw_reciprocity_check(self):
        for key, value in self.vdw_pw.items():
            r_key = (key[1], key[0])
            if r_key in self.vdw_pw:
                r_value = self.vdw_pw[r_key]
                if abs(value - r_value) > 0.001:
                    print("vdw(%s<->%s = %.3f <-> %.3f" % (key[0], key[1], value, r_value))
            else:
                print("vdw(%s->%s = %.3f but other way not reported" % (key[0], key[1], value))

        print("VDW two sides checked.")
        return

class CONNECT_param:
    def __init__(self, value_str):
        fields = value_str.split(",")
        self.orbital = fields[0].strip()
        self.connected = [x.strip().strip("\"") for x in fields[1:]]


class RADIUS_param:
    def __init__(self, value_str):
        fields = value_str.split(",")
        self.r_bound = float(fields[0])
        self.r_vdw = float(fields[1])
        self.e_vdw = float(fields[2])


class ENV:
    def __init__(self):
        self.runprm = {}
        self.param = {}
        self.load_runprm()
        #self.print_runprm()
        self.load_ftpl()

    def load_runprm(self):
        filename = "run.prm"
        lines = open(filename).readlines()
        for line in lines:
            entry_str = line.strip().split("#")[0]
            fields = entry_str.split()
            if len(fields) > 1:
                key_str = fields[-1]
                if key_str[0] == "(" and key_str[-1] == ")":
                    key = key_str.strip("()").strip()
                    value = fields[0]
                self.runprm[key] = value

    def print_runprm(self):
        for key, value in self.runprm.items():
            print("%s:%s" % (key, value))

    def load_ftpl(self):
        ftpldir = self.runprm["MCCE_HOME"]+"/param"
        cwd = os.getcwd()
        os.chdir(ftpldir)

        files = glob.glob("*.ftpl")
        files.sort()

        for fname in files:
            lines = open(fname).readlines()
            for line in lines:
                end = line.find("#")
                line = line[:end]
                fields = line.split(":")
                if len(fields) != 2:
                    continue

                key_string = fields[0].strip()
                keys = key_string.split(",")
                key1 = keys[0].strip().strip("\"")
                if len(keys) > 1:
                    key2 = keys[1].strip().strip("\"")
                else:
                    key2 = ""
                if len(keys) > 2:
                    key3 = keys[2].strip().strip("\"")
                else:
                    key3 = ""

                value_string = fields[1].strip()

                # Connectivity records
                if key1 == "CONNECT":
                    self.param[(key1,key2,key3)] = CONNECT_param(value_string)
                # VDW parameters
                elif key1 == "RADIUS":
                    self.param[(key1, key2, key3)] = RADIUS_param(value_string)

        os.chdir(cwd)

    def print_param(self):
        for key, value in self.param.items():
            key1, key2, key3 = key
            if key1 == "CONNECT":
                print("%s:%s, %s" % (key, value.orbital, value.connected))
            elif key1 == "RADIUS":
                print("%s: %6.3f, %6.3f %6.3f" % (key, value.r_bound, value.r_vdw, value.e_vdw))


def vdw_conf(conf1, conf2, verbose=False):
    vdw = 0.0
    for atom1 in conf1.atom:
        for atom2 in conf2.atom:
            vdw_a2a = vdw_atom(atom1, atom2)
            vdw += vdw_a2a
            if verbose:
                print("%s -> %s: %.3f" % (atom1.atomID, atom2.atomID, vdw_a2a))
    if vdw >= VDW_UPLIMIT:
        vdw = 999.0
    return vdw

def vdw_atom(atom1, atom2):
    # A good post: https://mattermodeling.stackexchange.com/questions/4845/how-to-create-a-lookup-table-of-%CF%B5-and-%CF%83-values-for-lennard-jones-potentials
    # Parameter source: http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
    #   CHARM36-jul2022.ff/ffnonbonded.itp
    #   σ and ε values are in nm and kJ/mol in this file
    # Lorentz-Berthelot combining rules:
    #   σij = 0.5*(σi + σj)
    #   ϵij = sqrt(ϵi * ϵj)
    # p_lj = ϵij[(σij/r)^12 - 2(σij/r)^6]
    # σij is the distance where LJ potential reaches minimum: -ϵij
    # r is the atom distance

    # Using the equation in vdw.c. Need to work on new parameter set.
    d2 = ddvv(atom1.xyz, atom2.xyz)

    if atom2 not in atom1.connect12 and atom2 not in atom1.connect13:
        if d2 > VDW_CUTOFF_FAR2:
            p_lj = 0.0
        elif d2 < VDW_CUTOFF_NEAR2:
            p_lj = 999.0
        else:
            r1 = atom1.r_vdw
            e1 = atom1.e_vdw
            r2 = atom2.r_vdw
            e2 = atom2.e_vdw
            if atom2 in atom1.connect14:
                scale = VDW_SCALE14
            else:
                scale = 1.0

            sig_min = r1 + r2
            eps = math.sqrt(e1 * e2)

            sig_d2 = sig_min * sig_min / d2
            sig_d6 = sig_d2 * sig_d2 * sig_d2
            sig_d12 = sig_d6 * sig_d6

            p_lj = eps * sig_d12 - 2. * eps * sig_d6
    else:
        p_lj = 0.0

    return p_lj


def vdw_by_conf_pair(protein, confID1, confID2):
    # detailed vdw between conf and conf
    for res1 in protein.residue:
        for conf1 in res1.conf:
            if conf1.confID != confID1:
                continue
            for res2 in protein.residue:
                for conf2 in res2.conf:
                    if conf2.confID != confID2:
                        continue
                    vdw = vdw_conf(conf1, conf2, verbose=True)
                    print("%s - %s: %.3f" % (conf1.confID, conf2.confID, vdw))

    return

env = ENV()


if __name__ == "__main__":
    #env.print_param()

    pdbfile = "step2_out.pdb"
    protein = Protein()
    protein.loadpdb(pdbfile)
    protein.make_connect12()
    protein.make_connect13()
    protein.make_connect14()

    # protein.print_connect12()
    # protein.print_connect13()
    # protein.print_connect14()
    # protein.exportpdb("a.pdb")
    # protein.print_atom_structure()
    protein.calc_vdw()
    protein.connect_reciprocity_check()
    protein.vdw_reciprocity_check()

    #vdw_by_conf_pair(protein, "GLYBKA0006_000", "GLYBKA0006_000", 0.001)
    #vdw_by_conf_pair(protein, "ASPBKA0002_000", "NTG01A0001_001", 0.001)


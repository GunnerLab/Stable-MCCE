#!/usr/bin/env python

import argparse, sys, math

# covalent radius to decide a bond. bond length: r1+r2
radius = {" H": 0.25,
          " N": 0.65,
          " C": 0.70,
          " O": 0.60,
          " P": 1.00,
          " S": 1.00,
          "NA": 1.80,
          "CL": 1.00
          }

elebd_radius = {" N": 1.5,
                " H": 1.0,
                " C": 1.7,
                " O": 1.4,
                " P": 1.85,
                " S": 1.85,
                " X": 1.85
}

vdw_parm = {" C": (2.000, 0.150),
            " H": (1.000, 0.020),
            " O": (1.600, 0.200),
            " N": (1.750, 0.160),
            " S": (2.000, 0.200),
            " P": (2.000, 0.200),
            " X": (2.000, 0.173)
            }

sp_orbitals = [" C", " N", " O", " P", " S"]
spd_orbitals = ["FE"]

tolerance_scale = 1.3  # (r1+r2) * this scale gives the bond upper limit, value between 1.2 to 1.5 recommended


def vector_normalize(v):
    vn = ()
    d = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if d < 1.0e-20:
        vn = (0.0, 0.0, 0.0)
    else:
        vn = (v[0] / d, v[1] / d, v[2] / d)

    return vn


def avv(v1, v2):
    v1 = vector_normalize(v1)
    v2 = vector_normalize(v2)
    t = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    if t > 1.0:
        t = 1.0
    elif t < -1.0:
        t = -1.0

    return math.acos(t)


def dvv(v1, v2):
    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    d2 = dx * dx + dy * dy + dz * dz
    return math.sqrt(d2)


def vector_vminusv(v1, v2):
    z = (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])
    return z


class Atom():
    def __init__(self):
        self.name = ""
        self.element = ""
        self.orbital = ""
        self.xyz = ()
        self.connect = []
        return

    def loadline(self, line):
        self.name = line[12:16]
        self.element = line[76:78]
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        return


class Pdb2ftpl:
    def __init__(self, arguments):
        self.resname = []
        self.atoms = self.file2atoms(arguments.pdbfile[0])
        if len(self.resname) != 1:
            print("%d residue names detected. The input pdb file can only have one residue." % len(self.resname))
            sys.exit()

        if arguments.d:
            records = self.records_from_distance()
        else:
            records = self.records_from_file(arguments.pdbfile[0])
        if not records:
            records = self.records_from_distance()

        # find connected atoms
        self.connect_from_records(records)

        # find hybrid type from connected atoms and bond angle
        self.make_bond_orbital()

        return

    def file2atoms(self, fname):
        atoms = []
        lines = open(fname).readlines()
        for line in lines:
            if len(line) < 54:
                continue
            if line[:6] == "ATOM  " or line[:6] == "HETATM":
                res = line[17:20]
                if res not in self.resname:
                    self.resname.append(res)
                atom = Atom()
                atom.loadline(line)
                atoms.append(atom)
        return atoms

    def records_from_file(self, fname):
        connect = {}
        lines = open(fname).readlines()
        for line in lines:
            if len(line) < 11:
                continue
            if line[:6] == "CONECT":
                fields = line.split()
                key = int(fields[1]) - 1
                value = [int(x) - 1 for x in fields[2:]]
                connect[key] = value
        return connect

    def records_from_distance(self):
        connect = {}
        n = len(self.atoms)
        for i in range(n - 1):
            for j in range(i + 1, n):
                d = dvv(self.atoms[i].xyz, self.atoms[j].xyz)
                bond_distance = radius[self.atoms[i].element] + radius[self.atoms[j].element]
                if d < bond_distance * tolerance_scale:
                    if i in connect:
                        connect[i].append(j)
                    else:
                        connect[i] = [j]
                    if j in connect:
                        connect[j].append(i)
                    else:
                        connect[j] = [i]

        return connect

    def connect_from_records(self, records):
        for i in range(len(self.atoms)):
            atom = self.atoms[i]
            if i in records.keys():
                # print(i, records[i])
                atom.connect = [self.atoms[k] for k in records[i]]
            else:
                atom.connect = []
        return

    def make_bond_orbital(self):
        for atom in self.atoms:
            if atom.element == " H":
                atom.orbital = "s"
            elif atom.element in sp_orbitals:  # sp, sp2, sp3
                if len(atom.connect) == 4:  # sp3
                    atom.orbital = "sp3"
                elif len(atom.connect) == 3 or len(atom.connect) == 2:
                    # bond angle: sp if 180, sp2 if 120, sp3 if 109
                    v1 = vector_vminusv(atom.connect[0].xyz, atom.xyz)
                    v2 = vector_vminusv(atom.connect[1].xyz, atom.xyz)
                    alpha = avv(v1, v2) / math.pi * 180
                    # print(atom.connect[0].name, atom.name, atom.connect[1].name, alpha)
                    if 104 < alpha < 115:
                        atom.orbital = "sp3"
                    elif 115 < alpha < 150:
                        atom.orbital = "sp2"
                    elif 150 < alpha < 180.1:
                        atom.orbital = "sp"
                    else:
                        print("%s - %s - %s bond angle = %.3f, can not interpret" % (
                        atom.connect[0].name, atom.name, atom.connect[1].name, alpha))
                        sys.exit()
                elif len(atom.connect) == 1:
                    # doesn't matter in the sense of geometry, but O on CH3-(CO)-CH3 is sp2 instead of sp3.
                    atom.orbital = "sp3"
                elif len(atom.connect) == 0:
                    atom.orbital = "ion"
                else:
                    atom.orbital = "udf"

        return

    def print_conflist(self):
        print("# Conformer definition")
        print("CONFLIST, %s: %s01" % (self.resname[0], self.resname[0]))

    def print_connect(self):
        print("# ATOM name and bonds")
        for atom in self.atoms:
            connected_atoms = ",".join(["\"%s\"" % x.name for x in atom.connect])
            print("CONNECT, \"%s\", %s01: %4s, %s" % (atom.name, self.resname[0], atom.orbital, connected_atoms))

    def print_charge(self):
        print("# ATOM charges")
        for atom in self.atoms:
            print("CHARGE, %s01, \"%s\": to_be_filled" % (self.resname[0], atom.name))

    def print_radius(self):
        print("# Atom radius, dielelctric boundary radius, VDW radius, and energy well depth")
        for atom in self.atoms:
            if atom.element in elebd_radius:
                rbd = elebd_radius[atom.element]
            else:
                rbd = elebd_radius[" X"]

            if atom.element in vdw_parm:
                rvdw, well = vdw_parm[atom.element]
            else:
                rvdw, well = vdw_parm[" X"]

            print("RADIUS, %s01, \"%s\": %6.3f, %6.3f, %6.3f" % (self.resname[0], atom.name, rbd, rvdw, well))

    def print_conformer(self):
        print("# Conformer parameters that appear in head3.lst: ne, Em0, nH, pKa0, rxn")
        print("CONFORMER, %s01:  Em0=0.0, pKa0=0.00, ne=0, nH=0, rxn02= to_be_filled, rxn04= to_be_filled, rxn08= to_be_filled" % (self.resname[0]))


if __name__ == "__main__":
    # Get the command arguments
    helpmsg = "Create a ftpl template file from a cofactor PDB file. The atoms in the input files are considered as one molecule."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-d", default=False, help="Ignore CONNECT, use distance to determine bond", action="store_true")
    parser.add_argument("pdbfile", metavar="pdbfile", nargs=1)
    args = parser.parse_args()

    ftpl = Pdb2ftpl(args)
    ftpl.print_conflist()
    print()
    ftpl.print_connect()
    print()
    ftpl.print_charge()
    print()
    ftpl.print_radius()
    print()
    ftpl.print_conformer()
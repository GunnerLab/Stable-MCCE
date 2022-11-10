#!/usr/bin/env python
# center and randomly rotate protein

import string

import numpy as np

from geom import *
import sys

class ATOM():
    def load_pdb(self, line):
        self.name = line[12:16]
        self.resname = line[17:20]
        self.chainid = line[21]
        self.seq = int(line[22:26].strip())
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def atoms2pdb(atoms):
    atom_count = 1
    pdblines = []
    for atom in atoms:
        pdblines.append("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f\n" % (atom_count, atom.name, atom.resname, atom.chainid,
                                                                atom.seq, atom.xyz[0], atom.xyz[1], atom.xyz[2]))
        if atom_count < 99999:
            atom_count += 1
    return pdblines


def center_protein(atoms):
    x_lo = x_up = atoms[0].xyz[0]
    y_lo = y_up = atoms[1].xyz[1]
    z_lo = z_up = atoms[2].xyz[2]
    for atom in atoms[1:]:
        if x_lo > atom.xyz[0]: x_lo = atom.xyz[0]
        if x_up < atom.xyz[0]: x_up = atom.xyz[0]
        if y_lo > atom.xyz[1]: y_lo = atom.xyz[1]
        if y_up < atom.xyz[1]: y_up = atom.xyz[1]
        if z_lo > atom.xyz[2]: z_lo = atom.xyz[2]
        if z_up < atom.xyz[2]: z_up = atom.xyz[2]

    c = ((x_lo+x_up)/2, (y_lo+y_up)/2, (z_lo+z_up)/2)
    op = OPERATION()
    op.move((-c[0], -c[1], -c[2]))

    for atom in atoms:
        atom.xyz = geom_apply(op, atom.xyz)

    return

if __name__ == "__main__":
    mol_atoms = []

    print("Reading input structure from %s ... " % sys.argv[1], end="", flush=True)
    lines=open(sys.argv[1]).readlines()
    for line in lines:
        if line[:6] == "ATOM  " or line[:6] == "HETATM":
            atom = ATOM()
            atom.load_pdb(line)
            mol_atoms.append(atom)
    print("Done")

    print("Move protein center to (0, 0, 0) ... ", end="", flush=True)
    center_protein(mol_atoms)
    print("Done")


    output_lines = atoms2pdb(mol_atoms)

    fname_str = sys.argv[1]
    fname = fname_str.split(".")[0]
    open("%s_center.pdb" % fname, "w").writelines(output_lines)
    print("Structure moved to origin is saved in %s_center.pdb" % fname)

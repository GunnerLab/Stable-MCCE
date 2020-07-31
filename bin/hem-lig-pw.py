#!/usr/bin/env python
"""
Reset hem to ligand pairwise interaction to 0.
If the the ligand pair is in the list, reset their pairwise interactions in opp files tp 0.
"""

import argparse
import math
import os
import sys

opp_dir = "energies"
head3lst ="head3.lst"
step2_out = "step2_out.pdb"

BOND_threshold = 2.7

ligands = ["HIL", "MEL"]
receptors = ["HEB", "HEC", "HEM"]


class Atom:
    def __init__(self):
        self.icount = 0
        self.name = ""
        self.element = ""  # two-letter code for element name
        self.resname = ""
        self.chainid = ""
        self.seqnum = 0
        self.icode = ""
        self.xyz = ()
        self.confname = ""
        return

    def loadline(self, line):
        self.icount = int(line[6:11])
        self.name = line[12:16]
        self.resname = line[17:20]
        self.chainid = line[21]
        self.seqnum = int(line[22:26])
        self.icode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.confname = "%3s%2s%9s" % (self.resname, line[80:82], line[21:30])

        if len(self.name.strip()) < 4:
            self.element = self.name[:2]
        else:
            self.element = " H"

        return



def ddvv(v1, v2):
    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    return dx*dx + dy*dy + dz*dz

def shortest_d(res1_atoms, res2_atoms):
    ddmin = 1.0E10
    for atom1 in res1_atoms:
        if atom1.element == " H":
            continue
        for atom2 in res2_atoms:
            if atom2.element == " H":
                continue
            dd = ddvv(atom1.xyz, atom2.xyz)
            if ddmin > dd:
                ddmin = dd
    return math.sqrt(ddmin)


def match_conf(atom, conformers):
    matched = []
    for conf in conformers:
        if (conf[:3], conf[5], int(conf[6:10]), conf[10]) == atom.resid:
            matched.append(conf)
    return matched

def set0(fname, conformer):
    if os.path.exists(fname):
        opplines = open(fname).readlines()
        newlines = []
        for line in opplines:
            fields = line.split()
            if len(fields) > 3:
                if fields[1] == conformer:
                    newline = "%s  +0.000   0.000   0.000   0.000\n" % (line[:21])
                    newlines.append(newline)
                else:
                    newlines.append(line)
            else:
                newlines.append(line)
        open(fname, "w").writelines(newlines)

    return

def reset_pw():
    conformers = []
    lines = open(head3lst).readlines()
    lines.pop(0)
    for line in lines:
        fields = line.split()
        if len(fields) < 3:
            continue
        conf = fields[1]
        conformers.append(conf)

    ligand_confs = []
    for conf in conformers:
        if conf[:3] in ligands:
            ligand_confs.append(conf)

    receptor_confs = []
    for conf in conformers:
        if conf[:3] in receptors:
            receptor_confs.append(conf)

    atoms = []
    lines = open(step2_out).readlines()
    for line in lines:
        if line[:6] == "ATOM  " or line[:6] == "HETATM":
            atom = Atom()
            atom.loadline(line)
            atoms.append(atom)

    for ligand_conf in ligand_confs:
        ligand_conf_atoms = []
        for atom in atoms:
            if ligand_conf == atom.confname:
                ligand_conf_atoms.append(atom)
        for receptor_conf in receptor_confs:
            receptor_conf_atoms = []
            for atom in atoms:
                if receptor_conf == atom.confname:
                    receptor_conf_atoms.append(atom)
            if shortest_d(ligand_conf_atoms, receptor_conf_atoms) < BOND_threshold:
                set0("energies/%s.opp" % ligand_conf, receptor_conf)
                set0("energies/%s.opp" % receptor_conf, ligand_conf)

    return


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Set pairwise interaction between ligands to be 0. The ligands are matched by both name and distance."
    parser = argparse.ArgumentParser(description=helpmsg)
    args = parser.parse_args()
    reset_pw()
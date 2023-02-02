#!/usr/bin/env python
"""
Compute vdw pairwise and write back to opp files
"""

import os
from pdbio import *

class ConfSelf:
    def __init__(self, line):
        fields = line.split()
        self.iconf = int(fields[0])
        self.conformer = fields[1]
        self.fl = fields[2]
        self.occ = float(fields[3])
        self.crg = float(fields[4])
        self.em0 = int(fields[5])
        self.pka0 = float(fields[6])
        self.ne = int(fields[7])
        self.nh = int(fields[8])
        self.vdw0 = float(fields[9])
        self.vdw1 = float(fields[10])
        self.tors = float(fields[11])
        self.epol = float(fields[12])
        self.dsolv = float(fields[13])
        self.extra = float(fields[14])
        self.history = fields[15]
        self.state = fields[16]

def update_opp(protein):
    # read head3.lst
    lines = open("head3.lst").readlines()
    for line in lines[1:]:
        conf_self = ConfSelf(line)




    for res in protein.residue:
        for conf in res.conf[1:]:
            fname = "energies/" + conf.confID + ".opp"
            if os.path.isfile(fname):
                lines = open(fname).readlines()

    return

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
    # protein.print_atom_structure()
    protein.calc_vdw()
    # protein.connect_reciprocity_check()
    # protein.vdw_reciprocity_check()
    update_opp(protein)

#!/usr/bin/env python
"""
Compute vdw pairwise and write back to opp files
"""

from pdbio import *

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
    protein.connect_reciprocity_check()
    protein.vdw_reciprocity_check()

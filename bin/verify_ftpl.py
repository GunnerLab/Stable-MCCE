#!/usr/bin/env python
"""
This program verifies one ftpl file over the ligand pdb file, and writes out the mcce pdb file.
"""

import sys
from pymccelib import *
import logging

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO, format='%(levelname)-s: %(message)s')

    if len(sys.argv) < 3:
        print("verify_ftpl.py tplfile pdbfile")
        print("Verify a ftpl over a pdb file, writes out step1_out.pdb")
        sys.exit()

    ftpl = sys.argv[1]
    pdb = sys.argv[2]
    env = Env()
    env.load_ftpl(ftpl)
    prot = Protein()

    # Verify if pdb can be loaded correctly

    # Verify VDW are present

    # Verify all atoms are connected within a residue

    # Verify rotate parameters exist for single bond

    # Verify torsion paramaters exist

    # Verify CONFORMER parameters are complete
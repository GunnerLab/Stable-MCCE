#!/usr/bin/env python
"""
This program analyze results.

Input:
 * fort.38
 * head3.lst

Output:
 * sum_crg.out
 * pK.out

Usage examples:

1. Run step 5 with default
    step5.py

2. Run step 5 with entropy correction
    step5.py --xts

3. Run step 5 at defined pH/Eh or mid-point
    step5.py -p 7
    step5.py -p m
"""

import argparse
import numpy as np

PH2KCAL = 1.364

class Titration:
    def __init__(self, args):
        self.titration_type = "PH"
        self.xts = args.xts
        self.mfe_point = args.p
        self.titration_points = []
        conformers = self.load_confs()
        self.residues = self.group_residues(conformers)

    def load_confs(self):
        fort38_lines = open("fort.38").readlines()
        head3_lines = open("head3.lst").readlines()
        line = fort38_lines.pop(0)
        fields = line.split()
        self.titration_type = fields[0].upper()
        self.titration_points = [float(x) for x in fields[1:]]
        head3_lines.pop(0)

    def group_residues(self, conformers):

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Run mcce step 5, generate net charge, fit titration curve, and do energy analysis on each ionizable residue."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("--xts", default=False, help="Enable entropy correction, default is false", action="store_true")
    parser.add_argument("-p", metavar="titration point", default="m", help="pH or Eh value, or \'m\' for midpoint")
    args = parser.parse_args()

    titration = Titration(args)
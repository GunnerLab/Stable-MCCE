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
import sys
import numpy as np

PH2KCAL = 1.364

class Conformer:
    def __init__(self):
        self.name = ""
        self.occ = []
        self.crg = 0.0
        self.em0 = 0.0
        self.pka0 = 0.0
        self.ne = 0
        self.nh = 0
        self.vdw0 = 0.0
        self.vdw1 = 0.0
        self.tors = 0.0
        self.epol = 0.0
        self.dsolv = 0.0
        self.extra = 0.0
        return

class Residue:
    def __init__(self):
        self.resid = ""
        self.charge_type = "0"
        self.conformers = []
        self.netcrg = []
        self.ground_state = []
        self.charged_state = []

        return

class Titration:
    def __init__(self, args):
        self.titration_type = "PH"
        self.xts = args.xts
        self.mfe_point = args.p
        self.titration_points = []
        conformers = self.load_confs()
        self.residues = self.group_residues(conformers)

        # get sum_crg
        npoints = len(self.titration_points)

        for res in self.residues:
            net_crg = [0.0 for i in range(npoints)]
            for conf in res.conformers:
                net_crg = [net_crg[i]+conf.crg*conf.occ[i] for i in range(npoints)]
            res.netcrg = net_crg



        return


    def load_confs(self):
        conformers = []

        fort38_lines = open("fort.38").readlines()
        head3_lines = open("head3.lst").readlines()
        line = fort38_lines.pop(0)
        fields = line.split()
        self.titration_type = fields[0].upper()
        self.titration_points = [float(x) for x in fields[1:]]
        head3_lines.pop(0)

        while fort38_lines and head3_lines:
            fort38_data = fort38_lines.pop(0)
            head3_data = head3_lines.pop(0)
            conf = Conformer()

            fields = fort38_data.split()
            conf.name = fields[0]
            conf.occ = [float(x) for x in fields[1:]]

            fields = head3_data.split()
            if conf.name != fields[1]:
                print("Confomer %s in fort.38 does not match %s in head3.lst" % (conf.name, fields[1]))
                sys.exit()

            conf.crg = float(fields[4])
            conf.em0 = float(fields[5])
            conf.pka0 = float(fields[6])
            conf.ne = int(fields[7])
            conf.nh = int(fields[8])
            conf.vdw0 = float(fields[9])
            conf.vdw1 = float(fields[10])
            conf.tors = float(fields[11])
            conf.epol = float(fields[12])
            conf.dsolv = float(fields[13])
            conf.extra = float(fields[14])
            conformers.append(conf)

        return conformers


    def group_residues(self, conformers):
        res_ids = []
        residues = {}
        for conf in conformers:
            resid = conf.name[:3]+conf.name[5:11]
            if resid in res_ids:
                res = residues[resid]
            else:
                res = Residue()
                res.resid = resid
                residues[resid] = res
                res_ids.append(resid)
            res.conformers.append(conf)


        return [residues[x] for x in res_ids]


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Run mcce step 5, generate net charge, fit titration curve, and do energy analysis on each ionizable residue."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("--xts", default=False, help="Enable entropy correction, default is false", action="store_true")
    parser.add_argument("-p", metavar="titration point", default="m", help="pH or Eh value, or \'m\' for midpoint")
    args = parser.parse_args()

    titration = Titration(args)
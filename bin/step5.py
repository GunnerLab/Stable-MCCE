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
from scipy.optimize import curve_fit

from sanity import *

PH2KCAL = 1.364
fname_sumcrg = "sum_crg2.out"
fname_pkout = "pK2.out"

def sigmoid(x, x0, k):
    e = np.exp(k * (x - x0))
    y = -e / (1.0 + e)
    return y


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
        self.state_flag = "0"
        self.conformers = []
        self.netcrg = []
        self.ground_state = []
        self.charged_state = []

        return

class Titration:
    def __init__(self, args):
        self.titration_type = "pH"
        self.titration_points = []
        self.residues = self.group_residues(self.load_confs())

        self.sum_crg()

        self.fitpka()

        return


    def sum_crg(self):
        # get sum_crg and charge originated from 1) proton, 2) electron - redox, 3) ion
        npoints = len(self.titration_points)
        total_crg = [0 for i in range(npoints)]
        total_proton_crg = [0 for i in range(npoints)]
        total_redox_crg = [0 for i in range(npoints)]
        total_ion_crg = [0 for i in range(npoints)]
        for res in self.residues:
            net_crg = [0.0 for i in range(npoints)]
            for conf in res.conformers:
                net_crg = [net_crg[i]+conf.crg*conf.occ[i] for i in range(npoints)]
                total_crg = [total_crg[i] + conf.crg*conf.occ[i] for i in range(npoints)]
                total_proton_crg = [total_proton_crg[i] + conf.nh*conf.occ[i] for i in range(npoints)]
                total_redox_crg = [total_redox_crg[i] - conf.ne*conf.occ[i] for i in range(npoints)]
                total_ion_crg = [total_ion_crg[i] + (conf.crg-conf.nh+conf.ne)*conf.occ[i] for i in range(npoints)]
            res.netcrg = net_crg

        # prepare for print
        lines = []
        headline = " %3s       %s\n" % (self.titration_type, " ".join(["%5.2f" % x for x in self.titration_points]))
        lines.append(headline)
        for res in self.residues:
            has_negative = False
            has_positive = False
            for conf in res.conformers:
                if conf.crg < -0.999:
                    has_negative = True
                if conf.crg > 0.999:
                    has_positive = True

            if has_positive and has_negative:
                state_flag = "*"
            elif has_negative:
                state_flag = "-"
            elif has_positive:
                state_flag = "+"
            else:
                state_flag = "0"

            res.state_flag = state_flag
            if state_flag != "0":    # print charged residues only
                res_str = res.resid[:3] + state_flag + res.resid[3:]
                lines.append("%s %s\n" % (res_str, " ".join(["%5.2f" % x for x in res.netcrg])))
        lines.append("---------\n")
        lines.append("Net_charge %s\n" % (" ".join(["%5.2f" % x for x in total_crg])))
        lines.append("Proton_crg %s\n" % (" ".join(["%5.2f" % x for x in total_proton_crg])))
        lines.append("Redox_crg  %s\n" % (" ".join(["%5.2f" % x for x in total_redox_crg])))
        lines.append("Ion_charge %s\n" % (" ".join(["%5.2f" % x for x in total_ion_crg])))
        open(fname_sumcrg, "w").writelines(lines)
        return

    def fitpka(self):
        titration_type = self.titration_type
        npoints = len(self.titration_points)
        if npoints < 2:
            print("Single point titration detected, can not fit titration curve.")
            sys.exit()

        titration_delta = abs(self.titration_points[1] - self.titration_points[0])
        titration_lb = self.titration_points[0]
        titration_ub = self.titration_points[-1]

        pkout = ["  %s             pKa/Em  n(slope) 1000*chi2\n" % titration_type]
        xvalue = np.array(self.titration_points)
        x0 = xvalue[0]
        k0 = xvalue[1] - xvalue[0]

        xdata = np.array([float(i) for i in range(npoints)])
        for res in self.residues:
            if res.state_flag != "0":
                res_str = res.resid[:3] + res.state_flag + res.resid[3:]

                # ydata of sigmoid function is always positive [0, 1] but our net charge is [-1, 0]
                ydata = np.array(res.netcrg)
                # print(xdata)
                # print("Before:",ydata)
                if ydata.mean() > 0:
                    ydata += -1.0  # acid[0, -1], base[+1, 0] -> [0, -1]
                if ydata[-1] - ydata[0] > 0.0:
                    ydata = -ydata - 1
                    # print("After:", ydata)
                msg = ""
                try:
                    (popt, pcov) = curve_fit(sigmoid, xdata, ydata, bounds=(0, [npoints, 4]))
                    # (popt, pcov) = curve_fit(sigmoid, xdata, ydata)
                except RuntimeError:
                    msg = "Titration out of range"
                except ValueError:
                    msg = "Input value not valid"

                #print(res_str, popt)
                if popt[0] < 0.001 or popt[0] > npoints - 1.001:  # x from 0 to 14 as 15 points
                    msg = "Titration out of range"

                if msg:
                    pkout.append("%s         %s\n" % (res_str, msg))
                else:
                    chi_squared = np.sum([(sigmoid(xdata, *popt) - ydata) ** 2])
                    midpoint = popt[0] * k0 + x0
                    if titration_type.upper() == "PH":
                        nslope = 0.4342 * popt[1] / titration_delta
                    elif titration_type.upper() == "EM":
                        nslope = 0.4342 * popt[1] / titration_delta * 58.0
                    elif titration_type.upper() == "CH" or titration_type.upper() == "EXTRA":
                        nslope = 0.4342 * popt[1] / titration_delta * PH2KCAL
                    else:
                        print("Why am I here?")
                    pka_str = "%9.3f %9.3f %9.3f" % (midpoint, nslope, chi_squared * 1000.0)
                    pkout.append("%s    %s\n" % (res_str, pka_str))

            open(fname_pkout, "w").writelines(pkout)

        return

    def load_confs(self):
        conformers = []

        fort38_lines = open("fort.38").readlines()
        head3_lines = open("head3.lst").readlines()
        line = fort38_lines.pop(0)
        fields = line.split()
        self.titration_type = fields[0]
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
    helpmsg = "Run mcce step 5, generate net charge, fit titration curve."
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", metavar="rule_file", default="", help='''User defined rule file. File format example:
# Residue_name: lower_bound, upper_bound
ASP-: 1, 6
HIS+: 6, 9
LYS+: 9, 13
ARG+: 10, 14
HEM+: 100, 400       
''')
    args = parser.parse_args()

    print("Compute charge of residues and titration curve ...", end = " ")
    titration = Titration(args)
    print("Done.")

    print("Checking abnormal ionizations ...")
    sanity_check(args)
    print("Done.")
#!/usr/bin/python

import sys
import numpy as np
from scipy.optimize import curve_fit


def sigmoid(x, x0, k):
    e = np.exp(k*(x-x0))
    y = -e/(1.0+e)
    return y


class Fort38:
    def __init__(self, fname):
        lines = open(fname).readlines()
        headline = lines.pop(0).strip()
        fields = headline.split()
        self.titration_type = fields[0].strip()
        self.titration_points = np.array(fields[1:], dtype=float)
        self.titration_step = self.titration_points[1] - self.titration_points[0]
        self.titration_lb = self.titration_points[0]
        self.titration_ub = self.titration_points[-1]
        self.occ_table = []
        self.conf_name = []
        for line in lines:
            fields = line.split()
            self.conf_name.append(fields[0])
            self.occ_table.append([float(x) for x in fields[1:]])
        return

    def load_charge(self, head3lst):
        occ_table = np.array(self.occ_table, dtype=np.float)
        crg_state = np.array(head3lst.conf_crg, dtype=np.float)
        self.crg_table = occ_table * crg_state[:, np.newaxis]
        return

class Head3List:
    def __init__(self, fname):
        lines = open(fname).readlines()
        lines.pop(0)
        self.conf_crg = []
        self.conf_name = []
        for line in lines:
            fields = line.split()
            self.conf_name.append(fields[1])
            self.conf_crg.append(float(fields[4]))
        return


class ResidueClass:
    def __init__(self):
        self.crg_states = []      # possible charge states
        self.netcharge = []       # net charge at various pH
        self.state_flag = ""      # * for multiple charges, + for positive and - for negative
        return


def diff(a, b):
    only_in_a = set(a) - set(b)
    only_in_b = set(b) - set(a)
    return only_in_a, only_in_b


if __name__ == "__main__":

    # read fort.38 and head3.lst
    fort38 = Fort38("fort.38")
    head3lst = Head3List("head3.lst")

    # verify two files have natching conformers
    only_in_fort38, only_in_head3lst = diff(fort38.conf_name, head3lst.conf_name)
    if only_in_fort38:
        print "These conformers only exist in fort.38: %s" % ",".join(only_in_fort38)
    if only_in_head3lst:
        print "These conformers only exist in head3.lst: %s" % ",".join(only_in_head3lst)
    if only_in_fort38 or only_in_head3lst:
        sys.exit()

    # apply charge from head3.lst to fort.38 occupancy table
    fort38.load_charge(head3lst)

    # group into residues
    residueDB = {}
    residue_names = []
    for i in range(len(fort38.conf_name)):
        # NTR01A0299_001 -> NTRA0299_
        resid = (fort38.conf_name[i][:3], fort38.conf_name[i][5:11])
        if resid not in residue_names:
            residue_names.append(resid)
            res = ResidueClass()
            res.netcharge = np.linspace(0.0, 0.0, len(fort38.titration_points))
            residueDB[resid] = res

        residueDB[resid].crg_states.append(head3lst.conf_crg[i])
        residueDB[resid].netcharge += fort38.crg_table[i]


    # charged residues, obtain a list charged residues, and their state flag
    charged_residue_names = []
    for resid in residue_names:
        has_negative = False
        has_positive = False
        for state in residueDB[resid].crg_states:
            if state < -0.999:
                has_negative = True
            if state > 0.999:
                has_positive = True

        if has_positive and has_negative:
            state_flag = "*"
            charged_residue_names.append(resid)
        elif has_negative:
            state_flag = "-"
            charged_residue_names.append(resid)
        elif has_positive:
            state_flag = "+"
            charged_residue_names.append(resid)
        else:
            state_flag =""

        residueDB[resid].state_flag = state_flag


    # write out sum_crg.out, put lines in a list, then write out in one command
    sumcrg = ["  %s          %s\n" % (fort38.titration_type, "".join([" %5.1f" % x for x in fort38.titration_points]))]
    for resid in charged_residue_names:
        charges = residueDB[resid].netcharge
        res_str = "%s%s%s" % (resid[0], residueDB[resid].state_flag, resid[1])
        crg_str = " ".join(["%5.2f" % x for x in residueDB[resid].netcharge])
        line = "%s     %s\n" % (res_str, crg_str)
        sumcrg.append(line)

    open("sum_crg2.out", "w").writelines(sumcrg)


    # Introduction on titration curve and sigmoid function
    # Curve fitting, how do we get the documentation

    # ===================================================================
    # Session 3 Curve fitting
    pkout = ["  %s             pKa/Em  n(slope) 1000*chi2\n" % fort38.titration_type]
    xdata = fort38.titration_points
    for resid in charged_residue_names:
        res_str = "%s%s%s" % (resid[0], residueDB[resid].state_flag, resid[1])
        # ydata of sigmoid function is always positive [0, 1] but our net charge can be
        ydata = residueDB[resid].netcharge
        #print xdata
        if ydata.mean() > 0:
            ydata += -1.0  # acid[0, -1], base[+1, 0] -> [0, -1]
        # print ydata.mean()

        # Initial guess of midpoint and slope
        ix = min(range(len(xdata)), key=lambda i: abs(ydata[i]+0.5))
        midpoint_guess = xdata[ix]
        n_guess = 1/0.4343   # default PH titration
        if fort38.titration_type.upper() == "EH":
            n_guess *= 58.0
        guess = (midpoint_guess, n_guess)

        msg = ""
        try:
            (popt, pcov) = curve_fit(sigmoid, xdata, ydata, guess)
        except RuntimeError:
            msg = "Titration out of range"
        except ValueError:
            msg = "Input value not valid"

        if popt[0] < fort38.titration_lb or popt[0] > fort38.titration_ub:
            msg = "Titration out of range"


        if msg:
            pkout.append("%s         %s\n" % (res_str, msg))
        else:
            chi_squared = np.sum([(sigmoid(xdata, *popt) - ydata)**2])
            if fort38.titration_type.upper() == "PH":
                nslope = 0.4342*popt[1]/fort38.titration_step
            elif fort38.titration_type.upper() == "EM":
                nslope = 0.4342*popt[1]/fort38.titration_step/58.0
            pka_str = "%9.3f %9.3f %9.3f" % (popt[0], 0.4342*popt[1], chi_squared*1000.0)
            pkout.append("%s    %s\n" % (res_str, pka_str))

        open("pK2.out", "w").writelines(pkout)
#!/usr/bin/env python2

import sys
import numpy as np

PH2KCAL = 1.364

# read fort.38
fort38_lines = open("fort.38").readlines()

# read head3.lst
head3_lines = open("head3.lst").readlines()

headline = fort38_lines.pop(0)
head3_lines.pop(0).strip()

# remove DM conformers
#fort38_lines = [x for x in fort38_lines if x[3:5] != "DM"]
#head3_lines = [x for x in head3_lines if x[9:11] != "DM"]

conformer = []  # 1D array
occ_table = []  # 2D array, here is the rows array, each column will be occ at one pH
crg_state = []  # 1D array, same size as conformer
while fort38_lines and head3_lines:
    fort38_data = fort38_lines.pop(0)
    fields = fort38_data.split()
    confname = fields[0]
    occ = fields[1:]

    head3_data = head3_lines.pop(0)
    fields = head3_data.split()

    if confname != fields[1]:
        print "Confomer %s in fort.38 does not match %s in head3.lst" % (confname, fields[1])
        sys.exit()

    crg = fields[4]

    conformer.append(confname)
    occ_table.append(occ)
    crg_state.append(crg)

# multiply occ_table with crg_state
nm_occ_table = np.array(occ_table, dtype=np.float)
nm_crg_state = np.array(crg_state, dtype=np.float)

nm_crg_table = nm_occ_table * nm_crg_state[:, np.newaxis]

# ===================================================================
# Session 2 Group the conformers into residues
# Reorganize the data, make it easier to reader and handle
# What we have so far:
# headline: the head line of fort.38
# conformer = a list of conformer names
# nm_crg_state = a np list of conformer charge states
# nm_crg_table = a 2D matrix of conformer net charge at titration pHs

# define a class for residue so that residue information can be stored at one place


class ResidueClass:
    def __init__(self):
        self.crg_states = []      # possible charge states
        self.netcharge = []       # net charge at various pH
        self.state_flag = ""      # * for multiple charges, + for positive and - for negative
        return


# group into residues
residues = {}
residue_names = []
for i in range(len(conformer)):
    # NTR01A0299_001 -> NTRA0299_
    resname = (conformer[i][:3], conformer[i][5:11])
    if resname not in residue_names:
        residue_names.append(resname)
        res = ResidueClass()
        res.netcharge = [0.0 for x in range(len(nm_crg_table[i]))]
        res.netcharge = np.linspace(0.0, 0.0, len(nm_crg_table[i]))
        residues[resname] = res

    residues[resname].crg_states.append(nm_crg_state[i])
    residues[resname].netcharge += nm_crg_table[i]


# charged residues, obtain a list charged residues, and their state flag
charged_residue_names = []
for resname in residue_names:
    has_negative = False
    has_positive = False
    for state in residues[resname].crg_states:
        if state < -0.999:
            has_negative = True
        if state > 0.999:
            has_positive = True

    if has_positive and has_negative:
        state_flag = "*"
        charged_residue_names.append(resname)
    elif has_negative:
        state_flag = "-"
        charged_residue_names.append(resname)
    elif has_positive:
        state_flag = "+"
        charged_residue_names.append(resname)
    else:
        state_flag =""

    residues[resname].state_flag = state_flag


# write out sum_crg.out, put lines in a list, then write out in one command
sumcrg = [headline]
for resname in charged_residue_names:
    charges = residues[resname].netcharge
    res_str = "%s%s%s" % (resname[0], residues[resname].state_flag, resname[1])
    crg_str = " ".join(["%5.2f" % x for x in residues[resname].netcharge])
    line = "%s     %s\n" % (res_str, crg_str)
    sumcrg.append(line)

open("sum_crg2.out", "w").writelines(sumcrg)


# Introduction on titration curve and sigmoid function
# Curve fitting, how do we get the documentation

# ===================================================================
# Session 3 Curve fitting

from scipy.optimize import curve_fit


def sigmoid(x, x0, k):
    e = np.exp(k*(x-x0))
    y = -e/(1.0+e)
    return y

fields = headline.split()
titration_type = fields[0].strip()
titration_delta = abs(float(fields[2]) - float(fields[1]))
titration_lb = float(fields[1])
titration_ub = float(fields[-1])


pkout = ["  %s             pKa/Em  n(slope) 1000*chi2\n" % titration_type]
xvalue = np.array([float(x) for x in headline[14:].split()])
x0 = xvalue[0]
k0 = xvalue[1] - xvalue[0]

xdata = np.array([float(i) for i in range(len(headline[14:].split()))])
for resname in charged_residue_names:
    res_str = "%s%s%s" % (resname[0], residues[resname].state_flag, resname[1])
    # ydata of sigmoid function is always positive [0, 1] but our net charge can be
    ydata = residues[resname].netcharge
    #print xdata
    #print "Before:",ydata
    if ydata.mean() > 0:
        ydata += -1.0  # acid[0, -1], base[+1, 0] -> [0, -1]
    if ydata[-1] - ydata[0] > 0.0:
        ydata =  -ydata - 1   
    #print "After:", ydata
    msg = ""
    try:
        (popt, pcov) = curve_fit(sigmoid, xdata, ydata, bounds=(0,[len(headline[14:].split()), 4]))
        #(popt, pcov) = curve_fit(sigmoid, xdata, ydata)
    except RuntimeError:
        msg = "Titration out of range"
    except ValueError:
        msg = "Input value not valid"

    if popt[0] < 0.001 or popt[0] > len(headline[14:].split())-0.001:
        msg = "Titration out of range"


    if msg:
        pkout.append("%s         %s\n" % (res_str, msg))
    else:
        chi_squared = np.sum([(sigmoid(xdata, *popt) - ydata)**2])
        midpoint = popt[0] * k0 + x0
        if titration_type.upper() == "PH":
            nslope = 0.4342*popt[1]/titration_delta
        elif titration_type.upper() == "EM":
            nslope = 0.4342*popt[1]/titration_delta*58.0
        elif titration_type.upper() == "CH" or titration_type.upper() == "EXTRA":
            nslope = 0.4342*popt[1]/titration_delta*PH2KCAL
        else:
            print("Why am I here?")
        pka_str = "%9.3f %9.3f %9.3f" % (midpoint, nslope, chi_squared*1000.0)
        pkout.append("%s    %s\n" % (res_str, pka_str))

    open("pK2.out", "w").writelines(pkout)

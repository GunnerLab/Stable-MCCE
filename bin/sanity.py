#!/usr/bin/env python
"""
This program runs sanity check on sum_crg.out.

Input:
 * sum_crg.out
 * run.prm

Output:
 * screen message

Usage examples:

1. Run sanity check with default
    sanity.py

2. Run sanity check with additional user defined rules
    sanity.py -f rulefile

Rule file format:

Residue_name: pH/Eh=(lower_bound, upper_bound)
ASP: ph=(0, 6)
HIS: ph=(6, 9)
LYS: ph=(9, 13)
ARG: ph=(10, 14)
HEM: eh=(100, 400)
"""

import argparse

rules = {"ASP-": (2, 6),
         "ACY-": (2, 6),
         "ARG+": (10, 15),
         "BCL-": (35, 45),
         "CTR-": (1, 5),
         "CYS-": (7, 11),
         "GLU-": (2, 6),
         "HIS+": (5, 11),
         "LYS+": (8, 13),
         "NTR+": (5, 9),
         "TYR-": (9, 12)
         }

def update_rules(fname):
    lines = open(fname).readlines()
    for line in lines:
        entry = line.split("#")[0]
        fields = entry.split(":")
        if len(fields) == 2:
            key = fields[0].strip()
            lower, upper = fields[1].split(",")
            rules[key] = (type, float(lower), float(upper))
    return

def sanity_check(args):
    # read user rules
    if args.f:
        update_rules(args.f)

    # go over sum_crg.out
    lines = open("sum_crg.out").readlines()
    headline = lines.pop(0)
    fields = headline.split()
    titration_type = fields[0].strip()
    titration_points = [float(x) for x in fields[1:]]
    # get initial pH from run.prm, it could be the pH for Eh titration
    ph0 = get_ph1()

    n_points = len(titration_points)
    if titration_type.upper() == "PH":
        phs = titration_points
    elif titration_type.upper() == "EH":
        phs = [ph0]

    for line in lines:
        fields = line.split()
        if len(fields) < n_points + 1:
            # reached the end
            break
        resid = fields[0]
        resname = resid[:4]
        if resname in rules:
            charges = [float(x) for x in fields[1:]]
            for i in range(len(titration_points)):
                if titration_type.upper() == "PH":
                    ph = phs[i]
                elif titration_type.upper() == "EH":
                    ph = phs[0]

                charge = charges[i]
                if abs(charge) < 0.001:
                    charge = 0.0

                # We have resname, pH, and charge. The question is if this charge under this pH is reasonable?
                msg = reasonable_charge(resname, ph, charge)
                if msg:
                    print("%s: charge %5.2f expected %s at pH=%4.1f" % (resid, charge, msg, ph))
                    break
    return

def reasonable_charge(resname, ph, charge):
    lower, upper = rules[resname]
    msg = ""
    # sane range
    # + (low ph)       --->               - (high ph)
    # Charge     <lower    in_range       >upper
    # Acid -     0 <-> -0.5  any           -0.5 <-> -1
    # Base +     1 <-> 0.5   any            0.5 <-> 0

    if ph <= lower:  # left
        if resname[-1] == "-":    # left side acid => can not be less than -0.5
            if charge < -0.5:
                msg = "near  0.0"   # "[-0.5,  0.0]"
        else:                     # left side base => can not be less than 0.5
            if charge < 0.5:
                msg = "near  1.0"   # "[ 0.5,  1.0]"
    elif ph >= upper:  # right
        if resname[-1] == "-":    # right side acid => can not be over -0.5
            if charge > -0.5:
                msg = "near -1.0"   # "[-1.0, -0.5]"
        else:                     # right side base => can not be over 0.5
            if charge > 0.5:
                msg = "near  0.0"   # "[ 0.0,  0.5]"

    return msg




def get_ph1():
    # get ph1 from run.prm
    lines = open("run.prm").readlines()
    ph = 7.0
    for line in lines:
        fields = line.strip().split()
        if "(TITR_PH0)" in fields:
            ph = float(fields[0])
    return ph

if __name__ == "__main__":
    # Get the command arguments
    helpmsg = "Run sum_crg.out sanity check. report unusual charge at titration point."
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

    sanity_check(args)

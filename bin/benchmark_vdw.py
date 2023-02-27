#!/usr/bin/env python

"""
Compare vdw values in vdw_pw.py created head3.lst and opp files with old mcce values
"""

import sys
import os

head3 = "head3.lst"
head3_old = "head3.lst_bak"
energies = "energies"
energies_old = "energies_bak"

if __name__ == "__main__":
    lines_a = open(head3).readlines()
    lines_b = open(head3_old).readlines()
    if len(lines_a) != len(lines_b):
        print("The number of lines in %s and %s do not match" % (head3, head3_old))
        sys.exit()

    vdw0_lines = ["Confname, old_vdw0, new_vdw0, difference\n"]
    vdw1_lines = ["Confname, old_vdw1, new_vdw, difference\n"]
    opp_lines = ["Confname1, Confname2, old_vdw, new_vdw, difference\n"]
    for iline in range(1, len(lines_a)):
        line = lines_a[iline]
        fields = line.split()
        confname_a = fields[1]
        vdw0_a = fields[9]
        vdw1_a = fields[10]

        line = lines_b[iline]
        fields = line.split()
        confname_b = fields[1]
        vdw0_b = fields[9]
        vdw1_b = fields[10]

        if confname_a == confname_b:
            vdw0_lines.append("%s, %7s, %7s, %10.3f\n" % (confname_a, vdw0_b, vdw0_a, float(vdw0_a)-float(vdw0_b)))
            vdw1_lines.append("%s, %7s, %7s, %10.3f\n" % (confname_a, vdw1_b, vdw1_a, float(vdw1_a)-float(vdw1_b)))
        else:
            print("Mismatch at line %d: %s <-> %s" % (iline, confname_b, confname_a))
            sys.exit()

        # prepare vdw dictionary from old opp files
        fname = "%s/%s.opp" % (energies_old, confname_a)
        vdw_old = {}
        if os.path.exists(fname):
            lines = open(fname).readlines()
            for line in lines:
                fields = line.split()
                confname_b = fields[1]
                vdw_value = fields[3]
                vdw_old[confname_b] = vdw_value

        # put old and new values together
        fname = "%s/%s.opp" % (energies, confname_a)
        if os.path.exists(fname):
            lines = open(fname).readlines()
            for line in lines:
                fields = line.split()
                confname_b = fields[1]
                vdw_value = fields[3]
                if confname_b in vdw_old:
                    vdw_value_old = vdw_old[confname_b]
                else:
                    print("skipping %s %s" % (confname_a, confname_b))
                    continue
                opp_lines.append("%s, %s, %7s, %7s, %10.3f\n" % (confname_a, confname_b, vdw_value_old, vdw_value, float(vdw_value)-float(vdw_value_old)))



    open("vdw0.csv", "w").writelines(vdw0_lines)
    open("vdw1.csv", "w").writelines(vdw1_lines)
    open("vdw.csv", "w").writelines(opp_lines)


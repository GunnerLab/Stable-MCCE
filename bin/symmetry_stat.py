#!/usr/bin/env python
"""Check step3 opp file symmetry"""
from os import path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import style


CUTOFF = 0.01

class Pairwise:
    def __init__(self):
        head3lst = "head3.lst"
        energies_folder = "energies"
        lines = open(head3lst).readlines()
        lines.pop(0)
        self.confnames = [ x.split()[1] for x in lines if len(x) > 80]
        self.nconf = len(self.confnames)
        self.ele = [[0.0 for i in range(self.nconf)] for j in range(self.nconf)]
        self.vdw = [[0.0 for i in range(self.nconf)] for j in range(self.nconf)]

        for (i, conf) in enumerate(self.confnames):
            fname = energies_folder+"/"+conf+".opp"
            if path.exists(fname):
                lines = open(fname).readlines()
                for line in lines:
                    if len(line) > 40:
                        fields = line.split()
                        conf2 = fields[1]
                        j = self.confnames.index(conf2)
                        self.ele[i][j] = float(fields[2])
                        self.vdw[i][j] = float(fields[3])
        return


def show_diff(pw):
    for i in range(pw.nconf-1):
        for j in range(i+1, pw.nconf):
            if abs(pw.ele[i][j]-pw.ele[j][i]) > CUTOFF:
                print("ELE %s <-> %s: %.3f <-> %.3f" % (pw.confnames[i], pw.confnames[j], pw.ele[i][j],
                                                            pw.ele[j][i]))
            if abs(pw.vdw[i][j]-pw.vdw[j][i]) > CUTOFF:
                print("VDW %s <-> %s: %.3f <-> %.3f" % (pw.confnames[i], pw.confnames[j], pw.vdw[i][j],
                                                           pw.vdw[j][i]))
    return


if __name__ == "__main__":
    pw = Pairwise()
    show_diff(pw)

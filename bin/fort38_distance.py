#!/usr/bin/env python

import argparse
import numpy as np

class Fort38:
    def __init__(self, fname):
        self.type = ""
        self.trange = []
        self.confnames = []
        self.confocc = []
        self.resnames = []
        self.residues = []
        self.load_fort38(fname)

    def load_fort38(self, fname):
        lines = open(fname).readlines()
        line = lines.pop(0)
        fields = line.split()
        self.type = fields[0]
        self.trange = [float(x) for x in fields[1:]]
        for line in lines:
            fields = line.split()
            self.confnames.append(fields[0])
            self.confocc.append([float(x) for x in fields[1:]])

        # group conformers into residues
        old_confname = self.confnames[0]
        old_resname = old_confname[:3]+old_confname[5:11]
        res = [0]
        for iconf in range(1, len(self.confnames)):
            resname = self.confnames[iconf][:3]+self.confnames[iconf][5:11]
            if old_resname == resname:
                res.append(iconf)
            else:
                self.residues.append(res)
                old_resname = resname
                res = [iconf]
        # add the last one
        self.residues.append(res)

        for res in self.residues:
            resname = self.confnames[res[0]][:3]+self.confnames[res[0]][5:11]
            self.resnames.append(resname)

        return


def bhata_distance(prob1, prob2):
    d_max = 10.0   # Max possible value set to this
    p1 = np.array((prob1)) / sum(prob1)
    p2 = np.array((prob2)) / sum(prob2)
    if len(p1) != len(p2):
        d = d_max
    else:
        bc = sum(np.sqrt(p1 * p2))
    #    print(bc, np.exp(-d_max))
        if bc <= np.exp(-d_max):
            d = d_max
        else:
            d = -np.log(bc)

    if d <= 0:
        d = 0
    return d



def compare_fort38(f1, f2):
    fort38_1 = Fort38(f1)
    fort38_2 = Fort38(f2)

    btd = []
    if fort38_1.resnames == fort38_2.resnames:
        for ires in range(len(fort38_1.residues)):
            res = fort38_1.residues[ires]
            p1_all = [fort38_1.confocc[iconf] for iconf in res]
            p2_all = [fort38_2.confocc[iconf] for iconf in res]
            btd_this_residue = []
            for ipoint in range(len(p1_all[0])):
                d = bhata_distance([x1[ipoint] for x1 in p1_all], [x2[ipoint] for x2 in p2_all])
                #print(fort38_1.resnames[ires], d)
                btd_this_residue.append(d)
            btd.append(btd_this_residue)
    else:
        print("Two fort.38 files don't have same residues")

    # print the distance table
    print(" %-8s %s" % (fort38_1.type, " ".join(["%5.1f" % x for x in fort38_1.trange])))
    for ires in range(len(fort38_1.residues)):
        print("%s %s" % (fort38_1.resnames[ires], " ".join(["%5.0f" % (x*1000) for x in btd[ires]])))

#    print(btd)

    return

if __name__ == "__main__":
    # Get the command arguments
    helpmsg = "Compare two fort.38 files and calculate the Bhattachayya distance at residue level."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("files", metavar="fort.38_file1 fort.38_file2", nargs=2)
    args = parser.parse_args()

    compare_fort38(args.files[0], args.files[1])

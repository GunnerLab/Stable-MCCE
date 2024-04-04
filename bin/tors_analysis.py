#!/usr/bin/env python

# comparison of VDW0, VDW1, TORS, DSOLV in head3.lst
old_folder = "../4lzt-old/"
new_folder = "./"

class SELF_E:
    def __init__(self):
        self.vdw0 = 0.0
        self.vdw1 = 0.0
        self.tors = 0.0
        self.epol = 0.0
        self.dsolv = 0.0

    def loadline(self, line):
        fields = line.split()
        confID = ""
        if len(fields) == 17:
            self.vdw0 = float(fields[9])
            self.vdw1 = float(fields[10])
            self.tors = float(fields[11])
            self.epol = float(fields[12])
            self.dsolv = float(fields[13])
            confID = fields[1]
        return confID

def load_head3(folder):
    fname = "%s/head3.lst" % folder

    lines = open(fname).readlines()

    head3 = {}
    for line in lines:
        self_e = SELF_E()
        confID = self_e.loadline(line)

        if confID:
            head3[confID] = self_e

    return head3


if __name__ == "__main__":
    new_head3 = load_head3("./")
    old_head3 = load_head3("../4lzt-old")

    confIDs = old_head3.keys()

    lines = ["Conformer, vdw0_old, vdw0_new, vdw1_old, vdw1_new\n"]
    for confID in confIDs:
        line = "%s, %.3f,  %.3f, %.3f, %.3f\n" % (confID, old_head3[confID].vdw0, new_head3[confID].vdw0, old_head3[confID].vdw1, new_head3[confID].vdw1)
        lines.append(line)
    open("vdw01.csv", "w").writelines(lines)

    lines = ["Conformer, tors_old, tors_new\n"]
    for confID in confIDs:
        line = "%s, %.3f,  %.3f\n" % (confID, old_head3[confID].tors, new_head3[confID].tors)
        lines.append(line)
    open("tors.csv", "w").writelines(lines)

    lines = ["Conformer, dsolv_old, dsolv_new\n"]
    for confID in confIDs:
        line = "%s, %.3f,  %.3f\n" % (confID, old_head3[confID].dsolv, new_head3[confID].dsolv)
        lines.append(line)
    open("dsolv.csv", "w").writelines(lines)

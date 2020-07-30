#!/usr/bin/python

import math
from os import listdir
from os.path import isfile, join

from pdbio import ddvv

VDW_CUTOFF_NEAR = 1
VDW_ELIMIT_NEAR = 999
VDW_CUTOFF_FAR = 10
VDW_ELIMIT_FAR = 0

cutoff_near2 = VDW_CUTOFF_NEAR * VDW_CUTOFF_NEAR
cutoff_far2 = VDW_CUTOFF_FAR  * VDW_CUTOFF_FAR

print_cutoff = 0.000001

class VDWPARM:
    def __init__(self):
        self.vdwparam = {}
        self.scaling = {}
        self.factor_14lj = 0.5 # will be modified in initdb()
        return

    def load(self, files):
        for file in files:
            lines = open(file).readlines()
            for line in lines:
                fields = line.split()
                if len(fields) > 0:
                    if fields[0] == "VDWAMBER":
                        Ctype = fields[1]
                        pair = fields[2]
                        if (Ctype, pair) in self.vdwparam:
                            print("WARNING: %s of %s entry duplicates detected." % (Ctype, pair))
                        self.vdwparam[(Ctype, pair)] = float(fields[3])
                    elif fields[0] == "SCALING":
                        Vtype = fields[1]
                        self.scaling[Vtype] = float(fields[2])
        return

vdwdb = VDWPARM()

def initdb():
    mcce_home = ""
    extra = ""
    factor_14lj = 0.5
    lines = open("run.prm").readlines()
    for line in lines:
        fields = line.split()
        if len(fields) >= 2:
            key = fields[-1].strip()
            if key[0] == "(" and key[-1] == ")":
                key = key.strip("()").strip()
                if key == "MCCE_HOME":
                    mcce_home = fields[0]
                elif key == "EXTRA":
                    extra = fields[0]
                elif key == "EPSILON_PROT":
                    epsilon_prot = fields[0]
                elif key == "FACTOR_14LJ":
                    factor_14lj = float(fields[0])

    # get a list of all parameter files
    if epsilon_prot == "4.0":
        parmdir = mcce_home+"/param04"
    elif epsilon_prot == "8.0":
        parmdir = mcce_home+"/param08"

    tplfiles = [join(parmdir, f) for f in listdir(parmdir) if f.endswith(".tpl")]
    tplfiles.append(extra)

    # pass vdw parameters to db
    vdwdb.load(tplfiles)
    vdwdb.factor_14lj = factor_14lj
    return


def vdw(atom1, atom2):
    v1 = atom1.xyz
    v2 = atom2.xyz
    d2 = ddvv(v1, v2)

    if d2 > cutoff_far2:
        e = VDW_ELIMIT_FAR
    elif d2 < cutoff_near2:
        #print v1, v2
        e = VDW_ELIMIT_NEAR
    else:
        pair = "%c-%c" % (atom1.name[1], atom2.name[1])
        if ("C12", pair) in vdwdb.vdwparam:
            C12 = vdwdb.vdwparam[("C12", pair)]
        else:
            C12 = vdwdb.vdwparam[("C12", "X-X")]
        if ("C6", pair) in vdwdb.vdwparam:
            C6 = vdwdb.vdwparam[("C6", pair)]
        else:
            C6 = vdwdb.vdwparam[("C6", "X-X")]

        d6 = d2*d2*d2
        d12 = d6*d6
        e = C12/d12 - C6/d6

    if abs(e) > print_cutoff:
        d = math.sqrt(d2)
        if atom1.confID == atom2.confID:
            mark = "*"
        else:
            mark = ""
        print("Atom %s <-> %s: %.3f  (d = %.3f) %s" % (atom1.atomID, atom2.atomID, e, d, mark))

    return e


def vdw_conf(conf1, conf2):
    e = 0.0
    for atom1 in conf1.atom:
        for atom2 in conf2.atom:
            if atom1 == atom2 or atom2 in atom1.connect12 or atom2 in atom1.connect13:
                continue
            elif atom2 in atom1.connect14:
                e += vdw(atom1, atom2) * vdwdb.factor_14lj
            else:
                e += vdw(atom1, atom2)

    if conf1 == conf2:
        e = 0.5 * e

    return e


if __name__ == "__main__":
    initdb()
    print vdwdb.vdwparam
    #print vdwdb.scaling

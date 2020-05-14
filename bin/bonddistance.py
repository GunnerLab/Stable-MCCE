#!/usr/bin/env python

import sys
import math

radius = {" H": 0.25,
          " N": 0.65,
          " C": 0.70,
          " O": 0.60,
          " P": 1.00,
          "NA": 1.80,
          "CL": 1.00
}
class Atom:
    def __init__(self):
        self.name = ""
        self.element = ""
        self.xyz = ()


def dvv(v1, v2):
    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    d2 = dx*dx + dy*dy + dz*dz
    return math.sqrt(d2)

fname = sys.argv[1]

lines = open(fname).readlines()

atoms = []
connect = {}

for line in lines:
    if len(line) < 11:
        continue
    if line[:6] == "ATOM  " or line[:6] == "HETATM":
        atom =Atom()
        atom.name = line[12:16]
        atom.element = line[76:78]
        atom.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        atoms.append(atom)
    elif line[:6] == "CONECT":
        fields = line.split()
        key = int(fields[1]) - 1
        value = [int(x) - 1 for x in fields[2:]]
        connect[key] = value

n = len(atoms)
for i in range(n-1):
    for j in range(i+1, n):
        d = dvv(atoms[i].xyz, atoms[j].xyz)
        if j in connect[i]:
            connected = "True"
        else:
            connected = "False"
        bond_distance = radius[atoms[i].element] + radius[atoms[j].element]
        if d < bond_distance * 1.6:
            predicted = "True"
        else:
            predicted = "False"
        print("%s - %s  %8.3f %5s %5s" % (atoms[i].name, atoms[j].name, d, predicted, connected))

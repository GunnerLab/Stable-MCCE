#!/usr/bin/env python
import sys

HATOMS = ["HG", "HD", "HE", "HH"]

lines = open(sys.argv[1]).readlines()
for line in lines:
    if line[:6] == "ATOM  " or line[:6] == "HETATM":
        if line[17:20] == "WAT":
            continue
        if line[13] == "H":
            continue
        if line[12:14] in HATOMS:
            continue
        print(line.strip("\n"))

#!/usr/bin/env python

import re
from energy import *

if __name__ == "__main__":
    ph = 0.0
    eh = 0.0

    if len(sys.argv) > 1:
        cutoff = float(sys.argv[1])
    else:
        cutoff = PW_PRINT_CUT

    lines = open("states").readlines()
    line = lines.pop(0)
    fields = line.strip().split(",")
    mc_parm = {}
    for field in fields:
        key, value = field.split("=")
        mc_parm[key.strip()] = float(value)

    T = mc_parm["T"]
    ph = mc_parm["ph"]
    eh = mc_parm["eh"]



    for line in lines:
        state = [int(x) for x in re.findall(r"[\w']+", line)]
        analyze_state_energy(state, ph=ph, eh=eh, cutoff=cutoff)

#    for conf in head3lst:
#        conf.printme()

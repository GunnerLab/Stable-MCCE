#!/usr/bin/env python

import argparse

fname_fort38 = "fort.38"
fname_head3lst = "head3.lst"
fname_pkout = "pK2.out"

ph2Kcal = 1.364
mev2Kcal = 0.0235
Kcal2kT = 1.688

# Normally a residue is divided into two groups based on total charge. If more than two charge states are detected or
# a mfe analysis is needed for groups with the same charge, you can specify the grouping rule here.
#              RES  : ground  excited
Special_res = {"_CU": (["+1"], ["+2"]),
               "UbQ": (["01"], ["-1"])}


class Resdiue:
    def __init__(self):
        return


def load_runprm(env):
    lines = open("run.prm").readlines()
    for line in lines:
        fields = line.split("#", 1)[0].strip().split()
        if len(fields) > 1 and fields[-1].find("(") >= 0 and fields[-1].find(")") > fields[-1].find("("):
            key = fields[-1].strip("()").strip()
            env[key] = fields[0]
    return


def load_scaling(env):
    lines = open(env["EXTRA"]).readlines()
    for line in lines:
        fields = line.split()
        if len(fields) > 2 and fields[0] == "SCALING":
            if fields[1] == "ELE":
                env["scale_ele"] = float(fields[2])
            elif fields[1] == "VDW":
                env["scale_vdw"] = float(fields[2])
            elif fields[1] == "VDW0":
                env["scale_vdw0"] = float(fields[2])
            elif fields[1] == "VDW1":
                env["scale_vdw1"] = float(fields[2])
            elif fields[1] == "TORS":
                env["scale_tors"] = float(fields[2])
            elif fields[1] == "DSOLV":
                env["scale_dsolv"] = float(fields[2])
    return


def mfe(args):
    # prepare mfe environment
    env = {"scale_ele": 1.0,
           "scale_vdw": 1.0,
           "scale_vdw0": 1.0,
           "scale_vdw1": 1.0,
           "scale_tors": 1.0,
           "scale_dsolv": 1.0,
           "mfe_res_name": args.residue[0],
           "mfe_pwcutoff": args.c}

    load_runprm(env)
    if env["EXTRA"]:
        load_scaling(env)

    if args.x == "t":
        env["mfe_xts"] = True
    elif args.x == "f":
        env["mfe_xts"] = False
    else:   # determine by run.prm
        if "MONTE_TSX" in env and env["MONTE_TSX"] == "t":
            env["mfe_xts"] = True
        else:
            env["mfe_xts"] = False

    # mfe_ph
    if args.p == "m":  # determined by pKa
        

    return


if __name__ == "__main__":
    # Get the command arguments
    helpmsg = "Calculate mean field energy ionization energy on ionazable residue at a specific pH/eH."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-p", metavar="pH/Eh", default="m",
                        help="pH or Eh mfe analysis is carried out, or \'m\' for midpoint")
    parser.add_argument("-x", metavar="TS_correction", default="r",
                        help="f: False, t: True, or r: determined by run.prm (default)")
    parser.add_argument("-c", metavar="cutoff", default="-0.01", help="pairwise cutoff in reporting", type=float)
    parser.add_argument("residue", metavar="residue", help="the residue name as in pK.out", nargs=1)

    args = parser.parse_args()

    mfe(args)

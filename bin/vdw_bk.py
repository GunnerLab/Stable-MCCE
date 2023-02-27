#!/usr/bin/env python

import argparse
from pdbio import *


def vdw1_breakdown(protein, confID, cutoff=-0.001):
    conf1 = None
    Found = False
    for res in protein.residue:
        if Found:
            break
        for conf in res.conf[1:]:
            if Found:
                break
            if conf.confID == confID:
                Found = True
                conf1 = conf
    if not Found:
        return

    vdw1 = 0.0
    for res in protein.residue:
        conf2 = res.conf[0]
        vdw = vdw_conf(conf1, conf2)
        vdw1 += vdw
        if abs(vdw) > cutoff:
            print("%s %7.3f" % (conf2.confID, vdw))
    print("Total          %7.3f" % vdw1)

    return

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Compute vdw1 breakdown at conformer level."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-c", metavar="cutoff", default="-0.001", help="Cutoff value of displaying conf vdw pairwise")
    parser.add_argument("conf", metavar="confID", nargs=1, help="Conformer ID as in head3.lst")
    args = parser.parse_args()

    conformer = args.conf[0]
    cutoff = float(args.c)

    pdbfile = "step2_out.pdb"
    protein = Protein()
    protein.loadpdb(pdbfile)
    protein.make_connect12()
    protein.make_connect13()
    protein.make_connect14()

    vdw1_breakdown(protein, conformer, cutoff)



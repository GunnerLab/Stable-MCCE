#!/usr/bin/env python
"""
Compute vdw pairwise and write back to opp files
"""

import argparse
from pdbio import *


def vdw_by_conf_pair(protein, confID1, confID2, cutoff, verbose):
    # detailed vdw between conf and conf
    found = False
    for res1 in protein.residue:
        if found:
            break
        for conf1 in res1.conf:
            if found:
                break
            if conf1.confID != confID1:
                continue
            for res2 in protein.residue:
                if found:
                    break
                for conf2 in res2.conf:
                    if conf2.confID != confID2:
                        continue
                    vdw = vdw_conf(conf1, conf2, verbose=True, cutoff=cutoff, display=verbose)
                    print("%s - %s: %.3f" % (conf1.confID, conf2.confID, vdw))

                    if confID1 == confID2:
                        print("! Self interaction, VDW between conformers is half of the atom vdw sum.")
                    found = True
                    break
    return

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Compute detailed conformer to conformer vdw."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-c", metavar="cutoff", default="-0.001", help="Cutoff value of displaying atom to atom vdw")
    parser.add_argument("-v", default=False, help="Turn on verbose mode, displaying more details", action="store_true")
    parser.add_argument("confs", metavar="confID", nargs=2, help="Conformer ID as in head3.lst, two IDs required.")
    args = parser.parse_args()

    conformers = [x.strip() for x in args.confs]
    cutoff = float(args.c)

    pdbfile = "step2_out.pdb"
    env.load_runprm()
    env.load_ftpl()
    protein = Protein()
    protein.loadpdb(pdbfile)
    protein.make_connect12()
    protein.make_connect13()
    protein.make_connect14()

    vdw_by_conf_pair(protein, conformers[0], conformers[1], cutoff, args.v)



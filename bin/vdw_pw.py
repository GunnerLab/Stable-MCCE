#!/usr/bin/env python
"""
Compute vdw pairwise and write back to opp files
"""

import os
from pdbio import *
import time
import argparse


class ConfSelf:
    def __init__(self, line):
        fields = line.split()
        self.iconf = int(fields[0])
        self.confID = fields[1]
        self.fl = fields[2]
        self.occ = float(fields[3])
        self.crg = float(fields[4])
        self.em0 = int(fields[5])
        self.pka0 = float(fields[6])
        self.ne = int(fields[7])
        self.nh = int(fields[8])
        self.vdw0 = float(fields[9])
        self.vdw1 = float(fields[10])
        self.tors = float(fields[11])
        self.epol = float(fields[12])
        self.dsolv = float(fields[13])
        self.extra = float(fields[14])
        self.history = fields[15]
        self.state = fields[16]


def update_opp(protein, verbose=False):
    # read head3.lst
    lines = open("head3.lst").readlines()
    conflist = []
    iconf_dict = {}
    for line in lines[1:]:
        conf_self = ConfSelf(line)
        conflist.append(conf_self)
        iconf_dict[conf_self.confID] = conf_self.iconf

    # get the iconf to print in opp
    for res in protein.residue:
        for conf in res.conf[1:]:
            conf.iconf = iconf_dict[conf.confID]

    for res in protein.residue:
        resID1 = res.resID
        for conf in res.conf[1:]:
            write_opp = False
            fname = "energies/" + conf.confID + ".opp"
            opp_pw = [[i, "", 0, 0, 0, 0, ""] for i in range(len(conflist)+1)]    # store retrieved pw in an array indexed by iconf instead of dict
            if os.path.isfile(fname):
                opplines = open(fname).readlines()
                for oppline in opplines:
                    oppline = oppline.strip()
                    if len(oppline) < 54:
                        oppline = oppline + "   "
                    fields = [int(oppline[:5]), oppline[6:20], float(oppline[20:29]), float(oppline[29:37]), float(oppline[37:45]), float(oppline[45:53]), oppline[53:]]
                    iconf = fields[0]
                    opp_pw[iconf] = fields
                write_opp = True
            else:  # check the biggest vdw and decide if a new opp file is needed
                # positive_interaction = np.max(protein.vdw_pw[conf.i])
                # negative_interaction = np.min(protein.vdw_pw[conf.i])
                row = protein.vdw_pw[conf.i].data
                positive_interaction = max(row[0])
                negative_interaction = min(row[0])
                if abs(positive_interaction) > 0.01 or abs(negative_interaction) > 0.01:
                    write_opp = True

            if write_opp:   # write new opp files
                if verbose:
                    print("   opp - %s ..." % conf.confID)

                #print("writing %s" % fname)
                new_opplines = []
                for res2 in protein.residue:
                    for conf2 in res2.conf[1:]:
                        resID2 = "%3s%4s%c" % (conf2.confID[:3], conf2.confID[6:10], conf2.confID[5])
                        if resID1 == resID2:
                            continue

                        # pw_key = (conf.confID, conf2.confID)
                        # if pw_key in protein.vdw_pw:
                        #     new_pw = protein.vdw_pw[pw_key]
                        # else:
                        #     new_pw = 0.0
                        new_pw = protein.vdw_pw[conf.i, conf2.i]

                        newline = ""
                        if conf2.confID in opp_pw:
                            if abs(float(opp_pw[conf2.confID][2])) > 0.001 or abs(float(opp_pw[conf2.confID][3])) > 0.001:
                                newline = "%05d %s %8.3f%8.3f%8.3f%8.3f %s\n" % (conf2.iconf,
                                                                               conf2.confID,
                                                                               float(opp_pw[conf2.confID][2]),
                                                                               new_pw,
                                                                               float(opp_pw[conf2.confID][4]),
                                                                               float(opp_pw[conf2.confID][5]),
                                                                               opp_pw[conf2.confID][6])
                            elif abs(new_pw) > 0.001:
                                newline = "%05d %s %8.3f%8.3f%8.3f%8.3f +\n" % (conf2.iconf,
                                                                               conf2.confID,
                                                                               0.0,
                                                                               new_pw,
                                                                               0.0,
                                                                               0.0)

                        if newline:
                            new_opplines.append(newline)

                    open(fname, "w").writelines(new_opplines)


    # update head3.lst
    new_lines = ["iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    history\n"]

    # match conformer by confID
    conf_dict = {}
    for res in protein.residue:
        for conf in res.conf:
            conf_dict[conf.confID] = conf

    for conf in conflist:
        if conf.confID in conf_dict:    # DM conformer vdw1 not in dict
            vdw0 = conf_dict[conf.confID].vdw0
        else:
            vdw0 = conf.vdw0

        if conf.confID in conf_dict:    # DM conformer vdw1 not in dict
            vdw1 = conf_dict[conf.confID].vdw1
        else:
            vdw1 = conf.vdw1
        newline = "%05d %s %s %4.2f %6.3f %5d %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10s %c\n" % (conf.iconf,
                                          conf.confID,
                                          conf.fl,
                                          conf.occ,
                                          conf.crg,
                                          conf.em0,
                                          conf.pka0,
                                          conf.ne,
                                          conf.nh,
                                          vdw0,
                                          vdw1,
                                          conf.tors,
                                          conf.epol,
                                          conf.dsolv,
                                          conf.extra,
                                          conf.history,
                                          conf.state)

        new_lines.append(newline)

    open("head3.lst", "w").writelines(new_lines)

    return


if __name__ == "__main__":
    # Get the command arguments
    helpmsg = "Compute detailed conformer to conformer vdw."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-v", default=False, help="Turn on verbose mode to show progress", action="store_true")
    args = parser.parse_args()

    start_time = time.time()

    print("Loading paramter and structure ...")
    current_time = time.time()
    pdbfile = "step2_out.pdb"
    env.load_runprm()
    env.load_ftpl()
    protein = Protein()
    protein.loadpdb(pdbfile)
    elapsed = time.time() - current_time
    print("Done, elapsed time %.3f seconds." % elapsed)


    print("Making atom connectivity ...")
    current_time = time.time()
    protein.make_connect12()
    protein.make_connect13()
    protein.make_connect14()
    elapsed = time.time() - current_time
    print("Done, elapsed time %.3f seconds." % elapsed)


    print("Calculating vdw ...")
    current_time = time.time()
    protein.calc_vdw(verbose=args.v)
    elapsed = time.time() - current_time
    print("Done, elapsed time %.3f seconds." % elapsed)

    print("Write out head3.lst and opps ...")
    current_time = time.time()
    update_opp(protein, verbose=args.v)
    elapsed = time.time() - current_time
    print("Done, elapsed time %.3f seconds." % elapsed)
    elapsed = time.time() - start_time
    print("Total elapsed time %.3f seconds." % elapsed)

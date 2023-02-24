#!/usr/bin/env python
"""
Compute vdw pairwise and write back to opp files
"""

import os
from pdbio import *

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

def update_opp(protein):
    # read head3.lst
    lines = open("head3.lst").readlines()
    conflist = []
    for line in lines[1:]:
        conf_self = ConfSelf(line)
        conflist.append(conf_self)

    for res in protein.residue:
        resID1 = res.resID
        for conf in res.conf[1:]:
            write_opp = False
            fname = "energies/" + conf.confID + ".opp"
            opp_pw = {}
            if os.path.isfile(fname):
                opplines = open(fname).readlines()
                for oppline in opplines:
                    fields = oppline.split()
                    if len(fields) < 7:
                        fields.append(" ")
                    id = fields[1]
                    opp_pw[id] = fields
                write_opp = True
            else:  # check the biggest vdw and decide if a new opp file is needed
                for key in protein.vdw_pw.keys():
                    if conf.confID in key:
                        write_opp = True
                        break
            if write_opp:   # write new opp files
                #print("writing %s" % fname)
                new_opplines = []
                for conf2 in conflist:
                    resID2 = "%3s%4s%c" % (conf2.confID[:3], conf2.confID[6:10], conf2.confID[5])
                    if resID1 == resID2:
                        continue
                    pw_key = (conf.confID, conf2.confID)
                    if pw_key in protein.vdw_pw:
                        new_pw = protein.vdw_pw[pw_key]
                    else:
                        new_pw = 0.0
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
    for conf in conflist:
        key = (conf.confID, conf.confID)
        if key in protein.vdw_pw:
            vdw0 = protein.vdw_pw[key]
        else:
            vdw0 = 0.0

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
                                          conf.vdw1,
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
    pdbfile = "step2_out.pdb"
    protein = Protein()
    protein.loadpdb(pdbfile)
    protein.make_connect12()
    protein.make_connect13()
    protein.make_connect14()
    protein.calc_vdw()
    update_opp(protein)

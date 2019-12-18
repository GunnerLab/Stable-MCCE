#!/usr/bin/python
"""
 Calculate vdw energy table.
 This is a separate calculation at atom level. The result is not used by MCCE. However the conformer level vdw terms
 should be consistent with those in opp files.
"""

import sys

import libvdw
from pdbio import *

if __name__ == "__main__":
    # find location of parameters
    libvdw.initdb()
    #print vdwdb.vdwparam
    #print vdwdb.scaling

    pdbfile = "step2_out.pdb"
    protein = Protein()
    protein.loadpdb(pdbfile)
    protein.make_connect12()
    protein.make_connect13()
    protein.make_connect14()

    if len(sys.argv) < 2:
        print "vdw.py conformerID cutoff"
        print "Example: vdw.py GLU-1A0060_005 0.00"
        print "Note: conformerID is in head3.lst."
        sys.exit()

    resName = sys.argv[1][:3]
    confSeq = sys.argv[1][5:]
    if len(sys.argv) > 2:
        libvdw.print_cutoff = float(sys.argv[2])

    chainID = confSeq[0]
    resSeq = int(confSeq[1:5])
    confNum = int(confSeq[-3:])
    resID = "%3s%04d%c" % (resName, resSeq, chainID)
    confID = "%3s%04d%c%03d" % (resName, resSeq, chainID, confNum)

    found_res = False
    found_conf = False
    for res in protein.residue:
        if found_res:
            break
        if res.resID == resID:
            found_res = True
            for conf in res.conf:
                if found_conf:
                    break
                if confID == conf.confID:
                    found_conf = True
                    vdw0 = 0.0  # intra
                    vdw1 = 0.0  # to backbone
                    for res2 in protein.residue:
                        vdwt = libvdw.vdw_conf(conf, res2.conf[0])
                        vdw1 += vdwt
                        if abs(vdwt) > libvdw.print_cutoff:
                            print "Backbone(Accumulative): %s -> %s: %.3f" % (conf.confID, res2.conf[0].confID, vdw1)
                        for conf2 in res2.conf[1:]:
                            if conf2 == conf:  # Intra
                                vdw0 = libvdw.vdw_conf(conf, conf2)
                                print "Intra: %s -> %s: %.3f *" % (conf.confID, conf2.confID, vdw0)
                            elif res == res2:  # same residue to other conformers
                                vdwt = 0.0
                            else:
                                vdwt = libvdw.vdw_conf(conf, conf2)
                                if abs(vdwt) > libvdw.print_cutoff:
                                    print "Pairwise: %s -> %s: %.3f" % (conf.confID, conf2.confID, vdwt)

    print "%s: vdw0=%.3f, vdw1=%.3f" % (sys.argv[1], vdw0, vdw1)

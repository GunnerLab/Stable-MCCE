#!/usr/bin/env python

import argparse
import numpy as np

class Fort38:
    def __init__(self, fname):
        self.type = ""
        self.trange = []
        self.confnames = []
        self.confocc = {}   # confname and pH as the key, occupancy in fort.38 as value
        self.load_fort38(fname)

    def load_fort38(self, fname):
        lines = open(fname).readlines()
        line = lines.pop(0)
        fields = line.split()
        self.type = fields[0]
        self.trange = ["%.1f" % float(x) for x in fields[1:]]
        for line in lines:
            fields = line.split()
            confname = fields[0]
            for tpoint in range(len(fields[1:])):
                pH = self.trange[tpoint]
                self.confocc[(confname, pH)] = float(fields[tpoint+1])
            self.confnames.append(confname)

        return


def bhata_distance(prob1, prob2):
    d_max = 10.0   # Max possible value set to this
    p1 = np.array((prob1)) / sum(prob1)
    p2 = np.array((prob2)) / sum(prob2)
    if len(p1) != len(p2):
        d = d_max
    else:
        bc = sum(np.sqrt(p1 * p2))
    #    print(bc, np.exp(-d_max))
        if bc <= np.exp(-d_max):
            d = d_max
        else:
            d = -np.log(bc)

    if d <= 0:
        d = 0
    return d


def merge_lists(list1, list2):
    """
    Return a merged list while preserving order.
    """

    ipos = 0
    list_merged = []
    for x in list1:
        if x not in list_merged:
            if x in list2:
                xpos = list2.index(x)
                list_merged += list2[ipos:xpos]
                ipos = xpos + 1
            list_merged.append(x)

    return list_merged


def compare_fort38(f1, f2):
    fort38_1 = Fort38(f1)
    fort38_2 = Fort38(f2)

    merged_confnames = merge_lists(fort38_1.confnames, fort38_2.confnames)
    merged_trange = merge_lists(fort38_1.trange, fort38_2.trange)

    # group conformers to residues
    resnames = []
    conformers_in_res = []
    res_confs = {}
    old_res = merged_confnames[0][:3]+merged_confnames[0][5:11]
    for confname in merged_confnames:
        res = confname[:3]+confname[5:11]
        if old_res != res:  # new residue
            resnames.append(old_res)
            res_confs[old_res] = conformers_in_res
            conformers_in_res = [confname]
            old_res = res
        else:  # same residue
            conformers_in_res.append(confname)
    # last one
    resnames.append(res)
    res_confs[res] = conformers_in_res


    # for res in resnames:
    #     print(res)
    #     for conf in res_confs[res]:
    #         print("   %s" % conf)

    btd = []
    n_tpoints = len(merged_trange)
    for res in resnames:
        btd_this_residue = []
        confs = res_confs[res]
        resconf_match = True
        for p in merged_trange:
            if p in fort38_1.trange and p in fort38_2.trange:
                p1 = []
                p2 = []
                for conf in confs:
                    key = (conf, p)
                    if key in fort38_1.confocc and key in fort38_2.confocc:
                        p1.append(fort38_1.confocc[key])
                        p2.append(fort38_2.confocc[key])
                    else:
                        btd_this_residue.append("*")
                        resconf_match = False
                        break
                #print(p1, p2)
                if resconf_match:
                    d = bhata_distance(p1, p2)
                    btd_this_residue.append("%5.0f" % (d*1000))
            else:
                btd_this_residue.append("*")
        #print(btd_this_residue)
        btd.append(btd_this_residue)

    # print the distance table
    print(" %-8s %s" % (fort38_1.type, " ".join(["%5s" % x for x in merged_trange])))
    for ires in range(len(resnames)):
        print("%s %s" % (resnames[ires], " ".join(["%5s" % x for x in btd[ires]])))

    # if fort38_1.resnames == fort38_2.resnames:
    #     for ires in range(len(fort38_1.residues)):
    #         res = fort38_1.residues[ires]
    #         p1_all = [fort38_1.confocc[iconf] for iconf in res]
    #         p2_all = [fort38_2.confocc[iconf] for iconf in res]
    #         btd_this_residue = []
    #         for ipoint in range(len(p1_all[0])):
    #             d = bhata_distance([x1[ipoint] for x1 in p1_all], [x2[ipoint] for x2 in p2_all])
    #             #print(fort38_1.resnames[ires], d)
    #             btd_this_residue.append(d)
    #         btd.append(btd_this_residue)
    # else:
    #     print("Two fort.38 files don't have same residues")
    #
    # # print the distance table
    # print(" %-8s %s" % (fort38_1.type, " ".join(["%5.1f" % x for x in fort38_1.trange])))
    # for ires in range(len(fort38_1.residues)):
    #     print("%s %s" % (fort38_1.resnames[ires], " ".join(["%5.0f" % (x*1000) for x in btd[ires]])))

#    print(btd)

    return

if __name__ == "__main__":
    # Get the command arguments
    helpmsg = """Compare two fort.38 files and calculate the Bhattachayya distance at residue level. 
The number is reported as 1000 * Bhattacharyya distance. 
The residues in two files don't have to be the same. 
For the same residue, conformers have to be the same."""

    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("files", metavar="fort.38_file1 fort.38_file2", nargs=2)
    args = parser.parse_args()

    compare_fort38(args.files[0], args.files[1])

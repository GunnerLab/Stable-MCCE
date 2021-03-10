#!/usr/bin/env python

import sys
import math
import numpy as np

ph2Kcal = 1.364
Kcal2kT = 1.688


class Microstate:
    def __init__(self, state, E, count):
        self.state = state
        self.E = E
        self.count = count


class Conformer:
    def __init__(self):
        self.iconf = 0
        self.ires = 0
        self.confid = ""
        self.resid = ""
        self.occ = 0.0
        self.crg = 0.0

    def load_from_head3lst(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3]+self.confid[5:11]
        self.crg = float(fields[4])

class MSout:
    def __init__(self, fname):
        self.T = 273.15
        self.pH = 7.0
        self.Eh = 0.0
        self.N_ms = 0
        self.N_uniq = 0
        self.lowest_E = 0.0
        self.highest_E = 0.0
        self.average_E = 0.0
        self.fixed_iconfs = []
        self.fixed_crg = 0.0
        self.fixed_ne = 0.0
        self.fixed_nh = 0.0
        self.free_residues = []   # free residues, referred by conformer indices
        self.iconf2ires = {}      # from conformer index to free residue index
        self.microstates = {}
        self.conformers = []
        self.load_msout(fname)

    def load_msout(self, fname):
        lines = open(fname).readlines()

        # Get a valid line
        while True:
            line = lines.pop(0).strip()
            if len(line) > 0 and line[0] != "#":
                break

        fields = line.split(",")
        for field in fields:
            key, value = field.split(":")
            key = key.strip().upper()
            value = float(value)
            if key == "T":
                self.T = value
            elif key == "PH":
                self.pH = value
            elif key == "EH":
                self.Eh = value

        # second line, confirm this is from Monte Carlo sampleing
        while True:
            line = lines.pop(0).strip()
            if len(line) > 0 and line[0] != "#":
                break

        key, value = line.split(":")
        if key.strip() != "METHOD" or value.strip() != "MONTERUNS":
            print("This file %s is not a valid microstate file" % fname)
            sys.exit(-1)

        # Third line, fixed conformer indicies
        while True:
            line = lines.pop(0).strip()
            if len(line) > 0 and line[0] != "#":
                break

        _, iconfs = line.split(":")
        self.fixed_iconfs = [int(i) for i in iconfs.split()]

        # 4th line, free residues
        while True:
            line = lines.pop(0).strip()
            if len(line) > 0 and line[0] != "#":
                break

        _, residues_str = line.split(":")
        residues = residues_str.split(";")
        self.free_residues = []
        for f in residues:
            if f.strip():
                self.free_residues.append([int(i) for i in f.split()])
        for i_res in range(len(self.free_residues)):
            for iconf in self.free_residues[i_res]:
                self.iconf2ires[iconf] = i_res

        # find the next MC record
        found_mc = False
        newmc = False
        self.N_ms = 0

        for line in lines:
            if line.find("MC:") == 0:   # ms starts
                found_mc = True
                newmc = True
                continue
            elif newmc:
                f1, f2 = line.split(":")
                current_state = [int(c) for c in f2.split()]
                newmc = False
                continue
            elif found_mc:
                fields = line.split(",")
                if len(fields) >= 3:
                    state_e = float(fields[0])
                    count = int(fields[1])
                    flipped = [int(c) for c in fields[2].split()]

                    for ic in flipped:
                        ir = self.iconf2ires[ic]
                        current_state[ir] = ic

                    ms = Microstate(list(current_state), state_e, count)
                    key = ",".join(["%d" % i for i in ms.state])
                    if key in self.microstates:
                        self.microstates[key].count += ms.count
                    else:
                        self.microstates[key] = ms

        # find N_ms, lowerst, highst, averge E
        self.N_ms = 0
        E_sum = 0.0
        self.lowest_E = next(iter(self.microstates.values())).E
        self.highest_E = next(iter(self.microstates.values())).E
        msvalues = self.microstates.values()
        self.N_uniq = len(msvalues)
        for ms in msvalues:
            self.N_ms += ms.count
            E_sum += ms.E * ms.count
            if self.lowest_E > ms.E:
                self.lowest_E = ms.E
            if self.highest_E < ms.E:
                self.highest_E = ms.E
        self.average_E = E_sum / self.N_ms


def groupms_byenergy(microstates, ticks):
    """
    This function takes in a list of microstates and a list of energy numbers (N values), divide the microstates into N
    bands by using the energy number as lower boundaries. The list of energy will be sorted from small to large.
    """
    N = len(ticks)
    ticks.sort()
    ticks.append(1.0e100)    # add a big number as the rightest-most boundary
    resulted_bands = [[] for i in range(N)]

    for ms in microstates:
        it = -1
        for itick in range(N):
            if ticks[itick] <= ms.E < ticks[itick+1]:
                it = itick
                break
        if it >= 0:
            resulted_bands[it].append(ms)

    return resulted_bands


def groupms_byiconf(microstates, iconfs):
    """
    This function takes in a list of microstates and a list of conformer indicies, divide microstates into two groups:
    the first one is those contain one of the given conformers, the second one is those contain none of the listed conformers.
    """
    ingroup = []
    outgroup = []
    for ms in microstates:
        contain = False
        for ic in iconfs:
            if ic in ms.state:
                ingroup.append(ms)
                contain = True
                break
        if not contain:
            outgroup.append(ms)

    return ingroup, outgroup


def groupms_byconfid(microstates, confids):
    """
    Group conformers by tge conformer IDs. IDs are in a list and ID is considered as a match as long as it is a
    substring of the conformer name. The selected microstates must have all conformers and returned in the first group,
    and the rest are in the second group.
    """
    ingroup = []
    outgroup = []
    for ms in microstates:
        contain = True
        names = [conformers[ic].confid for ic in ms.state]
        for confid in confids:
            innames = False
            for name in names:
                if confid in name:
                    innames = True
                    break
            contain = contain and innames
        if contain:
            ingroup.append(ms)
        else:
            outgroup.append(ms)

    return ingroup, outgroup




def ms_energy_stat(microstates):
    """
    Given a list of microstates, find the lowest energy, average energy, and highest energy
    """
    ms = next(iter(microstates))
    lowerst_E = highest_E = ms.E
    N_ms = 0
    total_E = 0.0
    for ms in microstates:
        if lowerst_E > ms.E:
            lowerst_E = ms.E
        elif highest_E < ms.E:
            highest_E = ms.E
        N_ms += ms.count
        total_E += ms.E*ms.count

    average_E = total_E/N_ms

    return lowerst_E, average_E, highest_E


def ms_convert2occ(microstates):
    """
    Given a list of microstates, convert to conformer occupancy of conformers appeared at least once in the microstates.
    """
    occurance = {}  # occurance of conformer, as a dictionary
    occ = {}
    N_ms = 0
    for ms in microstates:
        N_ms += ms.count
        for ic in ms.state:
            if ic in occurance:
                occurance[ic] += ms.count
            else:
                occurance[ic] = ms.count

    for key in occurance.keys():
        occ[key] = occurance[key]/N_ms

    return occ



def ms_counts(microstates):
    """
    Calculate total counts of microstates
    """
    N_ms = 0
    for ms in microstates:
        N_ms += ms.count

    return N_ms


def ms_charge(ms):
    "Compute microstate charge"
    crg = 0.0
    for ic in ms.state:
        crg += conformers[ic].crg
    return crg


def ms_convert2sumcrg(microstates, free_res):
    """
    Given a list of microstates, convert to net charge of each free residue.
    """
    iconf2ires = {}
    for i_res in range(len(free_res)):
        for iconf in free_res[i_res]:
            iconf2ires[iconf] = i_res

    charges_total = [0.0 for i in range(len(free_res))]
    N_ms = 0
    for ms in microstates:
        N_ms += ms.count
        for ic in ms.state:
            ir = iconf2ires[ic]
            charges_total[ir] += conformers[ic].crg * ms.count

    charges = [x/N_ms for x in charges_total]

    return charges



def read_conformers():
    conformers = []
    lines = open("head3.lst").readlines()
    lines.pop(0)
    for line in lines:
        conf = Conformer()
        conf.load_from_head3lst(line)
        conformers.append(conf)

    return conformers


def e2occ(energies):
    "Given a list of energy values in unit Kacl/mol, calculate the occupancy by Boltzmann Distribution."
    e = np.array(energies)
    e = e - min(e)
    Pi_raw = np.exp(-Kcal2kT*e)
    Pi_sum = sum(Pi_raw)
    Pi_norm = Pi_raw/Pi_sum

    return Pi_norm


def bhata_distance(prob1, prob2):
    d_max = 10000.0   # Max possible value set to this
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

    return d


def whatchanged_conf(msgroup1, msgroup2):
    "Given two group of microstates, calculate what changed at conformer level."
    occ1 = ms_convert2occ(msgroup1)
    occ2 = ms_convert2occ(msgroup2)

    all_keys = set(occ1.keys())
    all_keys |= set(occ2.keys())

    all_keys = list(all_keys)
    all_keys.sort()
    diff_occ = {}
    for key in all_keys:
        if key in occ1:
            p1 = occ1[key]
        else:
            p1 = 0.0
        if key in occ2:
            p2 = occ2[key]
        else:
            p2 = 0.0
        diff_occ[key] = p2 - p1

    return diff_occ


def whatchanged_res(msgroup1, msgroup2, free_res):
    "Return a list of Bhatachaya Distance of free residues."
    occ1 = ms_convert2occ(msgroup1)
    occ2 = ms_convert2occ(msgroup2)

    bhd = []
    for res in free_res:
        p1 = []
        p2 = []
        for ic in res:
            if ic in occ1:
                p1.append(occ1[ic])
            else:
                p1.append(0.0)
            if ic in occ2:
                p2.append(occ2[ic])
            else:
                p2.append(0.0)
        bhd.append(bhata_distance(p1, p2))

    return bhd

conformers = read_conformers()



if __name__ == "__main__":
    msout = MSout("ms_out/pH4eH0ms.txt")
    # e_step = (msout.highest_E - msout.lowest_E)/20
    # ticks = [msout.lowest_E + e_step*(i) for i in range(20)]
    # ms_in_bands = groupms_byenergy(msout.microstates.values(), ticks)
    # print([len(band) for band in ms_in_bands])
#     netural, charged = groupms_byiconf(msout.microstates.values(), [12, 13, 14, 15])
#     l_E, a_E, h_E = ms_energy_stat(msout.microstates.values())
#     print(l_E, a_E, h_E)

    # charge over energy bands
    # e_step = (msout.highest_E - msout.lowest_E) / 20
    # ticks = [msout.lowest_E + e_step*(i+1) for i in range(19)]
    # ms_in_bands = groupms_byenergy(msout.microstates.values(), ticks)
    # for band in ms_in_bands:
    #     band_total_crg = 0.0
    #     for ms in band:
    #         band_total_crg += ms_charge(ms)
    #     print(band_total_crg/ms_counts(band))

    # netural, charged = groupms_byiconf(msout.microstates.values(), [12, 13, 14, 15])
    # diff_occ = whatchanged_conf(netural, charged)
    # for key in diff_occ.keys():
    #     print("%3d, %s: %6.3f" % (key, conformers[key].confid, diff_occ[key]))

    # diff_bhd = whatchanged_res(netural, charged, msout.free_residues)
    # for ir in range(len(msout.free_residues)):
    #     print("%s: %6.4f" % (conformers[msout.free_residues[ir][0]].resid, diff_bhd[ir]))
    # charges = ms_convert2sumcrg(msout.microstates.values(), msout.free_residues)
    # for ir in range(len(msout.free_residues)):
    #     print("%s: %6.4f" % (conformers[msout.free_residues[ir][0]].resid, charges[ir]))
    microstates = list(msout.microstates.values())
    glu35_charged, _ = groupms_byconfid(microstates, ["GLU-1A0035"])
    print(len(microstates))
    print(len(glu35_charged))
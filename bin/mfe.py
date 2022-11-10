#!/usr/bin/env python

import argparse
import sys
import os
import math

fname_fort38 = "fort.38"
fname_head3lst = "head3.lst"
fname_pkout = "pK.out"
fname_sumcrg = "sum_crg.out"

ph2Kcal = 1.364
mev2Kcal = 0.0235
Kcal2kT = 1.688

# Normally a residue is divided into two groups based on total charge. If more than two charge states are detected or
# a mfe analysis is needed for groups with the same charge, you can specify the grouping rule here.
#              RES  : ground  excited
Special_res = {"_CU": (["+1"], ["+2"]),
               "UbQ": (["01"], ["-1"]),
               "HOH": (["01"], ["DM"]),
	       "T4Y": (["01"], ["+1", "+2"])}

class E_IONIZE:
    def __init__(self):
        self.mfe = []
        return

class Conformer:
    def __init__(self):
        self.name = ""
        self.fl = ""
        self.occ = []
        self.crg = 0.0
        self.em0 = 0.0
        self.pka0 = 0.0
        self.ne = 0
        self.nh = 0
        self.vdw0 = 0.0
        self.vdw1 = 0.0
        self.tors = 0.0
        self.epol = 0.0
        self.dsolv = 0.0
        self.extra = 0.0
        return


class Residue:
    def __init__(self):
        self.resid = ""
        self.conformers = []
        return


class Titration:
    def __init__(self):
        self.range = []
        self.type = ""
        self.conformers = []
        self.residues = []
        return

    def read_confs(self, env):
        scale_ele = env["scale_ele"]
        scale_vdw = env["scale_vdw"]
        scale_vdw0 = env["scale_vdw0"]
        scale_vdw1 = env["scale_vdw1"]
        scale_tors = env["scale_tors"]
        scale_dsolv = env["scale_dsolv"]

        fort38_lines = open(fname_fort38).readlines()
        head3_lines = open(fname_head3lst).readlines()
        line = fort38_lines.pop(0)
        fields = line.split()
        self.type = fields[0]
        self.range = [float(x) for x in fields[1:]]
        head3_lines.pop(0)

        while fort38_lines and head3_lines:
            fort38_data = fort38_lines.pop(0)
            head3_data = head3_lines.pop(0)
            conf = Conformer()

            fields = fort38_data.split()
            conf.name = fields[0]
            conf.occ = [float(x) for x in fields[1:]]

            fields = head3_data.split()
            if conf.name != fields[1]:
                print("Confomer %s in fort.38 does not match %s in head3.lst" % (conf.name, fields[1]))
                sys.exit()

            conf.fl = fields[2]
            conf.hocc = float(fields[3])  # occ in head3.lst
            conf.crg = float(fields[4])
            conf.em0 = float(fields[5])
            conf.pka0 = float(fields[6])
            conf.ne = int(fields[7])
            conf.nh = int(fields[8])
            conf.vdw0 = float(fields[9]) * scale_vdw0
            conf.vdw1 = float(fields[10]) * scale_vdw1
            conf.tors = float(fields[11]) * scale_tors
            conf.epol = float(fields[12]) * scale_ele
            conf.dsolv = float(fields[13]) * scale_dsolv
            conf.extra = float(fields[14])
            conf.self = conf.vdw0 + conf.vdw1 + conf.tors + conf.epol + conf.dsolv + conf.extra

            self.conformers.append(conf)

        return

    def group_residues(self):
        res_ids = []
        residues = {}
        for conf in self.conformers:
            resid = conf.name[:3] + conf.name[5:11]
            if resid in res_ids:
                res = residues[resid]
            else:
                res = Residue()
                res.resid = resid
                residues[resid] = res
                res_ids.append(resid)
            res.conformers.append(conf)

        self.residues = [residues[x] for x in res_ids]
        return


def get_mfe(env, residue, mfe_ph, mfe_xts, mfe_pwcut, pKa):
    scale_ele = env["scale_ele"]
    scale_vdw = env["scale_vdw"]
    scale_vdw0 = env["scale_vdw0"]
    scale_vdw1 = env["scale_vdw1"]
    scale_tors = env["scale_tors"]
    scale_dsolv = env["scale_dsolv"]

    this_resid = residue[:3]+residue[4:]

    ph1 = 7.0
    if "TITR_PH0" in env:
        ph1 = float(env["TITR_PH0"])

    eh1 = 0.0
    if "TITR_EH0" in env:
        eh1 = float(env["TITR_EH0"])

    report_lines = []

    titration = Titration()
    titration.read_confs(env)
    titration.group_residues()

    mfe_conformers = []
    resid = residue[:3]+residue[4:]
    found = False
    for res in titration.residues:
        if resid ==(res.resid):
            mfe_conformers = res.conformers
            found = True
            break

    if not found:
        print("Residue %s not found in head3.lst" % resid)
        sys.exit()

    # divide in to two states:
    ground_conformers = []
    exited_conformers = []
    if residue[:3] in Special_res:  # group by type name if defined
        ground_typeid = Special_res[residue[:3]][0]
        exited_typeid = Special_res[residue[:3]][1]
        for conf in mfe_conformers:
            if conf.name[3:5] in ground_typeid:
                ground_conformers.append(conf)
            elif conf.name[3:5] in exited_typeid:
                exited_conformers.append(conf)
    else:   # group conformers by charge if not defined
        for conf in mfe_conformers:
            if abs(conf.crg) < 0.01:
                ground_conformers.append(conf)
            else:
                exited_conformers.append(conf)

    for conf in mfe_conformers:
        conf.pHeffect = [0.0 for i in range(len(titration.range))]
        conf.Eheffect = [0.0 for i in range(len(titration.range))]
        conf.res_mfe = [[0.0 for i in range(len(titration.range))] for x in titration.residues]
        conf.mfe_total = [0.0 for i in range(len(titration.range))]
        conf.E_total = [0.0 for i in range(len(titration.range))]

        # pairwise energy table
        pairwise = {}
        opp_path = "energies/" + conf.name + ".opp"
        if os.path.exists(opp_path):
            lines = open(opp_path).readlines()
            for line in lines:
                fields = line.split()
                if len(fields) >= 3:
                    # check the vdw clash
                    if fields[3].find("999.000") >= 0:
                        pairwise[fields[1]] = 999.000
                    else:
                        if len(fields[3]) > 8:
                            fields[3] = fields[3][:-8]
                        pairwise[fields[1]] = float(fields[2]) * scale_ele + float(fields[3]) * scale_vdw

        # make conformer mfe
        for i in range(len(titration.range)):

            point = titration.range[i]
            conf_mfe = [0.0 for x in titration.range]

            # pH effect in Kcal/mol
            if titration.type.upper() == 'PH':
                conf.pHeffect[i] = (point - conf.pka0) * conf.nh * ph2Kcal
            else:
                conf.pHeffect[i] = (ph1 - conf.pka0) * conf.nh * ph2Kcal

            # Eh effect in Kcal/mol
            if titration.type.upper() == 'EH':
                conf.Eheffect[i] = (point - conf.em0) * conf.ne * mev2Kcal
            else:
                conf.Eheffect[i] = (eh1 - conf.em0) * conf.ne * mev2Kcal


            for j in range(len(titration.residues)):
                res = titration.residues[j]
                if res.resid == this_resid:
                    mfe = 0.0
                else:
                    mfe = 0.0
                    for conf2 in res.conformers:
                        if conf2.name in pairwise:
                            mfe += pairwise[conf2.name] * conf2.occ[i]
                # This mfe is at 1 titration point, from one residue

                conf.res_mfe[j][i] = mfe
                conf.mfe_total[i] += mfe

            # update conformer E_total
            conf.E_total[i] = conf.mfe_total[i] \
                                   + conf.pHeffect[i] \
                                   + conf.Eheffect[i] \
                                   + conf.self

    # Ground state recovered occ
    # Reference energy: lowest E
    Eref = ground_conformers[0].E_total[0]  # reference E (lowest of this res)
    for conf in ground_conformers:
        for i in range(len(titration.range)):
            if Eref > conf.E_total[i]:
                Eref = conf.E_total[i]

    # Calculate mfe occupancy of each conformer
    SigmaE = [0.0 for i in range(len(titration.range))]
    for i in range(len(titration.range)):
        Ei = [math.exp(-(conformer.E_total[i] - Eref) * Kcal2kT) for conformer in ground_conformers]
        SigmaE[i] = sum(Ei)

    for conformer in ground_conformers:
        conformer.rocc = [0.0 for x in titration.range]
        for i in range(len(titration.range)):
            if conformer.fl.upper() == 'T':
                # print conformer.hocc
                if conformer.hocc < 0.001:
                    conformer.rocc[i] = 0.0
                elif conformer.hocc > 0.999:
                    conformer.rocc[i] = 1.00
                else:
                    print()
                    "Error: partial occupancy was assigned in head3.lst"
                    sys.exit()
            else:
                conformer.rocc[i] = math.exp(-(conformer.E_total[i] - Eref) * Kcal2kT) / SigmaE[i]  # recovered occ

    # Exited state recovered occ
    # Reference energy: lowest E
    Eref = exited_conformers[0].E_total[0]  # reference E (lowest of this res)
    for conf in exited_conformers:
        for i in range(len(titration.range)):
            if Eref > conf.E_total[i]:
                Eref = conf.E_total[i]

    # Calculate mfe occupancy of each conformer
    SigmaE = [0.0 for i in range(len(titration.range))]
    for i in range(len(titration.range)):
        Ei = [math.exp(-(conformer.E_total[i] - Eref) * Kcal2kT) for conformer in exited_conformers]
        SigmaE[i] = sum(Ei)

    for conformer in exited_conformers:
        conformer.rocc = [0.0 for x in titration.range]
        for i in range(len(titration.range)):
            if conformer.fl.upper() == 'T':
                # print conformer.hocc
                if conformer.hocc < 0.001:
                    conformer.rocc[i] = 0.0
                elif conformer.hocc > 0.999:
                    conformer.rocc[i] = 1.00
                else:
                    print()
                    "Error: partial occupancy was assigned in head3.lst"
                    sys.exit()
            else:
                conformer.rocc[i] = math.exp(-(conformer.E_total[i] - Eref) * Kcal2kT) / SigmaE[i]  # recovered occ

    SigmaOcc = [0.0 for i in range(len(titration.range))]
    for conformer in ground_conformers:
        for i in range(len(titration.range)): SigmaOcc[i] += conformer.rocc[i]
    for conformer in ground_conformers:
        conformer.nocc = []
        for i in range(len(titration.range)):
            if SigmaOcc[i] < 1.0E-25 and conformer.rocc[i] < 1.0E-25:
                conformer.nocc.append(1.0)
            else:
                conformer.nocc.append(conformer.rocc[i] / SigmaOcc[i])

    SigmaOcc = [0.0 for i in range(len(titration.range))]
    for conformer in exited_conformers:
        for i in range(len(titration.range)): SigmaOcc[i] += conformer.rocc[i]
    for conformer in exited_conformers:
        conformer.nocc = []
        for i in range(len(titration.range)):
            if SigmaOcc[i] < 1.0E-25 and conformer.rocc[i] < 1.0E-25:
                conformer.nocc.append(1.0)
            else:
                conformer.nocc.append(conformer.rocc[i] / SigmaOcc[i])

    # energy terms of ground state
    ground_state = E_IONIZE()
    ground_state.vdw0 = [0.0 for x in titration.range]
    ground_state.vdw1 = [0.0 for x in titration.range]
    ground_state.tors = [0.0 for x in titration.range]
    ground_state.epol = [0.0 for x in titration.range]
    ground_state.dsolv = [0.0 for x in titration.range]
    ground_state.extra = [0.0 for x in titration.range]
    ground_state.pHeffect = [0.0 for x in titration.range]
    ground_state.Eheffect = [0.0 for x in titration.range]
    ground_state.mfe_total = [0.0 for x in titration.range]
    ground_state.res_mfe = [[0.0 for x in titration.range] for x in titration.residues]
    ground_state.E_total = [0.0 for x in titration.range]
    ground_state.TS = [0.0 for x in titration.range]
    for conformer in ground_conformers:
        for i in range(len(titration.range)):
            ground_state.vdw0[i] += conformer.nocc[i] * conformer.vdw0
            ground_state.vdw1[i] += conformer.nocc[i] * conformer.vdw1
            ground_state.tors[i] += conformer.nocc[i] * conformer.tors
            ground_state.epol[i] += conformer.nocc[i] * conformer.epol
            ground_state.dsolv[i] += conformer.nocc[i] * conformer.dsolv
            ground_state.extra[i] += conformer.nocc[i] * conformer.extra
            ground_state.pHeffect[i] += conformer.nocc[i] * conformer.pHeffect[i]
            ground_state.Eheffect[i] += conformer.nocc[i] * conformer.Eheffect[i]
            ground_state.mfe_total[i] += conformer.nocc[i] * conformer.mfe_total[i]
            ground_state.E_total[i] += conformer.nocc[i] * conformer.E_total[i]
            if conformer.nocc[i] > 0.000001:
                ground_state.TS[i] += -conformer.nocc[i] * math.log(conformer.nocc[i]) / Kcal2kT
        for j in range(len(titration.residues)):
            for i in range(len(titration.range)):
                ground_state.res_mfe[j][i] += conformer.nocc[i] * conformer.res_mfe[j][i]

    if mfe_xts:
        ground_state.G = [ground_state.E_total[i] - ground_state.TS[i] for i in range(len(titration.range))]
    else:
        ground_state.G = [ground_state.E_total[i] for i in range(len(titration.range))]

    # energy terms of charged state
    charged_state = E_IONIZE()
    charged_state.vdw0 = [0.0 for x in titration.range]
    charged_state.vdw1 = [0.0 for x in titration.range]
    charged_state.tors = [0.0 for x in titration.range]
    charged_state.epol = [0.0 for x in titration.range]
    charged_state.dsolv = [0.0 for x in titration.range]
    charged_state.extra = [0.0 for x in titration.range]
    charged_state.pHeffect = [0.0 for x in titration.range]
    charged_state.Eheffect = [0.0 for x in titration.range]
    charged_state.res_mfe = [[0.0 for x in titration.range] for x in titration.residues]
    charged_state.mfe_total = [0.0 for x in titration.range]
    charged_state.E_total = [0.0 for x in titration.range]
    charged_state.TS = [0.0 for x in titration.range]
    for conformer in exited_conformers:
        for i in range(len(titration.range)):
            charged_state.vdw0[i] += conformer.nocc[i] * conformer.vdw0
            charged_state.vdw1[i] += conformer.nocc[i] * conformer.vdw1
            charged_state.tors[i] += conformer.nocc[i] * conformer.tors
            charged_state.epol[i] += conformer.nocc[i] * conformer.epol
            charged_state.dsolv[i] += conformer.nocc[i] * conformer.dsolv
            charged_state.extra[i] += conformer.nocc[i] * conformer.extra
            charged_state.pHeffect[i] += conformer.nocc[i] * conformer.pHeffect[i]
            charged_state.Eheffect[i] += conformer.nocc[i] * conformer.Eheffect[i]
            charged_state.mfe_total[i] += conformer.nocc[i] * conformer.mfe_total[i]
            charged_state.E_total[i] += conformer.nocc[i] * conformer.E_total[i]
            if conformer.nocc[i] > 0.000001:
                charged_state.TS[i] += -conformer.nocc[i] * math.log(conformer.nocc[i]) / Kcal2kT
        for j in range(len(titration.residues)):
            for i in range(len(titration.range)):
                charged_state.res_mfe[j][i] += conformer.nocc[i] * conformer.res_mfe[j][i]

    if mfe_xts:
        charged_state.G = [charged_state.E_total[i] - charged_state.TS[i] for i in range(len(titration.range))]
    else:
        charged_state.G = [charged_state.E_total[i] for i in range(len(titration.range))]

    dG = E_IONIZE()
    dG.ground_state = ground_state
    dG.charged_state = charged_state
    dG.ground_confs = ground_conformers
    dG.charged_confs = exited_conformers
    dG.resID = this_resid

    # decide which two columns will be used to get residue mfe
    t_low_found = t_high_found = 0
    i_low = 0
    i_high = 0
    for i in range(len(titration.range)):
        if mfe_ph > titration.range[i] - 0.001:
            i_low = i
            t_low_found = 1

    for i in range(len(titration.range)):
        if mfe_ph < titration.range[i] + 0.001:
            i_high = i
            t_high_found = 1
            break

    if not t_low_found or not t_high_found:
        print("titration_point out off range")
        sys.exit()

    # get to this pH/Eh
    dG_point = E_IONIZE()

    if i_low == i_high:  # one point
        dG_point.vdw0 = charged_state.vdw0[i_low] - ground_state.vdw0[i_low]
        dG_point.vdw1 = charged_state.vdw1[i_low] - ground_state.vdw1[i_low]
        dG_point.tors = charged_state.tors[i_low] - ground_state.tors[i_low]
        dG_point.epol = charged_state.epol[i_low] - ground_state.epol[i_low]
        dG_point.dsolv = charged_state.dsolv[i_low] - ground_state.dsolv[i_low]
        dG_point.extra = charged_state.extra[i_low] - ground_state.extra[i_low]
        dG_point.pHeffect = charged_state.pHeffect[i_low] - ground_state.pHeffect[i_low]
        dG_point.Eheffect = charged_state.Eheffect[i_low] - ground_state.Eheffect[i_low]
        dG_point.mfe_total = charged_state.mfe_total[i_low] - ground_state.mfe_total[i_low]
        dG_point.E_total = charged_state.E_total[i_low] - ground_state.E_total[i_low]
        # entropy should be used only when MC is carried out
        dG_point.TS = charged_state.TS[i_low] - ground_state.TS[i_low]

        for i in range(len(charged_state.res_mfe)):
            dG_point.mfe.append(charged_state.res_mfe[i][i_low] - ground_state.res_mfe[i][i_low])
        dG_point.G = charged_state.G[i_low] - ground_state.G[i_low]

    else:  # scale average of two points
        k = (mfe_ph - titration.range[i_low]) / (titration.range[i_high] - titration.range[i_low])
        dG_point.vdw0 = (1 - k) * (charged_state.vdw0[i_low] - ground_state.vdw0[i_low]) \
                        + k * (charged_state.vdw0[i_high] - ground_state.vdw0[i_high])
        dG_point.vdw1 = (1 - k) * (charged_state.vdw1[i_low] - ground_state.vdw1[i_low]) \
                        + k * (charged_state.vdw1[i_high] - ground_state.vdw1[i_high])
        dG_point.tors = (1 - k) * (charged_state.tors[i_low] - ground_state.tors[i_low]) \
                        + k * (charged_state.tors[i_high] - ground_state.tors[i_high])
        dG_point.epol = (1 - k) * (charged_state.epol[i_low] - ground_state.epol[i_low]) \
                        + k * (charged_state.epol[i_high] - ground_state.epol[i_high])
        dG_point.dsolv = (1 - k) * (charged_state.dsolv[i_low] - ground_state.dsolv[i_low]) \
                         + k * (charged_state.dsolv[i_high] - ground_state.dsolv[i_high])
        dG_point.extra = (1 - k) * (charged_state.extra[i_low] - ground_state.extra[i_low]) \
                         + k * (charged_state.extra[i_high] - ground_state.extra[i_high])
        dG_point.pHeffect = (1 - k) * (charged_state.pHeffect[i_low] - ground_state.pHeffect[i_low]) \
                            + k * (charged_state.pHeffect[i_high] - ground_state.pHeffect[i_high])
        dG_point.Eheffect = (1 - k) * (charged_state.Eheffect[i_low] - ground_state.Eheffect[i_low]) \
                            + k * (charged_state.Eheffect[i_high] - ground_state.Eheffect[i_high])
        dG_point.mfe_total = (1 - k) * (charged_state.mfe_total[i_low] - ground_state.mfe_total[i_low]) \
                             + k * (charged_state.mfe_total[i_high] - ground_state.mfe_total[i_high])
        dG_point.E_total = (1 - k) * (charged_state.E_total[i_low] - ground_state.E_total[i_low]) \
                           + k * (charged_state.E_total[i_high] - ground_state.E_total[i_high])
        dG_point.TS = (1-k)*(charged_state.TS[i_low] - ground_state.TS[i_low])\
                        + k*(charged_state.TS[i_high] - ground_state.TS[i_high])

        for i in range(len(charged_state.res_mfe)):
            dG_point.mfe.append((1 - k) * (charged_state.res_mfe[i][i_low] - ground_state.res_mfe[i][i_low]) \
                                + k * (charged_state.res_mfe[i][i_high] - ground_state.res_mfe[i][i_high]))

        dG_point.G = (1 - k) * (charged_state.G[i_low] - ground_state.G[i_low]) \
                     + k * (charged_state.G[i_high] - ground_state.G[i_high])


    # net charge of residues at mfe_ph
    res_crg = {}
    lines = open(fname_sumcrg).readlines()
    lines.pop(0)
    for line in lines:
        fields = line.split()
        if len(fields) > 1:
            key = fields[0][:3]+fields[0][4:]
            res_crg[key] = [float(x) for x in fields[1:]]

    for res in titration.residues:
        if res.resid in res_crg:
            res.crg = res_crg[res.resid]
        else:
            res.crg = [0.0 for x in titration.range]

        if i_low == i_high:
            res.point_crg = res.crg[i_low]
        else:
            res.point_crg = (1-k)*res.crg[i_low] + k*res.crg[i_high]

    report_lines.append("Residue %s pKa/Em=%s\n" % (residue, pKa))
    report_lines.append("=================================\n")
    report_lines.append("Terms          pH     meV    Kcal\n")
    report_lines.append("---------------------------------\n")
    report_lines.append("vdw0     %8.2f%8.2f%8.2f\n" % (dG_point.vdw0 / ph2Kcal, dG_point.vdw0 / mev2Kcal, dG_point.vdw0))
    report_lines.append("vdw1     %8.2f%8.2f%8.2f\n" % (dG_point.vdw1 / ph2Kcal, dG_point.vdw1 / mev2Kcal, dG_point.vdw1))
    report_lines.append("tors     %8.2f%8.2f%8.2f\n" % (dG_point.tors / ph2Kcal, dG_point.tors / mev2Kcal, dG_point.tors))
    report_lines.append("ebkb     %8.2f%8.2f%8.2f\n" % (dG_point.epol / ph2Kcal, dG_point.epol / mev2Kcal, dG_point.epol))
    report_lines.append("dsol     %8.2f%8.2f%8.2f\n" % (dG_point.dsolv / ph2Kcal, dG_point.dsolv / mev2Kcal, dG_point.dsolv))
    report_lines.append("offset   %8.2f%8.2f%8.2f\n" % (dG_point.extra / ph2Kcal, dG_point.extra / mev2Kcal, dG_point.extra))
    report_lines.append("pH&pK0   %8.2f%8.2f%8.2f\n" % (dG_point.pHeffect / ph2Kcal, dG_point.pHeffect / mev2Kcal, dG_point.pHeffect))
    report_lines.append("Eh&Em0   %8.2f%8.2f%8.2f\n" % (dG_point.Eheffect / ph2Kcal, dG_point.Eheffect / mev2Kcal, dG_point.Eheffect))
    if mfe_xts:
        report_lines.append("-TS      %8.2f%8.2f%8.2f\n" % (-dG_point.TS/ph2Kcal, -dG_point.TS/mev2Kcal, -dG_point.TS))
    else:
        report_lines.append("-TS      %8.2f%8.2f%8.2f\n" % (0.0 / ph2Kcal, 0.0 / mev2Kcal, 0.0))
    report_lines.append("residues %8.2f%8.2f%8.2f\n" % (dG_point.mfe_total / ph2Kcal, dG_point.mfe_total / mev2Kcal, dG_point.mfe_total))
    report_lines.append("*********************************\n")
    report_lines.append("TOTAL    %8.2f%8.2f%8.2f%8.9s\n" % (dG_point.G / ph2Kcal, dG_point.G / mev2Kcal,
                                       dG_point.G, "  sum_crg"))
    report_lines.append("*********************************\n")
    for i in range(len(dG_point.mfe)):
        if abs(dG_point.mfe[i] / ph2Kcal) > mfe_pwcut:
            ion_state = titration.residues[i].point_crg
            report_lines.append("%-9s%8.2f%8.2f%8.2f%8.2f\n" % (titration.residues[i].resid, dG_point.mfe[i] / ph2Kcal,
                                          dG_point.mfe[i] / mev2Kcal, dG_point.mfe[i], ion_state))
    report_lines.append("=================================\n")


    return report_lines


def get_residue_pKa(res_name):
    pka = "Residue %s not found in file %s" % (res_name, fname_pkout)
    lines = open(fname_pkout).readlines()
    for line in lines:
        fields = line.split()

        if fields[0] == res_name:
            try:
                pka = float(fields[1])
                break
            except ValueError:
                pka = "Titration of residue %s out of range" % res_name

    return pka


def load_runprm(env):
    lines = open("run.prm").readlines()
    for line in lines:
        fields = line.split("#", 1)[0].strip().split()
        if len(fields) > 1 and fields[-1].find("(") >= 0 and fields[-1].find(")") > fields[-1].find("("):
            key = fields[-1].strip("()").strip()
            env[key] = fields[0]

    if "MONTE_TSX" not in env:
        env["MONTE_TSX"] = "f"

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

    mfe_xts = True
    if args.x == "t":
        mfe_xts = True
    elif args.x == "f":
        mfe_xts = False
    elif env["MONTE_TSX"].upper() == "T":
        mfe_xts = False    # Did entropy correction in Monte Carlo already, TS term should not be in mfe
    else:
        mfe_xts = True

    # mfe_ph
    pKa = get_residue_pKa(args.residue[0])
    if args.p == "m":  # determined by pKa
        mfe_ph = pKa
    else:
        mfe_ph = float(args.p)

    if isinstance(mfe_ph, float):
        mfe_pwcut = float(args.c)

        # By now we have residue name in args.residue[0], mfe_ph, mfe_xts and pw cutoff in arges.c
        report = get_mfe(env=env, residue=args.residue[0], mfe_ph=mfe_ph, mfe_xts=mfe_xts, mfe_pwcut=mfe_pwcut, pKa = pKa)
    else:
        report = ["%s\n" % mfe_ph]

    return report


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

    report = mfe(args)
    sys.stdout.writelines(report)


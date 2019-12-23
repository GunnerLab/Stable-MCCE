#!/usr/bin/env python
""" PyMCCE data structure """
import sys
import os.path
import numpy as np

ROOMT = 298.15
PH2KCAL = 1.364

float_values = ["EPSILON_PROT", "TITR_PH0", "TITR_PHD", "TITR_EH0", "TITR_EHD", "CLASH_DISTANCE",
                "BIG_PAIRWISE", "MONTE_T", "MONTE_REDUCE", "EXTRAE", "SCALING"]
int_values = ["TITR_STEPS", "MONTE_RUNS", "MONTE_TRACE", "MONTE_NITER", "MONTE_NEQ",
              "MONTE_NSTART", "MONTE_FLIPS"]


class Env:
    def __init__(self):
        # Hard define values
        self.version = "PyMCCE 0.1"
        self.fn_runprm = "run.prm"
        self.fn_conflist1 = "head1.lst"
        self.fn_conflist2 = "head2.lst"
        self.fn_conflist3 = "head3.lst"
        self.energy_table = "energies"
        # run.prm parameters key:value
        self.tpl = {}
        # tpl parameters (key1, key2, key3):value
        self.runprm = {}

        # load parameters
        self.load_runprm()
        self.read_extra()
        return

    def load_runprm(self):
        lines = open(self.fn_runprm).readlines()
        # Sample line: "t        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)"
        for line in lines:
            line = line.strip()
            line = line.split("#")[0]  # This cuts off everything after #
            left_p = line.rfind("(")
            right_p = line.rfind(")")

            if left_p > 0 and right_p > left_p + 1:
                key = line[left_p + 1:right_p]
                fields = line[:left_p].split()
                if len(fields) < 1:
                    continue
                else:
                    value = fields[0]
            else:
                continue

            self.set(key, value)
        return

    def set(self, key, value):
        # Known non-string value are converted, otherwise string presumed
        if key in float_values:
            self.runprm[key] = float(value)
        elif key in int_values:
            self.runprm[key] = int(value)
        else:
            self.runprm[key] = value
        return

    def print_runprm(self):
        for key in self.runprm.keys():
            print("%-25s:%s" % (key, str(self.runprm[key])))
        return

    def load_tpl(self, fname):
        """Load a tpl file."""
        lines = open(fname).readlines()
        for line in lines:
            line = line.split("#")[0]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            keys = [x.strip().strip("\"") for x in keys]
            keys = [x for x in keys if x]
            keys = tuple(keys)

            value_string = fields[1].strip()
            if keys[0] in float_values:
                self.tpl[keys] = float(value_string)
            elif keys[0] in int_values:
                self.tpl[keys] = int(value_string)
            else:
                self.tpl[keys] = value_string

        return

    def print_tpl(self):
        for key in self.tpl.keys():
            print("%-25s:%s" % (key, str(self.tpl[key])))
        return

    def read_extra(self):
        """Read extra.tpl."""
        self.load_tpl(self.runprm["EXTRA"])
        default_values_keys = [("SCALING", "VDW0"),
                               ("SCALING", "VDW1"),
                               ("SCALING", "VDW"),
                               ("SCALING", "TORS"),
                               ("SCALING", "ELE"),
                               ("SCALING", "DSOLV")]
        for element in default_values_keys:
            if element not in self.tpl:
                self.tpl[element] = 1.0

        return

    def print_scaling(self):
        """Print scaling factors."""
        # print self.param
        print("   Scaling factors:")
        print("   VDW0  = %.3f" % self.tpl[("SCALING", "VDW0")])
        print("   VDW1  = %.3f" % self.tpl[("SCALING", "VDW1")])
        print("   VDW   = %.3f" % self.tpl[("SCALING", "VDW")])
        print("   TORS  = %.3f" % self.tpl[("SCALING", "TORS")])
        print("   ELE   = %.3f" % self.tpl[("SCALING", "ELE")])
        print("   DSOLV = %.3f" % self.tpl[("SCALING", "DSOLV")])
        print("   Done\n")
        return


class Conformer:
    """ Conformer structure """

    def __init__(self, fields):
        # from head3.lst
        self.iConf = int(fields[0])
        self.confname = fields[1]
        self.flag = fields[2].lower()
        self.occ = float(fields[3])
        self.crg = float(fields[4])
        self.em0 = float(fields[5])
        self.pk0 = float(fields[6])
        self.ne = int(fields[7])
        self.nh = int(fields[8])
        self.vdw0 = float(fields[9]) * env.tpl[("SCALING", "VDW0")]
        self.vdw1 = float(fields[10]) * env.tpl[("SCALING", "VDW1")]
        self.tors = float(fields[11]) * env.tpl[("SCALING", "TORS")]
        self.epol = float(fields[12]) * env.tpl[("SCALING", "ELE")]
        self.dsolv = float(fields[13]) * env.tpl[("SCALING", "DSOLV")]
        self.extra = float(fields[14])
        self.history = fields[15]

        # from MC entropy sampling
        self.entropy = 0.0  # TS in kcal/mol, will be calculated at entropy sampling

        # needed by MC process
        self.E_self = 0.0  # self energy in head3.lst
        self.E_self_mfe = 0.0  # self energy including pairwise contribution from fixed residues
        self.counter = 0  # MC counters
        self.mc_occ = 0.0  # MC occ
        self.acc_occ = []  # a list of past mc_occ until rest, history is needed to test convergence
        self.occ_at_points = []  # average occ at titration points
        return

    def reset(self):
        """ Reset to head3lstoriginal """
        # from MC entropy sampling
        self.entropy = 0.0  # TS in kcal/mol, will be calculated at entropy sampling

        # needed by MC process
        self.E_self = 0.0  # self energy in head3.lst
        self.E_self_mfe = 0.0  # self energy including pairwise contribution from fixed residues
        self.counter = 0  # MC counters
        self.mc_occ = 0.0  # MC occ
        self.acc_occ = []  # a list of past mc_occ until rest, history is needed to test convergence
        self.occ_at_points = []  # average occ at titration points
        return

    def printme(self):
        print("%05d %s %c %4.2f %6.3f %5d %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %s" % (self.iConf,
                                                                                                         self.confname,
                                                                                                         self.flag,
                                                                                                         self.occ,
                                                                                                         self.crg,
                                                                                                         self.em0,
                                                                                                         self.pk0,
                                                                                                         self.ne,
                                                                                                         self.nh,
                                                                                                         self.vdw0,
                                                                                                         self.vdw1,
                                                                                                         self.tors,
                                                                                                         self.epol,
                                                                                                         self.dsolv,
                                                                                                         self.extra,
                                                                                                         self.entropy,
                                                                                                         self.history))


def load_head3lst():
    conformers = []
    fname = env.fn_conflist3
    lines = open(fname).readlines()
    lines.pop(0)
    for line in lines:
        fields = line.split()
        if len(fields) >= 16:
            conf = Conformer(fields)
            conformers.append(conf)

    confnames = [x.confname for x in conformers]
    for name in confnames:
        if len(name) != 14:
            print("%s is not a conformer name.")
            sys.exit()
        occurrence = confnames.count(name)
        if occurrence > 1:
            print("Conformer %s occurred %d times" % occurrence)
            sys.exit()

    return conformers

def load_pairwise():
    folder = env.energy_table
    n_conf = len(head3lst)
    confnames = [x.confname for x in head3lst]
    scale_ele = env.tpl[("SCALING", "ELE")]
    scale_vdw = env.tpl[("SCALING", "VDW")]
    pw = np.zeros(shape=(n_conf, n_conf))
    for ic in range(n_conf):
        conf = head3lst[ic]
        oppfile = "%s/%s.opp" % (folder, conf.confname)
        if os.path.isfile(oppfile):
            lines = open(oppfile)
            for line in lines:
                fields = line.split()
                if len(fields) < 6:
                    continue
                confname = fields[1]
                jc = confnames.index(confname)
                if jc < 0:
                    print("      Warning: %s in file %s is not a conformer" % (confname, oppfile))
                    continue
                ele = float(fields[2])
                vdw = float(fields[3])
                pw[ic][jc] = ele * scale_ele + vdw * scale_vdw

    # average the opposite sides
    for ic in range(n_conf-1):
        for jc in range(ic+1, n_conf):
            #if abs(pw[ic][jc] - pw[jc][ic]) > 0.000001:
            #    print("%s %.3f <-> %s %.3f" % (confnames[ic], pw[ic][jc], confnames[jc], pw[jc][ic]))
            averaged_pw = (pw[ic][jc] + pw[jc][ic]) * 0.5
            pw[ic][jc] = pw[jc][ic] = averaged_pw

    return pw

env = Env()
head3lst = load_head3lst()
pairwise = load_pairwise()


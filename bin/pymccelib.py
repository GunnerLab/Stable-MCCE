#!/usr/bin/env python

import logging

class Env:
    def __init__(self):
        # hard code values
        self.tpl = {}
        self.atomnames = {}
        self.radius = {}
        self.charge = {}
        return

    def load_ftpl(self, fname):
        """Load a ftpl file"""

        # default string values, or defined as below
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        logging.info("Loading ftpl file %s" % fname)
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

            # Make an atom list in the natural order of CONNECT record.
            if keys[0] == "CONNECT":
                atom = keys[1]
                conf = keys[2]
                if conf in self.atomnames:
                    self.atomnames[conf].append(atom)
                else:
                    self.atomnames[conf] = [atom]
        return

class Atom:
    def __init__(self):
        self.atomname = ""
        self.confname = ""
        self.resname = ""
        self.on = False
        self.iatom = "0"
        return

class Conformer:
    def __init__(self):
        self.confname = ""
        self.resname = ""
        self.atoms = []
        return

class Residue:
    def __init__(self):
        self.resname = []
        self.conformers = []
        return

class Protein:
    """Protein structure"""
    def __init__(self):
        self.residues = []
        return


    def pdb2mcce(self, pdb):
        """Convert pdb to mcce pdb"""
        atom_exceptions = [" H2 ", " OXT", " HXT"]
        mccelines = []
        lines = [x for x in open(pdb).readlines() if x[:6] == "ATOM  " or x[:6] == "HETATM"]

        icount = 0
        previous_resid = ()
        possible_confs = []
        for line in lines:
            # pdb line
            atomname = line[12:16]
            resname = line[17:20]
            chainid = line[21]
            seqnum = int(line[22:26])
            icode = line[26]
            xyz = line[30:54]

            current_resid = (resname, chainid, seqnum, icode)
            # mcce line, need to add conf_number, radius, charge, conf_type, conf_history
            if current_resid != previous_resid:
                possible_confs = [x.strip() for x in env.tpl[("CONFLIST", resname)].split(",")]
                logging.info("Identified a new residue %s: %s" % (resname, ", ".join(possible_confs)))
                previous_resid = current_resid
            Found = False
            for confname in possible_confs:
                if atomname in env.atomnames[confname]:
                    conf_type = confname[3:5]
                    conf_number = possible_confs.index(confname)
                    cname = confname
                    Found = True
                    break
            if not Found:
                # this atom is not found in all conformers
                if atomname not in atom_exceptions:
                    print("Atom \"%s\" in pdb file %s can not be assigned to any conformer" % (atomname, pdb))
                continue

            key = ("RADIUS", cname, atomname)
            if key in env.tpl:
                radius_str = env.tpl[key]
                rad, _, _ = radius_str.split(",")
                rad = float(rad)
            else:
                rad = 0.0

            key = ("CHARGE", cname, atomname)
            if key in env.tpl:
                charge_str = env.tpl[key]
                crg = float(charge_str)
            else:
                crg = 0.0

            conf_history = "________"
            newline = "ATOM  %5d %4s %s %c%4d%c%03d%s%8.3f    %8.3f      %s%s\n" % \
                      (icount, atomname, resname, chainid, seqnum, icode, conf_number, xyz, rad, crg, conf_type, conf_history)
            mccelines.append(newline)
            icount += 1

        return mccelines


env = Env()

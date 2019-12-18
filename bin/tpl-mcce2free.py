#!/usr/bin/env python
"""Convert mcce format tpl to free format."""

# bug, single atom residue doesn't have CONNECT in mcce tpl, but should have a CONNECT record in free tpl

import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)-s: %(message)s')


class Paramfile:
    def __init__(self, fname):
        self.mccedb = {}  # parameter database in mcce format
        self.extra_records = ["EXTRA", "SCALING", "VDWAMBER", "RADCOVAL", "BOND_ANG", "RELAX",
                              "TORSION"]  # records that can be directly translated
        self.tplout = []
        self.load_file(fname)
        self.vdw_complete()
        return

    def load_file(self, fname):
        lines = open(fname).readlines()

        extralines = []
        comment_lines = []
        for line in lines:
            if line[0] == "#":
                comment_lines.append(line)

            end = line.find("#")
            line = line[:end]
            if len(line) < 20:
                continue
            key1 = line[:9].strip()
            key2 = line[9:15].strip()
            key3 = line[15:19]
            value = line[20:]

            # Direct translate: This part directly translate mcce to free format
            if key1 in self.extra_records:
                line = "%s, %s, %s: %s\n" % (key1, key2, key3, value)
                extralines.append(line)
            else:
                self.mccedb[(key1, key2, key3)] = value

        # Collect all conformer names from the read file
        conformers = []
        for k in self.mccedb.keys():
            if k[0] == "CONFLIST":
                conformers += self.mccedb[k].split()
        logging.debug("# Detected these conformers: [%s]" % ', '.join(map(str, conformers)))

        # check consistency between ATOMNAME and IATOM
        for conf in conformers:
            if self.atom_consistency(conf):  # pased
                logging.debug("# Consistency test passed for ATOM records of conformer %s." % conf)
            else:
                logging.debug("# There are discrepancies in ATOM records of conformer %s shown above." % conf)

        # Make conflist
        # Include old comment lines
        self.tplout.append(">>>START of original comments, this file was converted from old format\n")
        for line in comment_lines:
            self.tplout.append(line)
        self.tplout.append("<<<END of original comments\n")

        self.tplout.append("\n# Values of the same key are appended and separated by \",\"\n")
        residue_names = [x[:3] for x in conformers]
        residues = list(set(residue_names))
        for residue in residues:
            line = "CONFLIST, %s: " % residue
            conflist = []
            for conf in conformers:
                if conf[:3] == residue:
                    conflist.append(conf)
            line += ", ".join(conflist)
            line += "\n"
            self.tplout.append(line)

        # Make atom records
        self.tplout.append("\n# Atom definition\n")
        for conf in conformers:
            self.tplout += self.make_atom(conf)

        # Make charge records
        self.tplout.append("\n# Atom charges\n")
        for conf in conformers:
            self.tplout += self.make_charge(conf)

        # Make radius records
        self.tplout.append("\n# Atom radius, dielelctric boundary radius, VDW radius, and energy well depth\n")
        for conf in conformers:
            self.tplout += self.make_radius(conf)

        # Make conformer parameters
        self.tplout.append("\n# Conformer parameters that appear in head3.lst: ne, Em0, nH, pKa0, rxn\n")
        self.tplout += self.make_confparm(conformers)

        # Make rotatable bonds
        self.tplout.append("\n# Rotatable bonds. The atoms extended in the bond direction will all be rotated.\n")
        lines = []
        for key in self.mccedb.keys():
            if key[0] == "ROTAMER":
                residue = key[1]
                value = self.mccedb[key]
                atom1 = "\"%s\"" % value[:4]
                atom2 = "\"%s\"" % value[5:9]
                bond = "%s - %s" % (atom1, atom2)
                line = "ROTATE, %s: %s\n" % (residue, bond)
                lines.append(line)
        self.tplout += lines

        # Direct translate: This part directly translate mcce to free format
        lines = []
        for key in self.mccedb.keys():
            if key[0] in self.extra_records:
                line = "%s, %s: %s\n" % (key[0], key[1], self.mccedb[key])
                lines.append(line)
        self.tplout += lines

        self.tplout += extralines

        return()

    def atom_consistency(self, conf):
        passed = False
        natom = int(self.mccedb["NATOM", conf, "    "])
        if natom == 0:
            passed = True
        else:
            for i in range(natom):
                try:
                    key = ("ATOMNAME", conf, "%4d" % i)
                    atomname = "{:<4}".format(self.mccedb[key][:4])
                except:
                    logging.debug("# Error in fetching number %d atom. Check ATOMNAME record of conformer %s" % (i, conf))
                    return passed
                try:
                    key = ("IATOM", conf, atomname)
                    iatom = int(self.mccedb[key].strip())
                except:
                    logging.debug("# Error in finding index for atom \"%s\" of conformer %s" % (atomname, conf))
                    return passed
                if iatom == i:
                    passed = True

        return passed

    def make_atom(self, conf):
        natom = int(self.mccedb["NATOM", conf, "    "])
        lines = []
        for i in range(natom):
            key = ("ATOMNAME", conf, "%4d" % i)
            atomname = "{:<4}".format(self.mccedb[key][:4])
            key = ("CONNECT", conf, atomname)
            connect = self.mccedb[key].rstrip()
            orbital_type = connect[:9].strip()
            if len(connect[10:]) < 1:
                nconnected = 0
            else:
                nconnected = int(len(connect[10:])/10)+1
            connected_atoms = []
            for j in range(1, nconnected+1):
                if connect[j*10:j*10+5].strip() == "LIG":
                    catomname = '%4s' % ("{:<4}".format(connect[j*10+5: j*10+9]))
                else:
                    serial_str = connect[j*10:j*10+5]
                    serial = int(serial_str)
                    catomname = '%4s' % ("{:<4}".format(connect[j*10+5: j*10+9]))
                    if serial != 0:
                        catomname = " ?  "
                connected_atoms.append(catomname)
            quoted = ['"%s"' % x for x in connected_atoms]
            str_value = ", ".join(quoted).rstrip(",")
            line = "CONNECT, \"%s\", %s: %s, %s\n" % (atomname, conf, orbital_type, str_value)
            lines.append(line)

        return lines

    def make_charge(self, conf):
        natom = int(self.mccedb["NATOM", conf, "    "])
        lines = []
        for i in range(natom):
            key = ("ATOMNAME", conf, "%4d" % i)
            atomname = "{:<4}".format(self.mccedb[key][:4])
            key = ("CHARGE", conf, atomname)
            if key in self.mccedb:
                charge = float(self.mccedb[key])
            else:
                charge = 0.0
            line = "CHARGE, %s, \"%4s\": %6.3f\n" % (conf, atomname, charge)
            lines.append(line)

        return lines

    def make_radius(self, conf):
        natom = int(self.mccedb["NATOM", conf, "    "])
        lines = []
        for i in range(natom):
            key = ("ATOMNAME", conf, "%4d" % i)
            atomname = "{:<4}".format(self.mccedb[key][:4])
            key = ("RADIUS", conf[:3], atomname)
            if key in self.mccedb:
                radius = float(self.mccedb[key])
            else:
                radius = 0.0
            line = "RADIUS, %s, \"%4s\": %6.3f, to_be_filled, to_be_filled\n" % (conf, atomname, radius)
            lines.append(line)
            #lines = list(set(lines))

        return lines

    def make_confparm(self, conformers):
        lines = []
        for conf in conformers:
            if conf[-2:] == "BK" or conf[-2:] == "DM": continue
            key = ("PROTON", conf, "    ")
            nH = int(self.mccedb[key].strip())
            key = ("ELECTRON", conf, "    ")
            ne = int(self.mccedb[key].strip())
            key = ("PKA", conf, "    ")
            pKa0 = float(self.mccedb[key].strip())
            key = ("EM", conf, "    ")
            Em0 = float(self.mccedb[key].strip())
            key = ("RXN", conf, "    ")
            rxn = float(self.mccedb[key].strip())

            line = "CONFORMER, %s: Em0=%6.1f, pKa0=%6.2f, ne=%2d, nH=%2d, rxn=%7.3f\n" % (conf, Em0, pKa0, ne, nH, rxn)
            lines.append(line)
        return lines

    def vdw_complete(self):
        newlines = []
        for line in self.tplout:
            i = line.find("#")
            nline = line[:i]
            fields = nline.split(":")
            if len(fields) != 2:
                newlines.append(line)
                continue
            keys = fields[0].strip().split(",")
            values = fields[1].strip().split(",")
            record_name = keys[0].strip().upper()
            if record_name == "RADIUS":
                nline = self.complete_radius(keys, values)
                newlines.append(nline)
            else:
                newlines.append(line)

        self.tplout = newlines
        return

    def complete_radius(self, keys, values):
        vdw_parm = {"C": (2.000, 0.150),
                    "H": (1.000, 0.020),
                    "O": (1.600, 0.200),
                    "N": (1.750, 0.160),
                    "S": (2.000, 0.200),
                    "X": (2.000, 0.173)
                    }
        atomname = keys[2].strip().strip("\"")
        if atomname[0] == "H": # H atom
            element = "H"
        else:
            element = atomname[1]
        if not (element in vdw_parm.keys()):
            element = "X"
        vdw = vdw_parm[element]
        line = "RADIUS, %s, %s: %s, %7.3f, %7.3f\n" % (keys[1], keys[2], values[0], vdw[0], vdw[1])
        return line


if __name__ == "__main__":
    filename = sys.argv[1]
    file_param = Paramfile(filename)
    sys.stdout.writelines(file_param.tplout)

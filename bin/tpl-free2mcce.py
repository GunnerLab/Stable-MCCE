#!/usr/bin/env python
"""
This program converts ftpl files in param ftpl directory and extra.tpl into a single mcce tpl file.
Any missing parameters will be loaded from a default tpl file "default.ftpl"
"""

import os
import glob
import shutil
import sys
import logging



peptide = {" N  ": (-1, " C  "),
           " C  ": (1,  " N  ")}

disulfur = ["CYD01"]

def residues_in_db(freedb):
    residues = []
    for key in freedb.keys():
        if key[0] == "CONFLIST" and key[1] not in residues:
            residues.append(key[1])
    return residues


def create_connect(key, value):
    atom = key[1]
    conf = key[2]
    fields = value.strip().split(",")
    orbital = fields[0]
    connected = [x.strip().strip("\"") for x in fields[1:]]

    nvalue = []
    for x in connected:
        if len(x) != 4:
            print("ATOM \"%s\" is not 4 characters" % x)
            sys.exit()

        if x == " ?  ":
            if atom in peptide:
                offset = "%-3d" % peptide[atom][0]
                connected_atom = peptide[atom][1]
            else:
                offset = "LIG"
                connected_atom = " ?  "
        else:
            offset = "0  "
            connected_atom = x
        if conf in disulfur and x == " SG ":
            offset = "LIG"
            connected_atom = " SG "


        nvalue.append("%s  %s" % (offset, connected_atom))

    nvalue_str = " ".join(nvalue)
    line = "CONNECT  %5s %4s %-5s     %s\n" % (conf, atom, orbital, nvalue_str)

    return line

def convert(filename, epsilon):
    freedb = {}  # parameter database in free format


    lines = open(filename).readlines()

    natoms = {}
    iatom = []
    atomname = []
    i_counter = {}
    connect = []
    em0_records = []
    pka0_records = []
    ne_records = []
    nh_records = []
    rxn_records = []
    radius_records = []
    charge_records = []
    rotamer_records = []

    for line in lines:
        end = line.find("#")
        line = line[:end]
        fields = line.split(":")
        if len(fields) != 2:
            continue

        key_string = fields[0].strip()
        keys = key_string.split(",")
        key1 = keys[0].strip().strip("\"")
        if len(keys) > 1:
            key2 = keys[1].strip().strip("\"")
        else:
            key2 = ""
        if len(keys) > 2:
            key3 = keys[2].strip().strip("\"")
        else:
            key3 = ""

        value_string = fields[1].strip()

        # make CONNECT on the go so CONNECT keeps the same order
        if key1 == "CONNECT":
            atom = key2
            conf = key3

            # CONNECT
            nline = create_connect((key1, key2, key3), value_string)
            connect.append(nline)

            # IATOM
            if conf in natoms:
                natoms[conf] += 1
            else:
                natoms[conf] = 1
            nline = "IATOM    %5s %4s %-3d\n" % (conf, atom, natoms[conf]-1)
            iatom.append(nline)

            # ATOMNAME
            nline = "ATOMNAME %5s  %3d %s\n" % (conf, natoms[conf]-1, atom)
            atomname.append(nline)

        elif key1 == "CONFORMER":
            conf = key2
            fields = [x.strip() for x in value_string.split(",")]
            for field in fields:
                subfields = [x.strip() for x in field.split("=")]
                if subfields[0].upper() == "EM0":
                    em0 = float(subfields[1])
                elif subfields[0].upper() == "PKA0":
                    pka0 = float(subfields[1])
                elif subfields[0].upper() == "NE":
                    ne = int(subfields[1])
                elif subfields[0].upper() == "NH":
                    nh = int(subfields[1])
                elif subfields[0].upper() == "RXN%02d" % int(epsilon):
                    rxn = float(subfields[1])
            em0_records.append( "EM       %5s      %-.2f\n" % (conf, em0))
            pka0_records.append("PKA      %5s      %-.2f\n" % (conf, pka0))
            ne_records.append("ELECTRON %5s      %-2d\n" % (conf, ne))
            nh_records.append("PROTON   %5s      %-2d\n" % (conf, nh))
            rxn_records.append("RXN      %5s      %-.2f\n" % (conf, rxn))

        elif key1 == "RADIUS":
            conf = key2
            atom = key3
            fields = value_string.split(",")
            ele_radius = float(fields[0])

            record = "RADIUS   %3s   %4s %4.2f\n" % (conf[:3], atom, ele_radius)
            if record not in radius_records:
                radius_records.append(record)

        elif key1 == "CHARGE":
            conf = key2
            atom = key3
            charge = float(value_string)

            if abs(charge) > 0.001:
                record = "CHARGE   %5s %4s %6.2f\n" % (conf, atom, charge)
                charge_records.append(record)

        key = (key1, key2, key3)
        if key in freedb:
            freedb[key] = freedb[key] + " , " + value_string
        else:
            freedb[key] = value_string

    # check if this database contains only one residue
    residues = residues_in_db(freedb)
    if len(residues) > 1:
        print("Multiple residues detected: %s" % ",".join(residues))
        sys.exit()

    residue = residues[0]

    nlines = []

    # CONFLIST
    fields = freedb["CONFLIST", residue, ""].split(",")
    conformers = []
    for field in fields:
        conf = field.strip()
        if len(conf) >= 1:
            conformers.append("%5s" % conf)
    key = "CONFLIST %3s        " % (residue)
    value = " ".join(conformers)
    line = "%s%s\n" %(key, value)
    nlines.append(line)

    natom_records = []
    for conf in conformers:
        if conf in natoms:
            counter = natoms[conf]
        else:
            counter = 0
        natom_records.append("NATOM    %5s      %-3d\n" % (conf, counter))

    # make heavy atom connect table for rotamer definition
    connect_table = {}
    for key in freedb.keys():
        if key[0] != "CONNECT":
            continue
        if key[2][:3] != residue:
            continue
        atom = key[1]
        if atom[1] == "H" or atom[0] == "H":
            continue
        fields = freedb[key].strip(",").split(",")
        connected = [x.strip().strip("\"") for x in fields[1:]]
        connected = [x for x in connected if x[0] != "H" and x[1] != "H" and x[1] != "?"]
        connect_table[atom] = connected

    # make rotamer records
    key = ("ROTATE", residue, "")
    counter = 0
    if key in freedb:
        fields = [x.strip() for x in freedb[key].split(",")]
        for field in fields:

            subfields = field.split("-")
            atom0 = subfields[0].strip().strip("\"")
            atom1 = subfields[1].strip().strip("\"")
            chained = [atom0, atom1]
            for atom in chained:  # Caution: loop over a increasing list
                if atom == atom0:
                    continue
                for x in connect_table[atom]:
                    if x not in chained:
                        chained.append(x)
            line = "ROTAMER  %3s %3d    %4s-%s %s\n" % (residue, counter, atom0, atom1, " ".join(chained[2:]))
            rotamer_records.append(line)
            counter += 1


    # atom records
    nlines.append("\n")
    nlines += natom_records
    nlines.append("\n")
    nlines += iatom
    nlines.append("\n")
    nlines += atomname
    nlines.append("\n")
    nlines += connect

    # conformer
    nlines.append("\n")
    nlines += nh_records
    nlines += pka0_records
    nlines += ne_records
    nlines += em0_records
    nlines += rxn_records

    # radius
    nlines.append("\n")
    nlines += radius_records

    # charge
    nlines.append("\n")
    nlines += charge_records

    # rotamer
    nlines.append("\n")
    nlines += rotamer_records

    return nlines

if __name__ == "__main__":
        # get param path, epsilon
    lines = open("run.prm").readlines()
    epsilon = "4.0"
    mcce_home = ""

    for line in lines:
        fields = line.split()
        if len(fields) >= 2:
            key = fields[-1].strip()
            value = fields[0].strip()
            if key == "(EPSILON_PROT)":
                epsilon = value
            elif key == "(MCCE_HOME)":
                mcce_home = value

    ftpldir = mcce_home+"/param"

    targetdir = "./param"
    if not os.path.exists(targetdir):
        os.mkdir(targetdir)


    # copy 00always_needed.tpl
    shutil.copyfile(ftpldir+"/00always_needed.tpl", targetdir+"/00always_needed.tpl")

    # Load all ftpl files
    cwd = os.getcwd()
    os.chdir(ftpldir)

    files = glob.glob("*.ftpl")
    files.sort()
    logging.info("   Loading ftpl files from %s" % ftpldir)

    tpllines = []
    for fname in files:
        tpllines += convert(fname, epsilon)

    os.chdir(cwd)

    open(targetdir+"/mcce.tpl", "w").writelines(tpllines)

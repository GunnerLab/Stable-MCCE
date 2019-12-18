#!/usr/bin/python
"""Sub routines to input and write out MCCE PDB files."""

import math

# bond distance threshold
BONDDISTANCE = 1.65
CUTOFF2 = BONDDISTANCE*BONDDISTANCE

def ddvv(xyz1, xyz2):
    """Distance squared between two vectors."""
    dx=xyz1[0]-xyz2[0]
    dy=xyz1[1]-xyz2[1]
    dz=xyz1[2]-xyz2[2]
    return dx*dx+dy*dy+dz*dz

class Atom:
    def __init__(self):
        self.serial = 0
        self.name = "  X "
        self.altLoc = " "
        self.resName = "UNK"
        self.chainID = "A"
        self.resSeq = 0
        self.iCode = " "
        self.confNum = 0
        self.atomID = ""
        self.confID = ""
        self.resID = ""
        self.xyz = (0.0, 0.0, 0.0)
        self.connect12 = []
        self.connect13 = []
        self.connect14 = []
        return

    def loadline(self, line):
        self.serial = int(line[6:11])
        self.name = line[12:16]
        self.altLoc = line[16]
        self.resName = line[17:20]
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.confNum = int(line[27:30])
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.atomID = "%4s%3s%04d%c%03d" % (self.name, self.resName, self.resSeq, self.chainID, self.confNum)
        self.confID = "%3s%04d%c%03d" % (self.resName, self.resSeq, self.chainID, self.confNum)
        self.resID = "%3s%04d%c" % (self.resName, self.resSeq, self.chainID)
        return

    def writeline(self):
        line = "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f\n" % \
               (self.serial, self.name, self.altLoc, self.resName, self.chainID, \
                self.resSeq, self.iCode, self.xyz[0], self.xyz[1], \
                self.xyz[2])
        return line


class Conformer:
    def __init__(self):
        self.confID = ""
        self.resID = ""
        self.atom = []
        return

class Residue:
    def __init__(self):
        self.resID = ""
        self.conf = []
        return

class Protein:
    def __init__(self):
        self.residue = []
        return

    def loadpdb(self, fname):
        rawlines = open(fname).readlines()
        lines = [x.strip("\n") for x in rawlines if x[:6] == "ATOM  " or x[:6] == "HETATM"]
        for line in lines:
            atom = Atom()
            atom.loadline(line)
            # reverse search if this atom belongs to an existing conformer and residue
            found_conf = False
            found_res = False
            for i_res in range(len(self.residue) - 1, -1, -1):
                insert_res = i_res
                for i_conf in range(len(self.residue[i_res].conf) - 1, -1, -1):
                    if self.residue[i_res].conf[i_conf].confID == atom.confID:
                        found_conf = True
                        self.residue[i_res].conf[i_conf].atom.append(atom)
                        break
                if found_conf:
                    break
                elif self.residue[i_res].resID == atom.resID:  # not in conf but residue ID matches
                    conf = Conformer()
                    conf.confID = atom.confID
                    conf.resID = atom.resID
                    conf.atom.append(atom)
                    self.residue[i_res].conf.append(conf)
                    found_res = True

            if not found_conf and not found_res:  # new residue
                conf = Conformer()
                conf.confID = atom.confID
                conf.resID = atom.resID
                conf.atom.append(atom)
                res = Residue()
                res.resID = conf.resID
                res.conf.append(conf)
                self.residue.append(res)

        # Insert an empty conformer for cofactors that do not have backbone
        for res in self.residue:
            if res.conf[0].confID[-3:] != "000":
                conf = Conformer()
                conf.confID = "%s000" % res.resID
                conf.resID = res.resID
                res.conf.insert(0, conf)

        return

    def make_connect12(self):
        # make connect table
        for res in self.residue:
            for conf in res.conf:
                # with backbone
                for atom in conf.atom:
                    for atom2 in res.conf[0].atom:
                        if atom != atom2:
                            if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                if atom2 not in atom.connect12:
                                    atom.connect12.append(atom2)
                    # with self
                    for atom2 in conf.atom:
                        if atom != atom2:
                            if ddvv(atom.xyz, atom2.xyz) < CUTOFF2:
                                if atom2 not in atom.connect12:
                                    atom.connect12.append(atom2)
        return

    def print_connect12(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print atom.atomID
                    for atom2 in atom.connect12:
                        print "   -> %s" % atom2.atomID
        return

    def make_connect13(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    for atom2 in atom.connect12:
                        if atom2 == atom: continue
                        for atom3 in atom2.connect12:
                            if (atom3 != atom) and (atom3 not in atom.connect12) and (atom3 not in atom.connect13):
                                atom.connect13.append(atom3)
        return

    def make_connect14(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    for atom3 in atom.connect13:
                        if atom3 == atom: continue
                        for atom4 in atom3.connect12:
                            if (atom4 != atom) \
                                    and (atom4 not in atom.connect12) \
                                    and (atom4 not in atom.connect13) \
                                    and (atom4 not in atom.connect14):
                                atom.connect14.append(atom4)
        return

    def print_connect13(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print atom.atomID
                    for atom2 in atom.connect13:
                        print "   -X-> %s" % atom2.atomID
        return

    def print_connect14(self):
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    print atom.atomID
                    for atom2 in atom.connect14:
                        print "   -X-X-> %s" % atom2.atomID
        return

    def exportpdb(self, fname):
        lines = []
        for res in self.residue:
            for conf in res.conf:
                for atom in conf.atom:
                    lines.append(atom.writeline())
        open(fname, "w").writelines(lines)
        return

    def print_atom_structure(self):
        for res in self.residue:
            print "Residue %s" % res.resID
            for conf in res.conf:
                print "-->Conformer %s" % conf.confID
                for atom in conf.atom:
                    print "---->Atom %s" % atom.atomID
        return

if __name__ == "__main__":
    pdbfile = "step2_out.pdb"
    protein = Protein()
    protein.loadpdb(pdbfile)
    protein.make_connect12()
    protein.make_connect13()
    protein.make_connect14()

    #protein.print_connect14()
    #protein.exportpdb("a.pdb")
    #protein.print_atom_structure()

#!/usr/bin/env python

import argparse

DEFAULT_RAD = 1.8
probe_rad = 1.4

radius = {" H": 1.2,
          " C": 1.7,
          " N": 1.55,
          " O": 1.52,
          " F": 1.47,
          " P": 1.8,
          " S": 1.8,
          "CL": 1.75,
          "CU": 1.4,
          " B": 1.92,
          "AL": 1.84,
          "NA": 2.27,
          "MG": 1.73,
          "SI": 2.1,
          "CA": 2.31,
          " K": 2.75,
          "FE": 1.63,
          "ZN": 1.39,
          "BR": 1.85
}


class Atom:
    def __init__(self):
        self.icount = 0
        self.name = ""
        self.element = ""  # two-letter code for element name
        self.resname = ""
        self.chainid = ""
        self.seqnum = 0
        self.icode = ""
        self.xyz = ()
        self.resid = ()  # (resname, chainid, seqnum, icode)
        self.rad = DEFAULT_RAD
        return

    def loadline(self, line):
        self.icount = int(line[6:11])
        self.name = line[12:16]
        self.resname = line[17:20]
        self.chainid = line[21]
        self.seqnum = int(line[22:26])
        self.icode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

        if len(self.name.strip()) < 4:
            self.element = self.name[:2]
        else:
            self.element = " H"

        self.resid = (self.resname, self.chainid, self.seqnum, self.icode)
        if self.element in radius:
            self.rad = radius[self.element]

        return

class Protein:
    def __init__(self):
        self.atoms = []
        self.residues = []
        return

    def loadpdb(self, fname):
        lines = [x for x in open(fname).readlines() if x[:6] == "ATOM  " or x[:6] == "HETATM"]
        for line in lines:
            atom = Atom()
            atom.loadline(line)
        return

    def group_residues(self):
        residue_ids = []    # this is a index in right sequence
        residue = {}        # this is residue records indexed by resid
        for atom in self.atoms:
            if atom.resid in residue_ids:
            # if atom.resid in residue_ids:    # faster?
                residue[atom.resid].atoms.append(atom)
            else:


        return

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Strip off exposed cofactors like water and ions based on solvent accessible surface area."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-c", metavar="RES", nargs="*", help="Specify cofactor names to strip off, default is HOH.")
    parser.add_argument("-s", metavar="exposure", help="Fraction exposure threshold to be cut. Default is 0.05.", default=0.05, type=float)
    parser.add_argument("-f", metavar="inputfile", help="Input file name.", required=True)
    parser.add_argument("-o", metavar="outputfile", help="Output file name, default is inputfile name with extension .stripped.")
    args = parser.parse_args()

    prot = Protein()
    prot.loadpdb(args.f)
    prot.group_residues()
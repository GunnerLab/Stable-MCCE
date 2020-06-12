#!/usr/bin/env pypy3

import argparse
import sys
import time
import math
import os


DEFAULT_RAD = 1.8
probe_rad = 1.4
area_k = 4.0*math.pi

BOX_SIZE = 2.3 + probe_rad    # roughly = max atom radius + probe radius


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


def fibonacci_sphere(samples):

    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return points

n_points = 36    # 122
point_preset = fibonacci_sphere(n_points)


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
        self.rad_ext = DEFAULT_RAD + probe_rad
        self.crg = 0.0
        self.sas = 0.0
        self.ibox = ()
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
            self.rad_ext = self.rad + probe_rad
        return

    def print_me(self):
        line = "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%8.3f    %8.3f\n" % (self.icount,
                                                                                          self.name,
                                                                                          self.resname,
                                                                                          self.chainid,
                                                                                          self.seqnum,
                                                                                          self.icode,
                                                                                          self.xyz[0],
                                                                                          self.xyz[1],
                                                                                          self.xyz[2],
                                                                                          self.rad,
                                                                                          self.crg)
        return line


class Residue:
    def __init__(self, resid):
        self.resid = resid
        self.atoms = []
        self.sas = 0.0
        self.sas_fraction = 0.0
        self.check = False
        return

    def print_me(self):
        lines = []
        for atom in self.atoms:
            lines.append(atom.print_me())
        return lines

class Protein:
    def __init__(self):
        self.atoms = []
        self.residues = []
        self.lower_bound = [1E10, 1E10, 1E10]
        self.upper_bound = [-1E10, -1E10, -1E10]
        self.boxes = {}
        return

    def loadpdb(self, fname):
        lines = [x for x in open(fname).readlines() if x[:6] == "ATOM  " or x[:6] == "HETATM"]
        for line in lines:
            atom = Atom()
            atom.loadline(line)
            self.atoms.append(atom)
        return

    def group_residues(self, cofactors):
        cur_resid = ()
        new_res = None
        self.residues = []
        for atom in self.atoms:
            # ignore residues we don't care
            if atom.resname not in cofactors:
                continue
            if atom.resid == cur_resid:
                new_res.atoms.append(atom)
            else:
                if new_res:
                    self.residues.append(new_res)
                new_res = Residue(atom.resid)
                new_res.atoms.append(atom)
                cur_resid = atom.resid
        # last one
        self.residues.append(new_res)

        point_preset = fibonacci_sphere(122)
        for res in self.residues:
            res.max_exposed = 0.0
            for atom in new_res.atoms:
                res.max_exposed += atomacc_in_res(atom, res.atoms, point_preset)

        return

    def print_me(self):
        for res in self.residues:
            print(res.resid)
            sys.stdout.writelines(res.print_me())

    def writepdb(self, fname):
        lines = []
        for atom in self.atoms:
            line = atom.print_me()
            lines.append(line)
        open(fname, "w").writelines(lines)
        return

    def writesas(self, fname):
        lines = []
        for res in self.residues:
            line = "%s %8.3f %8.3f\n" % (str(res.resid), res.sas, res.sas_fraction)
            lines.append(line)
        open(fname, "w").writelines(lines)
        return

    def make_regions(self):
        # make a 10x10x10 boxes so that we only need to consider atoms within neighboring boxes when compute sas
        for atom in self.atoms:
            if self.lower_bound[0] > atom.xyz[0]:
                self.lower_bound[0] = atom.xyz[0]
            if self.lower_bound[1] > atom.xyz[1]:
                self.lower_bound[1] = atom.xyz[1]
            if self.lower_bound[2] > atom.xyz[2]:
                self.lower_bound[2] = atom.xyz[2]
            if self.upper_bound[0] < atom.xyz[0]:
                self.upper_bound[0] = atom.xyz[0]
            if self.upper_bound[1] < atom.xyz[1]:
                self.upper_bound[1] = atom.xyz[1]
            if self.upper_bound[2] < atom.xyz[2]:
                self.upper_bound[2] = atom.xyz[2]
        # This boundary encloses all atoms. There is no need to expand boundary.

        print("      Box dimension (%.2f, %.2f, %.2f) - (%.2f, %.2f, %.2f)" % (self.lower_bound[0],
                                                                             self.lower_bound[1],
                                                                             self.lower_bound[2],
                                                                             self.upper_bound[0],
                                                                             self.upper_bound[1],
                                                                             self.upper_bound[2]))
        ncells_x = math.ceil((self.upper_bound[0] - self.lower_bound[0])/BOX_SIZE)
        ncells_y = math.ceil((self.upper_bound[1] - self.lower_bound[1])/BOX_SIZE)
        ncells_z = math.ceil((self.upper_bound[2] - self.lower_bound[2])/BOX_SIZE)
        print("      Number of boxes = %d x %d x %d = %d" % (ncells_x, ncells_y, ncells_z, ncells_x*ncells_y*ncells_z), end=",")

        self.boxes = {}
        for atom in self.atoms:
            ix = math.floor((atom.xyz[0] - self.lower_bound[0])/BOX_SIZE)
            iy = math.floor((atom.xyz[1] - self.lower_bound[1])/BOX_SIZE)
            iz = math.floor((atom.xyz[2] - self.lower_bound[2])/BOX_SIZE)
            ibox = (ix, iy, iz)
            atom.ibox = ibox
            if ibox in self.boxes:
                self.boxes[ibox].append(atom)
            else:
                self.boxes[ibox] = [atom]

        return

    def atom_sas(self, point_preset):
        for res in self.residues:
            for atom in res.atoms:
                n_points = len(point_preset)
                counter = n_points
                for p_raw in point_preset:
                    point = (p_raw[0] * atom.rad_ext + atom.xyz[0],
                             p_raw[1] * atom.rad_ext + atom.xyz[1],
                             p_raw[2] * atom.rad_ext + atom.xyz[2])

                    ibox = (math.floor((point[0] - self.lower_bound[0])/BOX_SIZE),
                            math.floor((point[1] - self.lower_bound[1])/BOX_SIZE),
                            math.floor((point[2] - self.lower_bound[2])/BOX_SIZE))

                    #print("---", ibox)
                    buried = False
                    for ix in range(ibox[0] - 1, ibox[0] + 2):
                        if buried:
                            break
                        for iy in range(ibox[1] - 1, ibox[1] + 2):
                            if buried:
                                break
                            for iz in range(ibox[2] - 1, ibox[2] + 2):
                                if buried:
                                    break
                                neighbor_box = (ix, iy, iz)
                                #print(neighbor_box)
                                if neighbor_box in self.boxes:
                                    for atom2 in self.boxes[neighbor_box]:
                                        if atom2 != atom:
                                            dx = point[0] - atom2.xyz[0]
                                            dy = point[1] - atom2.xyz[1]
                                            dz = point[2] - atom2.xyz[2]
                                            dd = dx*dx + dy*dy + dz*dz
                                            if dd < atom2.rad_ext * atom2.rad_ext:
                                                buried = True
                                                counter -= 1
                                                break

                #print(counter)
                atom.sas = area_k * atom.rad_ext * atom.rad_ext * counter / n_points
                #print(atom.name, atom.sas)

    def res_sas(self):
        for res in self.residues:
            total_exposed = 0.0
            for atom in res.atoms:
                total_exposed += atom.sas
            res.sas = total_exposed
            res.sas_fraction = total_exposed / res.max_exposed

        return


    def bulk_remove(self, args):
        cofactors = args.c
        box_allcofactor = {}
        for box in self.boxes.keys():
            box_allcofactor[box] = True
            for atom in self.boxes[box]:
                if atom.resname not in cofactors:
                    box_allcofactor[box] = False
                    break

        for res in self.residues:          
            touch_protein = False
            for atom in res.atoms:
                if touch_protein:
                    break
                ibox = atom.ibox
                for ix in range(ibox[0] - 2, ibox[0] + 3):
                    if touch_protein:
                        break
                    for iy in range(ibox[1] - 2, ibox[1] + 3):
                        if touch_protein:
                            break
                        for iz in range(ibox[2] - 2, ibox[2] + 3):
                            if touch_protein:
                                break
                            neighbor_box = (ix, iy, iz)
                            if neighbor_box in self.boxes:
                                if not box_allcofactor[neighbor_box]:
                                    touch_protein = True
                                    break
            if not touch_protein:
                res.sas_fraction = 1.0

        print("\n      Delete bulk cofactors ...", end=" ")
        remove_atom = set()
        remove_res = set()
        cutoff = float(args.s)
        for res in prot.residues:
            if res.sas_fraction > cutoff:
                remove_atom.update(res.atoms)
                remove_res.add(res)
        prot.atoms = list(set(prot.atoms) - remove_atom)
        prot.residues = list(set(prot.residues) - remove_res)
        n_stripped = len(remove_res)
        print("%d cofactors deleted" % n_stripped)

        return


def atomacc_in_res(atom, all_atoms, point_preset):
    r_extended = probe_rad + atom.rad
    n_points = len(point_preset)
    counter = n_points
    for p_raw in point_preset:
        point = (p_raw[0] * r_extended + atom.xyz[0],
                 p_raw[1] * r_extended + atom.xyz[1],
                 p_raw[2] * r_extended + atom.xyz[2])

        for atom2 in all_atoms:
            if atom2 != atom:
                dx = point[0] - atom2.xyz[0]
                dy = point[1] - atom2.xyz[1]
                dz = point[2] - atom2.xyz[2]
                dd = dx * dx + dy * dy + dz * dz
                if dd < atom2.rad_ext * atom2.rad_ext:
                    counter -= 1
                    break
    sas = area_k * atom.rad_ext * atom.rad_ext * counter / n_points

    return sas


def strip_surface(prot, cutoff, point_preset, ncycle = 10):
    n_stripped = 1

    current_cycle = 1
    while n_stripped and current_cycle <= ncycle:
        timeD = time.time()

        print("      Total atoms: %d; processing cofactors: %d ..." % (len(prot.atoms), len(prot.residues)))
        prot.make_regions()
        timeC = time.time()
        print("takes %.3f seconds" % (timeC-timeD))

        print("      Compute atom sas ...", end=" ")
        prot.atom_sas(point_preset)
        timeD = time.time()
        print("takes %.3f seconds" % (timeD-timeC))

        print("      Compute residue sas ...", end=" ")
        prot.res_sas()
        timeC = time.time()
        print("      takes in %.3f seconds" % (timeC-timeD))

        print("      Delete surface cofactors cycle %d ..." % current_cycle, end=" ")
        remove_atom = set()
        remove_res = set()
        for res in prot.residues:
            if res.sas_fraction > cutoff:
                remove_atom.update(res.atoms)
                remove_res.add(res)
        prot.atoms = list(set(prot.atoms) - remove_atom)
        prot.residues = list(set(prot.residues) - remove_res)
        n_stripped = len(remove_res)
        print("%d\n" % n_stripped)
        current_cycle += 1
    return


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Strip off exposed cofactors like water and ions based on solvent accessible surface area."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-c", metavar="RES", nargs="+", default=["HOH"], help="Specify cofactor names to strip off, default is HOH.")
    parser.add_argument("-s", metavar="exposure", help="Fraction exposure threshold to be cut. Default is 0.05.", default=0.05, type=float)
    parser.add_argument("-f", metavar="inputfile", help="Input file name.", required=True)
    parser.add_argument("-o", metavar="outputfile", help="Output file name, default is inputfile name with extension .stripped.")
    args = parser.parse_args()

    prot = Protein()

    timeLast = time.time()
    timeStart = timeLast
    print("Start at time 0.000 ...")

    print("   Read in structure ...", end=" ")
    prot.loadpdb(args.f)
    timeNow = time.time()
    print("takes %.3f seconds" % (timeNow-timeLast))
    timeLast = timeNow

    print("   Group into residues ...", end=" ")
    prot.group_residues(args.c)
    timeNow = time.time()
    print("takes %.3f seconds\n" % (timeNow-timeLast))
    timeLast = timeNow


    print("   Bulk strip off surface %s ..." % " ".join(args.c))
    print("      Total atom: %d; processed cofactors: %d" % (len(prot.atoms), len(prot.residues)))
    prot.make_regions()
    prot.bulk_remove(args)
    timeNow = time.time()
    print("   Bulk strip off process takes %.3f seconds\n" % (timeNow-timeLast))
    timeLast = timeNow

    n_steps = [12, 48, 100]
    for n in n_steps:
        print("   Layered strip off surface %s at precision +-%.3f ..." % (" ".join(args.c), 1.0/n))
        point_preset = fibonacci_sphere(n)

        if n == n_steps[-1]:
            ncycles = 9999
            cutoff = float(args.s)
        else:
            ncycles = 12
            cutoff = float(args.s) + 1.0 / n

        strip_surface(prot, cutoff, point_preset, ncycle=ncycles)
        timeNow = time.time()
        print("   takes %.3f seconds\n" % (timeNow - timeLast))
        timeLast = timeNow


    if args.o:
        fname = args.o
    else:
        filename, extension = os.path.splitext(args.f)
        if extension == ".pdb":
            fname = "%s-stripped.pdb" % filename
        else:
            fname = "%s-stripped.pdb" % args.f
    prot.writepdb(fname)
    filename, extension = os.path.splitext(args.f)
    prot.writesas("%s.acc" % filename)
    print("Output structure file: %s; accessibility: %s" % (fname, "%s.acc" % filename))

    timeEnd = time.time()
    print("Total time: %.3f" % (timeEnd - timeStart))

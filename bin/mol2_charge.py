#!/usr/bin/env python

import argparse
import sys


class Atom:
    def __init__(self, name):
        self.name = name
        self.element = ""
        self.serial = 0
        self.direct_neighbors = []
        self.all_neighbors = []  # [[step 1 neighbors], [step 2 neighbors], ...]
        self.fingerprint = []  # [atom, {step 1 neighbors}, {step 2 neighbors, ...]
        self.charge = 0.0

    def print_me(self):
        print("atom:\"%s\", element:\"%s\", %3d connect:(%s)" % (
            self.name, self.element, self.serial, ",".join([x.name for x in self.direct_neighbors])))


def load_ftplatoms(ftpl, conf):
    atoms = []
    lines = [x for x in open(ftpl).readlines() if x[:7] == "CONNECT"]
    connect = {}
    atom_by_name = {}
    for line in lines:
        key, value = line.split(":")
        _, atom_name, confname = [x.strip().strip('"') for x in key.split(",")]
        if confname == conf:
            atom = Atom(atom_name)
            atoms.append(atom)
            fields = [x.strip().strip('"') for x in value.split(",")]
            connect[atom.name] = fields[1:]
            atom_by_name[atom.name] = atom

    # immediate connect, neighbors
    counter = 1
    for atom in atoms:
        if len(atom.name.strip()) == 4 and atom.name[0] == "H":
            atom.element = " H"
        else:
            atom.element = atom.name[:2]

        atom.serial = counter
        counter += 1

        atom.direct_neighbors = [atom_by_name[x] for x in connect[atom.name]]

    # extended connect
    extend_connect(atoms)

    return atoms


def load_mol2atoms(mol2):
    atoms = []
    atom_by_serial = {}

    start = False
    lines = open(mol2).readlines()
    for line in lines:
        if line.strip() == "@<TRIPOS>ATOM":
            start = True
            continue
        if line.strip() == "@<TRIPOS>BOND":
            break

        if start:
            fields = line.split()
            if len(fields) > 8:

                atom = Atom(fields[1])
                if len(atom.name.strip()) == 4 and atom.name[0] == "H":
                    atom.element = " H"
                else:
                    atom.element = " " + atom.name[0]
                atom.serial = int(fields[0])
                # print(atom.name, atom.element)

                atom.charge = float(fields[8])
                atoms.append(atom)
                atom_by_serial[atom.serial] = atom

    start = False
    for line in lines:
        if line.strip() == "@<TRIPOS>BOND":
            start = True
            continue
        if line.strip() == "@<TRIPOS>SUBSTRUCTURE":
            break

        if start:
            fields = line.split()
            ia = int(fields[1])
            ib = int(fields[2])
            atom_by_serial[ia].direct_neighbors.append(atom_by_serial[ib])
            atom_by_serial[ib].direct_neighbors.append(atom_by_serial[ia])

    # extended connect
    extend_connect(atoms)

    return atoms


def extend_connect(atoms):
    for atom in atoms:
        all_neighbors = [atom.direct_neighbors]
        pool = [atom] + atom.direct_neighbors
        i_neighbor_group = 0

        while i_neighbor_group < len(all_neighbors):
            next_group = []
            for atom1 in all_neighbors[i_neighbor_group]:
                for atom2 in atom1.direct_neighbors:
                    if atom2 not in pool:
                        next_group.append(atom2)
                        pool.append(atom2)
            if next_group:
                all_neighbors.append(next_group)

            i_neighbor_group += 1

        atom.all_neighbors = all_neighbors

        atom.fingerprint = [atom.element]
        for xgroup in all_neighbors:
            elements = [x.element for x in xgroup]
            elements.sort()
            atom.fingerprint.append(elements)

        # print(atom.fingerprint)
    return


def assign_charge(ftpl, mol2, conf, show=False):
    charge_records = []
    ftpl2mol2 = {}   # key is ftpl atom name, value is mol2 atom

    # convert conformer in ftpl file to an atom list
    ftpl_atoms = load_ftplatoms(ftpl, conf)

    # convert mol2 to an atom list
    mol2_atoms = load_mol2atoms(mol2)

    # match atoms to mol2 atom list,
    atom1_list = [x for x in ftpl_atoms]
    atom2_list = [x for x in mol2_atoms]
    for atom1 in atom1_list:
        for atom2 in atom2_list:
            if atom1.fingerprint == atom2.fingerprint:
                ftpl2mol2[atom1.name] = atom2
                atom2_list.remove(atom2)
                break

    # assign charge
    if ftpl2mol2:
        for atom in ftpl_atoms:
            atom.charge = ftpl2mol2[atom.name].charge
            charge_records.append("CHARGE, %s, \"%s\": %7.2f\n" % (conf, atom.name, atom.charge))

    if show:
        if ftpl2mol2:
            print("=================")
            print(" ftpl  ---  mol2")
            print("=================")
            for atom in ftpl_atoms:
                print("\"%4s\" --- \"%s\"" % (atom.name, ftpl2mol2[atom.name].name))
            print("=================")

    return charge_records


if __name__ == "__main__":
    # Get the command arguments
    helpmsg = "Apply mol2 charge to a conformer in ftpl template. The mol2 file doesn't need to have the matching " \
              "atom names, as long as the element names are correct."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-m", metavar="mol2", help="mol2 file with charge", required=True)
    parser.add_argument("-c", metavar="conf", help="Conformer ID in conformer list to be matched.", required=True)
    parser.add_argument("--show", default=False, help="Show atom name matching, default false.", action='store_true')
    parser.add_argument("ftpl", metavar="ftpl", nargs=1)
    args = parser.parse_args()

    mol2_file = args.m
    ftpl_file = args.ftpl[0]

    charge_records = assign_charge(ftpl=ftpl_file, mol2=mol2_file, conf=args.c, show=args.show)
    if charge_records:
        sys.stdout.writelines(charge_records)
    else:
        print("Atoms don't match, the best match found:")

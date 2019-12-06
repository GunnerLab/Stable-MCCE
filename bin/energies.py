#!/usr/bin/env python
"""
This program set up and run MCCE step 3 in multiple threads.

Usage examples:

1. Run step 3 with 6 threads
    energies.py -p 6 [-e mcce]

2. Run step 3 to recreate head3.lst
    energies.py -r [-e mcce]

3. Run partial conformers. conformer 1 to 100
    energies.py -c 1, 100 [-e mcce]

-p number: run with number of processes
-r: refresh opp files vdw and ele, and head3.lst
-c start, end: run conformer from start to end

Assuming the underlying c code will:
1. Write only opp files that are calculated (so that multi-threads won't overwrite each other). Dummies and
un-calculated conformers don't have an opp file.
2. When writing opp file, no previous run is loaded.
3. head3.lst is generated at the end of step 3.


Mechanism:
 * mcce is modified to run part of the conformers, as controlled by run.prm.
 * if an opp file exists, a conformer is considered to be already computed and therefore loaded.
 * no matter how many conformers are calculated, the averaging and vdw terms are always recalculated.
 * by modifying run.prm, we can run partial step 3.
 * a dummy run after all partial runs will update the final head3.lst and vdw terms.
"""

import sys, argparse, shutil
import logging
import threading
import subprocess


class Prot:
    def __init__(self):
        self.nconf = self.loadpdb("step2_out.pdb")
        self.runprm = open("run.prm").readlines()
        shutil.copy("run.prm", ".run.prm.bak")
        return

    def loadpdb(self, fname):
        lines = open(fname).readlines()
        conformers = []
        for line in lines:
            # residue name, type and chain/sequence/icode/confNumber
            if line[80:82] != "BK":
                uniqID = line[17:20] + line[80:82] + line[21:30]
                if uniqID not in conformers:
                    conformers.append(uniqID)

        return len(conformers)

    def runprm_set(self, key, value):
        for i in range(len(self.runprm)):
            line = self.runprm[i]
            fields = line.split()
            if len(fields) > 1 and fields[-1] == "(" + key + ")":
                self.runprm[i] = "%-8s %s\n" % (value, line.split(" ", 1)[1].strip())
                break

        return

    def runprm_write(self):
        open("run.prm", "w").writelines(self.runprm)
        return

def mcce_thread(mcce, name):
    logging.info("Thread %s: staring", name)
    subprocess.run(mcce)
    logging.info("Thread %s: finishing", name)

    return


if __name__ == "__main__":
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")
    helpmsg = "Run mcce step 3, energy calculations, with multiple threads."

    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-p", metavar="processes", default=1, help="run mcce with number of processes, default to 1", type=int)
    parser.add_argument("-r", default=False, help="refresh opp files and head3.lst without running delphi", action="store_true")
    parser.add_argument("-c", metavar=('start', 'end'), default=[1, 9999], nargs=2, help="starting and ending "
                                                                                         "conformer, default to 1 and 9999", type=int)
    parser.add_argument("-e", metavar="/path/to/mcce", default="mcce", help="mcce executable location, default to \"mcce\"")
    args = parser.parse_args()

    prot = Prot()


    if args.r:
        prot.runprm_set("PBE_START", "1")
        prot.runprm_set("PBE_END", "-1")
        prot.runprm_write()
        logging.info("Create a single thread to run this job.")
        x  = threading.Thread(target=mcce_thread, args=(args.e, 1))
        x.start()
    elif
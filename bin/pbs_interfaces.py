#!/usr/bin/env python
import logging
import subprocess
import math
import struct


class PBS_DELPHI:
    def __init__(self):
        # consider loading these parameters from a file
        self.exe = "delphi"   # This has to be made available by the execution environment
        self.radius_probe = 1.4
        self.grids_per_ang = 2.0
        self.epsilon_prot = 4.0  # default value, be overwritten by run_options
        self.epsilon_solv = 80.0
        self.ionrad = 2.0
        self.salt = 0.15
        self.grids_delphi = 65

        return
    

    def depth(self, bound):
        # determine delphi focusing depth
        x_min = x_max = bound.single_bnd_xyzrcp[0].x
        y_min = y_max = bound.single_bnd_xyzrcp[0].y
        z_min = z_max = bound.single_bnd_xyzrcp[0].z
        for p in bound.single_bnd_xyzrcp[1:]:
            if x_min > p.x: x_min = p.x
            if x_max < p.x: x_max = p.x
            if y_min > p.y: y_min = p.y
            if y_max < p.y: y_max = p.y
            if z_min > p.z: z_min = p.z
            if z_max < p.z: z_max = p.z

        dx = x_max - x_min
        dy = y_max - y_min
        dz = z_max - z_min
        dm = max(dx, dy, dz)
        dm += self.radius_probe * 2 + 3.4  # expand the largest dimension by the probe radius and safety

        scale = self.grids_per_ang/(self.grids_delphi/(2*dm))  # scale is a multiplier on grid_per_ang required to reach the target resolution
        if scale <= 1.0:
            depth = 1
        else:
            depth = math.ceil(math.log(scale)/math.log(2.0)) + 1

        return depth


    def write_fort15(self, xyzrcp):
        header = "ATOM      1  O   LYS A   1    "
        with open("fort.15", "w") as fh:
            for p in xyzrcp:
                fh.write("%-30s%8.3f%8.3f%8.3f\n" % (header, p.x, p.y, p.z))
        return

    def write_fort13(self, xyzrcp):
        struct_fmt = '=ifffffi'
        with open("fort.13", "wb") as fh:
            for p in xyzrcp:
                record_unf = struct.pack(struct_fmt, 20, p.x, p.y, p.z, p.r, p.c, 20)
                fh.write(record_unf)
        return


    def run(self, bound, run_options):
        """PBE solver interface for delphi. 
        It will generate site p in both boundary conditions 
        and return rxn in single boundary condition.
        """

        # snippets to check the input and environment
        # Current working directory
        # cwd = os.getcwd()
        # print(cwd)
        # What are in bound
        #print(vars(bound))

        depth = self.depth(bound)


        # single side chain boundary condition
        # fort.13
        # The first run starts with fort.13 as dielectric boundary, the following runs will be focusing runs, using the phi
        # as input
        #
        self.write_fort13(bound.single_bnd_xyzrcp)
        self.write_fort15(bound.single_bnd_xyzrcp)


        # fort.27
        center = [0.0, 0.0, 0.0]
        weight = 0.0
        for p in bound.single_bnd_xyzrcp:
            w = abs(p.c)
            if w > 0.00001:
                center[0] += p.x * w
                center[1] += p.y * w
                center[2] += p.z * w
                weight += w

        if weight > 0.000001:
            center = [c/(weight+0.000001) for c in center]
        else:
            logging.error("PB solver shouldn't run a conformer has no charged atom.")
        with open("fort.27", "w") as fh:
            fh.write("ATOM  %5d  C   CEN  %04d    %8.3f%8.3f%8.3f\n" % (1, 1, center[0], center[1], center[2]))

        # fort.10
        self.epsilon_prot = run_options.d
        with open("fort.10", "w") as fh:
            fh.write("gsize=%d\n" % self.grids_delphi)
            fh.write("scale=%.2f\n" % self.grids_per_ang)
            fh.write("in(unpdb,file=\"fort.13\")\n")
            fh.write("indi=%.1f\n" % self.epsilon_prot)
            fh.write("exdi=%.1f\n" % self.epsilon_solv)
            fh.write("ionrad=%.1f\n" % self.ionrad)
            fh.write("salt=%.2f\n" % self.salt)
            fh.write("bndcon=2\n")
            fh.write("center(777, 777, 0)\n")
            fh.write("out(frc,file=\"run01.frc\")\n")
            fh.write("out(phi,file=\"run01.phi\")\n")
            fh.write("site(a,c,p)\n")
            fh.write("energy(g,an,sol)\n")   # g for grid energy, sol for corrected rxn

        # 1st delphi run
        result = subprocess.run([self.exe], capture_output=True)
        # print("STDOUT ===========")
        # print(result.stdout)
        delphi_log = [result.stdout]   # delphi logs of all depths will be stored in this list

        if result.stderr:
            print("delphi encountered error ===========")
            print(result.stderr)
            print("====================================")

        # subsequent delphi runs
        for i in range(1, depth):
            with open("fort.10", "w") as fh:
                fh.write("gsize=%d\n" % self.grids_delphi)
                fh.write("scale=%.2f\n" % self.grids_per_ang)
                fh.write("in(unpdb,file=\"fort.13\")\n")
                fh.write("in(phi,file=\"run%02d.phi\")\n" % i)
                fh.write("indi=%.1f\n" % self.epsilon_prot)
                fh.write("exdi=%.1f\n" % self.epsilon_solv)
                fh.write("ionrad=%.1f\n" % self.ionrad)
                fh.write("salt=%.2f\n" % self.salt)
                fh.write("bndcon=3\n")
                fh.write("center(777, 777, 0)\n")
                fh.write("out(frc,file=\"run%02d.frc\")\n" % (i+1))
                fh.write("out(phi,file=\"run%02d.phi\")\n" % (i+1))
                fh.write("site(a,c,p)\n")
                fh.write("energy(g,an,sol)\n")  # g for grid energy, sol for corrected rxn

            result = subprocess.run([self.exe], capture_output=True)
            delphi_log.append(result.stdout)


        # collect results from the log
        try:
            lines = open("run01.frc", "r").readlines()
        except OSError:
            logging.error("Could not open Delphi output file run01.frc.")






        # multi side chain boundary condition

        rxn = 0


        return rxn
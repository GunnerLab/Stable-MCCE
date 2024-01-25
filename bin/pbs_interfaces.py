#!/usr/bin/env python

import os
import math
import struct

class PBS_DELPHI:
    def __init__(self):
        # consider loading these parameters from a file
        self.exe = "delphi"   # This has to be made available by the execution environment
        self.radius_probe = 1.4
        self.grid_per_ang = 2.0
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

        scale = self.grid_per_ang/(self.grids_delphi/(2*dm))  # scale is a multiplier on grid_per_ang required to reach the target resolution
        if scale <= 1.0:
            depth = 1
        else:
            depth = math.ceil(math.log(scale)/math.log(2.0)) + 1

        return depth


    def run(self, bound):
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
        
        struct_fmt = '=ifffffi'
        with open("fort.13", "wb") as fh:
            for p in bound.single_bnd_xyzrcp:
                record_unf = struct.pack(struct_fmt, 20, p.x, p.y, p.z, p.r, p.c, 20)
                fh.write(record_unf)
        




        # multi side chain boundary condition

        rxn = 0


        return rxn
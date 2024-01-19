#!/usr/bin/env python

import os

def pbs_delphi(bound):
    """PBE solver interface for delphi. 
    It will generate site p in both boundary conditions 
    and return rxn in single boundary condition.
    """
    exe = "delphi"   # This has to be made available by the execution environment

    # snippets to check the input and environment
    # Current working directory
    # cwd = os.getcwd()
    # print(cwd)
    # What are in bound
    print(vars(bound))

    # determine delphi focusing depth
    x_min = x_max = bound.single_bnd_xyzrcp[0].x


    # single side chain boundary condition
    

    # multi side chain boundary condition

    rxn = 0


    return rxn
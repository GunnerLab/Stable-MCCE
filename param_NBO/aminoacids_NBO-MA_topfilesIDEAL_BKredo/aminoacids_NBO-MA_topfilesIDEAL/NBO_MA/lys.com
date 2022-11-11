%nprocshared=24
%mem=6GB
%chk=lys.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

1 1
 N(PDBName=N,ResName=LYS,ResNum=1_A)                  -2.82400000   -0.09600000    1.10100000
 C(PDBName=CA,ResName=LYS,ResNum=1_A)                 -1.53500000    0.03600000    0.43300000
 C(PDBName=C,ResName=LYS,ResNum=1_A)                  -1.50100000    1.15900000   -0.59600000
 O(PDBName=O,ResName=LYS,ResNum=1_A)                  -0.57200000    1.94000000   -0.71400000
 C(PDBName=CB,ResName=LYS,ResNum=1_A)                 -1.15500000   -1.31900000   -0.20600000
 C(PDBName=CG,ResName=LYS,ResNum=1_A)                  0.23500000   -1.42100000   -0.86200000
 C(PDBName=CD,ResName=LYS,ResNum=1_A)                  1.43900000   -1.11400000    0.07400000
 C(PDBName=CE,ResName=LYS,ResNum=1_A)                  2.05400000    0.24300000   -0.25500000
 N(PDBName=NZ,ResName=LYS,ResNum=1_A)                  3.15600000    0.60600000    0.74600000
 H(PDBName=H1,ResName=LYS,ResNum=1_A)                 -3.56900000   -0.39500000    0.47500000
 H(PDBName=HA,ResName=LYS,ResNum=1_A)                 -0.79600000    0.29500000    1.20400000
 H(PDBName=HB1,ResName=LYS,ResNum=1_A)                -1.90700000   -1.57700000   -0.96400000
 H(PDBName=HG1,ResName=LYS,ResNum=1_A)                 0.28600000   -0.77500000   -1.74900000
 H(PDBName=HD1,ResName=LYS,ResNum=1_A)                 2.20200000   -1.89600000   -0.03300000
 H(PDBName=HE1,ResName=LYS,ResNum=1_A)                 1.32400000    1.05600000   -0.21600000
 H(PDBName=HB2,ResName=LYS,ResNum=1_A)                -1.25100000   -2.07000000    0.58800000
 H(PDBName=HG2,ResName=LYS,ResNum=1_A)                 0.34100000   -2.44300000   -1.24000000
 H(PDBName=HD2,ResName=LYS,ResNum=1_A)                 1.11000000   -1.13600000    1.12300000
 H(PDBName=HE2,ResName=LYS,ResNum=1_A)                 2.53600000    0.24700000   -1.23700000
 H(PDBName=H,ResName=LYS,ResNum=1_A)                  -2.39700000    1.25100000   -1.24500000
 H(PDBName=H2,ResName=LYS,ResNum=1_A)                 -3.11700000    0.75800000    1.57000000
 H                                                     3.96195566    0.03942042    0.57446854
 H                                                     2.82726043    0.44496076    1.67658942
 H                                                     3.39792889    1.57050596    0.64017628

 1 2 1.0 10 1.0 21 1.0
 2 3 1.0 5 1.0 11 1.0
 3 20 1.0 4 2.0
 4
 5 12 1.0 6 1.0 16 1.0
 6 13 1.0 17 1.0 7 1.0
 7 8 1.0 14 1.0 18 1.0
 8 19 1.0 15 1.0 9 1.0
 9 22 1.0 23 1.0 24 1.0
 10
 11
 12
 13
 14
 15
 16
 17
 18
 19
 20
 21
 22
 23
 24









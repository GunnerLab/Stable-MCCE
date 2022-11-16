%nprocshared=24
%mem=6GB
%chk=cysteine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=CYS,ResNum=1)                     0.10200000   -0.14300000    0.00000000
 C(PDBName=CA,ResName=CYS,ResNum=1)                    1.55800000   -0.02900000    0.00000000
 C(PDBName=C,ResName=CYS,ResNum=1)                     1.98600000    1.41900000    0.00000000
 O(PDBName=O,ResName=CYS,ResNum=1)                     1.11924313    2.33130647    0.00000000
 C(PDBName=CB,ResName=CYS,ResNum=1)                    2.11100000   -0.80700000   -1.21000000
 S(PDBName=SG,ResName=CYS,ResNum=1)                    3.91600000   -0.85600000   -1.16800000
 H(PDBName=H1,ResName=CYS,ResNum=1)                   -0.54000000    0.72500000    0.00000000
 H(PDBName=H2,ResName=CYS,ResNum=1)                   -0.39700000   -1.10100000    0.00000000
 H(PDBName=HA,ResName=CYS,ResNum=1)                    1.94200000   -0.48500000    0.93100000
 H(PDBName=HB2,ResName=CYS,ResNum=1)                   1.74300000   -1.84900000   -1.21800000
 H(PDBName=HB3,ResName=CYS,ResNum=1)                   1.77800000   -0.35700000   -2.16800000
 H(PDBName=HG,ResName=CYS,ResNum=1)                    4.12300000    0.03600000   -2.13400000
 H                                                     3.02629016    1.66939244    0.00000000

 1 2 1.0 7 1.0 8 1.0
 2 3 1.0 5 1.0 9 1.0
 3 4 2.0 13 1.0
 4
 5 6 1.0 10 1.0 11 1.0
 6 12 1.0
 7
 8
 9
 10
 11
 12
 13




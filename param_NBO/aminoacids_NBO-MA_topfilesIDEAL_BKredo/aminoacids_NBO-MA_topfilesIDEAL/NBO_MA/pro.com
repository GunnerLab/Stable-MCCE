%nprocshared=24
%mem=6GB
%chk=proline.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=PRO,ResNum=1)                     0.11200000    0.68400000    0.12400000
 C(PDBName=CA,ResName=PRO,ResNum=1)                    1.62200000    0.53000000   -0.03900000
 C(PDBName=C,ResName=PRO,ResNum=1)                     2.25500000    1.52800000   -1.00500000
 O(PDBName=O,ResName=PRO,ResNum=1)                     1.53723603    2.38210745   -1.58713908
 C(PDBName=CB,ResName=PRO,ResNum=1)                    1.85000000   -0.86900000   -0.59300000
 C(PDBName=CG,ResName=PRO,ResNum=1)                    0.59800000   -1.03900000   -1.47000000
 C(PDBName=CD,ResName=PRO,ResNum=1)                   -0.52300000   -0.42900000   -0.62200000
 H(PDBName=H,ResName=PRO,ResNum=1)                    -0.37500000    1.42300000    0.65200000
 H(PDBName=HA,ResName=PRO,ResNum=1)                    2.11200000    0.66600000    0.91400000
 H(PDBName=HB2,ResName=PRO,ResNum=1)                   2.73600000   -0.86100000   -1.25000000
 H(PDBName=HB3,ResName=PRO,ResNum=1)                   1.81800000   -1.59500000    0.22900000
 H(PDBName=HG2,ResName=PRO,ResNum=1)                   0.56100000   -0.36200000   -2.33900000
 H(PDBName=HG3,ResName=PRO,ResNum=1)                   0.48600000   -2.11700000   -1.64600000
 H(PDBName=HD2,ResName=PRO,ResNum=1)                  -1.31400000   -0.03000000   -1.27900000
 H(PDBName=HD3,ResName=PRO,ResNum=1)                  -0.88400000   -1.18000000    0.09800000
 H                                                     3.30903962    1.50136534   -1.18718417

 1 2 1.0 7 1.0 8 1.0
 2 3 1.0 5 1.0 9 1.0
 3 4 2.0 16 1.0
 4
 5 6 1.0 10 1.0 11 1.0
 6 7 1.0 12 1.0 13 1.0
 7 14 1.0 15 1.0
 8
 9
 10
 11
 12
 13
 14
 15
 16




%nprocshared=24
%mem=6GB
%chk=his.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=HIS,ResNum=1)                     0.10200000   -0.14300000    0.00000000
 C(PDBName=CA,ResName=HIS,ResNum=1)                    1.55800000   -0.02900000    0.00000000
 C(PDBName=C,ResName=HIS,ResNum=1)                     1.98600000    1.41900000    0.00000000
 O(PDBName=O,ResName=HIS,ResNum=1)                     1.11924313    2.33130647    0.00000000
 C(PDBName=CB,ResName=HIS,ResNum=1)                    2.04700000   -0.73200000   -1.27700000
 C(PDBName=CG,ResName=HIS,ResNum=1)                    2.95200000   -1.89200000   -0.97800000
 N(PDBName=ND1,ResName=HIS,ResNum=1)                   3.33500000   -2.32000000    0.29100000
 C(PDBName=CD2,ResName=HIS,ResNum=1)                   3.51700000   -2.67800000   -1.97200000
 C(PDBName=CE1,ResName=HIS,ResNum=1)                   4.11800000   -3.36300000   -0.04600000
 N(PDBName=NE2,ResName=HIS,ResNum=1)                   4.28200000   -3.64100000   -1.36800000
 H(PDBName=H1,ResName=HIS,ResNum=1)                   -0.54000000    0.72500000    0.00000000
 H(PDBName=H2,ResName=HIS,ResNum=1)                   -0.39700000   -1.10100000    0.00000000
 H(PDBName=HA,ResName=HIS,ResNum=1)                    1.92800000   -0.46800000    0.94100000
 H(PDBName=HB2,ResName=HIS,ResNum=1)                   1.20600000   -1.11300000   -1.89300000
 H(PDBName=HB3,ResName=HIS,ResNum=1)                   2.59800000   -0.03100000   -1.93500000
 H(PDBName=HD2,ResName=HIS,ResNum=1)                   3.36200000   -2.53000000   -3.03700000
 H(PDBName=HE2,ResName=HIS,ResNum=1)                   4.83700000   -4.38800000   -1.81000000
 H(PDBName=HE1,ResName=HIS,ResNum=1)                   4.57900000   -3.92300000    0.75800000
 H                                                     3.02629016    1.66939244    0.00000000

 1 2 1.0 11 1.0 12 1.0
 2 3 1.0 5 1.0 13 1.0
 3 4 2.0 19 1.0
 4
 5 6 1.0 14 1.0 15 1.0
 6 7 1.0 8 1.5
 7 9 1.5
 8 10 1.5 16 1.0
 9 10 1.5 18 1.0
 10 17 1.0
 11
 12
 13
 14
 15
 16
 17
 18
 19




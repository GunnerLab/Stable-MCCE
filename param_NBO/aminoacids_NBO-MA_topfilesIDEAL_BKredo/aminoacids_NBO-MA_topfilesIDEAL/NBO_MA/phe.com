%nprocshared=24
%mem=6GB
%chk=Phe.chk
#p b3lyp/6-31+g(d,p) geom=connectivity polar

Title Card Required

0 1
 N(PDBName=N,ResName=PHE,ResNum=1)                     0.11000000   -0.20000000    0.00000000
 C(PDBName=CA,ResName=PHE,ResNum=1)                    1.56100000   -0.03800000    0.00000000
 C(PDBName=C,ResName=PHE,ResNum=1)                     1.94300000    1.42300000    0.00000000
 O(PDBName=O,ResName=PHE,ResNum=1)                     1.04780103    2.30741470    0.00000000
 C(PDBName=CB,ResName=PHE,ResNum=1)                    2.17500000   -0.71900000   -1.26100000
 C(PDBName=CG,ResName=PHE,ResNum=1)                    1.89000000   -2.21800000   -1.43600000
 C(PDBName=CD1,ResName=PHE,ResNum=1)                   0.73800000   -2.62900000   -2.11700000
 C(PDBName=CD2,ResName=PHE,ResNum=1)                   2.73300000   -3.17800000   -0.87000000
 C(PDBName=CE1,ResName=PHE,ResNum=1)                   0.42400000   -3.98000000   -2.21400000
 C(PDBName=CE2,ResName=PHE,ResNum=1)                   2.42400000   -4.53100000   -0.97700000
 C(PDBName=CZ,ResName=PHE,ResNum=1)                    1.26800000   -4.93200000   -1.64600000
 H(PDBName=H1,ResName=PHE,ResNum=1)                   -0.56100000    0.64700000    0.00000000
 H(PDBName=H2,ResName=PHE,ResNum=1)                   -0.35800000   -1.17300000    0.00000000
 H(PDBName=HA,ResName=PHE,ResNum=1)                    1.97300000   -0.50000000    0.91700000
 H(PDBName=HB2,ResName=PHE,ResNum=1)                   1.84000000   -0.18200000   -2.17300000
 H(PDBName=HB3,ResName=PHE,ResNum=1)                   3.27200000   -0.56900000   -1.26200000
 H(PDBName=HD1,ResName=PHE,ResNum=1)                   0.06100000   -1.89400000   -2.53300000
 H(PDBName=HD2,ResName=PHE,ResNum=1)                   3.61900000   -2.87300000   -0.33100000
 H(PDBName=HE1,ResName=PHE,ResNum=1)                  -0.48300000   -4.28500000   -2.71500000
 H(PDBName=HE2,ResName=PHE,ResNum=1)                   3.07500000   -5.26900000   -0.53100000
 H(PDBName=HZ,ResName=PHE,ResNum=1)                    1.02000000   -5.98000000   -1.71800000
 H                                                     2.97484348    1.70619432    0.00000000

 1 2 1.0 12 1.0 13 1.0
 2 3 1.0 5 1.0 14 1.0
 3 4 2.0 22 1.0
 4
 5 6 1.0 15 1.0 16 1.0
 6 7 1.5 8 1.5
 7 9 1.5 17 1.0
 8 10 1.5 18 1.0
 9 11 1.5 19 1.0
 10 11 1.5 20 1.0
 11 21 1.0
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




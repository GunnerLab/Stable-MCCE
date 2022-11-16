%nprocshared=23
%mem=10GB
%chk=thr.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=THR,ResNum=1)                     0.08000000   -0.08300000    0.00000000
 C(PDBName=CA,ResName=THR,ResNum=1)                    1.53900000   -0.03300000    0.00000000
 C(PDBName=C,ResName=THR,ResNum=1)                     2.03200000    1.39400000    0.00000000
 O(PDBName=O,ResName=THR,ResNum=1)                     1.20739388    2.34457630    0.00000000
 C(PDBName=CB,ResName=THR,ResNum=1)                    2.11600000   -0.81800000   -1.22600000
 O(PDBName=OG1,ResName=THR,ResNum=1)                   1.74300000   -2.18800000   -1.16200000
 C(PDBName=CG2,ResName=THR,ResNum=1)                   3.65400000   -0.84200000   -1.36100000
 H(PDBName=H1,ResName=THR,ResNum=1)                   -0.52300000    0.81300000    0.00000000
 H(PDBName=H2,ResName=THR,ResNum=1)                   -0.46100000   -1.01800000    0.00000000
 H(PDBName=HA,ResName=THR,ResNum=1)                    1.90400000   -0.50800000    0.92900000
 H(PDBName=HB,ResName=THR,ResNum=1)                    1.69100000   -0.37800000   -2.15400000
 H(PDBName=HG21,ResName=THR,ResNum=1)                  3.98600000   -1.43300000   -2.23500000
 H(PDBName=HG22,ResName=THR,ResNum=1)                  4.08300000    0.16800000   -1.50200000
 H(PDBName=HG23,ResName=THR,ResNum=1)                  4.13800000   -1.27900000   -0.46700000
 H(PDBName=HG1,ResName=THR,ResNum=1)                   2.06600000   -2.59200000   -1.97300000
 H                                                     3.08255080    1.59708377    0.00000000

 1 2 1.0 8 1.0 9 1.0
 2 3 1.0 5 1.0 10 1.0
 3 4 2.0 16 1.0
 4
 5 6 1.0 7 1.0 11 1.0
 6 15 1.0
 7 12 1.0 13 1.0 14 1.0
 8
 9
 10
 11
 12
 13
 14
 15
 16




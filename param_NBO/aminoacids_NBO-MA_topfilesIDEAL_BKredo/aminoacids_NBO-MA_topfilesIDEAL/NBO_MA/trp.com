%nprocshared=24
%mem=6GB
%chk=tryptophan.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=TRP,ResNum=1)                     0.38300000    8.10900000    0.19000000
 C(PDBName=CA,ResName=TRP,ResNum=1)                    1.83900000    8.22000000    0.18800000
 C(PDBName=C,ResName=TRP,ResNum=1)                     2.27200000    9.66600000    0.15700000
 O(PDBName=O,ResName=TRP,ResNum=1)                     1.40843154   10.58111007    0.13715975
 C(PDBName=CB,ResName=TRP,ResNum=1)                    2.40300000    7.47500000   -1.05400000
 C(PDBName=CG,ResName=TRP,ResNum=1)                    2.19500000    5.95700000   -1.04100000
 C(PDBName=CD1,ResName=TRP,ResNum=1)                   1.14800000    5.26600000   -1.68600000
 C(PDBName=CD2,ResName=TRP,ResNum=1)                   2.90600000    5.00200000   -0.34600000
 N(PDBName=NE1,ResName=TRP,ResNum=1)                   1.19300000    3.88500000   -1.41300000
 C(PDBName=CE2,ResName=TRP,ResNum=1)                   2.29100000    3.74700000   -0.57800000
 C(PDBName=CE3,ResName=TRP,ResNum=1)                   4.03700000    5.10600000    0.50500000
 C(PDBName=CZ2,ResName=TRP,ResNum=1)                   2.81000000    2.58300000    0.03000000
 C(PDBName=CZ3,ResName=TRP,ResNum=1)                   4.53600000    3.94200000    1.09000000
 C(PDBName=CH2,ResName=TRP,ResNum=1)                   3.93300000    2.69800000    0.85500000
 H(PDBName=H1,ResName=TRP,ResNum=1)                   -0.21000000    8.91400000    0.17300000
 H(PDBName=H2,ResName=TRP,ResNum=1)                   -0.08100000    7.22300000    0.20900000
 H(PDBName=HA,ResName=TRP,ResNum=1)                    2.23800000    7.77100000    1.11700000
 H(PDBName=HB2,ResName=TRP,ResNum=1)                   1.96900000    7.89800000   -1.98200000
 H(PDBName=HB3,ResName=TRP,ResNum=1)                   3.49000000    7.66200000   -1.15100000
 H(PDBName=HD1,ResName=TRP,ResNum=1)                   0.38100000    5.74900000   -2.27700000
 H(PDBName=HE1,ResName=TRP,ResNum=1)                   0.53800000    3.15600000   -1.71600000
 H(PDBName=HE3,ResName=TRP,ResNum=1)                   4.50300000    6.06200000    0.69800000
 H(PDBName=HZ2,ResName=TRP,ResNum=1)                   2.34400000    1.62500000   -0.14200000
 H(PDBName=HZ3,ResName=TRP,ResNum=1)                   5.40100000    4.00200000    1.73400000
 H(PDBName=HH2,ResName=TRP,ResNum=1)                   4.34400000    1.81300000    1.31800000
 H                                                     3.31315745    9.91270862    0.15189946

 1 2 1.0 15 1.0 16 1.0
 2 3 1.0 5 1.0 17 1.0
 3 4 2.0 26 1.0
 4
 5 6 1.0 18 1.0 19 1.0
 6 7 1.5 8 2.0
 7 9 1.0 20 1.0
 8 10 1.5 11 1.5
 9 10 1.0 21 1.0
 10 12 1.5
 11 13 1.5 22 1.0
 12 14 1.5 23 1.0
 13 14 1.5 24 1.0
 14 25 1.0
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
 25
 26




%nprocshared=24
%mem=6GB
%chk=isoleucine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=ILE,ResNum=1)                     0.11600000   -0.21600000    0.00000000
 C(PDBName=CA,ResName=ILE,ResNum=1)                    1.56500000   -0.03600000    0.00000000
 C(PDBName=C,ResName=ILE,ResNum=1)                     1.92800000    1.43000000    0.00000000
 O(PDBName=O,ResName=ILE,ResNum=1)                     1.02137160    2.30269439    0.00000000
 C(PDBName=CB,ResName=ILE,ResNum=1)                    2.24200000   -0.81000000   -1.20300000
 C(PDBName=CG1,ResName=ILE,ResNum=1)                   1.85000000   -2.31900000   -1.31000000
 C(PDBName=CG2,ResName=ILE,ResNum=1)                   3.79900000   -0.71500000   -1.18500000
 C(PDBName=CD1,ResName=ILE,ResNum=1)                   2.22300000   -3.02000000   -2.63300000
 H(PDBName=H1,ResName=ILE,ResNum=1)                   -0.33900000   -1.19600000    0.00000000
 H(PDBName=H2,ResName=ILE,ResNum=1)                   -0.56500000    0.62200000    0.00000000
 H(PDBName=HA,ResName=ILE,ResNum=1)                    1.95400000   -0.45500000    0.94600000
 H(PDBName=HB,ResName=ILE,ResNum=1)                    1.88700000   -0.32100000   -2.13500000
 H(PDBName=HG12,ResName=ILE,ResNum=1)                  2.25900000   -2.88400000   -0.45000000
 H(PDBName=HG13,ResName=ILE,ResNum=1)                  0.75200000   -2.42300000   -1.21200000
 H(PDBName=HG21,ResName=ILE,ResNum=1)                  4.26100000   -1.20900000   -2.05800000
 H(PDBName=HG22,ResName=ILE,ResNum=1)                  4.17200000    0.32400000   -1.22900000
 H(PDBName=HG23,ResName=ILE,ResNum=1)                  4.23400000   -1.18200000   -0.28100000
 H(PDBName=HD11,ResName=ILE,ResNum=1)                  1.84600000   -4.05900000   -2.65500000
 H(PDBName=HD12,ResName=ILE,ResNum=1)                  1.79100000   -2.49900000   -3.50900000
 H(PDBName=HD13,ResName=ILE,ResNum=1)                  3.31500000   -3.07700000   -2.78800000
 H                                                     2.95607215    1.72659342    0.00000000

 1 2 1.0 9 1.0 10 1.0
 2 3 1.0 5 1.0 11 1.0
 3 4 2.0 21 1.0
 4
 5 6 1.0 7 1.0 12 1.0
 6 8 1.0 13 1.0 14 1.0
 7 15 1.0 16 1.0 17 1.0
 8 18 1.0 19 1.0 20 1.0
 9
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




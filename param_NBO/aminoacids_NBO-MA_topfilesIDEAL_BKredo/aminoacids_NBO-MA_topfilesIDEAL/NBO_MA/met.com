%nprocshared=24
%mem=6GB
%chk=methionine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=MET,ResNum=1)                     0.07100000   -0.21500000    0.00000000
 C(PDBName=CA,ResName=MET,ResNum=1)                    1.51900000   -0.02400000    0.00000000
 C(PDBName=C,ResName=MET,ResNum=1)                     1.87100000    1.44400000    0.00000000
 O(PDBName=O,ResName=MET,ResNum=1)                     0.95794587    2.30996924    0.00000000
 C(PDBName=CB,ResName=MET,ResNum=1)                    2.14600000   -0.74500000   -1.22500000
 C(PDBName=CG,ResName=MET,ResNum=1)                    1.97300000   -2.27800000   -1.26600000
 S(PDBName=SD,ResName=MET,ResNum=1)                    2.99600000   -2.97300000   -2.57500000
 C(PDBName=CE,ResName=MET,ResNum=1)                    2.28800000   -4.62500000   -2.61700000
 H(PDBName=H1,ResName=MET,ResNum=1)                   -0.61500000    0.61900000    0.00000000
 H(PDBName=H2,ResName=MET,ResNum=1)                   -0.37600000   -1.19800000    0.00000000
 H(PDBName=HA,ResName=MET,ResNum=1)                    1.93200000   -0.45700000    0.93000000
 H(PDBName=HB2,ResName=MET,ResNum=1)                   1.73900000   -0.31600000   -2.16300000
 H(PDBName=HB3,ResName=MET,ResNum=1)                   3.23100000   -0.52800000   -1.26200000
 H(PDBName=HG2,ResName=MET,ResNum=1)                   2.24800000   -2.74500000   -0.30200000
 H(PDBName=HG3,ResName=MET,ResNum=1)                   0.91300000   -2.53600000   -1.45300000
 H(PDBName=HE1,ResName=MET,ResNum=1)                   2.33700000   -5.10300000   -1.62200000
 H(PDBName=HE2,ResName=MET,ResNum=1)                   2.83300000   -5.26200000   -3.33500000
 H(PDBName=HE3,ResName=MET,ResNum=1)                   1.22900000   -4.58200000   -2.93000000
 H                                                     2.89685181    1.74818427    0.00000000

 1 2 1.0 9 1.0 10 1.0
 2 3 1.0 5 1.0 11 1.0
 3 4 2.0 19 1.0
 4
 5 6 1.0 12 1.0 13 1.0
 6 7 1.0 14 1.0 15 1.0
 7 8 1.0
 8 16 1.0 17 1.0 18 1.0
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




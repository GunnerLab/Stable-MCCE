%nprocshared=24
%mem=6GB
%chk=glutamic_acid.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

-1 1
 N(PDBName=N,ResName=GLU,ResNum=1)                     0.14500000   -0.22300000    0.00000000
 C(PDBName=CA,ResName=GLU,ResNum=1)                    1.59200000   -0.02700000    0.00000000
 C(PDBName=C,ResName=GLU,ResNum=1)                     1.93800000    1.44300000    0.00000000
 O(PDBName=O,ResName=GLU,ResNum=1)                     1.02134070    2.30515212    0.00000000
 C(PDBName=CB,ResName=GLU,ResNum=1)                    2.18700000   -0.74600000   -1.24200000
 C(PDBName=CG,ResName=GLU,ResNum=1)                    1.79600000   -0.18800000   -2.65100000
 C(PDBName=CD,ResName=GLU,ResNum=1)                    2.28600000   -0.93800000   -3.89200000
 O(PDBName=OE1,ResName=GLU,ResNum=1)                   2.00400000   -0.59400000   -5.03200000
 O(PDBName=OE2,ResName=GLU,ResNum=1)                   3.06300000   -2.02000000   -3.60700000
 H(PDBName=H1,ResName=GLU,ResNum=1)                   -0.54500000    0.60700000    0.00000000
 H(PDBName=H2,ResName=GLU,ResNum=1)                   -0.29800000   -1.20800000    0.00000000
 H(PDBName=HA,ResName=GLU,ResNum=1)                    2.01000000   -0.47200000    0.92200000
 H(PDBName=HB2,ResName=GLU,ResNum=1)                   3.29100000   -0.74300000   -1.16000000
 H(PDBName=HB3,ResName=GLU,ResNum=1)                   1.91600000   -1.81900000   -1.19400000
 H(PDBName=HG2,ResName=GLU,ResNum=1)                   0.69500000   -0.12200000   -2.73200000
 H(PDBName=HG3,ResName=GLU,ResNum=1)                   2.13800000    0.85900000   -2.75400000
 H                                                     2.96257371    1.75146184    0.00000000

 1 2 1.0 10 1.0 11 1.0
 2 3 1.0 5 1.0 12 1.0
 3 4 2.0 17 1.0
 4
 5 6 1.0 13 1.0 14 1.0
 6 7 1.0 15 1.0 16 1.0
 7 8 2.0 9 1.0
 8
 9
 10
 11
 12
 13
 14
 15
 16
 17



%nprocshared=24
%mem=6GB
%chk=valine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=VAL,ResNum=1)                     0.09500000   -0.18100000    0.00000000
 C(PDBName=CA,ResName=VAL,ResNum=1)                    1.54700000   -0.02600000    0.00000000
 C(PDBName=C,ResName=VAL,ResNum=1)                     1.93500000    1.43400000    0.00000000
 O(PDBName=O,ResName=VAL,ResNum=1)                     1.04335474    2.32199734    0.00000000
 C(PDBName=CB,ResName=VAL,ResNum=1)                    2.16600000   -0.77200000   -1.24800000
 C(PDBName=CG1,ResName=VAL,ResNum=1)                   1.97700000   -2.30800000   -1.21700000
 C(PDBName=CG2,ResName=VAL,ResNum=1)                   3.68300000   -0.55100000   -1.50100000
 H(PDBName=H1,ResName=VAL,ResNum=1)                   -0.57100000    0.66900000    0.00000000
 H(PDBName=H2,ResName=VAL,ResNum=1)                   -0.37700000   -1.15200000    0.00000000
 H(PDBName=HA,ResName=VAL,ResNum=1)                    1.95200000   -0.46800000    0.92800000
 H(PDBName=HB,ResName=VAL,ResNum=1)                    1.63800000   -0.39500000   -2.15000000
 H(PDBName=HG11,ResName=VAL,ResNum=1)                  2.52900000   -2.78200000   -0.38400000
 H(PDBName=HG12,ResName=VAL,ResNum=1)                  0.91400000   -2.59000000   -1.10700000
 H(PDBName=HG13,ResName=VAL,ResNum=1)                  2.32300000   -2.78300000   -2.15500000
 H(PDBName=HG21,ResName=VAL,ResNum=1)                  4.30700000   -0.87600000   -0.64700000
 H(PDBName=HG22,ResName=VAL,ResNum=1)                  4.04300000   -1.09400000   -2.39600000
 H(PDBName=HG23,ResName=VAL,ResNum=1)                  3.92800000    0.51100000   -1.69600000
 H                                                     2.96797080    1.71305435    0.00000000

 1 2 1.0 8 1.0 9 1.0
 2 3 1.0 5 1.0 10 1.0
 3 4 2.0 18 1.0
 4
 5 6 1.0 7 1.0 11 1.0
 6 12 1.0 13 1.0 14 1.0
 7 15 1.0 16 1.0 17 1.0
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
 18




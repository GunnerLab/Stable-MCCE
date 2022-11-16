%nprocshared=24
%mem=6GB
%chk=arginine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity polar

Title Card Required

1 1
 N(PDBName=N,ResName=ARG,ResNum=1)                     0.73500000    2.21900000    1.38900000
 C(PDBName=CA,ResName=ARG,ResNum=1)                    2.18900000    2.28500000    1.27400000
 C(PDBName=C,ResName=ARG,ResNum=1)                     2.60000000    3.01000000    0.01400000
 O(PDBName=O,ResName=ARG,ResNum=1)                     1.72248421    3.45951118   -0.76796949
 C(PDBName=CB,ResName=ARG,ResNum=1)                    2.76000000    0.83800000    1.30800000
 C(PDBName=CG,ResName=ARG,ResNum=1)                    4.30200000    0.73100000    1.47300000
 C(PDBName=CD,ResName=ARG,ResNum=1)                    4.78600000   -0.72900000    1.46000000
 N(PDBName=NE,ResName=ARG,ResNum=1)                    6.26700000   -0.77000000    1.64200000
 C(PDBName=CZ,ResName=ARG,ResNum=1)                    7.00700000   -1.87400000    1.67300000
 N(PDBName=NH1,ResName=ARG,ResNum=1)                   6.51800000   -3.07400000    1.53800000
 N(PDBName=NH2,ResName=ARG,ResNum=1)                   8.28500000   -1.75200000    1.84600000
 H(PDBName=H1,ResName=ARG,ResNum=1)                    0.11500000    2.62700000    0.68200000
 H(PDBName=H2,ResName=ARG,ResNum=1)                    0.27300000    1.76300000    2.18300000
 H(PDBName=HA,ResName=ARG,ResNum=1)                    2.57500000    2.86400000    2.13400000
 H(PDBName=HB2,ResName=ARG,ResNum=1)                   2.27700000    0.27600000    2.13200000
 H(PDBName=HB3,ResName=ARG,ResNum=1)                   2.43800000    0.30700000    0.38800000
 H(PDBName=HG2,ResName=ARG,ResNum=1)                   4.81200000    1.29200000    0.66200000
 H(PDBName=HG3,ResName=ARG,ResNum=1)                   4.60400000    1.23100000    2.41600000
 H(PDBName=HD2,ResName=ARG,ResNum=1)                   4.26900000   -1.29600000    2.26400000
 H(PDBName=HD3,ResName=ARG,ResNum=1)                   4.49000000   -1.20500000    0.49900000
 H(PDBName=HE,ResName=ARG,ResNum=1)                    6.80600000    0.09400000    1.76100000
 H(PDBName=HH11,ResName=ARG,ResNum=1)                  7.15400000   -3.87200000    1.57300000
 H(PDBName=HH12,ResName=ARG,ResNum=1)                  5.50500000   -3.10000000    1.40300000
 H(PDBName=HH21,ResName=ARG,ResNum=1)                  8.64900000   -0.80300000    1.94900000
 H(PDBName=HH22,ResName=ARG,ResNum=1)                  8.83700000   -2.61200000    1.86800000
 H                                                     3.63724729    3.14129828   -0.21354959

 1 2 1.0 12 1.0 13 1.0
 2 3 1.0 5 1.0 14 1.0
 3 4 2.0 26 1.0
 4
 5 6 1.0 15 1.0 16 1.0
 6 7 1.0 17 1.0 18 1.0
 7 8 1.0 19 1.0 20 1.0
 8 9 1.5 21 1.0
 9 10 2.0 11 2.0
 10 22 1.0 23 1.0
 11 24 1.0 25 1.0
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
 23
 24
 25
 26




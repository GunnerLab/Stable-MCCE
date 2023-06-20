%nprocshared=24
%mem=6GB
%chk=serine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=SER,ResNum=1)                     0.04300000   -0.05600000    0.00000000
 C(PDBName=CA,ResName=SER,ResNum=1)                    1.50200000   -0.02000000    0.00000000
 C(PDBName=C,ResName=SER,ResNum=1)                     2.00900000    1.40200000    0.00000000
 O(PDBName=O,ResName=SER,ResNum=1)                     1.19379383    2.36064981    0.00000000
 C(PDBName=CB,ResName=SER,ResNum=1)                    2.07300000   -0.82900000   -1.19200000
 O(PDBName=OG,ResName=SER,ResNum=1)                    1.82400000   -0.20700000   -2.45700000
 H(PDBName=H1,ResName=SER,ResNum=1)                   -0.55200000    0.84500000    0.00000000
 H(PDBName=H2,ResName=SER,ResNum=1)                   -0.50700000   -0.98600000    0.00000000
 H(PDBName=HA,ResName=SER,ResNum=1)                    1.85400000   -0.48600000    0.93900000
 H(PDBName=HB2,ResName=SER,ResNum=1)                   3.16800000   -0.94700000   -1.07200000
 H(PDBName=HB3,ResName=SER,ResNum=1)                   1.67500000   -1.86200000   -1.20300000
 H(PDBName=HG,ResName=SER,ResNum=1)                    0.87100000   -0.17800000   -2.58000000
 H                                                     3.06149956    1.59472953    0.00000000

 1 2 1.0 7 1.0 8 1.0
 2 3 1.0 5 1.0 9 1.0
 3 4 2.0 13 1.0
 4
 5 6 1.0 10 1.0 11 1.0
 6 12 1.0
 7
 8
 9
 10
 11
 12
 13




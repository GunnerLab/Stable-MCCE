%chk=glycine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=GLY,ResNum=1)                     0.16200000   -0.20200000    0.00000000
 C(PDBName=CA,ResName=GLY,ResNum=1)                    1.61200000   -0.03100000    0.00000000
 C(PDBName=C,ResName=GLY,ResNum=1)                     1.98500000    1.43200000    0.00000000
 O(PDBName=O,ResName=GLY,ResNum=1)                     1.08442076    2.31093549    0.00000000
 H(PDBName=H1,ResName=GLY,ResNum=1)                   -0.29900000   -1.17900000    0.00000000
 H(PDBName=H2,ResName=GLY,ResNum=1)                   -0.51300000    0.64100000    0.00000000
 H(PDBName=HA2,ResName=GLY,ResNum=1)                   2.05500000   -0.51900000    0.88800000
 H(PDBName=HA3,ResName=GLY,ResNum=1)                   2.05500000   -0.51900000   -0.88800000
 H                                                     3.01509614    1.72148564    0.00000000

 1 2 1.0 5 1.0 6 1.0
 2 3 1.0 7 1.0 8 1.0
 3 4 2.0 9 1.0
 4
 5
 6
 7
 8
 9




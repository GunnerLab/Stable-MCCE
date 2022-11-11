%nprocshared=24
%mem=6GB
%chk=glutamine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=GLN,ResNum=1)                     0.08900000   -0.15300000    0.00000000
 C(PDBName=CA,ResName=GLN,ResNum=1)                    1.54400000   -0.03100000    0.00000000
 C(PDBName=C,ResName=GLN,ResNum=1)                     1.96500000    1.41900000    0.00000000
 O(PDBName=O,ResName=GLN,ResNum=1)                     1.09385457    2.32711685    0.00000000
 C(PDBName=CB,ResName=GLN,ResNum=1)                    2.10300000   -0.76900000   -1.23600000
 C(PDBName=CG,ResName=GLN,ResNum=1)                    3.66100000   -0.78800000   -1.38400000
 C(PDBName=CD,ResName=GLN,ResNum=1)                    4.28600000   -1.45000000   -2.61800000
 O(PDBName=OE1,ResName=GLN,ResNum=1)                   5.49500000   -1.45000000   -2.78900000
 N(PDBName=NE2,ResName=GLN,ResNum=1)                   3.51900000   -2.00500000   -3.52200000
 H(PDBName=H1,ResName=GLN,ResNum=1)                   -0.48600000    0.61500000    0.00000000
 H(PDBName=H2,ResName=GLN,ResNum=1)                   -0.35000000   -1.00700000    0.00000000
 H(PDBName=HA,ResName=GLN,ResNum=1)                    1.93900000   -0.49500000    0.92300000
 H(PDBName=HB2,ResName=GLN,ResNum=1)                   1.73900000   -1.81600000   -1.22300000
 H(PDBName=HB3,ResName=GLN,ResNum=1)                   1.66000000   -0.31800000   -2.14800000
 H(PDBName=HG2,ResName=GLN,ResNum=1)                   4.05700000    0.24400000   -1.39500000
 H(PDBName=HG3,ResName=GLN,ResNum=1)                   4.12200000   -1.24900000   -0.49100000
 H(PDBName=HE21,ResName=GLN,ResNum=1)                  2.51700000   -1.91100000   -3.35400000
 H(PDBName=HE22,ResName=GLN,ResNum=1)                  4.01100000   -2.36500000   -4.34200000
 H                                                     3.00407082    1.67440523    0.00000000

 1 2 1.0 10 1.0 11 1.0
 2 3 1.0 5 1.0 12 1.0
 3 4 2.0 19 1.0
 4
 5 6 1.0 13 1.0 14 1.0
 6 7 1.0 15 1.0 16 1.0
 7 8 2.0 9 2.0
 8
 9 17 1.0 18 1.0
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




%nprocshared=24
%mem=6GB
%chk=tyrosine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=TYR,ResNum=1)                     0.11500000   -0.10500000    0.00000000
 C(PDBName=CA,ResName=TYR,ResNum=1)                    1.57400000   -0.04700000    0.00000000
 C(PDBName=C,ResName=TYR,ResNum=1)                     2.05900000    1.38300000    0.00000000
 O(PDBName=O,ResName=TYR,ResNum=1)                     1.22903031    2.32889686    0.00000000
 C(PDBName=CB,ResName=TYR,ResNum=1)                    2.15700000   -0.75300000   -1.25900000
 C(PDBName=CG,ResName=TYR,ResNum=1)                    1.85600000   -2.24800000   -1.41200000
 C(PDBName=CD1,ResName=TYR,ResNum=1)                   2.56700000   -3.19800000   -0.67300000
 C(PDBName=CD2,ResName=TYR,ResNum=1)                   0.85400000   -2.66800000   -2.29200000
 C(PDBName=CE1,ResName=TYR,ResNum=1)                   2.27300000   -4.55300000   -0.81000000
 C(PDBName=CE2,ResName=TYR,ResNum=1)                   0.56200000   -4.02200000   -2.42800000
 C(PDBName=CZ,ResName=TYR,ResNum=1)                    1.27200000   -4.96300000   -1.68600000
 O(PDBName=OH,ResName=TYR,ResNum=1)                    0.98500000   -6.29300000   -1.81500000
 H(PDBName=H1,ResName=TYR,ResNum=1)                   -0.49300000    0.78800000    0.00000000
 H(PDBName=H2,ResName=TYR,ResNum=1)                   -0.42100000   -1.04200000    0.00000000
 H(PDBName=HA,ResName=TYR,ResNum=1)                    1.95200000   -0.53400000    0.91800000
 H(PDBName=HB2,ResName=TYR,ResNum=1)                   1.82400000   -0.21300000   -2.17000000
 H(PDBName=HB3,ResName=TYR,ResNum=1)                   3.25800000   -0.63200000   -1.27700000
 H(PDBName=HD1,ResName=TYR,ResNum=1)                   3.34200000   -2.88800000    0.01400000
 H(PDBName=HD2,ResName=TYR,ResNum=1)                   0.29000000   -1.94200000   -2.86300000
 H(PDBName=HE1,ResName=TYR,ResNum=1)                   2.81600000   -5.28700000   -0.23400000
 H(PDBName=HE2,ResName=TYR,ResNum=1)                  -0.21900000   -4.33200000   -3.10600000
 H(PDBName=HH,ResName=TYR,ResNum=1)                    0.27600000   -6.38900000   -2.45300000
 H                                                     3.10838529    1.59202277    0.00000000

 1 2 1.0 13 1.0 14 1.0
 2 3 1.0 5 1.0 15 1.0
 3 4 2.0 23 1.0
 4
 5 6 1.0 16 1.0 17 1.0
 6 7 1.5 8 1.5
 7 9 1.5 18 1.0
 8 10 1.5 19 1.0
 9 11 1.5 20 1.0
 10 11 1.5 21 1.0
 11 12 1.0
 12 22 1.0
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




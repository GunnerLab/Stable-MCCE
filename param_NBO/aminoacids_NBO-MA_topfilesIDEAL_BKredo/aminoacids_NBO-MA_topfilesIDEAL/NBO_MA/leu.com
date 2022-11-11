%nprocshared=24
%mem=6GB
%chk=leucine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=LEU,ResNum=1_A)                  -2.06200000   -1.47200000    0.33700000
 C(PDBName=CA,ResName=LEU,ResNum=1_A)                 -1.14400000   -0.38700000    0.26600000
 C(PDBName=C,ResName=LEU,ResNum=1_A)                  -1.80200000    0.70700000   -0.53400000
 O(PDBName=O,ResName=LEU,ResNum=1_A)                  -1.66800000    1.90900000   -0.18700000
 C(PDBName=CB,ResName=LEU,ResNum=1_A)                  0.24800000   -0.82700000   -0.32600000
 C(PDBName=CG,ResName=LEU,ResNum=1_A)                  1.44000000    0.20600000   -0.14300000
 C(PDBName=CD1,ResName=LEU,ResNum=1_A)                 2.66900000   -0.22800000   -0.90600000
 C(PDBName=CD2,ResName=LEU,ResNum=1_A)                 1.80100000    0.45200000    1.33900000
 H(PDBName=H,ResName=LEU,ResNum=1_A)                  -2.38400000    0.46200000   -1.39800000
 H(PDBName=HA,ResName=LEU,ResNum=1_A)                 -0.94200000    0.00600000    1.27300000
 H(PDBName=HB1,ResName=LEU,ResNum=1_A)                 0.54200000   -1.74200000    0.19500000
 H(PDBName=HB2,ResName=LEU,ResNum=1_A)                 0.10500000   -1.14100000   -1.36100000
 H(PDBName=HG,ResName=LEU,ResNum=1_A)                  1.08400000    1.16100000   -0.53400000
 H(PDBName=HD11,ResName=LEU,ResNum=1_A)                2.41100000   -0.61000000   -1.89900000
 H(PDBName=HD12,ResName=LEU,ResNum=1_A)                3.46200000    0.52300000   -0.92600000
 H(PDBName=HD13,ResName=LEU,ResNum=1_A)                3.10500000   -1.09800000   -0.41000000
 H(PDBName=HD21,ResName=LEU,ResNum=1_A)                2.57300000    1.21100000    1.42700000
 H(PDBName=HD22,ResName=LEU,ResNum=1_A)                1.06700000    0.83800000    2.05200000
 H(PDBName=HD23,ResName=LEU,ResNum=1_A)                2.25700000   -0.41600000    1.82100000
 H(PDBName=H1,ResName=LEU,ResNum=1_A)                 -1.85900000   -2.33100000   -0.13300000
 H(PDBName=H2,ResName=LEU,ResNum=1_A)                 -2.91100000   -1.37600000    0.85800000

 1 2 1.0 20 1.0 21 1.0
 2 3 1.0 5 1.0 10 1.0
 3 4 2.0 9 1.0
 4
 5 6 1.0 11 1.0 12 1.0
 6 7 1.0 8 1.0 13 1.0
 7 14 1.0 15 1.0 16 1.0
 8 17 1.0 18 1.0 19 1.0
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




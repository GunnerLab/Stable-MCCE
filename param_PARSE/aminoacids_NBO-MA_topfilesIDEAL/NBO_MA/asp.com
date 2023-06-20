%nprocshared=24
%mem=16GB
%chk=aspartic_acid.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

-1 1
 N(PDBName=N,ResName=ASP,ResNum=1)                     1.88400000    0.02000000    1.37200000
 C(PDBName=CA,ResName=ASP,ResNum=1)                    3.32900000    0.36800000    1.33600000
 C(PDBName=C,ResName=ASP,ResNum=1)                     3.74100000    1.88300000    1.36000000
 O(PDBName=O,ResName=ASP,ResNum=1)                     2.85766500    2.77377968    1.45900202
 C(PDBName=CB,ResName=ASP,ResNum=1)                    3.87400000   -0.37600000    0.10200000
 C(PDBName=CG,ResName=ASP,ResNum=1)                    5.41600000   -0.54800000    0.03100000
 O(PDBName=OD1,ResName=ASP,ResNum=1)                   6.04000000   -0.34100000    1.07500000
 O(PDBName=OD2,ResName=ASP,ResNum=1)                   5.83700000   -0.88700000   -1.07600000
 H(PDBName=H1,ResName=ASP,ResNum=1)                    1.31700000    0.84900000    1.47800000
 H(PDBName=H2,ResName=ASP,ResNum=1)                    1.67600000   -0.57000000    0.57900000
 H(PDBName=HA,ResName=ASP,ResNum=1)                    3.78200000   -0.08300000    2.24100000
 H(PDBName=HB2,ResName=ASP,ResNum=1)                   3.48000000   -1.40900000    0.09200000
 H(PDBName=HB3,ResName=ASP,ResNum=1)                   3.52300000    0.08000000   -0.83800000
 H                                                     4.77284053    2.15796330    1.29217451

 1 2 1.0 9 1.0 10 1.0
 2 3 1.0 5 1.0 11 1.0
 3 4 2.0 14 1.0
 4
 5 6 1.0 12 1.0 13 1.0
 6 7 2.0 8 2.0
 7
 8
 9
 10
 11
 12
 13
 14




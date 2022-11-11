%nprocshared=24
%mem=6GB
%chk=alanine.chk
#p b3lyp/6-31+g(d,p) geom=connectivity pop=nbo

Title Card Required

0 1
 N(PDBName=N,ResName=ALA,ResNum=1)                     0.03900000   -0.02800000    0.00000000
 C(PDBName=CA,ResName=ALA,ResNum=1)                    1.49900000   -0.04300000    0.00000000
 C(PDBName=C,ResName=ALA,ResNum=1)                     2.05500000    1.36100000    0.00000000
 O(PDBName=O,ResName=ALA,ResNum=1)                     1.27347445    2.34723412    0.01142959
 C(PDBName=CB,ResName=ALA,ResNum=1)                    1.95600000   -0.86600000   -1.21700000
 H(PDBName=H1,ResName=ALA,ResNum=1)                   -0.52400000    0.89400000    0.00000000
 H(PDBName=H2,ResName=ALA,ResNum=1)                   -0.54300000   -0.93800000    0.00000000
 H(PDBName=HA,ResName=ALA,ResNum=1)                    1.84700000   -0.53400000    0.92800000
 H(PDBName=HB1,ResName=ALA,ResNum=1)                   3.05800000   -0.93900000   -1.27400000
 H(PDBName=HB2,ResName=ALA,ResNum=1)                   1.57100000   -1.90300000   -1.18100000
 H(PDBName=HB3,ResName=ALA,ResNum=1)                   1.61000000   -0.42500000   -2.17200000
 H                                                     3.11348528    1.51725130   -0.00971842

 1 2 1.0 6 1.0 7 1.0
 2 3 1.0 5 1.0 8 1.0
 3 4 2.0 12 1.0
 4
 5 9 1.0 10 1.0 11 1.0
 6
 7
 8
 9
 10
 11
 12




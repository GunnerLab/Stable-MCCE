#2,3-dimethoxy, 5,6-methyl-Benzoquinone (Q-)  by Agnes 7.9.08      
  
CONFLIST 556        556BK 556-1 556DM

NATOM    556DM      0
NATOM    556BK      0
NATOM    556-1      26

IATOM    556-1  C6  0
IATOM    556-1  C5  1
IATOM    556-1  CM5 2
IATOM    556-1  O3  3
IATOM    556-1  C4  4
IATOM    556-1  C1  5
IATOM    556-1  C2  6
IATOM    556-1  C3  7
IATOM    556-1  CM3 8
IATOM    556-1  O1  9
IATOM    556-1  O4  10
IATOM    556-1  CM2 11
IATOM    556-1 2H1  12
IATOM    556-1 2H2  13
IATOM    556-1 2H3  14
IATOM    556-1 3H1  15
IATOM    556-1 3H2  16
IATOM    556-1 3H3  17
IATOM    556-1  O2  18
IATOM    556-1  CM6 19
IATOM    556-1 5H1  20
IATOM    556-1 5H2  21
IATOM    556-1 5H3  22
IATOM    556-1 6H1  23
IATOM    556-1 6H2  24
IATOM    556-1 6H3  25


ATOMNAME 556-1    0  C6
ATOMNAME 556-1    1  C5
ATOMNAME 556-1    2  CM5
ATOMNAME 556-1    3  O3
ATOMNAME 556-1    4  C4
ATOMNAME 556-1    5  C1
ATOMNAME 556-1    6  C2
ATOMNAME 556-1    7  C3
ATOMNAME 556-1    8  CM3
ATOMNAME 556-1    9  O1
ATOMNAME 556-1   10  O4
ATOMNAME 556-1   11  CM2
ATOMNAME 556-1   12 2H1
ATOMNAME 556-1   13 2H2
ATOMNAME 556-1   14 2H3
ATOMNAME 556-1   15 3H1
ATOMNAME 556-1   16 3H2
ATOMNAME 556-1   17 3H3
ATOMNAME 556-1   18  O2
ATOMNAME 556-1   19  CM6
ATOMNAME 556-1   20 5H1
ATOMNAME 556-1   21 5H2
ATOMNAME 556-1   22 5H3
ATOMNAME 556-1   23 6H1
ATOMNAME 556-1   24 6H2
ATOMNAME 556-1   25 6H3


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   556-1      0
PKA      556-1      0.0
ELECTRON 556-1      1
EM       556-1      -75
RXN      556-1      -13.7

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

CONNECT  556-1  C6  sp2        0    C5   0    C1   0    CM6
CONNECT  556-1  C5  sp2        0    C6   0    C4   0    CM5
CONNECT  556-1  CM5 sp3        0    C5   0   5H1   0   5H2   0   5H3
CONNECT  556-1  O3  sp2        0    C3   0    CM3
CONNECT  556-1  C4  sp2        0    C5   0    C3   0    O4
CONNECT  556-1  C1  sp2        0    C6   0    C2   0    O1
CONNECT  556-1  C2  sp2        0    C1   0    C3   0    O2
CONNECT  556-1  C3  sp2        0    C4   0    C2   0    O3
CONNECT  556-1  CM3 sp3        0    O3   0   3H1   0   3H2   0   3H3
CONNECT  556-1  CM2 sp3        0    O2   0   2H1   0   2H2   0   2H3
CONNECT  556-1  O1  s          0    C1
CONNECT  556-1  O4  s          0    C4
CONNECT  556-1 2H1  s          0    CM2
CONNECT  556-1 2H2  s          0    CM2
CONNECT  556-1 2H3  s          0    CM2
CONNECT  556-1 3H1  s          0    CM3
CONNECT  556-1 3H2  s          0    CM3
CONNECT  556-1 3H3  s          0    CM3
CONNECT  556-1  O2  sp2        0    C2   0    CM2
CONNECT  556-1  CM6 sp3        0    C6   0   6H1   0   6H2   0   6H3
CONNECT  556-1 5H1  s          0    C5
CONNECT  556-1 5H2  s          0    C5
CONNECT  556-1 5H3  s          0    C5
CONNECT  556-1 6H1  s          0    CM6
CONNECT  556-1 6H2  s          0    CM6
CONNECT  556-1 6H3  s          0    CM6

#de-protonated---------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   556    C6  1.70
RADIUS   556    C5  1.70
RADIUS   556    CM5 1.70
RADIUS   556    O3  1.40
RADIUS   556    C4  1.70
RADIUS   556    C1  1.70
RADIUS   556    C2  1.70
RADIUS   556    C3  1.70
RADIUS   556    CM3 1.70
RADIUS   556    O1  1.40
RADIUS   556    O4  1.40
RADIUS   556    CM2 1.70
RADIUS   556   2H1  1.00
RADIUS   556   2H2  1.00
RADIUS   556   2H3  1.00
RADIUS   556   3H1  1.00
RADIUS   556   3H2  1.00
RADIUS   556   3H3  1.00
RADIUS   556    O2  1.40
RADIUS   556    CM6 1.70
RADIUS   556   5H1  1.00
RADIUS   556   5H2  1.00
RADIUS   556   5H3  1.00
RADIUS   556   6H1  1.00
RADIUS   556   6H2  1.00
RADIUS   556   6H3  1.00



#NEUTRAL------
#23456789A123456789B123456789C
 
#p ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=800)   Agnes 7.9.08
CHARGE   556-1  O1  -0.57
CHARGE   556-1  C1   0.32
CHARGE   556-1  C2   0.09
CHARGE   556-1  O2  -0.36
CHARGE   556-1  CM2  0.23
CHARGE   556-1  C3   0.08
CHARGE   556-1  CM3  0.24
CHARGE   556-1  C4   0.32
CHARGE   556-1  O4  -0.57
CHARGE   556-1  C5  -0.19 
CHARGE   556-1  CM5  0.30
CHARGE   556-1 2H1  -0.03 
CHARGE   556-1 2H2   0.00
CHARGE   556-1 2H3  -0.02 
CHARGE   556-1 3H1  -0.04
CHARGE   556-1 3H2  -0.01 
CHARGE   556-1 3H3  -0.02 
CHARGE   556-1  O3  -0.37
CHARGE   556-1 5H1  -0.08 
CHARGE   556-1 5H2  -0.08 
CHARGE   556-1 5H3  -0.08 
CHARGE   556-1  C6  -0.20 
CHARGE   556-1  CM6  0.30
CHARGE   556-1 6H1  -0.08
CHARGE   556-1 6H2  -0.08 
CHARGE   556-1 6H3  -0.08 
 

#ParaNam|Res  |Atom|Param/toggle
TRANS    556          t

#TORSION  556-1  CM2  O2   C2   C1   f        0.0         2    0.00
#TORSION  556-1  CM3  O3   C3   C4   f        0.0         2    0.00

TORSION  556-1  CM2  O2   C2   C3   f      0.0           2    180.00
TORSION  556-1  CM3  O3   C3   C2   f      0.0           2    180.00
TORSION  556-1  CM2  O2   C2   C1   f      0.215         6    180.00     -0.195        4    180.00    -0.427         3    180.00     4.208         2    180.00    -0.632         1    180.00
TORSION  556-1  CM3  O3   C3   C4   f      0.215         6    180.00     -0.195        4    180.00    -0.427         3    180.00     4.208         2    180.00    -0.632         1    180.00
TORSION  556-1  2H1  CM2  O2   C2   f     -0.180         3      0.00
TORSION  556-1  3H1  CM3  O3   C3   f     -0.180         3      0.00
TORSION  556-1  6H1  CM6  C6   C5   f     -0.180         3      0.00
TORSION  556-1  5H1  CM5  C5   C6   f     -0.180         3      0.00


#TORSION  556-1  C2   O2   CM2  2H1  f      0.180         3      0.00
#TORSION  556-1  C3   O3   CM3  3H1  f      0.180         3      0.00
#TORSION  556-1  C5   C6   CM6  6H1  f      0.180         3      0.00
#TORSION  556-1  C6   C5   CM5  5H1  f      0.180         3      0.00


#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  556   0     C2 - O2   CM2
ROTAMER  556   1     C3 - O3   CM3
#=========================================================================


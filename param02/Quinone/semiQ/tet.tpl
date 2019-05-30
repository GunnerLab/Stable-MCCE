#Tetramethyl-Benzoquinone  (Q-)      
#Agnes 7.9.08  
CONFLIST TET        TETBK TET-1 TETDM

NATOM    TETDM      0
NATOM    TETBK      0
NATOM    TET-1      24

IATOM    TET-1  C6  0
IATOM    TET-1  C5  1
IATOM    TET-1  CM6  2
IATOM    TET-1  CM3 3
IATOM    TET-1  C4  4
IATOM    TET-1  C1  5
IATOM    TET-1  C2  6
IATOM    TET-1  C3  7
IATOM    TET-1  CM5 8
IATOM    TET-1  O1  9
IATOM    TET-1  O4  10
IATOM    TET-1  CM2 11
IATOM    TET-1 2H1  12
IATOM    TET-1 2H2  13
IATOM    TET-1 2H3  14
IATOM    TET-1 5H1  15
IATOM    TET-1 5H2  16
IATOM    TET-1 5H3  17
IATOM    TET-1 3H1  18
IATOM    TET-1 3H2  19
IATOM    TET-1 3H3  20
IATOM    TET-1 6H1  21
IATOM    TET-1 6H2  22
IATOM    TET-1 6H3  23

ATOMNAME TET-1    0  C6
ATOMNAME TET-1    1  C5
ATOMNAME TET-1    2  CM6
ATOMNAME TET-1    3  CM3
ATOMNAME TET-1    4  C4
ATOMNAME TET-1    5  C1
ATOMNAME TET-1    6  C2
ATOMNAME TET-1    7  C3
ATOMNAME TET-1    8  CM5
ATOMNAME TET-1    9  O1
ATOMNAME TET-1   10  O4
ATOMNAME TET-1   11  CM2
ATOMNAME TET-1   12 2H1
ATOMNAME TET-1   13 2H2
ATOMNAME TET-1   14 2H3
ATOMNAME TET-1   15 5H1
ATOMNAME TET-1   16 5H2
ATOMNAME TET-1   17 5H3
ATOMNAME TET-1   18 3H1
ATOMNAME TET-1   19 3H2
ATOMNAME TET-1   20 3H3
ATOMNAME TET-1   21 6H1
ATOMNAME TET-1   22 6H2
ATOMNAME TET-1   23 6H3


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   TET-1      0
PKA      TET-1      0.0
ELECTRON TET-1      1
EM       TET-1      -255 
RXN      TET-1      -13.281

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

CONNECT  TET-1  C6  sp2        0    C5   0    C1   0    CM6
CONNECT  TET-1  C5  sp2        0    C6   0    C4   0    CM5
CONNECT  TET-1  CM6 sp3        0    C6   0   6H1   0   6H2   0   6H3
CONNECT  TET-1  CM3 sp3        0    C3   0   3H1   0   3H2   0   3H3
CONNECT  TET-1  C4  sp2        0    C5   0    C3   0    O4
CONNECT  TET-1  C1  sp2        0    C6   0    C2   0    O1
CONNECT  TET-1  C2  sp2        0    C1   0    C3   0    CM2
CONNECT  TET-1  C3  sp2        0    C4   0    C2   0    CM3
CONNECT  TET-1  CM5 sp3        0    C5   0   5H1   0   5H2   0   5H3
CONNECT  TET-1  CM2 sp3        0    C2   0   2H1   0   2H2   0   2H3
CONNECT  TET-1  O1  s          0    C1
CONNECT  TET-1  O4  s          0    C4
CONNECT  TET-1 2H1  s          0    CM2
CONNECT  TET-1 2H2  s          0    CM2
CONNECT  TET-1 2H3  s          0    CM2
CONNECT  TET-1 5H1  s          0    CM5
CONNECT  TET-1 5H2  s          0    CM5
CONNECT  TET-1 5H3  s          0    CM5
CONNECT  TET-1 3H1  s          0    CM3
CONNECT  TET-1 3H2  s          0    CM3
CONNECT  TET-1 3H3  s          0    CM3
CONNECT  TET-1 6H1  s          0    CM6
CONNECT  TET-1 6H2  s          0    CM6
CONNECT  TET-1 6H3  s          0    CM6


#de-protonated---------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   TET    C6  1.70
RADIUS   TET    C5  1.70
RADIUS   TET    CM6 1.70
RADIUS   TET    CM3 1.70
RADIUS   TET    C4  1.70
RADIUS   TET    C1  1.70
RADIUS   TET    C2  1.70
RADIUS   TET    C3  1.70
RADIUS   TET    O1  1.40
RADIUS   TET    O4  1.40
RADIUS   TET    CM5 1.70
RADIUS   TET    CM2 1.70
RADIUS   TET   2H1  1.00
RADIUS   TET   2H2  1.00
RADIUS   TET   2H3  1.00
RADIUS   TET   5H1  1.00
RADIUS   TET   5H2  1.00
RADIUS   TET   5H3  1.00
RADIUS   TET   3H1  1.00
RADIUS   TET   3H2  1.00
RADIUS   TET   3H3  1.00
RADIUS   TET   6H1  1.00
RADIUS   TET   6H2  1.00
RADIUS   TET   6H3  1.00


#NEUTRAL------
#23456789A123456789B123456789C

# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   Agnes 7.9.08
CHARGE   TET-1  O1  -0.60
CHARGE   TET-1  C1   0.37
CHARGE   TET-1  C4   0.37
CHARGE   TET-1  O4  -0.60
CHARGE   TET-1  C2  -0.16 
CHARGE   TET-1  CM2  0.20
CHARGE   TET-1  C3  -0.15 
CHARGE   TET-1  CM5  0.20
CHARGE   TET-1  C5  -0.16
CHARGE   TET-1  CM3  0.17
CHARGE   TET-1  C6  -0.15 
CHARGE   TET-1  CM6  0.17
CHARGE   TET-1 2H1  -0.06 
CHARGE   TET-1 2H2  -0.06 
CHARGE   TET-1 2H3  -0.06 
CHARGE   TET-1 5H1  -0.06 
CHARGE   TET-1 5H2  -0.06 
CHARGE   TET-1 5H3  -0.06 
CHARGE   TET-1 3H1  -0.05
CHARGE   TET-1 3H2  -0.05 
CHARGE   TET-1 3H3  -0.05 
CHARGE   TET-1 6H1  -0.05 
CHARGE   TET-1 6H2  -0.05 
CHARGE   TET-1 6H3  -0.05 
 
#TORSION  TET    HO2  O2   C2   C3   f      1.800         2    180.00
TORSION  TET-1  2H1  CM2  C2   C1   f     -0.180         3      0.00
TORSION  TET-1  3H1  CM3  C3   C4   f     -0.180         3      0.00
TORSION  TET-1  6H1  CM6  C6   C5   f     -0.180         3      0.00
TORSION  TET-1  5H1  CM5  C5   C6   f     -0.180         3      0.00



#ParaNam|Res  |Atom|Param/toggle
TRANS    TET          t


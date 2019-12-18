#Trimethyl-Benzoquinone (Q-)
#Agnes 7.11.08  
CONFLIST TRI        TRIBK TRI-1 TRIDM

NATOM    TRIDM      0
NATOM    TRIBK      0
NATOM    TRI-1      21

IATOM    TRI-1  C6  0
IATOM    TRI-1  C5  1
IATOM    TRI-1  H6  2
IATOM    TRI-1  CM3 3
IATOM    TRI-1  C4  4
IATOM    TRI-1  C1  5
IATOM    TRI-1  C2  6
IATOM    TRI-1  C3  7
IATOM    TRI-1  CM5 8
IATOM    TRI-1  O1  9
IATOM    TRI-1  O4  10
IATOM    TRI-1  CM2 11
IATOM    TRI-1 2H1  12
IATOM    TRI-1 2H2  13
IATOM    TRI-1 2H3  14
IATOM    TRI-1 5H1  15
IATOM    TRI-1 5H2  16
IATOM    TRI-1 5H3  17
IATOM    TRI-1 3H1  18
IATOM    TRI-1 3H2  19
IATOM    TRI-1 3H3  20


ATOMNAME TRI-1    0  C6
ATOMNAME TRI-1    1  C5
ATOMNAME TRI-1    2  H6
ATOMNAME TRI-1    3  CM3
ATOMNAME TRI-1    4  C4
ATOMNAME TRI-1    5  C1
ATOMNAME TRI-1    6  C2
ATOMNAME TRI-1    7  C3
ATOMNAME TRI-1    8  CM5
ATOMNAME TRI-1    9  O1
ATOMNAME TRI-1   10  O4
ATOMNAME TRI-1   11  CM2
ATOMNAME TRI-1   12 2H1
ATOMNAME TRI-1   13 2H2
ATOMNAME TRI-1   14 2H3
ATOMNAME TRI-1   15 5H1
ATOMNAME TRI-1   16 5H2
ATOMNAME TRI-1   17 5H3
ATOMNAME TRI-1   18 3H1
ATOMNAME TRI-1   19 3H2
ATOMNAME TRI-1   20 3H3

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   TRI-1      0
PKA      TRI-1      0.0
ELECTRON TRI-1      1
EM       TRI-1      -165
RXN      TRI-1      -13.8


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

CONNECT  TRI-1  C6  sp2        0    C5   0    C1   0    H6
CONNECT  TRI-1  C5  sp2        0    C6   0    C4   0    CM5
CONNECT  TRI-1  H6  s          0    C6
CONNECT  TRI-1  CM3 sp3        0    C3   0   3H1   0   3H2   0   3H3
CONNECT  TRI-1  C4  sp2        0    C5   0    C3   0    O4
CONNECT  TRI-1  C1  sp2        0    C6   0    C2   0    O1
CONNECT  TRI-1  C2  sp2        0    C1   0    C3   0    CM2
CONNECT  TRI-1  C3  sp2        0    C4   0    C2   0    CM3
CONNECT  TRI-1  CM5 sp3        0    C5   0   5H1   0   5H2   0   5H3
CONNECT  TRI-1  CM2 sp3        0    C2   0   2H1   0   2H2   0   2H3
CONNECT  TRI-1  O1  s          0    C1
CONNECT  TRI-1  O4  s          0    C4
CONNECT  TRI-1 2H1  s          0    CM2
CONNECT  TRI-1 2H2  s          0    CM2
CONNECT  TRI-1 2H3  s          0    CM2
CONNECT  TRI-1 5H1  s          0    CM5
CONNECT  TRI-1 5H2  s          0    CM5
CONNECT  TRI-1 5H3  s          0    CM5
CONNECT  TRI-1 3H1  s          0    CM3
CONNECT  TRI-1 3H2  s          0    CM3
CONNECT  TRI-1 3H3  s          0    CM3

#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   TRI    C6  1.70
RADIUS   TRI    C5  1.70
RADIUS   TRI    H6  1.00
RADIUS   TRI    CM3 1.70
RADIUS   TRI    C4  1.70
RADIUS   TRI    C1  1.70
RADIUS   TRI    C2  1.70
RADIUS   TRI    C3  1.70
RADIUS   TRI    O1  1.40
RADIUS   TRI    O4  1.40
RADIUS   TRI    CM5 1.70
RADIUS   TRI    CM2 1.70
RADIUS   TRI   2H1  1.00
RADIUS   TRI   2H2  1.00
RADIUS   TRI   2H3  1.00
RADIUS   TRI   5H1  1.00
RADIUS   TRI   5H2  1.00
RADIUS   TRI   5H3  1.00
RADIUS   TRI   3H1  1.00
RADIUS   TRI   3H2  1.00
RADIUS   TRI   3H3  1.00


#NEUTRAL------
#23456789A123456789B123456789C

# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   Agnes 7.10.08
CHARGE   TRI-1  O1  -0.70
CHARGE   TRI-1  C1   0.57
CHARGE   TRI-1  C4   0.36
CHARGE   TRI-1  O4  -0.63
CHARGE   TRI-1  C2  -0.17 
CHARGE   TRI-1  CM2  0.08
CHARGE   TRI-1  C3  -0.14 
CHARGE   TRI-1  CM5  0.11
CHARGE   TRI-1  C5   0.00
CHARGE   TRI-1  CM3  0.15
CHARGE   TRI-1  C6  -0.46 
CHARGE   TRI-1  H6   0.13
CHARGE   TRI-1 2H1  -0.01 
CHARGE   TRI-1 2H2  -0.03 
CHARGE   TRI-1 2H3  -0.03 
CHARGE   TRI-1 5H1  -0.04 
CHARGE   TRI-1 5H2  -0.04 
CHARGE   TRI-1 5H3  -0.04 
CHARGE   TRI-1 3H1  -0.04 
CHARGE   TRI-1 3H2  -0.04 
CHARGE   TRI-1 3H3  -0.05 
 
#TORSION  TRI    HO2  O2   C2   C3   f      1.800         2    180.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    TRI           t


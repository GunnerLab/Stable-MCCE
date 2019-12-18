#2,5-dimethyl-Benzoquinone (Q-)       
# Agnes 7.9.08  
CONFLIST 25D        25DBK 25D-1 25DDM

NATOM    25DDM      0
NATOM    25DBK      0
NATOM    25D-1      18

IATOM    25D-1  C6  0
IATOM    25D-1  C5  1
IATOM    25D-1  H6  2
IATOM    25D-1  H3  3
IATOM    25D-1  C4  4
IATOM    25D-1  C1  5
IATOM    25D-1  C2  6
IATOM    25D-1  C3  7
IATOM    25D-1  CM5 8
IATOM    25D-1  O1  9
IATOM    25D-1  O4  10
IATOM    25D-1  CM2 11
IATOM    25D-1 2H1  12
IATOM    25D-1 2H2  13
IATOM    25D-1 2H3  14
IATOM    25D-1 5H1  15
IATOM    25D-1 5H2  16
IATOM    25D-1 5H3  17


ATOMNAME 25D-1    0  C6
ATOMNAME 25D-1    1  C5
ATOMNAME 25D-1    2  H6
ATOMNAME 25D-1    3  H3
ATOMNAME 25D-1    4  C4
ATOMNAME 25D-1    5  C1
ATOMNAME 25D-1    6  C2
ATOMNAME 25D-1    7  C3
ATOMNAME 25D-1    8  CM5
ATOMNAME 25D-1    9  O1
ATOMNAME 25D-1   10  O4
ATOMNAME 25D-1   11  CM2
ATOMNAME 25D-1   12 2H1
ATOMNAME 25D-1   13 2H2
ATOMNAME 25D-1   14 2H3
ATOMNAME 25D-1   15 5H1
ATOMNAME 25D-1   16 5H2
ATOMNAME 25D-1   17 5H3



#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   25D-1      0
PKA      25D-1      0.0
ELECTRON 25D-1      1
EM       25D-1      -66
RXN      25D-1      -14.10

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|


CONNECT  25D-1  C6  sp2        0    C5   0    C1   0    H6
CONNECT  25D-1  C5  sp2        0    C6   0    C4   0    CM5
CONNECT  25D-1  H6  s          0    C6
CONNECT  25D-1  H3  s          0    C3
CONNECT  25D-1  C4  sp2        0    C5   0    C3   0    O4
CONNECT  25D-1  C1  sp2        0    C6   0    C2   0    O1
CONNECT  25D-1  C2  sp2        0    C1   0    C3   0    CM2
CONNECT  25D-1  C3  sp2        0    C4   0    C2   0    H3
CONNECT  25D-1  CM5 sp3        0    C5   0   5H1   0   5H2   0   5H3
CONNECT  25D-1  CM2 sp3        0    C2   0   2H1   0   2H2   0   2H3
CONNECT  25D-1  O1  s          0    C1
CONNECT  25D-1  O4  s          0    C4
CONNECT  25D-1 2H1  s          0    CM2
CONNECT  25D-1 2H2  s          0    CM2
CONNECT  25D-1 2H3  s          0    CM2
CONNECT  25D-1 5H1  s          0    CM5
CONNECT  25D-1 5H2  s          0    CM5
CONNECT  25D-1 5H3  s          0    CM5


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   25D    C6  1.70
RADIUS   25D    C5  1.70
RADIUS   25D    H6  1.00
RADIUS   25D    H3  1.00
RADIUS   25D    C4  1.70
RADIUS   25D    C1  1.70
RADIUS   25D    C2  1.70
RADIUS   25D    C3  1.70
RADIUS   25D    O1  1.40
RADIUS   25D    O4  1.40
RADIUS   25D    CM5 1.70
RADIUS   25D    CM2 1.70
RADIUS   25D   2H1  1.00
RADIUS   25D   2H2  1.00
RADIUS   25D   2H3  1.00
RADIUS   25D   5H1  1.00
RADIUS   25D   5H2  1.00
RADIUS   25D   5H3  1.00


#NEUTRAL------
#23456789A123456789B123456789C

# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   Agnes 7.9.08
CHARGE   25D-1  O1  -0.70
CHARGE   25D-1  C1   0.54
CHARGE   25D-1  C4   0.54
CHARGE   25D-1  O4  -0.70
CHARGE   25D-1  C2  -0.03 
CHARGE   25D-1  CM2  0.06
CHARGE   25D-1  C3  -0.42 
CHARGE   25D-1  CM5  0.06
CHARGE   25D-1  C5  -0.03 
CHARGE   25D-1  H3   0.13
CHARGE   25D-1  C6  -0.42 
CHARGE   25D-1  H6   0.13
CHARGE   25D-1 2H1  -0.03
CHARGE   25D-1 2H2  -0.02 
CHARGE   25D-1 2H3  -0.03 
CHARGE   25D-1 5H1  -0.02 
CHARGE   25D-1 5H2  -0.03
CHARGE   25D-1 5H3  -0.03

 
#TORSION  25D    HO2  O2   C2   C3   f      1.800         2    180.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    25D          t


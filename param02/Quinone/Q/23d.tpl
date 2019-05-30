#2,3-dimethyl-Benzoquinone        
#Gennady Khirich 7.31.06  
CONFLIST 23D        23DBK 23D01 23DDM

NATOM    23DDM      0
NATOM    23DBK      0
NATOM    23D01      18

IATOM    23D01  C6  0
IATOM    23D01  C5  1
IATOM    23D01  H6  2
IATOM    23D01  H5  3
IATOM    23D01  C4  4
IATOM    23D01  C1  5
IATOM    23D01  C2  6
IATOM    23D01  C3  7
IATOM    23D01  CM3 8
IATOM    23D01  O1  9
IATOM    23D01  O4  10
IATOM    23D01  CM2 11
IATOM    23D01 2H1  12
IATOM    23D01 2H2  13
IATOM    23D01 2H3  14
IATOM    23D01 3H1  15
IATOM    23D01 3H2  16
IATOM    23D01 3H3  17

ATOMNAME 23D01    0  C6 
ATOMNAME 23D01    1  C5 
ATOMNAME 23D01    2  H6 
ATOMNAME 23D01    3  H5 
ATOMNAME 23D01    4  C4 
ATOMNAME 23D01    5  C1 
ATOMNAME 23D01    6  C2 
ATOMNAME 23D01    7  C3 
ATOMNAME 23D01    8  CM3 
ATOMNAME 23D01    9  O1 
ATOMNAME 23D01   10  O4 
ATOMNAME 23D01   11  CM2 
ATOMNAME 23D01   12 2H1
ATOMNAME 23D01   13 2H2
ATOMNAME 23D01   14 2H3
ATOMNAME 23D01   15 3H1
ATOMNAME 23D01   16 3H2
ATOMNAME 23D01   17 3H3

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   23D01      0

PKA      23D01      0.0

ELECTRON 23D01      0

EM       23D01      0.0

RXN      23D01      -2.071


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  23D01  C6  sp2        0    C5   0    C1   0    H6
CONNECT  23D01  C5  sp2        0    C6   0    C4   0    H5
CONNECT  23D01  H6  s          0    C6
CONNECT  23D01  H5  s          0    C5
CONNECT  23D01  C4  sp2        0    C5   0    C3   0    O4
CONNECT  23D01  C1  sp2        0    C6   0    C2   0    O1
CONNECT  23D01  C2  sp2        0    C1   0    C3   0    CM2
CONNECT  23D01  C3  sp2        0    C4   0    C2   0    CM3
CONNECT  23D01  CM3 sp3        0    C3   0   3H1   0   3H2   0   3H3  
CONNECT  23D01  CM2 sp3        0    C2   0   2H1   0   2H2   0   2H3   
CONNECT  23D01  O1  s          0    C1
CONNECT  23D01  O4  s          0    C4
CONNECT  23D01 2H1  s          0    CM2
CONNECT  23D01 2H2  s          0    CM2
CONNECT  23D01 2H3  s          0    CM2
CONNECT  23D01 3H1  s          0    CM3
CONNECT  23D01 3H2  s          0    CM3
CONNECT  23D01 3H3  s          0    CM3

#de-protonated---------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   23D    C6  1.70
RADIUS   23D    C5  1.70
RADIUS   23D    H6  1.00
RADIUS   23D    H5  1.00
RADIUS   23D    C4  1.70
RADIUS   23D    C1  1.70
RADIUS   23D    C2  1.70
RADIUS   23D    C3  1.70
RADIUS   23D    O1  1.40
RADIUS   23D    O4  1.40
RADIUS   23D    CM3 1.70
RADIUS   23D    CM2 1.70
RADIUS   23D   2H1  1.00
RADIUS   23D   2H2  1.00
RADIUS   23D   2H3  1.00
RADIUS   23D   3H1  1.00
RADIUS   23D   3H2  1.00
RADIUS   23D   3H3  1.00


#NEUTRAL------
#23456789A123456789B123456789C
# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   gk 7.31.06
CHARGE   23D01  O1  -0.49
CHARGE   23D01  C1   0.58
CHARGE   23D01  C4   0.58
CHARGE   23D01  O4  -0.50
CHARGE   23D01  C2  -0.12
CHARGE   23D01  CM2  0.01
CHARGE   23D01  C3  -0.12
CHARGE   23D01  CM3  0.02
CHARGE   23D01  C5  -0.21
CHARGE   23D01  H5   0.14
CHARGE   23D01  C6  -0.21
CHARGE   23D01  H6   0.14
CHARGE   23D01 2H1   0.03
CHARGE   23D01 2H2   0.03
CHARGE   23D01 2H3   0.03
CHARGE   23D01 3H1   0.03
CHARGE   23D01 3H2   0.03
CHARGE   23D01 3H3   0.03
 

#ParaNam|Res  |Atom|Param/toggle
TRANS    23D          t


TORSION  23D01   C5  C4   C3   C2   f      0.0           2    180.00
TORSION  23D01   C1  C2   C3   C4   f      0.0           2    180.00
TORSION  23D01  3H1  CM3  C3   C4   f     -0.180         3      0.00
TORSION  23D01  2H1  CM2  C2   C1   f     -0.180         3      0.00



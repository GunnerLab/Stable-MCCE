#Benzoquinone        
#Gennady Khirich 7.31.06  
CONFLIST BZQ        BZQBK BZQ01 BZQDM

NATOM    BZQDM      0
NATOM    BZQBK      0
NATOM    BZQ01      12

IATOM    BZQ01  C6  0
IATOM    BZQ01  C5  1
IATOM    BZQ01  H6  2
IATOM    BZQ01  H5  3
IATOM    BZQ01  C4  4
IATOM    BZQ01  C1  5
IATOM    BZQ01  C2  6
IATOM    BZQ01  C3  7
IATOM    BZQ01  H3  8
IATOM    BZQ01  O1  9
IATOM    BZQ01  O4  10
IATOM    BZQ01  H2  11

ATOMNAME BZQ01    0  C6 
ATOMNAME BZQ01    1  C5 
ATOMNAME BZQ01    2  H6 
ATOMNAME BZQ01    3  H5 
ATOMNAME BZQ01    4  C4 
ATOMNAME BZQ01    5  C1 
ATOMNAME BZQ01    6  C2 
ATOMNAME BZQ01    7  C3 
ATOMNAME BZQ01    8  H3 
ATOMNAME BZQ01    9  O1 
ATOMNAME BZQ01   10  O4 
ATOMNAME BZQ01   11  H2 



#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   BZQ01      0

PKA      BZQ01      0.0

ELECTRON BZQ01      0

EM       BZQ01      0.0

RXN      BZQ01      -2.318


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  BZQ01  C5  sp2        0    C4   0    C6   0    H5
CONNECT  BZQ01  C6  sp2        0    C1   0    C5   0    H6
CONNECT  BZQ01  C4  sp2        0    C5   0    C3   0    O4
CONNECT  BZQ01  C1  sp2        0    C6   0    C2   0    O1
CONNECT  BZQ01  C2  sp2        0    C1   0    C3   0    H2
CONNECT  BZQ01  C3  sp2        0    C4   0    C2   0    H3
CONNECT  BZQ01  H3  s          0    C3
CONNECT  BZQ01  H2  s          0    C2
CONNECT  BZQ01  O1  s          0    C1
CONNECT  BZQ01  O4  s          0    C4
CONNECT  BZQ01  H6  s          0    C6
CONNECT  BZQ01  H5  s          0    C5


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   BZQ    C6  1.70
RADIUS   BZQ    C5  1.70
RADIUS   BZQ    H6  1.00
RADIUS   BZQ    H5  1.00
RADIUS   BZQ    C4  1.70
RADIUS   BZQ    C1  1.70
RADIUS   BZQ    C2  1.70
RADIUS   BZQ    C3  1.70
RADIUS   BZQ    O1  1.40
RADIUS   BZQ    O4  1.40
RADIUS   BZQ    H3  1.00
RADIUS   BZQ    H2  1.00


#NEUTRAL------
#23456789A123456789B123456789C
# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   gk 7.31.06
CHARGE   BZQ01  O1  -0.51
CHARGE   BZQ01  C1   0.65
CHARGE   BZQ01  C4   0.65
CHARGE   BZQ01  O4  -0.51
CHARGE   BZQ01  C2  -0.21
CHARGE   BZQ01  H2   0.14
CHARGE   BZQ01  C3  -0.21
CHARGE   BZQ01  H3   0.14
CHARGE   BZQ01  C5  -0.21
CHARGE   BZQ01  H5   0.14
CHARGE   BZQ01  C6  -0.21
CHARGE   BZQ01  H6   0.14

 
#TORSION  BZQ    HO2  O2   C2   C3   f      1.800         2    180.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    BZQ          t


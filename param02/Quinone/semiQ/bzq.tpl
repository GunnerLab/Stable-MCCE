#Benzoquinone (Q-)       
#Agnes 7.9.08  
CONFLIST BZQ        BZQBK BZQ-1 BZQDM

NATOM    BZQDM      0
NATOM    BZQBK      0
NATOM    BZQ-1      12

IATOM    BZQ-1  C6  0
IATOM    BZQ-1  C5  1
IATOM    BZQ-1  H6  2
IATOM    BZQ-1  H5  3
IATOM    BZQ-1  C4  4
IATOM    BZQ-1  C1  5
IATOM    BZQ-1  C2  6
IATOM    BZQ-1  C3  7
IATOM    BZQ-1  H3  8
IATOM    BZQ-1  O1  9
IATOM    BZQ-1  O4  10
IATOM    BZQ-1  H2  11


ATOMNAME BZQ-1    0  C6
ATOMNAME BZQ-1    1  C5
ATOMNAME BZQ-1    2  H6
ATOMNAME BZQ-1    3  H5
ATOMNAME BZQ-1    4  C4
ATOMNAME BZQ-1    5  C1
ATOMNAME BZQ-1    6  C2
ATOMNAME BZQ-1    7  C3
ATOMNAME BZQ-1    8  H3
ATOMNAME BZQ-1    9  O1
ATOMNAME BZQ-1   10  O4
ATOMNAME BZQ-1   11  H2



#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   BZQ-1      0
PKA      BZQ-1      0.0
ELECTRON BZQ-1      1
EM       BZQ-1      99 # 80 
RXN      BZQ-1      -14.793

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

CONNECT  BZQ-1  C5  sp2        0    C4   0    C6   0    H5
CONNECT  BZQ-1  C6  sp2        0    C1   0    C5   0    H6
CONNECT  BZQ-1  C4  sp2        0    C5   0    C3   0    O4
CONNECT  BZQ-1  C1  sp2        0    C6   0    C2   0    O1
CONNECT  BZQ-1  C2  sp2        0    C1   0    C3   0    H2
CONNECT  BZQ-1  C3  sp2        0    C4   0    C2   0    H3
CONNECT  BZQ-1  H3  s          0    C3
CONNECT  BZQ-1  H2  s          0    C2
CONNECT  BZQ-1  O1  s          0    C1
CONNECT  BZQ-1  O4  s          0    C4
CONNECT  BZQ-1  H6  s          0    C6
CONNECT  BZQ-1  H5  s          0    C5

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

 
#TORSION  BZQ    HO2  O2   C2   C3   f      1.800         2    180.00

# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   Agnes 7.31.08
CHARGE   BZQ-1  O1  -0.75
CHARGE   BZQ-1  C1   0.61
CHARGE   BZQ-1  C4   0.61
CHARGE   BZQ-1  O4  -0.75
CHARGE   BZQ-1  C2  -0.27 
CHARGE   BZQ-1  H2   0.09
CHARGE   BZQ-1  C3  -0.27
CHARGE   BZQ-1  H3   0.09
CHARGE   BZQ-1  C5  -0.27 
CHARGE   BZQ-1  H5   0.09
CHARGE   BZQ-1  C6  -0.27
CHARGE   BZQ-1  H6   0.09 

#ParaNam|Res  |Atom|Param/toggle
TRANS    BZQ          t


CONFLIST SO4        SO4BK SO4-2 SO4DM 

NATOM    SO4BK      0
NATOM    SO4-2      5
NATOM    SO4DM      0

IATOM    SO4-2  S   0
IATOM    SO4-2  O1  1
IATOM    SO4-2  O2  2
IATOM    SO4-2  O3  3
IATOM    SO4-2  O4  4

ATOMNAME SO4-2    0  S  
ATOMNAME SO4-2    1  O1 
ATOMNAME SO4-2    2  O2 
ATOMNAME SO4-2    3  O3 
ATOMNAME SO4-2    4  O4 


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
RXN      SO4-2      -62.2

#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  SO4-2  S   sp3        0    O1   0    O2   0    O3   0    O4
CONNECT  SO4-2  O1  sp3        0    S
CONNECT  SO4-2  O2  sp3        0    S
CONNECT  SO4-2  O3  sp3        0    S
CONNECT  SO4-2  O4  sp3        0    S

CHARGE   SO4-2  S     0.00
CHARGE   SO4-2  O1   -0.50
CHARGE   SO4-2  O2   -0.50
CHARGE   SO4-2  O3   -0.50
CHARGE   SO4-2  O4   -0.50

RADIUS   SO4    S     1.70
RADIUS   SO4    O     1.00


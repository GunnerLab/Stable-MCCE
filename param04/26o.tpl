#2,6-dimethoxy-Benzoquinone        
#Gennady Khirich 7.31.06  
CONFLIST 26O        26OBK 26O-1 26ODM

NATOM    26ODM      0
NATOM    26OBK      0
NATOM    26O-1      20

IATOM    26O-1  C6  0
IATOM    26O-1  C5  1
IATOM    26O-1  H5  2
IATOM    26O-1  H3  3
IATOM    26O-1  C4  4
IATOM    26O-1  C1  5
IATOM    26O-1  C2  6
IATOM    26O-1  C3  7
IATOM    26O-1  O1  8
IATOM    26O-1  O4  9
IATOM    26O-1  O6  10
IATOM    26O-1  CM6 11
IATOM    26O-1 2H1  12
IATOM    26O-1 2H2  13
IATOM    26O-1 2H3  14
IATOM    26O-1 6H1  15
IATOM    26O-1 6H2  16
IATOM    26O-1 6H3  17
IATOM    26O-1  O2  18
IATOM    26O-1  CM2 19

ATOMNAME 26O-1    0  C6 
ATOMNAME 26O-1    1  C5 
ATOMNAME 26O-1    2  H5 
ATOMNAME 26O-1    3  H3 
ATOMNAME 26O-1    4  C4 
ATOMNAME 26O-1    5  C1 
ATOMNAME 26O-1    6  C2 
ATOMNAME 26O-1    7  C3 
ATOMNAME 26O-1    8  O1 
ATOMNAME 26O-1    9  O4 
ATOMNAME 26O-1   10  O6 
ATOMNAME 26O-1   11  CM6 
ATOMNAME 26O-1   12 2H1
ATOMNAME 26O-1   13 2H2
ATOMNAME 26O-1   14 2H3
ATOMNAME 26O-1   15 6H1
ATOMNAME 26O-1   16 6H2
ATOMNAME 26O-1   17 6H3
ATOMNAME 26O-1   18  O2
ATOMNAME 26O-1   19  CM2


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   26O-1      0

PKA      26O-1      0.0

ELECTRON 26O-1      1

EM       26O-1      0.0

RXN      26O-1      -17.550


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  26O-1  C6  sp2        0    C5   0    C1   0    O6
CONNECT  26O-1  C5  sp2        0    C6   0    C4   0    H5
CONNECT  26O-1  H5  s          0    C5
CONNECT  26O-1  H3  s          0    C3
CONNECT  26O-1  C4  sp2        0    C5   0    C3   0    O4
CONNECT  26O-1  C1  sp2        0    C6   0    C2   0    O1
CONNECT  26O-1  C2  sp2        0    C1   0    C3   0    O2
CONNECT  26O-1  C3  sp2        0    C4   0    C2   0    H3
CONNECT  26O-1  CM6 sp3        0    O6   0   6H1   0   6H2   0   6H3  
CONNECT  26O-1  CM2 sp3        0    O2   0   2H1   0   2H2   0   2H3   
CONNECT  26O-1  O1  s          0    C1
CONNECT  26O-1  O4  s          0    C4
CONNECT  26O-1 2H1  s          0    CM2
CONNECT  26O-1 2H2  s          0    CM2
CONNECT  26O-1 2H3  s          0    CM2
CONNECT  26O-1 6H1  s          0    CM6
CONNECT  26O-1 6H2  s          0    CM6
CONNECT  26O-1 6H3  s          0    CM6
CONNECT  26O-1  O2  sp2        0    C2   0    CM2
CONNECT  26O-1  O6  sp2        0    C6   0    CM6


#de-protonated---------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   26O    C6  1.70
RADIUS   26O    C5  1.70
RADIUS   26O    H5  1.00
RADIUS   26O    H3  1.00
RADIUS   26O    C4  1.70
RADIUS   26O    C1  1.70
RADIUS   26O    C2  1.70
RADIUS   26O    C3  1.70
RADIUS   26O    O1  1.40
RADIUS   26O    O4  1.40
RADIUS   26O    CM6 1.70
RADIUS   26O    CM2 1.70
RADIUS   26O   2H1  1.00
RADIUS   26O   2H2  1.00
RADIUS   26O   2H3  1.00
RADIUS   26O   6H1  1.00
RADIUS   26O   6H2  1.00
RADIUS   26O   6H3  1.00
RADIUS   26O    O2  1.40
RADIUS   26O    O6  1.40

#-------|-----|----|----|----|----|----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
#        CONF  ATOM ATOM ATOM ATOM relx Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma
TORSION  26O    C1   C6   O6   CM6  f      -7.730        0      0.00     14.029        1      0.00      6.384        2    360.00      5.300        3    720.00      0.784        4    720.00      1.329        5    900.00

#NEUTRAL------
#23456789A123456789B123456789C
CHARGE   26O-1  O1  -0.60
CHARGE   26O-1  C1   0.27
CHARGE   26O-1  C4   0.61
CHARGE   26O-1  O4  -0.72
CHARGE   26O-1  C2   0.25
CHARGE   26O-1  CM2  0.20
CHARGE   26O-1  C3  -0.33
CHARGE   26O-1  CM6  0.20
CHARGE   26O-1  C5  -0.30
CHARGE   26O-1  H3   0.00
CHARGE   26O-1  C6   0.24
CHARGE   26O-1  H5   0.00
CHARGE   26O-1 2H1  -0.00
CHARGE   26O-1 2H2  -0.00
CHARGE   26O-1 2H3  -0.00
CHARGE   26O-1 6H1  -0.00
CHARGE   26O-1 6H2   0.00
CHARGE   26O-1 6H3  -0.00
CHARGE   26O-1  O2  -0.39
CHARGE   26O-1  O6  -0.43



#ParaNam|Res  |Atom|Param/toggle
TRANS    26O          t

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  26O   0     C2 - O2   CM2
ROTAMER  26O   1     C6 - O6   CM6
#=========================================================================


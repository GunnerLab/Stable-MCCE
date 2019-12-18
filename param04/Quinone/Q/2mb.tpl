#2-methyl-Benzoquinone        
#Gennady Khirich 7.31.06  
CONFLIST 2MB        2MBBK 2MB01 2MBDM

NATOM    2MBDM      0
NATOM    2MBBK      0
NATOM    2MB01      15

IATOM    2MB01  C6  0
IATOM    2MB01  C5  1
IATOM    2MB01  H6  2
IATOM    2MB01  H5  3
IATOM    2MB01  C4  4
IATOM    2MB01  C1  5
IATOM    2MB01  C2  6
IATOM    2MB01  C3  7
IATOM    2MB01  H3  8
IATOM    2MB01  O1  9
IATOM    2MB01  O4  10
IATOM    2MB01  CM2 11
IATOM    2MB01 2H1  12
IATOM    2MB01 2H2  13
IATOM    2MB01 2H3  14

ATOMNAME 2MB01    0  C6 
ATOMNAME 2MB01    1  C5 
ATOMNAME 2MB01    2  H6 
ATOMNAME 2MB01    3  H5 
ATOMNAME 2MB01    4  C4 
ATOMNAME 2MB01    5  C1 
ATOMNAME 2MB01    6  C2 
ATOMNAME 2MB01    7  C3 
ATOMNAME 2MB01    8  H3 
ATOMNAME 2MB01    9  O1 
ATOMNAME 2MB01   10  O4 
ATOMNAME 2MB01   11  CM2 
ATOMNAME 2MB01   12 2H1
ATOMNAME 2MB01   13 2H2
ATOMNAME 2MB01   14 2H3

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   2MB01      0

PKA      2MB01      0.0

ELECTRON 2MB01      0

EM       2MB01      0.0

RXN      2MB01      -2.293


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  2MB01  C6  sp2        0    C5   0    C1   0    H6
CONNECT  2MB01  C5  sp2        0    C6   0    C4   0    H5
CONNECT  2MB01  H6  s          0    C6
CONNECT  2MB01  H5  s          0    C5
CONNECT  2MB01  C4  sp2        0    C5   0    C3   0    O4
CONNECT  2MB01  C1  sp2        0    C6   0    C2   0    O1
CONNECT  2MB01  C2  sp2        0    C1   0    C3   0    CM2
CONNECT  2MB01  C3  sp2        0    C4   0    C2   0    H3
CONNECT  2MB01  H3  s          0    C3    
CONNECT  2MB01  CM2 sp3        0    C2   0   2H1   0   2H2   0   2H3   
CONNECT  2MB01  O1  s          0    C1
CONNECT  2MB01  O4  s          0    C4
CONNECT  2MB01 2H1  s          0    CM2
CONNECT  2MB01 2H2  s          0    CM2
CONNECT  2MB01 2H3  s          0    CM2

#de-protonated---------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   2MB    C6  1.70
RADIUS   2MB    C5  1.70
RADIUS   2MB    H6  1.00
RADIUS   2MB    H5  1.00
RADIUS   2MB    C4  1.70
RADIUS   2MB    C1  1.70
RADIUS   2MB    C2  1.70
RADIUS   2MB    C3  1.70
RADIUS   2MB    O1  1.40
RADIUS   2MB    O4  1.40
RADIUS   2MB    H3  1.00
RADIUS   2MB    CM2 1.70
RADIUS   2MB   2H1  1.00
RADIUS   2MB   2H2  1.00
RADIUS   2MB   2H3  1.00


#NEUTRAL------
#23456789A123456789B123456789C
# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   gk 7.31.06
CHARGE   2MB01  O1  -0.49
CHARGE   2MB01  C1   0.56
CHARGE   2MB01  C4   0.69
CHARGE   2MB01  O4  -0.53
CHARGE   2MB01  C2   0.05
CHARGE   2MB01  CM2 -0.10
CHARGE   2MB01  C3  -0.37
CHARGE   2MB01  H3   0.17
CHARGE   2MB01  C5  -0.23
CHARGE   2MB01  H5   0.15
CHARGE   2MB01  C6  -0.19
CHARGE   2MB01  H6   0.14
CHARGE   2MB01 2H1   0.05
CHARGE   2MB01 2H2   0.05
CHARGE   2MB01 2H3   0.05
 
#TORSION  2MB    HO2  O2   C2   C3   f      1.800         2    180.00

TORSION  2MB01  2H1  CM2  C2   C1   f     -0.180         3      0.00


#ParaNam|Res  |Atom|Param/toggle
TRANS    2MB          t


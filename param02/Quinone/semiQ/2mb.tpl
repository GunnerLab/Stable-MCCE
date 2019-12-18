#2-methyl-Benzoquinone (Q-)       
#Agnes 7.9.08  
CONFLIST 2MB        2MBBK 2MB-1 2MBDM

NATOM    2MBDM      0
NATOM    2MBBK      0
NATOM    2MB-1      15

IATOM    2MB-1  C6  0
IATOM    2MB-1  C5  1
IATOM    2MB-1  H6  2
IATOM    2MB-1  H5  3
IATOM    2MB-1  C4  4
IATOM    2MB-1  C1  5
IATOM    2MB-1  C2  6
IATOM    2MB-1  C3  7
IATOM    2MB-1  H3  8
IATOM    2MB-1  O1  9
IATOM    2MB-1  O4  10
IATOM    2MB-1  CM2 11
IATOM    2MB-1 2H1  12
IATOM    2MB-1 2H2  13
IATOM    2MB-1 2H3  14


ATOMNAME 2MB-1    0  C6
ATOMNAME 2MB-1    1  C5
ATOMNAME 2MB-1    2  H6
ATOMNAME 2MB-1    3  H5
ATOMNAME 2MB-1    4  C4
ATOMNAME 2MB-1    5  C1
ATOMNAME 2MB-1    6  C2
ATOMNAME 2MB-1    7  C3
ATOMNAME 2MB-1    8  H3
ATOMNAME 2MB-1    9  O1
ATOMNAME 2MB-1   10  O4
ATOMNAME 2MB-1   11  CM2
ATOMNAME 2MB-1   12 2H1
ATOMNAME 2MB-1   13 2H2
ATOMNAME 2MB-1   14 2H3

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   2MB-1      0
PKA      2MB-1      0.0
ELECTRON 2MB-1      1
EM       2MB-1      23
RXN      2MB-1      -14.86

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  2MB-1  C6  sp2        0    C5   0    C1   0    H6
CONNECT  2MB-1  C5  sp2        0    C6   0    C4   0    H5
CONNECT  2MB-1  H6  s          0    C6
CONNECT  2MB-1  H5  s          0    C5
CONNECT  2MB-1  C4  sp2        0    C5   0    C3   0    O4
CONNECT  2MB-1  C1  sp2        0    C6   0    C2   0    O1
CONNECT  2MB-1  C2  sp2        0    C1   0    C3   0    CM2
CONNECT  2MB-1  C3  sp2        0    C4   0    C2   0    H3
CONNECT  2MB-1  H3  s          0    C3
CONNECT  2MB-1  CM2 sp3        0    C2   0   2H1   0   2H2   0   2H3
CONNECT  2MB-1  O1  s          0    C1
CONNECT  2MB-1  O4  s          0    C4
CONNECT  2MB-1 2H1  s          0    CM2
CONNECT  2MB-1 2H2  s          0    CM2
CONNECT  2MB-1 2H3  s          0    CM2


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

# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   Agnes 7.9.08
CHARGE   2MB-1  O1  -0.70
CHARGE   2MB-1  C1   0.50
CHARGE   2MB-1  C4   0.64
CHARGE   2MB-1  O4  -0.74
CHARGE   2MB-1  C2   0.01
CHARGE   2MB-1  CM2 -0.02 
CHARGE   2MB-1  C3  -0.43
CHARGE   2MB-1  H3   0.13
CHARGE   2MB-1  C5  -0.29
CHARGE   2MB-1  H5   0.09
CHARGE   2MB-1  C6  -0.26
CHARGE   2MB-1  H6   0.09
CHARGE   2MB-1 2H1   0.00
CHARGE   2MB-1 2H2   0.00
CHARGE   2MB-1 2H3  -0.02 
 
#TORSION  2MB    HO2  O2   C2   C3   f      1.800         2    180.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    2MB          t


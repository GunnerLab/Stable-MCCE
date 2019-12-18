#2,6-dimethoxy-Benzoquinone        
#Gennady Khirich 7.31.06  
CONFLIST 26O        26OBK 26O01 26ODM

NATOM    26ODM      0
NATOM    26OBK      0
NATOM    26O01      20

IATOM    26O01  C6  0
IATOM    26O01  C5  1
IATOM    26O01  H5  2
IATOM    26O01  H3  3
IATOM    26O01  C4  4
IATOM    26O01  C1  5
IATOM    26O01  C2  6
IATOM    26O01  C3  7
IATOM    26O01  CM6 8
IATOM    26O01  O1  9
IATOM    26O01  O4  10
IATOM    26O01  CM2 11
IATOM    26O01 2H1  12
IATOM    26O01 2H2  13
IATOM    26O01 2H3  14
IATOM    26O01 6H1  15
IATOM    26O01 6H2  16
IATOM    26O01 6H3  17
IATOM    26O01  O2  18
IATOM    26O01  O6  19

ATOMNAME 26O01    0  C6 
ATOMNAME 26O01    1  C5 
ATOMNAME 26O01    2  H5 
ATOMNAME 26O01    3  H3 
ATOMNAME 26O01    4  C4 
ATOMNAME 26O01    5  C1 
ATOMNAME 26O01    6  C2 
ATOMNAME 26O01    7  C3 
ATOMNAME 26O01    8  CM6 
ATOMNAME 26O01    9  O1 
ATOMNAME 26O01   10  O4 
ATOMNAME 26O01   11  CM2 
ATOMNAME 26O01   12 2H1
ATOMNAME 26O01   13 2H2
ATOMNAME 26O01   14 2H3
ATOMNAME 26O01   15 6H1
ATOMNAME 26O01   16 6H2
ATOMNAME 26O01   17 6H3
ATOMNAME 26O01   18  O2
ATOMNAME 26O01   19  O6


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   26O01      0

PKA      26O01      0.0

ELECTRON 26O01      0

EM       26O01      0.0

RXN      26O01      -2.859


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  26O01  C6  sp2        0    C5   0    C1   0    O6
CONNECT  26O01  C5  sp2        0    C6   0    C4   0    H5
CONNECT  26O01  H5  s          0    C5
CONNECT  26O01  H3  s          0    C3
CONNECT  26O01  C4  sp2        0    C5   0    C3   0    O4
CONNECT  26O01  C1  sp2        0    C6   0    C2   0    O1
CONNECT  26O01  C2  sp2        0    C1   0    C3   0    O2
CONNECT  26O01  C3  sp2        0    C4   0    C2   0    H3
CONNECT  26O01  CM6 sp3        0    O6   0   6H1   0   6H2   0   6H3  
CONNECT  26O01  CM2 sp3        0    O2   0   2H1   0   2H2   0   2H3   
CONNECT  26O01  O1  s          0    C1
CONNECT  26O01  O4  s          0    C4
CONNECT  26O01 2H1  s          0    CM2
CONNECT  26O01 2H2  s          0    CM2
CONNECT  26O01 2H3  s          0    CM2
CONNECT  26O01 6H1  s          0    CM6
CONNECT  26O01 6H2  s          0    CM6
CONNECT  26O01 6H3  s          0    CM6
CONNECT  26O01  O2  sp2        0    C2   0    CM2
CONNECT  26O01  O6  sp2        0    C6   0    CM6


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



#NEUTRAL------
#23456789A123456789B123456789C
# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   gk 7.31.06
CHARGE   26O01  O1  -0.43
CHARGE   26O01  C1   0.42
CHARGE   26O01  C4   0.67
CHARGE   26O01  O4  -0.53
CHARGE   26O01  C2   0.26
CHARGE   26O01  CM2  0.19
CHARGE   26O01  C3  -0.44
CHARGE   26O01  CM6  0.08
CHARGE   26O01  C5  -0.42
CHARGE   26O01  H3   0.18
CHARGE   26O01  C6   0.24
CHARGE   26O01  H5   0.19
CHARGE   26O01 2H1   0.03
CHARGE   26O01 2H2   0.03
CHARGE   26O01 2H3   0.03
CHARGE   26O01 6H1   0.05
CHARGE   26O01 6H2   0.05
CHARGE   26O01 6H3   0.05
CHARGE   26O01  O2  -0.34
CHARGE   26O01  O6  -0.31
 
#TORSION  26O    HO2  O2   C2   C3   f      1.800         2    180.00
TORSION  26O01  CM2  O2   C2   C3   f      0.0           2    180.00
TORSION  26O01  CM2  O2   C2   C1   f      0.215         6    180.00     -0.195        4    180.00    -0.427         3    180.00     4.208         2    180.00    -0.632         1    180.00
TORSION  26O01  2H1  CM2  O2   C2   f     -0.180         3      0.00
TORSION  26O01  CM6  O6   C6   C5   f      0.0           2    180.00
TORSION  26O01  CM6  O6   C6   C1   f      0.215         6    180.00     -0.195        4    180.00    -0.427         3    180.00     4.208         2    180.00    -0.632         1    180.00
TORSION  26O01  6H1  CM6  O6   C6   f     -0.180         3      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    26O          t

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  26O   0     C2 - O2   CM2
ROTAMER  26O   1     C6 - O6   CM6
#=========================================================================




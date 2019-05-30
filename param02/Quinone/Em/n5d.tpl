#5,8-dimethyl,1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST N5D        N5DBK N5D01 N5D-1 N5DDM

NATOM    N5DDM      0
NATOM    N5DBK      0
NATOM    N5D01      24
NATOM    N5D-1      24

IATOM    N5D01  C1  0
IATOM    N5D01  O1  1
IATOM    N5D01  C2  2
IATOM    N5D01  H2  3
IATOM    N5D01  C3  4
IATOM    N5D01  H3  5
IATOM    N5D01  C4  6
IATOM    N5D01  O4  7
IATOM    N5D01  C5  8
IATOM    N5D01  C5M 9
IATOM    N5D01  C6  10
IATOM    N5D01  H6  11
IATOM    N5D01  C7  12
IATOM    N5D01  H7  13
IATOM    N5D01  C8  14
IATOM    N5D01  C8M 15
IATOM    N5D01  C9  16
IATOM    N5D01  C10 17
IATOM    N5D01 1H5M 18
IATOM    N5D01 2H5M 19
IATOM    N5D01 3H5M 20
IATOM    N5D01 1H8M 21
IATOM    N5D01 2H8M 22
IATOM    N5D01 3H8M 23
IATOM    N5D-1  C1  0
IATOM    N5D-1  O1  1
IATOM    N5D-1  C2  2
IATOM    N5D-1  H2  3
IATOM    N5D-1  C3  4
IATOM    N5D-1  H3  5
IATOM    N5D-1  C4  6
IATOM    N5D-1  O4  7
IATOM    N5D-1  C5  8
IATOM    N5D-1  C5M 9
IATOM    N5D-1  C6  10
IATOM    N5D-1  H6  11
IATOM    N5D-1  C7  12
IATOM    N5D-1  H7  13
IATOM    N5D-1  C8  14
IATOM    N5D-1  C8M 15
IATOM    N5D-1  C9  16
IATOM    N5D-1  C10 17
IATOM    N5D-1 1H5M 18
IATOM    N5D-1 2H5M 19
IATOM    N5D-1 3H5M 20
IATOM    N5D-1 1H8M 21
IATOM    N5D-1 2H8M 22
IATOM    N5D-1 3H8M 23


ATOMNAME N5D01    0  C1
ATOMNAME N5D01    1  O1 
ATOMNAME N5D01    2  C2 
ATOMNAME N5D01    3  H2 
ATOMNAME N5D01    4  C3 
ATOMNAME N5D01    5  H3 
ATOMNAME N5D01    6  C4 
ATOMNAME N5D01    7  O4 
ATOMNAME N5D01    8  C5 
ATOMNAME N5D01    9  C5M 
ATOMNAME N5D01   10  C6 
ATOMNAME N5D01   11  H6 
ATOMNAME N5D01   12  C7 
ATOMNAME N5D01   13  H7 
ATOMNAME N5D01   14  C8 
ATOMNAME N5D01   15  C8M
ATOMNAME N5D01   16  C9 
ATOMNAME N5D01   17  C10 
ATOMNAME N5D01   18 1H5M
ATOMNAME N5D01   19 2H5M
ATOMNAME N5D01   20 3H5M 
ATOMNAME N5D01   21 1H8M
ATOMNAME N5D01   22 2H8M
ATOMNAME N5D01   23 3H8M 
ATOMNAME N5D-1    0  C1
ATOMNAME N5D-1    1  O1
ATOMNAME N5D-1    2  C2
ATOMNAME N5D-1    3  H2
ATOMNAME N5D-1    4  C3
ATOMNAME N5D-1    5  H3
ATOMNAME N5D-1    6  C4
ATOMNAME N5D-1    7  O4
ATOMNAME N5D-1    8  C5
ATOMNAME N5D-1    9  C5M
ATOMNAME N5D-1   10  C6
ATOMNAME N5D-1   11  H6
ATOMNAME N5D-1   12  C7
ATOMNAME N5D-1   13  H7
ATOMNAME N5D-1   14  C8
ATOMNAME N5D-1   15  C8M
ATOMNAME N5D-1   16  C9
ATOMNAME N5D-1   17  C10
ATOMNAME N5D-1   18 1H5M
ATOMNAME N5D-1   19 2H5M
ATOMNAME N5D-1   20 3H5M
ATOMNAME N5D-1   21 1H8M
ATOMNAME N5D-1   22 2H8M
ATOMNAME N5D-1   23 3H8M

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   N5D01      0

PKA      N5D01      0.0

ELECTRON N5D01      0

EM       N5D01      0.0

RXN      N5D01      -3.605

PROTON   N5D-1      0

PKA      N5D-1      0.0

ELECTRON N5D-1      1

EM       N5D-1      -219

RXN      N5D-1      -15.965

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  N5D01  C1  sp2       0     O1  0     C2  0     C9
CONNECT  N5D01  O1  s         0     C1
CONNECT  N5D01  C2  sp2       0     H2  0     C1  0     C3
CONNECT  N5D01  H2  s         0     C2
CONNECT  N5D01  C3  sp2       0     C2  0     H3  0     C4
CONNECT  N5D01  H3  s         0     C3
CONNECT  N5D01  C4  sp2       0     C3  0     O4  0     C10
CONNECT  N5D01  O4  s         0     C4
CONNECT  N5D01  C5  sp2       0     C10 0     C5M 0     C6
CONNECT  N5D01  C5M sp3       0     C5  0    1H5M 0    2H5M 0    3H5M
CONNECT  N5D01  C6  sp2       0     C5  0     C7  0     H6
CONNECT  N5D01  H6  s         0     C6
CONNECT  N5D01  C7  sp2       0     C8  0     C6  0     H7
CONNECT  N5D01  H7  s         0     C7
CONNECT  N5D01  C8  sp2       0     C7  0     C9  0     C8M
CONNECT  N5D01  C8M sp3       0     C8  0    1H8M 0    2H8M 0    3H8M
CONNECT  N5D01  C9  sp2       0     C1  0     C8  0     C10
CONNECT  N5D01  C10 sp2       0     C4  0     C5  0     C9
CONNECT  N5D01 1H5M s         0     C5M
CONNECT  N5D01 2H5M s         0     C5M
CONNECT  N5D01 3H5M s         0     C5M
CONNECT  N5D01 1H8M s         0     C8M
CONNECT  N5D01 2H8M s         0     C8M
CONNECT  N5D01 3H8M s         0     C8M
CONNECT  N5D-1  C1  sp2       0     O1  0     C2  0     C9
CONNECT  N5D-1  O1  s         0     C1
CONNECT  N5D-1  C2  sp2       0     H2  0     C1  0     C3
CONNECT  N5D-1  H2  s         0     C2
CONNECT  N5D-1  C3  sp2       0     C2  0     H3  0     C4
CONNECT  N5D-1  H3  s         0     C3
CONNECT  N5D-1  C4  sp2       0     C3  0     O4  0     C10
CONNECT  N5D-1  O4  s         0     C4
CONNECT  N5D-1  C5  sp2       0     C10 0     C5M 0     C6
CONNECT  N5D-1  C5M sp3       0     C5  0    1H5M 0    2H5M 0    3H5M
CONNECT  N5D-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  N5D-1  H6  s         0     C6
CONNECT  N5D-1  C7  sp2       0     C8  0     C6  0     H7
CONNECT  N5D-1  H7  s         0     C7
CONNECT  N5D-1  C8  sp2       0     C7  0     C9  0     C8M
CONNECT  N5D-1  C8M sp3       0     C8  0    1H8M 0    2H8M 0    3H8M
CONNECT  N5D-1  C9  sp2       0     C1  0     C8  0     C10
CONNECT  N5D-1  C10 sp2       0     C4  0     C5  0     C9
CONNECT  N5D-1 1H5M s         0     C5M
CONNECT  N5D-1 2H5M s         0     C5M
CONNECT  N5D-1 3H5M s         0     C5M
CONNECT  N5D-1 1H8M s         0     C8M
CONNECT  N5D-1 2H8M s         0     C8M
CONNECT  N5D-1 3H8M s         0     C8M


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   N5D    C1  1.70
RADIUS   N5D    O1  1.40
RADIUS   N5D    C2  1.70
RADIUS   N5D    H2  1.00
RADIUS   N5D    C3  1.70
RADIUS   N5D    H3  1.00
RADIUS   N5D    C4  1.70
RADIUS   N5D    O4  1.40
RADIUS   N5D    C5  1.70
RADIUS   N5D    C5M 1.70
RADIUS   N5D    C6  1.70
RADIUS   N5D    H6  1.00
RADIUS   N5D    C7  1.70
RADIUS   N5D    H7  1.00
RADIUS   N5D    C8  1.70
RADIUS   N5D    C8M 1.70
RADIUS   N5D    C9  1.70
RADIUS   N5D    C10 1.70
RADIUS   N5D   1H5M 1.00
RADIUS   N5D   2H5M 1.00
RADIUS   N5D   3H5M 1.00
RADIUS   N5D   1H8M 1.00
RADIUS   N5D   2H8M 1.00
RADIUS   N5D   3H8M 1.00

#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600) Agnes 10/09
CHARGE   N5D01  C1   0.57
CHARGE   N5D01  O1  -0.50
CHARGE   N5D01  C2  -0.18
CHARGE   N5D01  H2   0.12
CHARGE   N5D01  C3  -0.18
CHARGE   N5D01  H3   0.13
CHARGE   N5D01  C4   0.58
CHARGE   N5D01  O4  -0.52
CHARGE   N5D01  C5   0.21
CHARGE   N5D01  C5M -0.20
CHARGE   N5D01  C6  -0.15
CHARGE   N5D01  H6   0.12
CHARGE   N5D01  C7  -0.17
CHARGE   N5D01  H7   0.12 
CHARGE   N5D01  C8   0.21
CHARGE   N5D01  C8M -0.24
CHARGE   N5D01  C9  -0.18
CHARGE   N5D01  C10 -0.18
CHARGE   N5D01 1H5M  0.04
CHARGE   N5D01 2H5M  0.04
CHARGE   N5D01 3H5M  0.12
CHARGE   N5D01 1H8M  0.09
CHARGE   N5D01 2H8M  0.09
CHARGE   N5D01 3H8M  0.06

CHARGE   N5D-1  C1   0.49
CHARGE   N5D-1  O1  -0.66
CHARGE   N5D-1  C2  -0.30
CHARGE   N5D-1  H2   0.12
CHARGE   N5D-1  C3  -0.27
CHARGE   N5D-1  H3   0.12
CHARGE   N5D-1  C4   0.48
CHARGE   N5D-1  O4  -0.68
CHARGE   N5D-1  C5   0.22
CHARGE   N5D-1  C5M -0.28
CHARGE   N5D-1  C6  -0.27
CHARGE   N5D-1  H6   0.12
CHARGE   N5D-1  C7  -0.25
CHARGE   N5D-1  H7   0.11
CHARGE   N5D-1  C8   0.21
CHARGE   N5D-1  C8M -0.37
CHARGE   N5D-1  C9  -0.14
CHARGE   N5D-1  C10 -0.15
CHARGE   N5D-1 1H5M  0.12
CHARGE   N5D-1 2H5M  0.12
CHARGE   N5D-1 3H5M  0.05
CHARGE   N5D-1 1H8M  0.03
CHARGE   N5D-1 2H8M  0.03
CHARGE   N5D-1 3H8M  0.13


#TORSION  N5D    O2   C2   C1   O1   f        0.0         1    180.00  
#TORSION  N5D    C8   C9   C1   O1   f        0.0         1    180.00  

#TORSION  N5D    C7   C6   C5   C10  f        0.0         1      0.00 
#TORSION  N5D    H5   C5   C10  C4   f        0.0         1      0.00
#TORSION  N5D    C6   C7   C8   C9   f        0.0         1      0.00

#TORSION  N5D    C7   C8   C9   C1   f        0.0         1      0.00 
#TORSION  N5D    C5   C6   C7   C8   f        0.0         1      0.00
#TORSION  N5D    H5   C5   C6   C7   f        0.0         1    180.00
#TORSION  N5D    C6   C5   C10  C9   f        0.0         1    180.00 
#TORSION  N5D    C6   C5   C10  C4   f        0.0         1      0.00 

#TORSION  N5D    H5   C5   C10  C4   f        0.0         1      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    N5D          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     N5D   0     C9 - C10- C1

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  N5D   0     C1 - C4   WHOLE_CONF
ROTAMER  N5D   1     C9 - C10  WHOLE_CONF


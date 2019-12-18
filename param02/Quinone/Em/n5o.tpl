#5-methoxy, 1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST N5O        N5OBK N5O01 N5O-1 N5ODM

NATOM    N5ODM      0
NATOM    N5OBK      0
NATOM    N5O01      22
NATOM    N5O-1      22

IATOM    N5O01  C1  0
IATOM    N5O01  O1  1
IATOM    N5O01  C2  2
IATOM    N5O01  H2  3
IATOM    N5O01  C3  4
IATOM    N5O01  H3  5
IATOM    N5O01  C4  6
IATOM    N5O01  O4  7
IATOM    N5O01  C5  8
IATOM    N5O01  C5M 9
IATOM    N5O01  C6  10
IATOM    N5O01  H6  11
IATOM    N5O01  C7  12
IATOM    N5O01  H7  13
IATOM    N5O01  C8  14
IATOM    N5O01  H8  15
IATOM    N5O01  C9  16
IATOM    N5O01  C10 17
IATOM    N5O01 1H5M 18
IATOM    N5O01 2H5M 19
IATOM    N5O01 3H5M 20
IATOM    N5O01  O5  21
IATOM    N5O-1  C1  0
IATOM    N5O-1  O1  1
IATOM    N5O-1  C2  2
IATOM    N5O-1  H2  3
IATOM    N5O-1  C3  4
IATOM    N5O-1  H3  5
IATOM    N5O-1  C4  6
IATOM    N5O-1  O4  7
IATOM    N5O-1  C5  8
IATOM    N5O-1  C5M 9
IATOM    N5O-1  C6  10
IATOM    N5O-1  H6  11
IATOM    N5O-1  C7  12
IATOM    N5O-1  H7  13
IATOM    N5O-1  C8  14
IATOM    N5O-1  H8  15
IATOM    N5O-1  C9  16
IATOM    N5O-1  C10 17
IATOM    N5O-1 1H5M 18
IATOM    N5O-1 2H5M 19
IATOM    N5O-1 3H5M 20
IATOM    N5O-1  O5  21


ATOMNAME N5O01    0  C1
ATOMNAME N5O01    1  O1 
ATOMNAME N5O01    2  C2 
ATOMNAME N5O01    3  H2 
ATOMNAME N5O01    4  C3 
ATOMNAME N5O01    5  H3 
ATOMNAME N5O01    6  C4 
ATOMNAME N5O01    7  O4 
ATOMNAME N5O01    8  C5 
ATOMNAME N5O01    9  C5M 
ATOMNAME N5O01   10  C6 
ATOMNAME N5O01   11  H6 
ATOMNAME N5O01   12  C7 
ATOMNAME N5O01   13  H7 
ATOMNAME N5O01   14  C8 
ATOMNAME N5O01   15  H8 
ATOMNAME N5O01   16  C9 
ATOMNAME N5O01   17  C10 
ATOMNAME N5O01   18 1H5M
ATOMNAME N5O01   19 2H5M
ATOMNAME N5O01   20 3H5M 
ATOMNAME N5O01   21  O5
ATOMNAME N5O-1    0  C1
ATOMNAME N5O-1    1  O1
ATOMNAME N5O-1    2  C2
ATOMNAME N5O-1    3  H2
ATOMNAME N5O-1    4  C3
ATOMNAME N5O-1    5  H3
ATOMNAME N5O-1    6  C4
ATOMNAME N5O-1    7  O4
ATOMNAME N5O-1    8  C5
ATOMNAME N5O-1    9  C5M
ATOMNAME N5O-1   10  C6
ATOMNAME N5O-1   11  H6
ATOMNAME N5O-1   12  C7
ATOMNAME N5O-1   13  H7
ATOMNAME N5O-1   14  C8
ATOMNAME N5O-1   15  H8
ATOMNAME N5O-1   16  C9
ATOMNAME N5O-1   17  C10
ATOMNAME N5O-1   18 1H5M
ATOMNAME N5O-1   19 2H5M
ATOMNAME N5O-1   20 3H5M
ATOMNAME N5O-1   21  O5

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   N5O01      0

PKA      N5O01      0.0

ELECTRON N5O01      0

EM       N5O01      0.0

RXN      N5O01      -5.838

PROTON   N5O-1      0

PKA      N5O-1      0.0

ELECTRON N5O-1      1

EM       N5O-1      -183

RXN      N5O-1      -18.442


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  N5O01  C1  sp2       0     O1  0     C2  0     C9
CONNECT  N5O01  O1  s         0     C1
CONNECT  N5O01  C2  sp2       0     H2  0     C1  0     C3
CONNECT  N5O01  H2  s         0     C2
CONNECT  N5O01  C3  sp2       0     C2  0     H3  0     C4
CONNECT  N5O01  H3  s         0     C3
CONNECT  N5O01  C4  sp2       0     C3  0     O4  0     C10
CONNECT  N5O01  O4  s         0     C4
CONNECT  N5O01  C5  sp2       0     C10 0     O5  0     C6
CONNECT  N5O01  C5M sp3       0     O5  0    1H5M 0    2H5M 0    3H5M
CONNECT  N5O01  C6  sp2       0     C5  0     C7  0     H6
CONNECT  N5O01  H6  s         0     C6
CONNECT  N5O01  C7  sp2       0     C8  0     C6  0     H7
CONNECT  N5O01  H7  s         0     C7
CONNECT  N5O01  C8  sp2       0     C7  0     C9  0     H8
CONNECT  N5O01  H8  s         0     C8
CONNECT  N5O01  C9  sp2       0     C1  0     C8  0     C10
CONNECT  N5O01  C10 sp2       0     C4  0     C5  0     C9
CONNECT  N5O01 1H5M s         0     C5M
CONNECT  N5O01 2H5M s         0     C5M
CONNECT  N5O01 3H5M s         0     C5M
CONNECT  N5O01  O5  sp2       0     C5M 0     C5

CONNECT  N5O-1  C1  sp2       0     O1  0     C2  0     C9
CONNECT  N5O-1  O1  s         0     C1
CONNECT  N5O-1  C2  sp2       0     H2  0     C1  0     C3
CONNECT  N5O-1  H2  s         0     C2
CONNECT  N5O-1  C3  sp2       0     C2  0     H3  0     C4
CONNECT  N5O-1  H3  s         0     C3
CONNECT  N5O-1  C4  sp2       0     C3  0     O4  0     C10
CONNECT  N5O-1  O4  s         0     C4
CONNECT  N5O-1  C5  sp2       0     C10 0     O5  0     C6
CONNECT  N5O-1  C5M sp3       0     O5  0    1H5M 0    2H5M 0    3H5M
CONNECT  N5O-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  N5O-1  H6  s         0     C6
CONNECT  N5O-1  C7  sp2       0     C8  0     C6  0     H7
CONNECT  N5O-1  H7  s         0     C7
CONNECT  N5O-1  C8  sp2       0     C7  0     C9  0     H8
CONNECT  N5O-1  H8  s         0     C8
CONNECT  N5O-1  C9  sp2       0     C1  0     C8  0     C10
CONNECT  N5O-1  C10 sp2       0     C4  0     C5  0     C9
CONNECT  N5O-1 1H5M s         0     C5M
CONNECT  N5O-1 2H5M s         0     C5M
CONNECT  N5O-1 3H5M s         0     C5M
CONNECT  N5O-1  O5  sp2       0     C5M 0     C5


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   N5O    C1  1.70
RADIUS   N5O    O1  1.40
RADIUS   N5O    C2  1.70
RADIUS   N5O    H2  1.00
RADIUS   N5O    C3  1.70
RADIUS   N5O    H3  1.00
RADIUS   N5O    C4  1.70
RADIUS   N5O    O4  1.40
RADIUS   N5O    C5  1.70
RADIUS   N5O    C5M 1.70
RADIUS   N5O    C6  1.70
RADIUS   N5O    H6  1.00
RADIUS   N5O    C7  1.70
RADIUS   N5O    H7  1.00
RADIUS   N5O    C8  1.70
RADIUS   N5O    H8  1.00
RADIUS   N5O    C9  1.70
RADIUS   N5O    C10 1.70
RADIUS   N5O   1H5M 1.00
RADIUS   N5O   2H5M 1.00
RADIUS   N5O   3H5M 1.00
RADIUS   N5O    O5  1.40

#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   N5O01  C1   0.56
CHARGE   N5O01  O1  -0.53
CHARGE   N5O01  C2  -0.18
CHARGE   N5O01  H2   0.15
CHARGE   N5O01  C3  -0.26
CHARGE   N5O01  H3   0.16
CHARGE   N5O01  C4   0.65
CHARGE   N5O01  O4  -0.53
CHARGE   N5O01  C5   0.56
CHARGE   N5O01  C5M  0.19
CHARGE   N5O01  C6  -0.28
CHARGE   N5O01  H6   0.16
CHARGE   N5O01  C7  -0.04
CHARGE   N5O01  H7   0.12
CHARGE   N5O01  C8  -0.14
CHARGE   N5O01  H8   0.14
CHARGE   N5O01  C9  -0.02
CHARGE   N5O01  C10 -0.30
CHARGE   N5O01 1H5M  0.00
CHARGE   N5O01 2H5M  0.06
CHARGE   N5O01 3H5M  0.05
CHARGE   N5O01  O5  -0.52

CHARGE   N5O-1  C1   0.41
CHARGE   N5O-1  O1  -0.66
CHARGE   N5O-1  C2  -0.24
CHARGE   N5O-1  H2   0.12
CHARGE   N5O-1  C3  -0.34
CHARGE   N5O-1  H3   0.13
CHARGE   N5O-1  C4   0.53
CHARGE   N5O-1  O4  -0.64
CHARGE   N5O-1  C5   0.51
CHARGE   N5O-1  C5M  0.32
CHARGE   N5O-1  C6  -0.28
CHARGE   N5O-1  H6   0.13
CHARGE   N5O-1  C7  -0.17
CHARGE   N5O-1  H7   0.10
CHARGE   N5O-1  C8  -0.15
CHARGE   N5O-1  H8   0.11
CHARGE   N5O-1  C9   0.04
CHARGE   N5O-1  C10 -0.28
CHARGE   N5O-1 1H5M -0.05
CHARGE   N5O-1 2H5M  0.01
CHARGE   N5O-1 3H5M  0.00
CHARGE   N5O-1  O5  -0.60

#TORSION  N5O    O2   C2   C1   O1   f        0.0         1    180.00  
#TORSION  N5O    C8   C9   C1   O1   f        0.0         1    180.00  

#TORSION  N5O    C7   C6   C5   C10  f        0.0         1      0.00 
#TORSION  N5O    H5   C5   C10  C4   f        0.0         1      0.00
#TORSION  N5O    C6   C7   C8   C9   f        0.0         1      0.00

#TORSION  N5O    C7   C8   C9   C1   f        0.0         1      0.00 
#TORSION  N5O    C5   C6   C7   C8   f        0.0         1      0.00
#TORSION  N5O    H5   C5   C6   C7   f        0.0         1    180.00
#TORSION  N5O    C6   C5   C10  C9   f        0.0         1    180.00 
#TORSION  N5O    C6   C5   C10  C4   f        0.0         1      0.00 

#TORSION  N5O    H5   C5   C10  C4   f        0.0         1      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    N5O          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
SPIN     N5O   0     C9 - C10- C1

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
#ROTAMER  N5O   0     C1 - C4   WHOLE_CONF
#ROTAMER  N5O   1     C9 - C10  WHOLE_CONF

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  N5O   0     C5 - O5   C5M
#=========================================================================


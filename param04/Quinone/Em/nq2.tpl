#1,2-Naphthoquinone        
#Agnes 10/09  
CONFLIST NQ2        NQ2BK NQ201 NQ2-1 NQ2DM

NATOM    NQ2DM      0
NATOM    NQ2BK      0
NATOM    NQ201      18
NATOM    NQ2-1      18

IATOM    NQ201  C1  0
IATOM    NQ201  O1  1
IATOM    NQ201  C2  2
IATOM    NQ201  O2  3
IATOM    NQ201  C3  4
IATOM    NQ201  H3  5
IATOM    NQ201  C4  6
IATOM    NQ201  H4  7
IATOM    NQ201  C5  8
IATOM    NQ201  H5  9
IATOM    NQ201  C6  10
IATOM    NQ201  H6  11
IATOM    NQ201  C7  12
IATOM    NQ201  H7  13
IATOM    NQ201  C8  14
IATOM    NQ201  H8  15
IATOM    NQ201  C9  16
IATOM    NQ201  C10 17
IATOM    NQ2-1  C1  0
IATOM    NQ2-1  O1  1
IATOM    NQ2-1  C2  2
IATOM    NQ2-1  O2  3
IATOM    NQ2-1  C3  4
IATOM    NQ2-1  H3  5
IATOM    NQ2-1  C4  6
IATOM    NQ2-1  H4  7
IATOM    NQ2-1  C5  8
IATOM    NQ2-1  H5  9
IATOM    NQ2-1  C6  10
IATOM    NQ2-1  H6  11
IATOM    NQ2-1  C7  12
IATOM    NQ2-1  H7  13
IATOM    NQ2-1  C8  14
IATOM    NQ2-1  H8  15
IATOM    NQ2-1  C9  16
IATOM    NQ2-1  C10 17


ATOMNAME NQ201    0  C1
ATOMNAME NQ201    1  O1 
ATOMNAME NQ201    2  C2 
ATOMNAME NQ201    3  O2 
ATOMNAME NQ201    4  C3 
ATOMNAME NQ201    5  H3 
ATOMNAME NQ201    6  C4 
ATOMNAME NQ201    7  H4 
ATOMNAME NQ201    8  C5 
ATOMNAME NQ201    9  H5 
ATOMNAME NQ201   10  C6 
ATOMNAME NQ201   11  H6 
ATOMNAME NQ201   12  C7 
ATOMNAME NQ201   13  H7 
ATOMNAME NQ201   14  C8 
ATOMNAME NQ201   15  H8 
ATOMNAME NQ201   16  C9 
ATOMNAME NQ201   17  C10 
ATOMNAME NQ2-1    0  C1
ATOMNAME NQ2-1    1  O1
ATOMNAME NQ2-1    2  C2
ATOMNAME NQ2-1    3  O2
ATOMNAME NQ2-1    4  C3
ATOMNAME NQ2-1    5  H3
ATOMNAME NQ2-1    6  C4
ATOMNAME NQ2-1    7  H4
ATOMNAME NQ2-1    8  C5
ATOMNAME NQ2-1    9  H5
ATOMNAME NQ2-1   10  C6
ATOMNAME NQ2-1   11  H6
ATOMNAME NQ2-1   12  C7
ATOMNAME NQ2-1   13  H7
ATOMNAME NQ2-1   14  C8
ATOMNAME NQ2-1   15  H8
ATOMNAME NQ2-1   16  C9
ATOMNAME NQ2-1   17  C10


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   NQ201      0

PKA      NQ201      0.0

ELECTRON NQ201      0

EM       NQ201      0.0

RXN      NQ201      -4.034

PROTON   NQ2-1      0

PKA      NQ2-1      0.0

ELECTRON NQ2-1      1

EM       NQ2-1      27

RXN      NQ2-1      -17.283

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  NQ201  C1  sp2       0     O1  0     C2  0     C9
CONNECT  NQ201  O1  s         0     C1
CONNECT  NQ201  C2  sp2       0     O2  0     C1  0     C3
CONNECT  NQ201  O2  s         0     C2
CONNECT  NQ201  C3  sp2       0     C2  0     H3  0     C4
CONNECT  NQ201  H3  s         0     C3
CONNECT  NQ201  C4  sp2       0     C3  0     H4  0     C10
CONNECT  NQ201  H4  s         0     C4
CONNECT  NQ201  C5  sp2       0     C10 0     H5  0     C6
CONNECT  NQ201  H5  s         0     C5
CONNECT  NQ201  C6  sp2       0     C5  0     C7  0     H6
CONNECT  NQ201  H6  s         0     C6
CONNECT  NQ201  C7  sp2       0     C8  0     C6  0     H7
CONNECT  NQ201  H7  s         0     C7
CONNECT  NQ201  C8  sp2       0     C7  0     C9  0     H8
CONNECT  NQ201  H8  s         0     C8
CONNECT  NQ201  C9  sp2       0     C1  0     C8  0     C10
CONNECT  NQ201  C10 sp2       0     C4  0     C5  0     C9


CONNECT  NQ2-1  C1  sp2       0     O1  0     C2  0     C9
CONNECT  NQ2-1  O1  s         0     C1
CONNECT  NQ2-1  C2  sp2       0     O2  0     C1  0     C3
CONNECT  NQ2-1  O2  s         0     C2
CONNECT  NQ2-1  C3  sp2       0     C2  0     H3  0     C4
CONNECT  NQ2-1  H3  s         0     C3
CONNECT  NQ2-1  C4  sp2       0     C3  0     H4  0     C10
CONNECT  NQ2-1  H4  s         0     C4
CONNECT  NQ2-1  C5  sp2       0     C10 0     H5  0     C6
CONNECT  NQ2-1  H5  s         0     C5
CONNECT  NQ2-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  NQ2-1  H6  s         0     C6
CONNECT  NQ2-1  C7  sp2       0     C8  0     C6  0     H7
CONNECT  NQ2-1  H7  s         0     C7
CONNECT  NQ2-1  C8  sp2       0     C7  0     C9  0     H8
CONNECT  NQ2-1  H8  s         0     C8
CONNECT  NQ2-1  C9  sp2       0     C1  0     C8  0     C10
CONNECT  NQ2-1  C10 sp2       0     C4  0     C5  0     C9


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   NQ2    C1  1.70
RADIUS   NQ2    O1  1.40
RADIUS   NQ2    C2  1.70
RADIUS   NQ2    O2  1.40
RADIUS   NQ2    C3  1.70
RADIUS   NQ2    H3  1.00
RADIUS   NQ2    C4  1.70
RADIUS   NQ2    H4  1.00
RADIUS   NQ2    C5  1.70
RADIUS   NQ2    H5  1.00
RADIUS   NQ2    C6  1.70
RADIUS   NQ2    H6  1.00
RADIUS   NQ2    C7  1.70
RADIUS   NQ2    H7  1.00
RADIUS   NQ2    C8  1.70
RADIUS   NQ2    H8  1.00
RADIUS   NQ2    C9  1.70
RADIUS   NQ2    C10 1.70


#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   NQ201  C1   0.42
CHARGE   NQ201  O1  -0.45
CHARGE   NQ201  C2   0.54
CHARGE   NQ201  O2  -0.49
CHARGE   NQ201  C3  -0.26
CHARGE   NQ201  H3   0.14
CHARGE   NQ201  C4  -0.17
CHARGE   NQ201  H4   0.14
CHARGE   NQ201  C5  -0.20
CHARGE   NQ201  H5   0.11
CHARGE   NQ201  C6   0.01
CHARGE   NQ201  H6   0.08
CHARGE   NQ201  C7  -0.14
CHARGE   NQ201  H7   0.10
CHARGE   NQ201  C8   0.01
CHARGE   NQ201  H8   0.11
CHARGE   NQ201  C9  -0.14
CHARGE   NQ201  C10  0.19

CHARGE   NQ2-1  C1   0.28
CHARGE   NQ2-1  O1  -0.60
CHARGE   NQ2-1  C2   0.43
CHARGE   NQ2-1  O2  -0.65
CHARGE   NQ2-1  C3  -0.26
CHARGE   NQ2-1  H3   0.13
CHARGE   NQ2-1  C4  -0.39
CHARGE   NQ2-1  H4   0.15
CHARGE   NQ2-1  C5  -0.28
CHARGE   NQ2-1  H5   0.10
CHARGE   NQ2-1  C6  -0.06
CHARGE   NQ2-1  H6   0.07
CHARGE   NQ2-1  C7  -0.27
CHARGE   NQ2-1  H7   0.11
CHARGE   NQ2-1  C8  -0.04
CHARGE   NQ2-1  H8   0.09
CHARGE   NQ2-1  C9  -0.06
CHARGE   NQ2-1  C10  0.25


TORSION  NQ2    O2   C2   C1   O1   f        0.0         1    180.00  
TORSION  NQ2    C8   C9   C1   O1   f        0.0         1    180.00  

#TORSION  NQ2    C7   C6   C5   C10  f        0.0         1      0.00 
#TORSION  NQ2    H5   C5   C10  C4   f        0.0         1      0.00
#TORSION  NQ2    C6   C7   C8   C9   f        0.0         1      0.00

#TORSION  NQ2    C7   C8   C9   C1   f        0.0         1      0.00 
#TORSION  NQ2    C5   C6   C7   C8   f        0.0         1      0.00
#TORSION  NQ2    H5   C5   C6   C7   f        0.0         1    180.00
#TORSION  NQ2    C6   C5   C10  C9   f        0.0         1    180.00 
#TORSION  NQ2    C6   C5   C10  C4   f        0.0         1      0.00 

#TORSION  NQ2    H5   C5   C10  C4   f        0.0         1      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    NQ2          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     NQ2   0     C9 - C10- C1

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  NQ2   0     C1 - C4   WHOLE_CONF
ROTAMER  NQ2   1     C9 - C10  WHOLE_CONF


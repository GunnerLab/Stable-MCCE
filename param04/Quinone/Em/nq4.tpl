#1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST NQ4        NQ4BK NQ401 NQ4-1 NQ4DM

NATOM    NQ4DM      0
NATOM    NQ4BK      0
NATOM    NQ401      18
NATOM    NQ4-1      18

IATOM    NQ401  C1  0
IATOM    NQ401  O1  1
IATOM    NQ401  C2  2
IATOM    NQ401  H2  3
IATOM    NQ401  C3  4
IATOM    NQ401  H3  5
IATOM    NQ401  C4  6
IATOM    NQ401  O4  7
IATOM    NQ401  C5  8
IATOM    NQ401  H5  9
IATOM    NQ401  C6  10
IATOM    NQ401  H6  11
IATOM    NQ401  C7  12
IATOM    NQ401  H7  13
IATOM    NQ401  C8  14
IATOM    NQ401  H8  15
IATOM    NQ401  C9  16
IATOM    NQ401  C10 17
IATOM    NQ4-1  C1  0
IATOM    NQ4-1  O1  1
IATOM    NQ4-1  C2  2
IATOM    NQ4-1  H2  3
IATOM    NQ4-1  C3  4
IATOM    NQ4-1  H3  5
IATOM    NQ4-1  C4  6
IATOM    NQ4-1  O4  7
IATOM    NQ4-1  C5  8
IATOM    NQ4-1  H5  9
IATOM    NQ4-1  C6  10
IATOM    NQ4-1  H6  11
IATOM    NQ4-1  C7  12
IATOM    NQ4-1  H7  13
IATOM    NQ4-1  C8  14
IATOM    NQ4-1  H8  15
IATOM    NQ4-1  C9  16
IATOM    NQ4-1  C10 17

ATOMNAME NQ401    0  C1
ATOMNAME NQ401    1  O1 
ATOMNAME NQ401    2  C2 
ATOMNAME NQ401    3  H2 
ATOMNAME NQ401    4  C3 
ATOMNAME NQ401    5  H3 
ATOMNAME NQ401    6  C4 
ATOMNAME NQ401    7  O4 
ATOMNAME NQ401    8  C5 
ATOMNAME NQ401    9  H5 
ATOMNAME NQ401   10  C6 
ATOMNAME NQ401   11  H6 
ATOMNAME NQ401   12  C7 
ATOMNAME NQ401   13  H7 
ATOMNAME NQ401   14  C8 
ATOMNAME NQ401   15  H8 
ATOMNAME NQ401   16  C9 
ATOMNAME NQ401   17  C10 
ATOMNAME NQ4-1    0  C1
ATOMNAME NQ4-1    1  O1
ATOMNAME NQ4-1    2  C2
ATOMNAME NQ4-1    3  H2
ATOMNAME NQ4-1    4  C3
ATOMNAME NQ4-1    5  H3
ATOMNAME NQ4-1    6  C4
ATOMNAME NQ4-1    7  O4
ATOMNAME NQ4-1    8  C5
ATOMNAME NQ4-1    9  H5
ATOMNAME NQ4-1   10  C6
ATOMNAME NQ4-1   11  H6
ATOMNAME NQ4-1   12  C7
ATOMNAME NQ4-1   13  H7
ATOMNAME NQ4-1   14  C8
ATOMNAME NQ4-1   15  H8
ATOMNAME NQ4-1   16  C9
ATOMNAME NQ4-1   17  C10



#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   NQ401      0

PKA      NQ401      0.0

ELECTRON NQ401      0

EM       NQ401      0.0

RXN      NQ401      -3.695

PROTON   NQ4-1      0

PKA      NQ4-1      0.0

ELECTRON NQ4-1      1

EM       NQ4-1      -83

RXN      NQ4-1      -16.502

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  NQ401  C1  sp2       0     O1  0     C2  0     C9
CONNECT  NQ401  O1  s         0     C1
CONNECT  NQ401  C2  sp2       0     H2  0     C1  0     C3
CONNECT  NQ401  H2  s         0     C2
CONNECT  NQ401  C3  sp2       0     C2  0     H3  0     C4
CONNECT  NQ401  H3  s         0     C3
CONNECT  NQ401  C4  sp2       0     C3  0     O4  0     C10
CONNECT  NQ401  O4  s         0     C4
CONNECT  NQ401  C5  sp2       0     C10 0     H5  0     C6
CONNECT  NQ401  H5  s         0     C5
CONNECT  NQ401  C6  sp2       0     C5  0     C7  0     H6
CONNECT  NQ401  H6  s         0     C6
CONNECT  NQ401  C7  sp2       0     C8  0     C6  0     H7
CONNECT  NQ401  H7  s         0     C7
CONNECT  NQ401  C8  sp2       0     C7  0     C9  0     H8
CONNECT  NQ401  H8  s         0     C8
CONNECT  NQ401  C9  sp2       0     C1  0     C8  0     C10
CONNECT  NQ401  C10 sp2       0     C4  0     C5  0     C9

CONNECT  NQ4-1  C1  sp2       0     O1  0     C2  0     C9
CONNECT  NQ4-1  O1  s         0     C1
CONNECT  NQ4-1  C2  sp2       0     H2  0     C1  0     C3
CONNECT  NQ4-1  H2  s         0     C2
CONNECT  NQ4-1  C3  sp2       0     C2  0     H3  0     C4
CONNECT  NQ4-1  H3  s         0     C3
CONNECT  NQ4-1  C4  sp2       0     C3  0     O4  0     C10
CONNECT  NQ4-1  O4  s         0     C4
CONNECT  NQ4-1  C5  sp2       0     C10 0     H5  0     C6
CONNECT  NQ4-1  H5  s         0     C5
CONNECT  NQ4-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  NQ4-1  H6  s         0     C6
CONNECT  NQ4-1  C7  sp2       0     C8  0     C6  0     H7
CONNECT  NQ4-1  H7  s         0     C7
CONNECT  NQ4-1  C8  sp2       0     C7  0     C9  0     H8
CONNECT  NQ4-1  H8  s         0     C8
CONNECT  NQ4-1  C9  sp2       0     C1  0     C8  0     C10
CONNECT  NQ4-1  C10 sp2       0     C4  0     C5  0     C9

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   NQ4    C1  1.70
RADIUS   NQ4    O1  1.40
RADIUS   NQ4    C2  1.70
RADIUS   NQ4    H2  1.00
RADIUS   NQ4    C3  1.70
RADIUS   NQ4    H3  1.00
RADIUS   NQ4    C4  1.70
RADIUS   NQ4    O4  1.40
RADIUS   NQ4    C5  1.70
RADIUS   NQ4    H5  1.00
RADIUS   NQ4    C6  1.70
RADIUS   NQ4    H6  1.00
RADIUS   NQ4    C7  1.70
RADIUS   NQ4    H7  1.00
RADIUS   NQ4    C8  1.70
RADIUS   NQ4    H8  1.00
RADIUS   NQ4    C9  1.70
RADIUS   NQ4    C10 1.70


#NEUTRAL------
#23456789A123456789B123456789C
# opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)   Agnes 10/09
CHARGE   NQ401  C1   0.55
CHARGE   NQ401  O1  -0.50
CHARGE   NQ401  C2  -0.18
CHARGE   NQ401  H2   0.13
CHARGE   NQ401  C3  -0.18
CHARGE   NQ401  H3   0.13
CHARGE   NQ401  C4   0.54
CHARGE   NQ401  O4  -0.50
CHARGE   NQ401  C5  -0.05
CHARGE   NQ401  H5   0.10
CHARGE   NQ401  C6  -0.06
CHARGE   NQ401  H6   0.09
CHARGE   NQ401  C7  -0.07
CHARGE   NQ401  H7   0.09
CHARGE   NQ401  C8  -0.04
CHARGE   NQ401  H8   0.10
CHARGE   NQ401  C9  -0.09
CHARGE   NQ401  C10 -0.06
CHARGE   NQ4-1  C1   0.40
CHARGE   NQ4-1  O1  -0.66
CHARGE   NQ4-1  C2  -0.28
CHARGE   NQ4-1  H2   0.12
CHARGE   NQ4-1  C3  -0.31
CHARGE   NQ4-1  H3   0.13
CHARGE   NQ4-1  C4   0.43
CHARGE   NQ4-1  O4  -0.66
CHARGE   NQ4-1  C5  -0.11
CHARGE   NQ4-1  H5   0.10
CHARGE   NQ4-1  C6  -0.15
CHARGE   NQ4-1  H6   0.08
CHARGE   NQ4-1  C7  -0.16
CHARGE   NQ4-1  H7   0.09
CHARGE   NQ4-1  C8  -0.13
CHARGE   NQ4-1  H8   0.10
CHARGE   NQ4-1  C9   0.03
CHARGE   NQ4-1  C10 -0.02

#TORSION  NQ4    O2   C2   C1   O1   f        0.0         1    180.00  
#TORSION  NQ4    C8   C9   C1   O1   f        0.0         1    180.00  

#TORSION  NQ4    C7   C6   C5   C10  f        0.0         1      0.00 
#TORSION  NQ4    H5   C5   C10  C4   f        0.0         1      0.00
#TORSION  NQ4    C6   C7   C8   C9   f        0.0         1      0.00

#TORSION  NQ4    C7   C8   C9   C1   f        0.0         1      0.00 
#TORSION  NQ4    C5   C6   C7   C8   f        0.0         1      0.00
#TORSION  NQ4    H5   C5   C6   C7   f        0.0         1    180.00
#TORSION  NQ4    C6   C5   C10  C9   f        0.0         1    180.00 
#TORSION  NQ4    C6   C5   C10  C4   f        0.0         1      0.00 

#TORSION  NQ4    H5   C5   C10  C4   f        0.0         1      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    NQ4          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     NQ4   0     C9 - C10- C1

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  NQ4   0     C1 - C4   WHOLE_CONF
ROTAMER  NQ4   1     C9 - C10  WHOLE_CONF

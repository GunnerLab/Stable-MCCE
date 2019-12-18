#5-methyl, 1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST N5M        N5MBK N5M-1 N5MDM

NATOM    N5MDM      0
NATOM    N5MBK      0
NATOM    N5M-1      21

IATOM    N5M-1  C1  0
IATOM    N5M-1  O1  1
IATOM    N5M-1  C2  2
IATOM    N5M-1  H2  3
IATOM    N5M-1  C3  4
IATOM    N5M-1  H3  5
IATOM    N5M-1  C4  6
IATOM    N5M-1  O4  7
IATOM    N5M-1  C5  8
IATOM    N5M-1  C5M 9
IATOM    N5M-1  C6  10
IATOM    N5M-1  H6  11
IATOM    N5M-1  C7  12
IATOM    N5M-1  H7  13
IATOM    N5M-1  C8  14
IATOM    N5M-1  H8  15
IATOM    N5M-1  C9  16
IATOM    N5M-1  C10 17
IATOM    N5M-1 1H5M 18
IATOM    N5M-1 2H5M 19
IATOM    N5M-1 3H5M 20


ATOMNAME N5M-1    0  C1
ATOMNAME N5M-1    1  O1 
ATOMNAME N5M-1    2  C2 
ATOMNAME N5M-1    3  H2 
ATOMNAME N5M-1    4  C3 
ATOMNAME N5M-1    5  H3 
ATOMNAME N5M-1    6  C4 
ATOMNAME N5M-1    7  O4 
ATOMNAME N5M-1    8  C5 
ATOMNAME N5M-1    9  C5M 
ATOMNAME N5M-1   10  C6 
ATOMNAME N5M-1   11  H6 
ATOMNAME N5M-1   12  C7 
ATOMNAME N5M-1   13  H7 
ATOMNAME N5M-1   14  C8 
ATOMNAME N5M-1   15  H8 
ATOMNAME N5M-1   16  C9 
ATOMNAME N5M-1   17  C10 
ATOMNAME N5M-1   18 1H5M
ATOMNAME N5M-1   19 2H5M
ATOMNAME N5M-1   20 3H5M 

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   N5M-1      0

PKA      N5M-1      0.0

ELECTRON N5M-1      1

EM       N5M-1      0.0

RXN      N5M-1      -16.665


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  N5M-1  C1  sp2       0     O1  0     C2  0     C9
CONNECT  N5M-1  O1  s         0     C1
CONNECT  N5M-1  C2  sp2       0     H2  0     C1  0     C3
CONNECT  N5M-1  H2  s         0     C2
CONNECT  N5M-1  C3  sp2       0     C2  0     H3  0     C4
CONNECT  N5M-1  H3  s         0     C3
CONNECT  N5M-1  C4  sp2       0     C3  0     O4  0     C10
CONNECT  N5M-1  O4  s         0     C4
CONNECT  N5M-1  C5  sp2       0     C10 0     C5M 0     C6
CONNECT  N5M-1  C5M sp3       0     C5  0    1H5M 0    2H5M 0    3H5M
CONNECT  N5M-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  N5M-1  H6  s         0     C6
CONNECT  N5M-1  C7  sp2       0     C8  0     C6  0     H7
CONNECT  N5M-1  H7  s         0     C7
CONNECT  N5M-1  C8  sp2       0     C7  0     C9  0     H8
CONNECT  N5M-1  H8  s         0     C8
CONNECT  N5M-1  C9  sp2       0     C1  0     C8  0     C10
CONNECT  N5M-1  C10 sp2       0     C4  0     C5  0     C9
CONNECT  N5M-1 1H5M s         0     C5M
CONNECT  N5M-1 2H5M s         0     C5M
CONNECT  N5M-1 3H5M s         0     C5M


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   N5M    C1  1.70
RADIUS   N5M    O1  1.40
RADIUS   N5M    C2  1.70
RADIUS   N5M    H2  1.00
RADIUS   N5M    C3  1.70
RADIUS   N5M    H3  1.00
RADIUS   N5M    C4  1.70
RADIUS   N5M    O4  1.40
RADIUS   N5M    C5  1.70
RADIUS   N5M    C5M 1.70
RADIUS   N5M    C6  1.70
RADIUS   N5M    H6  1.00
RADIUS   N5M    C7  1.70
RADIUS   N5M    H7  1.00
RADIUS   N5M    C8  1.70
RADIUS   N5M    H8  1.00
RADIUS   N5M    C9  1.70
RADIUS   N5M    C10 1.70
RADIUS   N5M   1H5M 1.00
RADIUS   N5M   2H5M 1.00
RADIUS   N5M   3H5M 1.00

#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity pop=chelpg scf(maxcycle=600)   Agnes 10/09
CHARGE   N5M-1  C1   0.38
CHARGE   N5M-1  O1  -0.66
CHARGE   N5M-1  C2  -0.23
CHARGE   N5M-1  H2   0.11
CHARGE   N5M-1  C3  -0.34
CHARGE   N5M-1  H3   0.13
CHARGE   N5M-1  C4   0.51
CHARGE   N5M-1  O4  -0.68
CHARGE   N5M-1  C5   0.27
CHARGE   N5M-1  C5M -0.30
CHARGE   N5M-1  C6  -0.29
CHARGE   N5M-1  H6   0.12
CHARGE   N5M-1  C7  -0.12
CHARGE   N5M-1  H7   0.08
CHARGE   N5M-1  C8  -0.18
CHARGE   N5M-1  H8   0.12
CHARGE   N5M-1  C9   0.08
CHARGE   N5M-1  C10 -0.22
CHARGE   N5M-1 1H5M  0.04
CHARGE   N5M-1 2H5M  0.04
CHARGE   N5M-1 3H5M  0.14



#TORSION  N5M    O2   C2   C1   O1   f        0.0         1    180.00  
#TORSION  N5M    C8   C9   C1   O1   f        0.0         1    180.00  

#TORSION  N5M    C7   C6   C5   C10  f        0.0         1      0.00 
#TORSION  N5M    H5   C5   C10  C4   f        0.0         1      0.00
#TORSION  N5M    C6   C7   C8   C9   f        0.0         1      0.00

#TORSION  N5M    C7   C8   C9   C1   f        0.0         1      0.00 
#TORSION  N5M    C5   C6   C7   C8   f        0.0         1      0.00
#TORSION  N5M    H5   C5   C6   C7   f        0.0         1    180.00
#TORSION  N5M    C6   C5   C10  C9   f        0.0         1    180.00 
#TORSION  N5M    C6   C5   C10  C4   f        0.0         1      0.00 

#TORSION  N5M    H5   C5   C10  C4   f        0.0         1      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    N5M          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     N5M   0     C9 - C10- C1
#SPIN     N5M   1     C8 - C5 - C10
#SPIN     N5M   2     C1 - C4 - C9

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  N5M   0     C1 - C4   WHOLE_CONF
ROTAMER  N5M   1     C9 - C10  WHOLE_CONF
ROTAMER  N5M   2     C2 - C7   WHOLE_CONF

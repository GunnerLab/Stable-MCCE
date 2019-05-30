#5-methyl, 1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST N5M        N5MBK N5M01 N5MDM

NATOM    N5MDM      0
NATOM    N5MBK      0
NATOM    N5M01      21

IATOM    N5M01  C1  0
IATOM    N5M01  O1  1
IATOM    N5M01  C2  2
IATOM    N5M01  H2  3
IATOM    N5M01  C3  4
IATOM    N5M01  H3  5
IATOM    N5M01  C4  6
IATOM    N5M01  O4  7
IATOM    N5M01  C5  8
IATOM    N5M01  C5M 9
IATOM    N5M01  C6  10
IATOM    N5M01  H6  11
IATOM    N5M01  C7  12
IATOM    N5M01  H7  13
IATOM    N5M01  C8  14
IATOM    N5M01  H8  15
IATOM    N5M01  C9  16
IATOM    N5M01  C10 17
IATOM    N5M01 1H5M 18
IATOM    N5M01 2H5M 19
IATOM    N5M01 3H5M 20


ATOMNAME N5M01    0  C1
ATOMNAME N5M01    1  O1 
ATOMNAME N5M01    2  C2 
ATOMNAME N5M01    3  H2 
ATOMNAME N5M01    4  C3 
ATOMNAME N5M01    5  H3 
ATOMNAME N5M01    6  C4 
ATOMNAME N5M01    7  O4 
ATOMNAME N5M01    8  C5 
ATOMNAME N5M01    9  C5M 
ATOMNAME N5M01   10  C6 
ATOMNAME N5M01   11  H6 
ATOMNAME N5M01   12  C7 
ATOMNAME N5M01   13  H7 
ATOMNAME N5M01   14  C8 
ATOMNAME N5M01   15  H8 
ATOMNAME N5M01   16  C9 
ATOMNAME N5M01   17  C10 
ATOMNAME N5M01   18 1H5M
ATOMNAME N5M01   19 2H5M
ATOMNAME N5M01   20 3H5M 

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   N5M01      0

PKA      N5M01      0.0

ELECTRON N5M01      0

EM       N5M01      0.0

RXN      N5M01      -3.789 


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  N5M01  C1  sp2       0     O1  0     C2  0     C9
CONNECT  N5M01  O1  s         0     C1
CONNECT  N5M01  C2  sp2       0     H2  0     C1  0     C3
CONNECT  N5M01  H2  s         0     C2
CONNECT  N5M01  C3  sp2       0     C2  0     H3  0     C4
CONNECT  N5M01  H3  s         0     C3
CONNECT  N5M01  C4  sp2       0     C3  0     O4  0     C10
CONNECT  N5M01  O4  s         0     C4
CONNECT  N5M01  C5  sp2       0     C10 0     C5M 0     C6
CONNECT  N5M01  C5M sp3       0     C5  0    1H5M 0    2H5M 0    3H5M
CONNECT  N5M01  C6  sp2       0     C5  0     C7  0     H6
CONNECT  N5M01  H6  s         0     C6
CONNECT  N5M01  C7  sp2       0     C8  0     C6  0     H7
CONNECT  N5M01  H7  s         0     C7
CONNECT  N5M01  C8  sp2       0     C7  0     C9  0     H8
CONNECT  N5M01  H8  s         0     C8
CONNECT  N5M01  C9  sp2       0     C1  0     C8  0     C10
CONNECT  N5M01  C10 sp2       0     C4  0     C5  0     C9
CONNECT  N5M01 1H5M s         0     C5M
CONNECT  N5M01 2H5M s         0     C5M
CONNECT  N5M01 3H5M s         0     C5M


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
# opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   N5M01  C1   0.54
CHARGE   N5M01  O1  -0.52
CHARGE   N5M01  C2  -0.18
CHARGE   N5M01  H2   0.15
CHARGE   N5M01  C3  -0.25
CHARGE   N5M01  H3   0.16
CHARGE   N5M01  C4   0.63
CHARGE   N5M01  O4  -0.53
CHARGE   N5M01  C5   0.31
CHARGE   N5M01  C5M -0.46
CHARGE   N5M01  C6  -0.22
CHARGE   N5M01  H6   0.11
CHARGE   N5M01  C7  -0.06
CHARGE   N5M01  H7   0.11
CHARGE   N5M01  C8  -0.10
CHARGE   N5M01  H8   0.12
CHARGE   N5M01  C9  -0.02
CHARGE   N5M01  C10 -0.23
CHARGE   N5M01 1H5M  0.15
CHARGE   N5M01 2H5M  0.14
CHARGE   N5M01 3H5M  0.15



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

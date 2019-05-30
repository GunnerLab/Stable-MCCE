#2,3-dimethyl, 1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST NDI        NDIBK NDI01 NDIDM

NATOM    NDIDM      0
NATOM    NDIBK      0
NATOM    NDI01      24

IATOM    NDI01  C1  0
IATOM    NDI01  O1  1
IATOM    NDI01  C2  2
IATOM    NDI01  C2M 3
IATOM    NDI01  C3  4
IATOM    NDI01  C3M 5
IATOM    NDI01  C4  6
IATOM    NDI01  O4  7
IATOM    NDI01  C5  8
IATOM    NDI01  H5  9
IATOM    NDI01  C6  10
IATOM    NDI01  H6  11
IATOM    NDI01  C7  12
IATOM    NDI01  H7  13
IATOM    NDI01  C8  14
IATOM    NDI01  H8  15
IATOM    NDI01  C9  16
IATOM    NDI01  C10 17
IATOM    NDI01 1H2M 18
IATOM    NDI01 2H2M 19
IATOM    NDI01 3H2M 20
IATOM    NDI01 1H3M 21
IATOM    NDI01 2H3M 22
IATOM    NDI01 3H3M 23



ATOMNAME NDI01    0  C1
ATOMNAME NDI01    1  O1 
ATOMNAME NDI01    2  C2 
ATOMNAME NDI01    3  C2M
ATOMNAME NDI01    4  C3 
ATOMNAME NDI01    5  C3M
ATOMNAME NDI01    6  C4 
ATOMNAME NDI01    7  O4 
ATOMNAME NDI01    8  C5 
ATOMNAME NDI01    9  H5 
ATOMNAME NDI01   10  C6 
ATOMNAME NDI01   11  H6 
ATOMNAME NDI01   12  C7 
ATOMNAME NDI01   13  H7 
ATOMNAME NDI01   14  C8 
ATOMNAME NDI01   15  H8 
ATOMNAME NDI01   16  C9 
ATOMNAME NDI01   17  C10 
ATOMNAME NDI01   18 1H2M
ATOMNAME NDI01   19 2H2M
ATOMNAME NDI01   20 3H2M
ATOMNAME NDI01   21 1H3M
ATOMNAME NDI01   22 2H3M
ATOMNAME NDI01   23 3H3M


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   NDI01      0

PKA      NDI01      0.0

ELECTRON NDI01      0

EM       NDI01      0.0

RXN      NDI01      -3.327


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  NDI01  C1  sp2       0     O1  0     C2  0     C9
CONNECT  NDI01  O1  s         0     C1
CONNECT  NDI01  C2  sp2       0     C2M 0     C1  0     C3
CONNECT  NDI01  C2M sp3       0     C2  0    1H2M 0    2H2M 0    3H2M
CONNECT  NDI01  C3  sp2       0     C2  0     C3M 0     C4
CONNECT  NDI01  C3M sp3       0     C3  0    1H3M 0    2H3M 0    3H3M
CONNECT  NDI01  C4  sp2       0     C3  0     O4  0     C10
CONNECT  NDI01  O4  s         0     C4
CONNECT  NDI01  C5  sp2       0     C10 0     H5  0     C6
CONNECT  NDI01  H5  s         0     C5
CONNECT  NDI01  C6  sp2       0     C5  0     C7  0     H6
CONNECT  NDI01  H6  s         0     C6
CONNECT  NDI01  C7  sp2       0     C8  0     C6  0     H7
CONNECT  NDI01  H7  s         0     C7
CONNECT  NDI01  C8  sp2       0     C7  0     C9  0     H8
CONNECT  NDI01  H8  s         0     C8
CONNECT  NDI01  C9  sp2       0     C1  0     C8  0     C10
CONNECT  NDI01  C10 sp2       0     C4  0     C5  0     C9
CONNECT  NDI01 1H2M s         0     C2M
CONNECT  NDI01 2H2M s         0     C2M
CONNECT  NDI01 3H2M s         0     C2M
CONNECT  NDI01 1H3M s         0     C3M
CONNECT  NDI01 2H3M s         0     C3M
CONNECT  NDI01 3H3M s         0     C3M



#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   NDI    C1  1.70
RADIUS   NDI    O1  1.40
RADIUS   NDI    C2  1.70
RADIUS   NDI    C2M 1.70
RADIUS   NDI    C3  1.70
RADIUS   NDI    C3M 1.70
RADIUS   NDI    C4  1.70
RADIUS   NDI    O4  1.40
RADIUS   NDI    C5  1.70
RADIUS   NDI    H5  1.00
RADIUS   NDI    C6  1.70
RADIUS   NDI    H6  1.00
RADIUS   NDI    C7  1.70
RADIUS   NDI    H7  1.00
RADIUS   NDI    C8  1.70
RADIUS   NDI    H8  1.00
RADIUS   NDI    C9  1.70
RADIUS   NDI    C10 1.70
RADIUS   NDI   1H2M 1.00
RADIUS   NDI   2H2M 1.00
RADIUS   NDI   3H2M 1.00
RADIUS   NDI   1H3M 1.00
RADIUS   NDI   2H3M 1.00
RADIUS   NDI   3H3M 1.00


#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   NDI01  C1   0.45
CHARGE   NDI01  O1  -0.47
CHARGE   NDI01  C2  -0.10
CHARGE   NDI01  C2M  0.04
CHARGE   NDI01  C3  -0.08
CHARGE   NDI01  C3M -0.09
CHARGE   NDI01  C4   0.50
CHARGE   NDI01  O4  -0.49
CHARGE   NDI01  C5  -0.01
CHARGE   NDI01  H5   0.09
CHARGE   NDI01  C6  -0.11
CHARGE   NDI01  H6   0.12
CHARGE   NDI01  C7  -0.12
CHARGE   NDI01  H7   0.12
CHARGE   NDI01  C8  -0.00
CHARGE   NDI01  H8   0.08
CHARGE   NDI01  C9  -0.07
CHARGE   NDI01  C10 -0.09
CHARGE   NDI01 1H2M  0.02
CHARGE   NDI01 2H2M  0.02
CHARGE   NDI01 3H2M  0.03
CHARGE   NDI01 1H3M  0.06
CHARGE   NDI01 2H3M  0.05
CHARGE   NDI01 3H3M  0.05




#TORSION  NDI    O2   C2   C1   O1   f        0.0         1    180.00  
#TORSION  NDI    C8   C9   C1   O1   f        0.0         1    180.00  

#TORSION  NDI    C7   C6   C5   C10  f        0.0         1      0.00 
#TORSION  NDI    H5   C5   C10  C4   f        0.0         1      0.00
#TORSION  NDI    C6   C7   C8   C9   f        0.0         1      0.00

#TORSION  NDI    C7   C8   C9   C1   f        0.0         1      0.00 
#TORSION  NDI    C5   C6   C7   C8   f        0.0         1      0.00
#TORSION  NDI    H5   C5   C6   C7   f        0.0         1    180.00
#TORSION  NDI    C6   C5   C10  C9   f        0.0         1    180.00 
#TORSION  NDI    C6   C5   C10  C4   f        0.0         1      0.00 

#TORSION  NDI    H5   C5   C10  C4   f        0.0         1      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    NDI          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     NDI   0     C9 - C10- C1
#SPIN     NDI   1     C8 - C5 - C10
#SPIN     NDI   2     C1 - C4 - C9

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  NDI   0     C1 - C4   WHOLE_CONF
ROTAMER  NDI   1     C9 - C10  WHOLE_CONF
ROTAMER  NDI   2     C2 - C7   WHOLE_CONF

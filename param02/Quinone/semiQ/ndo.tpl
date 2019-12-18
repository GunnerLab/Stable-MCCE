#2,3-dimethoxy, 1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST NDO        NDOBK NDO-1 NDODM

NATOM    NDODM      0
NATOM    NDOBK      0
NATOM    NDO-1      26

IATOM    NDO-1  C1  0
IATOM    NDO-1  O1  1
IATOM    NDO-1  C2  2
IATOM    NDO-1  C2M 3
IATOM    NDO-1  C3  4
IATOM    NDO-1  C3M 5
IATOM    NDO-1  C4  6
IATOM    NDO-1  O4  7
IATOM    NDO-1  C5  8
IATOM    NDO-1  H5  9
IATOM    NDO-1  C6  10
IATOM    NDO-1  H6  11
IATOM    NDO-1  C7  12
IATOM    NDO-1  H7  13
IATOM    NDO-1  C8  14
IATOM    NDO-1  H8  15
IATOM    NDO-1  C9  16
IATOM    NDO-1  C10 17
IATOM    NDO-1 1H2M 18
IATOM    NDO-1 2H2M 19
IATOM    NDO-1 3H2M 20
IATOM    NDO-1 1H3M 21
IATOM    NDO-1 2H3M 22
IATOM    NDO-1 3H3M 23
IATOM    NDO-1  O2  24
IATOM    NDO-1  O3  25


ATOMNAME NDO-1    0  C1
ATOMNAME NDO-1    1  O1 
ATOMNAME NDO-1    2  C2 
ATOMNAME NDO-1    3  C2M
ATOMNAME NDO-1    4  C3 
ATOMNAME NDO-1    5  C3M
ATOMNAME NDO-1    6  C4 
ATOMNAME NDO-1    7  O4 
ATOMNAME NDO-1    8  C5 
ATOMNAME NDO-1    9  H5 
ATOMNAME NDO-1   10  C6 
ATOMNAME NDO-1   11  H6 
ATOMNAME NDO-1   12  C7 
ATOMNAME NDO-1   13  H7 
ATOMNAME NDO-1   14  C8 
ATOMNAME NDO-1   15  H8 
ATOMNAME NDO-1   16  C9 
ATOMNAME NDO-1   17  C10 
ATOMNAME NDO-1   18 1H2M
ATOMNAME NDO-1   19 2H2M
ATOMNAME NDO-1   20 3H2M
ATOMNAME NDO-1   21 1H3M
ATOMNAME NDO-1   22 2H3M
ATOMNAME NDO-1   23 3H3M
ATOMNAME NDO-1   24  O2
ATOMNAME NDO-1   25  O3

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   NDO-1      0

PKA      NDO-1      0.0

ELECTRON NDO-1      1

EM       NDO-1      0.0

RXN      NDO-1      -15.736


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  NDO-1  C1  sp2       0     O1  0     C2  0     C9
CONNECT  NDO-1  O1  s         0     C1
CONNECT  NDO-1  C2  sp2       0     O2  0     C1  0     C3
CONNECT  NDO-1  C2M sp3       0     O2  0    1H2M 0    2H2M 0    3H2M
CONNECT  NDO-1  C3  sp2       0     C2  0     O3  0     C4
CONNECT  NDO-1  C3M sp3       0     O3  0    1H3M 0    2H3M 0    3H3M
CONNECT  NDO-1  C4  sp2       0     C3  0     O4  0     C10
CONNECT  NDO-1  O4  s         0     C4
CONNECT  NDO-1  C5  sp2       0     C10 0     H5  0     C6
CONNECT  NDO-1  H5  s         0     C5
CONNECT  NDO-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  NDO-1  H6  s         0     C6
CONNECT  NDO-1  C7  sp2       0     C8  0     C6  0     H7
CONNECT  NDO-1  H7  s         0     C7
CONNECT  NDO-1  C8  sp2       0     C7  0     C9  0     H8
CONNECT  NDO-1  H8  s         0     C8
CONNECT  NDO-1  C9  sp2       0     C1  0     C8  0     C10
CONNECT  NDO-1  C10 sp2       0     C4  0     C5  0     C9
CONNECT  NDO-1 1H2M s         0     C2M
CONNECT  NDO-1 2H2M s         0     C2M
CONNECT  NDO-1 3H2M s         0     C2M
CONNECT  NDO-1 1H3M s         0     C3M
CONNECT  NDO-1 2H3M s         0     C3M
CONNECT  NDO-1 3H3M s         0     C3M
CONNECT  NDO-1  O2  sp2       0     C2M 0     C2
CONNECT  NDO-1  O3  sp2       0     C3M 0     C3


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   NDO    C1  1.70
RADIUS   NDO    O1  1.40
RADIUS   NDO    C2  1.70
RADIUS   NDO    C2M 1.70
RADIUS   NDO    C3  1.70
RADIUS   NDO    C3M 1.70
RADIUS   NDO    C4  1.70
RADIUS   NDO    O4  1.40
RADIUS   NDO    C5  1.70
RADIUS   NDO    H5  1.00
RADIUS   NDO    C6  1.70
RADIUS   NDO    H6  1.00
RADIUS   NDO    C7  1.70
RADIUS   NDO    H7  1.00
RADIUS   NDO    C8  1.70
RADIUS   NDO    H8  1.00
RADIUS   NDO    C9  1.70
RADIUS   NDO    C10 1.70
RADIUS   NDO   1H2M 1.00
RADIUS   NDO   2H2M 1.00
RADIUS   NDO   3H2M 1.00
RADIUS   NDO   1H3M 1.00
RADIUS   NDO   2H3M 1.00
RADIUS   NDO   3H3M 1.00
RADIUS   NDO    O2  1.40
RADIUS   NDO    O3  1.40


#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity pop=chelpg scf(maxcycle=600)   Agnes 10/09
CHARGE   NDO-1  C1   0.19
CHARGE   NDO-1  O1  -0.53
CHARGE   NDO-1  C2   0.17
CHARGE   NDO-1  C2M  0.25
CHARGE   NDO-1  C3   0.12
CHARGE   NDO-1  C3M  0.25
CHARGE   NDO-1  C4   0.26
CHARGE   NDO-1  O4  -0.54
CHARGE   NDO-1  C5   0.02
CHARGE   NDO-1  H5   0.04
CHARGE   NDO-1  C6  -0.21
CHARGE   NDO-1  H6   0.10
CHARGE   NDO-1  C7  -0.17
CHARGE   NDO-1  H7   0.09
CHARGE   NDO-1  C8  -0.02
CHARGE   NDO-1  H8   0.05
CHARGE   NDO-1  C9  -0.03
CHARGE   NDO-1  C10 -0.09
CHARGE   NDO-1 1H2M  0.00
CHARGE   NDO-1 2H2M -0.01
CHARGE   NDO-1 3H2M -0.01
CHARGE   NDO-1 1H3M -0.01
CHARGE   NDO-1 2H3M -0.01
CHARGE   NDO-1 3H3M  0.00
CHARGE   NDO-1  O2  -0.46
CHARGE   NDO-1  O3  -0.45



#TORSION  NDO    O2   C2   C1   O1   f        0.0         1    180.00  
#TORSION  NDO    C8   C9   C1   O1   f        0.0         1    180.00  

#TORSION  NDO    C7   C6   C5   C10  f        0.0         1      0.00 
#TORSION  NDO    H5   C5   C10  C4   f        0.0         1      0.00
#TORSION  NDO    C6   C7   C8   C9   f        0.0         1      0.00

#TORSION  NDO    C7   C8   C9   C1   f        0.0         1      0.00 
#TORSION  NDO    C5   C6   C7   C8   f        0.0         1      0.00
#TORSION  NDO    H5   C5   C6   C7   f        0.0         1    180.00
#TORSION  NDO    C6   C5   C10  C9   f        0.0         1    180.00 
#TORSION  NDO    C6   C5   C10  C4   f        0.0         1      0.00 

#TORSION  NDO    H5   C5   C10  C4   f        0.0         1      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    NDO          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     NDO   0     C9 - C10- C1
#SPIN     NDO   1     C8 - C5 - C10
#SPIN     NDO   2     C1 - C4 - C9

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  NDO   0     C1 - C4   WHOLE_CONF
ROTAMER  NDO   1     C9 - C10  WHOLE_CONF
ROTAMER  NDO   2     C2 - C7   WHOLE_CONF

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  NDO   3     C2 - O2   C2M
ROTAMER  NDO   4     C3 - O3   C3M
#=========================================================================


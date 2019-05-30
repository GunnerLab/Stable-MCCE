#2,3-dimethoxy, 1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST NDO        NDOBK NDO01 NDODM

NATOM    NDODM      0
NATOM    NDOBK      0
NATOM    NDO01      26

IATOM    NDO01  C1  0
IATOM    NDO01  O1  1
IATOM    NDO01  C2  2
IATOM    NDO01  C2M 3
IATOM    NDO01  C3  4
IATOM    NDO01  C3M 5
IATOM    NDO01  C4  6
IATOM    NDO01  O4  7
IATOM    NDO01  C5  8
IATOM    NDO01  H5  9
IATOM    NDO01  C6  10
IATOM    NDO01  H6  11
IATOM    NDO01  C7  12
IATOM    NDO01  H7  13
IATOM    NDO01  C8  14
IATOM    NDO01  H8  15
IATOM    NDO01  C9  16
IATOM    NDO01  C10 17
IATOM    NDO01 1H2M 18
IATOM    NDO01 2H2M 19
IATOM    NDO01 3H2M 20
IATOM    NDO01 1H3M 21
IATOM    NDO01 2H3M 22
IATOM    NDO01 3H3M 23
IATOM    NDO01  O2  24
IATOM    NDO01  O3  25


ATOMNAME NDO01    0  C1
ATOMNAME NDO01    1  O1 
ATOMNAME NDO01    2  C2 
ATOMNAME NDO01    3  C2M
ATOMNAME NDO01    4  C3 
ATOMNAME NDO01    5  C3M
ATOMNAME NDO01    6  C4 
ATOMNAME NDO01    7  O4 
ATOMNAME NDO01    8  C5 
ATOMNAME NDO01    9  H5 
ATOMNAME NDO01   10  C6 
ATOMNAME NDO01   11  H6 
ATOMNAME NDO01   12  C7 
ATOMNAME NDO01   13  H7 
ATOMNAME NDO01   14  C8 
ATOMNAME NDO01   15  H8 
ATOMNAME NDO01   16  C9 
ATOMNAME NDO01   17  C10 
ATOMNAME NDO01   18 1H2M
ATOMNAME NDO01   19 2H2M
ATOMNAME NDO01   20 3H2M
ATOMNAME NDO01   21 1H3M
ATOMNAME NDO01   22 2H3M
ATOMNAME NDO01   23 3H3M
ATOMNAME NDO01   24  O2
ATOMNAME NDO01   25  O3

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   NDO01      0

PKA      NDO01      0.0

ELECTRON NDO01      0

EM       NDO01      0.0

RXN      NDO01      -4.190


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  NDO01  C1  sp2       0     O1  0     C2  0     C9
CONNECT  NDO01  O1  s         0     C1
CONNECT  NDO01  C2  sp2       0     O2  0     C1  0     C3
CONNECT  NDO01  C2M sp3       0     O2  0    1H2M 0    2H2M 0    3H2M
CONNECT  NDO01  C3  sp2       0     C2  0     O3  0     C4
CONNECT  NDO01  C3M sp3       0     O3  0    1H3M 0    2H3M 0    3H3M
CONNECT  NDO01  C4  sp2       0     C3  0     O4  0     C10
CONNECT  NDO01  O4  s         0     C4
CONNECT  NDO01  C5  sp2       0     C10 0     H5  0     C6
CONNECT  NDO01  H5  s         0     C5
CONNECT  NDO01  C6  sp2       0     C5  0     C7  0     H6
CONNECT  NDO01  H6  s         0     C6
CONNECT  NDO01  C7  sp2       0     C8  0     C6  0     H7
CONNECT  NDO01  H7  s         0     C7
CONNECT  NDO01  C8  sp2       0     C7  0     C9  0     H8
CONNECT  NDO01  H8  s         0     C8
CONNECT  NDO01  C9  sp2       0     C1  0     C8  0     C10
CONNECT  NDO01  C10 sp2       0     C4  0     C5  0     C9
CONNECT  NDO01 1H2M s         0     C2M
CONNECT  NDO01 2H2M s         0     C2M
CONNECT  NDO01 3H2M s         0     C2M
CONNECT  NDO01 1H3M s         0     C3M
CONNECT  NDO01 2H3M s         0     C3M
CONNECT  NDO01 3H3M s         0     C3M
CONNECT  NDO01  O2  sp2       0     C2M 0     C2
CONNECT  NDO01  O3  sp2       0     C3M 0     C3


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
# opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   NDO01  C1   0.40
CHARGE   NDO01  O1  -0.47
CHARGE   NDO01  C2   0.19
CHARGE   NDO01  C2M  0.08
CHARGE   NDO01  C3   0.05
CHARGE   NDO01  C3M  0.07
CHARGE   NDO01  C4   0.39
CHARGE   NDO01  O4  -0.44
CHARGE   NDO01  C5   0.02
CHARGE   NDO01  H5   0.07
CHARGE   NDO01  C6  -0.13
CHARGE   NDO01  H6   0.12
CHARGE   NDO01  C7  -0.13
CHARGE   NDO01  H7   0.12
CHARGE   NDO01  C8  -0.02
CHARGE   NDO01  H8   0.11
CHARGE   NDO01  C9  -0.05
CHARGE   NDO01  C10 -0.10
CHARGE   NDO01 1H2M  0.07
CHARGE   NDO01 2H2M  0.07
CHARGE   NDO01 3H2M  0.06
CHARGE   NDO01 1H3M  0.07
CHARGE   NDO01 2H3M  0.06 
CHARGE   NDO01 3H3M  0.07
CHARGE   NDO01  O2  -0.34
CHARGE   NDO01  O3  -0.34



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


#9,10-anthraquinone        
#Agnes 10/09  
CONFLIST AQ9        AQ9BK AQ901 AQ9DM

NATOM    AQ9DM      0
NATOM    AQ9BK      0
NATOM    AQ901      24

IATOM    AQ901  C1  0
IATOM    AQ901  H1  1
IATOM    AQ901  C2  2
IATOM    AQ901  H2  3
IATOM    AQ901  C3  4
IATOM    AQ901  H3  5
IATOM    AQ901  C4  6
IATOM    AQ901  H4  7
IATOM    AQ901  C5  8
IATOM    AQ901  H5  9
IATOM    AQ901  C6  10
IATOM    AQ901  H6  11
IATOM    AQ901  C7  12
IATOM    AQ901  H7  13
IATOM    AQ901  C8  14
IATOM    AQ901  H8  15
IATOM    AQ901  C9  16
IATOM    AQ901  C10 17
IATOM    AQ901  O9  18
IATOM    AQ901  O10 19
IATOM    AQ901  C11 20
IATOM    AQ901  C12 21
IATOM    AQ901  C13 22
IATOM    AQ901  C14 23


ATOMNAME AQ901    0  C1
ATOMNAME AQ901    1  H1 
ATOMNAME AQ901    2  C2 
ATOMNAME AQ901    3  H2 
ATOMNAME AQ901    4  C3 
ATOMNAME AQ901    5  H3 
ATOMNAME AQ901    6  C4 
ATOMNAME AQ901    7  H4 
ATOMNAME AQ901    8  C5 
ATOMNAME AQ901    9  H5 
ATOMNAME AQ901   10  C6 
ATOMNAME AQ901   11  H6 
ATOMNAME AQ901   12  C7 
ATOMNAME AQ901   13  H7 
ATOMNAME AQ901   14  C8 
ATOMNAME AQ901   15  H8 
ATOMNAME AQ901   16  C9 
ATOMNAME AQ901   17  C10 
ATOMNAME AQ901   18  O9
ATOMNAME AQ901   19  O10
ATOMNAME AQ901   20  C11
ATOMNAME AQ901   21  C12
ATOMNAME AQ901   22  C13 
ATOMNAME AQ901   23  C14


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   AQ901      0

PKA      AQ901      0.0

ELECTRON AQ901      0

EM       AQ901      0.0

RXN      AQ901      -3.533


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  AQ901  C1  sp2       0     H1  0     C2  0     C13
CONNECT  AQ901  H1  s         0     C1
CONNECT  AQ901  C2  sp2       0     H2  0     C1  0     C3
CONNECT  AQ901  H2  s         0     C2
CONNECT  AQ901  C3  sp2       0     C2  0     H3  0     C4
CONNECT  AQ901  H3  s         0     C3
CONNECT  AQ901  C4  sp2       0     C3  0     H4  0     C14
CONNECT  AQ901  H4  s         0     C4
CONNECT  AQ901  C5  sp2       0     C12 0     H5  0     C6
CONNECT  AQ901  H5  s         0     C5
CONNECT  AQ901  C6  sp2       0     C5  0     C7  0     H6
CONNECT  AQ901  H6  s         0     C6
CONNECT  AQ901  C7  sp2       0     C8  0     C6  0     H7
CONNECT  AQ901  H7  s         0     C7
CONNECT  AQ901  C8  sp2       0     C7  0     C11 0     H8
CONNECT  AQ901  H8  s         0     C8
CONNECT  AQ901  C9  sp2       0     C11 0     O9  0     C13
CONNECT  AQ901  C10 sp2       0     C14 0     O10 0     C12
CONNECT  AQ901  O9  s         0     C9
CONNECT  AQ901  O10 s         0     C10
CONNECT  AQ901  C11 sp2       0     C8  0     C9  0     C12
CONNECT  AQ901  C12 sp2       0     C11 0     C10 0     C5 
CONNECT  AQ901  C13 sp2       0     C1  0     C9  0     C14
CONNECT  AQ901  C14 sp2       0     C4  0     C10 0     C13



#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   AQ9    C1  1.70
RADIUS   AQ9    H1  1.00
RADIUS   AQ9    C2  1.70
RADIUS   AQ9    H2  1.00
RADIUS   AQ9    C3  1.70
RADIUS   AQ9    H3  1.00
RADIUS   AQ9    C4  1.70
RADIUS   AQ9    H4  1.00
RADIUS   AQ9    C5  1.70
RADIUS   AQ9    H5  1.00
RADIUS   AQ9    C6  1.70
RADIUS   AQ9    H6  1.00
RADIUS   AQ9    C7  1.70
RADIUS   AQ9    H7  1.00
RADIUS   AQ9    C8  1.70
RADIUS   AQ9    H8  1.00
RADIUS   AQ9    C9  1.70
RADIUS   AQ9    C10 1.70
RADIUS   AQ9    O9  1.40
RADIUS   AQ9    O10 1.40
RADIUS   AQ9    C11 1.70
RADIUS   AQ9    C12 1.70
RADIUS   AQ9    C13 1.70
RADIUS   AQ9    C14 1.70


#NEUTRAL------
#23456789A123456789B123456789C
# opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   AQ901  C1  -0.04
CHARGE   AQ901  H1   0.10
CHARGE   AQ901  C2  -0.10
CHARGE   AQ901  H2   0.12
CHARGE   AQ901  C3  -0.12
CHARGE   AQ901  H3   0.12
CHARGE   AQ901  C4  -0.02
CHARGE   AQ901  H4   0.10
CHARGE   AQ901  C5  -0.04
CHARGE   AQ901  H5   0.10
CHARGE   AQ901  C6  -0.10
CHARGE   AQ901  H6   0.10 
CHARGE   AQ901  C7  -0.12
CHARGE   AQ901  H7   0.12
CHARGE   AQ901  C8  -0.04
CHARGE   AQ901  H8   0.10
CHARGE   AQ901  C9   0.42
CHARGE   AQ901  C10  0.42
CHARGE   AQ901  O9  -0.48
CHARGE   AQ901  O10 -0.48
CHARGE   AQ901  C11 -0.05
CHARGE   AQ901  C12 -0.03
CHARGE   AQ901  C13 -0.03
CHARGE   AQ901  C14 -0.05






#ParaNam|Res  |Atom|Param/toggle
TRANS    AQ9          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     AQ9   0     C9 - C10- C1
#SPIN     AQ9   1     C8 - C5 - C10
#SPIN     AQ9   2     C1 - C4 - C9

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  AQ9   0     C9 - C10  WHOLE_CONF
ROTAMER  AQ9   1     C11- C13  WHOLE_CONF


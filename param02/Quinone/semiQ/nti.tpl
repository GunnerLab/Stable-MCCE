#2,3,5-trimethyl, 1,4-Naphthoquinone        
#Agnes 10/09  
CONFLIST NTI        NTIBK NTI-1 NTIDM

NATOM    NTIDM      0
NATOM    NTIBK      0
NATOM    NTI-1      27

IATOM    NTI-1  C1  0
IATOM    NTI-1  O1  1
IATOM    NTI-1  C2  2
IATOM    NTI-1  C2M 3
IATOM    NTI-1  C3  4
IATOM    NTI-1  C3M 5
IATOM    NTI-1  C4  6
IATOM    NTI-1  O4  7
IATOM    NTI-1  C5  8
IATOM    NTI-1  C5M 9
IATOM    NTI-1  C6  10
IATOM    NTI-1  H6  11
IATOM    NTI-1  C7  12
IATOM    NTI-1  H7  13
IATOM    NTI-1  C8  14
IATOM    NTI-1  H8  15
IATOM    NTI-1  C9  16
IATOM    NTI-1  C10 17
IATOM    NTI-1 1H2M 18
IATOM    NTI-1 2H2M 19
IATOM    NTI-1 3H2M 20
IATOM    NTI-1 1H3M 21
IATOM    NTI-1 2H3M 22
IATOM    NTI-1 3H3M 23
IATOM    NTI-1 1H5M 24
IATOM    NTI-1 2H5M 25
IATOM    NTI-1 3H5M 26



ATOMNAME NTI-1    0  C1
ATOMNAME NTI-1    1  O1 
ATOMNAME NTI-1    2  C2 
ATOMNAME NTI-1    3  C2M
ATOMNAME NTI-1    4  C3 
ATOMNAME NTI-1    5  C3M
ATOMNAME NTI-1    6  C4 
ATOMNAME NTI-1    7  O4 
ATOMNAME NTI-1    8  C5 
ATOMNAME NTI-1    9  C5M
ATOMNAME NTI-1   10  C6 
ATOMNAME NTI-1   11  H6 
ATOMNAME NTI-1   12  C7 
ATOMNAME NTI-1   13  H7 
ATOMNAME NTI-1   14  C8 
ATOMNAME NTI-1   15  H8 
ATOMNAME NTI-1   16  C9 
ATOMNAME NTI-1   17  C10 
ATOMNAME NTI-1   18 1H2M
ATOMNAME NTI-1   19 2H2M
ATOMNAME NTI-1   20 3H2M
ATOMNAME NTI-1   21 1H3M
ATOMNAME NTI-1   22 2H3M
ATOMNAME NTI-1   23 3H3M
ATOMNAME NTI-1   24 1H5M
ATOMNAME NTI-1   25 2H5M
ATOMNAME NTI-1   26 3H5M

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   NTI-1      0

PKA      NTI-1      0.0

ELECTRON NTI-1      1

EM       NTI-1      0.0

RXN      NTI-1      -15.373


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  NTI-1  C1  sp2       0     O1  0     C2  0     C9
CONNECT  NTI-1  O1  s         0     C1
CONNECT  NTI-1  C2  sp2       0     C2M 0     C1  0     C3
CONNECT  NTI-1  C2M sp3       0     C2  0    1H2M 0    2H2M 0    3H2M
CONNECT  NTI-1  C3  sp2       0     C2  0     C3M 0     C4
CONNECT  NTI-1  C3M sp3       0     C3  0    1H3M 0    2H3M 0    3H3M
CONNECT  NTI-1  C4  sp2       0     C3  0     O4  0     C10
CONNECT  NTI-1  O4  s         0     C4
CONNECT  NTI-1  C5  sp2       0     C10 0     C5M 0     C6
CONNECT  NTI-1  C5M sp3       0     C5  0    1H5M 0    2H5M 0    3H5M
CONNECT  NTI-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  NTI-1  H6  s         0     C6
CONNECT  NTI-1  C7  sp2       0     C8  0     C6  0     H7
CONNECT  NTI-1  H7  s         0     C7
CONNECT  NTI-1  C8  sp2       0     C7  0     C9  0     H8
CONNECT  NTI-1  H8  s         0     C8
CONNECT  NTI-1  C9  sp2       0     C1  0     C8  0     C10
CONNECT  NTI-1  C10 sp2       0     C4  0     C5  0     C9
CONNECT  NTI-1 1H2M s         0     C2M
CONNECT  NTI-1 2H2M s         0     C2M
CONNECT  NTI-1 3H2M s         0     C2M
CONNECT  NTI-1 1H3M s         0     C3M
CONNECT  NTI-1 2H3M s         0     C3M
CONNECT  NTI-1 3H3M s         0     C3M
CONNECT  NTI-1 1H5M s         0     C5M
CONNECT  NTI-1 2H5M s         0     C5M
CONNECT  NTI-1 3H5M s         0     C5M



#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   NTI    C1  1.70
RADIUS   NTI    O1  1.40
RADIUS   NTI    C2  1.70
RADIUS   NTI    C2M 1.70
RADIUS   NTI    C3  1.70
RADIUS   NTI    C3M 1.70
RADIUS   NTI    C4  1.70
RADIUS   NTI    O4  1.40
RADIUS   NTI    C5  1.70
RADIUS   NTI    C5M 1.70
RADIUS   NTI    C6  1.70
RADIUS   NTI    H6  1.00
RADIUS   NTI    C7  1.70
RADIUS   NTI    H7  1.00
RADIUS   NTI    C8  1.70
RADIUS   NTI    H8  1.00
RADIUS   NTI    C9  1.70
RADIUS   NTI    C10 1.70
RADIUS   NTI   1H2M 1.00
RADIUS   NTI   2H2M 1.00
RADIUS   NTI   3H2M 1.00
RADIUS   NTI   1H3M 1.00
RADIUS   NTI   2H3M 1.00
RADIUS   NTI   3H3M 1.00
RADIUS   NTI   1H5M 1.00
RADIUS   NTI   2H5M 1.00
RADIUS   NTI   3H5M 1.00


#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   NTI-1  C1   0.26
CHARGE   NTI-1  O1  -0.59
CHARGE   NTI-1  C2  -0.08
CHARGE   NTI-1  C2M  0.10
CHARGE   NTI-1  C3  -0.19
CHARGE   NTI-1  C3M  0.02
CHARGE   NTI-1  C4   0.43
CHARGE   NTI-1  O4  -0.61
CHARGE   NTI-1  C5   0.29
CHARGE   NTI-1  C5M -0.27
CHARGE   NTI-1  C6  -0.30
CHARGE   NTI-1  H6   0.12
CHARGE   NTI-1  C7  -0.16
CHARGE   NTI-1  H7   0.10
CHARGE   NTI-1  C8  -0.10
CHARGE   NTI-1  H8   0.08
CHARGE   NTI-1  C9   0.05
CHARGE   NTI-1  C10 -0.25
CHARGE   NTI-1 1H2M -0.02
CHARGE   NTI-1 2H2M -0.02
CHARGE   NTI-1 3H2M -0.02
CHARGE   NTI-1 1H3M  0.03
CHARGE   NTI-1 2H3M  0.03
CHARGE   NTI-1 3H3M -0.03
CHARGE   NTI-1 1H5M  0.03
CHARGE   NTI-1 2H5M  0.03
CHARGE   NTI-1 3H5M  0.07



#TORSION  NTI    O2   C2   C1   O1   f        0.0         1    180.00  
#TORSION  NTI    C8   C9   C1   O1   f        0.0         1    180.00  

#TORSION  NTI    C7   C6   C5   C10  f        0.0         1      0.00 
#TORSION  NTI    C5M  C5   C10  C4   f        0.0         1      0.00
#TORSION  NTI    C6   C7   C8   C9   f        0.0         1      0.00

#TORSION  NTI    C7   C8   C9   C1   f        0.0         1      0.00 
#TORSION  NTI    C5   C6   C7   C8   f        0.0         1      0.00
#TORSION  NTI    C5M  C5   C6   C7   f        0.0         1    180.00
#TORSION  NTI    C6   C5   C10  C9   f        0.0         1    180.00 
#TORSION  NTI    C6   C5   C10  C4   f        0.0         1      0.00 

#TORSION  NTI    C5M  C5   C10  C4   f        0.0         1      0.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    NTI          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     NTI   0     C9 - C10- C1
#SPIN     NTI   1     C8 - C5 - C10
#SPIN     NTI   2     C1 - C4 - C9

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  NTI   0     C1 - C4   WHOLE_CONF
ROTAMER  NTI   1     C9 - C10  WHOLE_CONF
ROTAMER  NTI   2     C2 - C7   WHOLE_CONF

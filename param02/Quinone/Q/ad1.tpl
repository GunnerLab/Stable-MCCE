#1,5-dimethoxy,9,10-anthraquinone        
#Agnes 10/09  
CONFLIST AD1        AD1BK AD101 AD1DM

NATOM    AD1DM      0
NATOM    AD1BK      0
NATOM    AD101      32

IATOM    AD101  C1  0
IATOM    AD101  C1M 1
IATOM    AD101  C2  2
IATOM    AD101  H2  3
IATOM    AD101  C3  4
IATOM    AD101  H3  5
IATOM    AD101  C4  6
IATOM    AD101  H4  7
IATOM    AD101  C5  8
IATOM    AD101  O5  9
IATOM    AD101  C6  10
IATOM    AD101  H6  11
IATOM    AD101  C7  12
IATOM    AD101  H7  13
IATOM    AD101  C8  14
IATOM    AD101  H8  15
IATOM    AD101  C9  16
IATOM    AD101  C10 17
IATOM    AD101  O9  18
IATOM    AD101  O10 19
IATOM    AD101  C11 20
IATOM    AD101  C12 21
IATOM    AD101  C13 22
IATOM    AD101  C14 23
IATOM    AD101 1H1M 24
IATOM    AD101 2H1M 25
IATOM    AD101 3H1M 26
IATOM    AD101  O1  27
IATOM    AD101 1H5M 28
IATOM    AD101 2H5M 29
IATOM    AD101 3H5M 30
IATOM    AD101  C5M 31

ATOMNAME AD101    0  C1
ATOMNAME AD101    1  C1M
ATOMNAME AD101    2  C2 
ATOMNAME AD101    3  H2 
ATOMNAME AD101    4  C3 
ATOMNAME AD101    5  H3 
ATOMNAME AD101    6  C4 
ATOMNAME AD101    7  H4 
ATOMNAME AD101    8  C5 
ATOMNAME AD101    9  O5 
ATOMNAME AD101   10  C6 
ATOMNAME AD101   11  H6 
ATOMNAME AD101   12  C7 
ATOMNAME AD101   13  H7 
ATOMNAME AD101   14  C8 
ATOMNAME AD101   15  H8 
ATOMNAME AD101   16  C9 
ATOMNAME AD101   17  C10 
ATOMNAME AD101   18  O9
ATOMNAME AD101   19  O10
ATOMNAME AD101   20  C11
ATOMNAME AD101   21  C12
ATOMNAME AD101   22  C13 
ATOMNAME AD101   23  C14
ATOMNAME AD101   24 1H1M
ATOMNAME AD101   25 2H1M
ATOMNAME AD101   26 3H1M
ATOMNAME AD101   27  O1
ATOMNAME AD101   28 1H5M
ATOMNAME AD101   29 2H5M
ATOMNAME AD101   30 3H5M
ATOMNAME AD101   31  C5M

#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   AD101      0

PKA      AD101      0.0

ELECTRON AD101      0

EM       AD101      0.0

RXN      AD101      -4.589


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  AD101  C1  sp2       0     O1  0     C2  0     C13
CONNECT  AD101  C1M sp3       0     O1  0    1H1M 0    2H1M 0    3H1M
CONNECT  AD101  C2  sp2       0     H2  0     C1  0     C3
CONNECT  AD101  H2  s         0     C2
CONNECT  AD101  C3  sp2       0     C2  0     H3  0     C4
CONNECT  AD101  H3  s         0     C3
CONNECT  AD101  C4  sp2       0     C3  0     H4  0     C14
CONNECT  AD101  H4  s         0     C4
CONNECT  AD101  C5  sp2       0     C12 0     O5  0     C6
CONNECT  AD101  O5  sp2       0     C5  0     C5M
CONNECT  AD101  C6  sp2       0     C5  0     C7  0     H6
CONNECT  AD101  H6  s         0     C6
CONNECT  AD101  C7  sp2       0     C8  0     C6  0     H7
CONNECT  AD101  H7  s         0     C7
CONNECT  AD101  C8  sp2       0     C7  0     C11 0     H8
CONNECT  AD101  H8  s         0     C8
CONNECT  AD101  C9  sp2       0     C11 0     O9  0     C13
CONNECT  AD101  C10 sp2       0     C14 0     O10 0     C12
CONNECT  AD101  O9  s         0     C9
CONNECT  AD101  O10 s         0     C10
CONNECT  AD101  C11 sp2       0     C8  0     C9  0     C12
CONNECT  AD101  C12 sp2       0     C11 0     C10 0     C5 
CONNECT  AD101  C13 sp2       0     C1  0     C9  0     C14
CONNECT  AD101  C14 sp2       0     C4  0     C10 0     C13
CONNECT  AD101 1H1M s         0     C1M
CONNECT  AD101 2H1M s         0     C1M
CONNECT  AD101 3H1M s         0     C1M
CONNECT  AD101  O1  sp2       0     C1M 0     C1
CONNECT  AD101 1H5M s         0     C5M
CONNECT  AD101 2H5M s         0     C5M
CONNECT  AD101 3H5M s         0     C5M
CONNECT  AD101  C5M sp3       0     O5  0    1H5M 0    2H5M 0    3H5M


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   AD1    C1  1.70
RADIUS   AD1    C1M 1.70
RADIUS   AD1    C2  1.70
RADIUS   AD1    H2  1.00
RADIUS   AD1    C3  1.70
RADIUS   AD1    H3  1.00
RADIUS   AD1    C4  1.70
RADIUS   AD1    H4  1.00
RADIUS   AD1    C5  1.70
RADIUS   AD1    O5  1.40
RADIUS   AD1    C6  1.70
RADIUS   AD1    H6  1.00
RADIUS   AD1    C7  1.70
RADIUS   AD1    H7  1.00
RADIUS   AD1    C8  1.70
RADIUS   AD1    H8  1.00
RADIUS   AD1    C9  1.70
RADIUS   AD1    C10 1.70
RADIUS   AD1    O9  1.40
RADIUS   AD1    O10 1.40
RADIUS   AD1    C11 1.70
RADIUS   AD1    C12 1.70
RADIUS   AD1    C13 1.70
RADIUS   AD1    C14 1.70
RADIUS   AD1   1H1M 1.00
RADIUS   AD1   2H1M 1.00
RADIUS   AD1   3H1M 1.00
RADIUS   AD1    O1  1.40
RADIUS   AD1   1H5M 1.00
RADIUS   AD1   2H5M 1.00
RADIUS   AD1   3H5M 1.00
RADIUS   AD1    C5M 1.70

#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   AD101  C1   0.48
CHARGE   AD101  C1M  0.01
CHARGE   AD101  C2  -0.27
CHARGE   AD101  H2   0.15
CHARGE   AD101  C3  -0.10
CHARGE   AD101  H3   0.13
CHARGE   AD101  C4  -0.12
CHARGE   AD101  H4   0.13
CHARGE   AD101  C5   0.48
CHARGE   AD101  O5  -0.38
CHARGE   AD101  C6  -0.27
CHARGE   AD101  H6   0.15
CHARGE   AD101  C7  -0.10
CHARGE   AD101  H7   0.13
CHARGE   AD101  C8  -0.12
CHARGE   AD101  H8   0.13  
CHARGE   AD101  C9   0.53
CHARGE   AD101  C10  0.52
CHARGE   AD101  O9  -0.52
CHARGE   AD101  O10 -0.52
CHARGE   AD101  C11 -0.00
CHARGE   AD101  C12 -0.26
CHARGE   AD101  C13 -0.27
CHARGE   AD101  C14  0.00
CHARGE   AD101 1H1M  0.06
CHARGE   AD101 2H1M  0.06
CHARGE   AD101 3H1M  0.11
CHARGE   AD101  O1  -0.38
CHARGE   AD101 1H5M  0.06
CHARGE   AD101 2H5M  0.06
CHARGE   AD101 3H5M  0.11
CHARGE   AD101  C5M  0.01

TORSION  AD1   1H1M  C1M  O1   C1   f        0.0         1       0.00
TORSION  AD1    C1M  O1   C1   C2   f        1.8         2     180.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    AD1          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     AD1   0     C9 - C10- C1
#SPIN     AD1   1     C8 - C5 - C10
#SPIN     AD1   2     C1 - C4 - C9


#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  AD1   0     C9 - C10  WHOLE_CONF
ROTAMER  AD1   1     C11- C13  WHOLE_CONF
#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  AD1   2     C1 - O1   C1M
ROTAMER  AD1   3     C5 - O5   C5M
#=========================================================================


#2,7-dimethoxy,9,10-anthraquinone        
#Agnes 10/09  
CONFLIST AD7        AD7BK AD7-1 AD7DM

NATOM    AD7DM      0
NATOM    AD7BK      0
NATOM    AD7-1      32

IATOM    AD7-1  C1  0
IATOM    AD7-1  H1  1
IATOM    AD7-1  C2  2
IATOM    AD7-1  C2M 3
IATOM    AD7-1  C3  4
IATOM    AD7-1  H3  5
IATOM    AD7-1  C4  6
IATOM    AD7-1  H4  7
IATOM    AD7-1  C5  8
IATOM    AD7-1  H5  9
IATOM    AD7-1  C6  10
IATOM    AD7-1  H6  11
IATOM    AD7-1  C7  12
IATOM    AD7-1  O7  13
IATOM    AD7-1  C8  14
IATOM    AD7-1  H8  15
IATOM    AD7-1  C9  16
IATOM    AD7-1  C10 17
IATOM    AD7-1  O9  18
IATOM    AD7-1  O10 19
IATOM    AD7-1  C11 20
IATOM    AD7-1  C12 21
IATOM    AD7-1  C13 22
IATOM    AD7-1  C14 23
IATOM    AD7-1 1H2M 24
IATOM    AD7-1 2H2M 25
IATOM    AD7-1 3H2M 26
IATOM    AD7-1  O2  27
IATOM    AD7-1  C7M 28    
IATOM    AD7-1 1H7M 29
IATOM    AD7-1 2H7M 30
IATOM    AD7-1 3H7M 31

ATOMNAME AD7-1    0  C1
ATOMNAME AD7-1    1  H1 
ATOMNAME AD7-1    2  C2 
ATOMNAME AD7-1    3  C2M
ATOMNAME AD7-1    4  C3 
ATOMNAME AD7-1    5  H3 
ATOMNAME AD7-1    6  C4 
ATOMNAME AD7-1    7  H4 
ATOMNAME AD7-1    8  C5 
ATOMNAME AD7-1    9  H5 
ATOMNAME AD7-1   10  C6 
ATOMNAME AD7-1   11  H6 
ATOMNAME AD7-1   12  C7 
ATOMNAME AD7-1   13  O7 
ATOMNAME AD7-1   14  C8 
ATOMNAME AD7-1   15  H8 
ATOMNAME AD7-1   16  C9 
ATOMNAME AD7-1   17  C10 
ATOMNAME AD7-1   18  O9
ATOMNAME AD7-1   19  O10
ATOMNAME AD7-1   20  C11
ATOMNAME AD7-1   21  C12
ATOMNAME AD7-1   22  C13 
ATOMNAME AD7-1   23  C14
ATOMNAME AD7-1   24 1H2M
ATOMNAME AD7-1   25 2H2M
ATOMNAME AD7-1   26 3H2M
ATOMNAME AD7-1   27  O2
ATOMNAME AD7-1   28  C7M
ATOMNAME AD7-1   29 1H7M
ATOMNAME AD7-1   30 2H7M
ATOMNAME AD7-1   31 3H7M



#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   AD7-1      0

PKA      AD7-1      0.0

ELECTRON AD7-1      1

EM       AD7-1      0.0

RXN      AD7-1      -17.221


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  AD7-1  C1  sp2       0     H1  0     C2  0     C13
CONNECT  AD7-1  H1  s         0     C1                                
CONNECT  AD7-1  C2  sp2       0     O2  0     C1  0     C3
CONNECT  AD7-1  C2M sp3       0     O2  0    1H2M 0    2H2M 0    3H2M
CONNECT  AD7-1  C3  sp2       0     C2  0     H3  0     C4
CONNECT  AD7-1  H3  s         0     C3
CONNECT  AD7-1  C4  sp2       0     C3  0     H4  0     C14
CONNECT  AD7-1  H4  s         0     C4
CONNECT  AD7-1  C5  sp2       0     C12 0     H5  0     C6
CONNECT  AD7-1  H5  s         0     C5
CONNECT  AD7-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  AD7-1  H6  s         0     C6
CONNECT  AD7-1  C7  sp2       0     C8  0     C6  0     O7
CONNECT  AD7-1  O7  sp2       0     C7  0     C7M
CONNECT  AD7-1  C8  sp2       0     C7  0     C11 0     H8
CONNECT  AD7-1  H8  s         0     C8
CONNECT  AD7-1  C9  sp2       0     C11 0     O9  0     C13
CONNECT  AD7-1  C10 sp2       0     C14 0     O10 0     C12
CONNECT  AD7-1  O9  s         0     C9
CONNECT  AD7-1  O10 s         0     C10
CONNECT  AD7-1  C11 sp2       0     C8  0     C9  0     C12
CONNECT  AD7-1  C12 sp2       0     C11 0     C10 0     C5 
CONNECT  AD7-1  C13 sp2       0     C1  0     C9  0     C14
CONNECT  AD7-1  C14 sp2       0     C4  0     C10 0     C13
CONNECT  AD7-1 1H2M s         0     C2M
CONNECT  AD7-1 2H2M s         0     C2M
CONNECT  AD7-1 3H2M s         0     C2M
CONNECT  AD7-1  O2  sp2       0     C2M 0     C2
CONNECT  AD7-1 1H7M s         0     C7M
CONNECT  AD7-1 2H7M s         0     C7M
CONNECT  AD7-1 3H7M s         0     C7M
CONNECT  AD7-1  C7M sp3       0     O7  0    1H7M 0    2H7M 0    3H7M

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   AD7    C1  1.70
RADIUS   AD7    H1  1.00
RADIUS   AD7    C2  1.70
RADIUS   AD7    C2M 1.70
RADIUS   AD7    C3  1.70
RADIUS   AD7    H3  1.00
RADIUS   AD7    C4  1.70
RADIUS   AD7    H4  1.00
RADIUS   AD7    C5  1.70
RADIUS   AD7    H5  1.00
RADIUS   AD7    C6  1.70
RADIUS   AD7    H6  1.00
RADIUS   AD7    C7  1.70
RADIUS   AD7    O7  1.40
RADIUS   AD7    C8  1.70
RADIUS   AD7    H8  1.00
RADIUS   AD7    C9  1.70
RADIUS   AD7    C10 1.70
RADIUS   AD7    O9  1.40
RADIUS   AD7    O10 1.40
RADIUS   AD7    C11 1.70
RADIUS   AD7    C12 1.70
RADIUS   AD7    C13 1.70
RADIUS   AD7    C14 1.70
RADIUS   AD7   1H2M 1.00
RADIUS   AD7   2H2M 1.00
RADIUS   AD7   3H2M 1.00
RADIUS   AD7    O2  1.40
RADIUS   AD7   1H7M 1.00
RADIUS   AD7   2H7M 1.00
RADIUS   AD7   3H7M 1.00
RADIUS   AD7    C7M 1.70


#NEUTRAL------
#23456789A123456789B123456789C
#opt ub3lyp/lanl2dz nosymm geom=connectivity pop=chelpg scf(maxcycle=600)  Agnes 10/09
CHARGE   AD7-1  C1  -0.30
CHARGE   AD7-1  H1   0.14
CHARGE   AD7-1  C2   0.45
CHARGE   AD7-1  C2M  0.21
CHARGE   AD7-1  C3  -0.35
CHARGE   AD7-1  H3   0.14
CHARGE   AD7-1  C4  -0.05
CHARGE   AD7-1  H4   0.10 
CHARGE   AD7-1  C5   0.02
CHARGE   AD7-1  H5   0.08
CHARGE   AD7-1  C6  -0.39
CHARGE   AD7-1  H6   0.14
CHARGE   AD7-1  C7   0.45
CHARGE   AD7-1  O7  -0.50
CHARGE   AD7-1  C8  -0.26
CHARGE   AD7-1  H8   0.14
CHARGE   AD7-1  C9   0.29
CHARGE   AD7-1  C10  0.37
CHARGE   AD7-1  O9  -0.61
CHARGE   AD7-1  O10 -0.62
CHARGE   AD7-1  C11 -0.01
CHARGE   AD7-1  C12 -0.12
CHARGE   AD7-1  C13  0.05
CHARGE   AD7-1  C14 -0.12
CHARGE   AD7-1 1H2M  0.01
CHARGE   AD7-1 2H2M  0.01
CHARGE   AD7-1 3H2M  0.03
CHARGE   AD7-1  O2  -0.52
CHARGE   AD7-1 1H7M  0.03
CHARGE   AD7-1 2H7M  0.03
CHARGE   AD7-1 3H7M  0.05
CHARGE   AD7-1  C7M  0.11




#ParaNam|Res  |Atom|Param/toggle
TRANS    AD7          t


#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     AD7   0     C9 - C10- C1
#SPIN     AD7   1     C8 - C5 - C10
#SPIN     AD7   2     C1 - C4 - C9

#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  AD7   0     C9 - C10  WHOLE_CONF
ROTAMER  AD7   1     C11- C13  WHOLE_CONF
#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  AD7   2     C2 - O2   C2M
ROTAMER  AD7   3     C7 - O7   C7M
#=========================================================================


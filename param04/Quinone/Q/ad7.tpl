#2,7-dimethoxy,9,10-anthraquinone        
#Agnes 10/09  
CONFLIST AD7        AD7BK AD701 AD7DM

NATOM    AD7DM      0
NATOM    AD7BK      0
NATOM    AD701      32

IATOM    AD701  C1  0
IATOM    AD701  H1  1
IATOM    AD701  C2  2
IATOM    AD701  C2M 3
IATOM    AD701  C3  4
IATOM    AD701  H3  5
IATOM    AD701  C4  6
IATOM    AD701  H4  7
IATOM    AD701  C5  8
IATOM    AD701  H5  9
IATOM    AD701  C6  10
IATOM    AD701  H6  11
IATOM    AD701  C7  12
IATOM    AD701  O7  13
IATOM    AD701  C8  14
IATOM    AD701  H8  15
IATOM    AD701  C9  16
IATOM    AD701  C10 17
IATOM    AD701  O9  18
IATOM    AD701  O10 19
IATOM    AD701  C11 20
IATOM    AD701  C12 21
IATOM    AD701  C13 22
IATOM    AD701  C14 23
IATOM    AD701 1H2M 24
IATOM    AD701 2H2M 25
IATOM    AD701 3H2M 26
IATOM    AD701  O2  27
IATOM    AD701  C7M 28    
IATOM    AD701 1H7M 29
IATOM    AD701 2H7M 30
IATOM    AD701 3H7M 31

ATOMNAME AD701    0  C1
ATOMNAME AD701    1  H1 
ATOMNAME AD701    2  C2 
ATOMNAME AD701    3  C2M
ATOMNAME AD701    4  C3 
ATOMNAME AD701    5  H3 
ATOMNAME AD701    6  C4 
ATOMNAME AD701    7  H4 
ATOMNAME AD701    8  C5 
ATOMNAME AD701    9  H5 
ATOMNAME AD701   10  C6 
ATOMNAME AD701   11  H6 
ATOMNAME AD701   12  C7 
ATOMNAME AD701   13  O7 
ATOMNAME AD701   14  C8 
ATOMNAME AD701   15  H8 
ATOMNAME AD701   16  C9 
ATOMNAME AD701   17  C10 
ATOMNAME AD701   18  O9
ATOMNAME AD701   19  O10
ATOMNAME AD701   20  C11
ATOMNAME AD701   21  C12
ATOMNAME AD701   22  C13 
ATOMNAME AD701   23  C14
ATOMNAME AD701   24 1H2M
ATOMNAME AD701   25 2H2M
ATOMNAME AD701   26 3H2M
ATOMNAME AD701   27  O2
ATOMNAME AD701   28  C7M
ATOMNAME AD701   29 1H7M
ATOMNAME AD701   30 2H7M
ATOMNAME AD701   31 3H7M



#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   AD701      0

PKA      AD701      0.0

ELECTRON AD701      0

EM       AD701      0.0

RXN      AD701      -5.723


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  AD701  C1  sp2       0     H1  0     C2  0     C13
CONNECT  AD701  H1  s         0     C1                                
CONNECT  AD701  C2  sp2       0     O2  0     C1  0     C3
CONNECT  AD701  C2M sp3       0     O2  0    1H2M 0    2H2M 0    3H2M
CONNECT  AD701  C3  sp2       0     C2  0     H3  0     C4
CONNECT  AD701  H3  s         0     C3
CONNECT  AD701  C4  sp2       0     C3  0     H4  0     C14
CONNECT  AD701  H4  s         0     C4
CONNECT  AD701  C5  sp2       0     C12 0     H5  0     C6
CONNECT  AD701  H5  s         0     C5
CONNECT  AD701  C6  sp2       0     C5  0     C7  0     H6
CONNECT  AD701  H6  s         0     C6
CONNECT  AD701  C7  sp2       0     C6  0     O7  0     C8  
CONNECT  AD701  O7  sp2       0     C7  0     C7M 
CONNECT  AD701  C8  sp2       0     C7  0     C11 0     H8
CONNECT  AD701  H8  s         0     C8
CONNECT  AD701  C9  sp2       0     C11 0     O9  0     C13
CONNECT  AD701  C10 sp2       0     C14 0     O10 0     C12
CONNECT  AD701  O9  s         0     C9
CONNECT  AD701  O10 s         0     C10
CONNECT  AD701  C11 sp2       0     C8  0     C9  0     C12
CONNECT  AD701  C12 sp2       0     C11 0     C10 0     C5 
CONNECT  AD701  C13 sp2       0     C1  0     C9  0     C14
CONNECT  AD701  C14 sp2       0     C4  0     C10 0     C13
CONNECT  AD701 1H2M s         0     C2M
CONNECT  AD701 2H2M s         0     C2M
CONNECT  AD701 3H2M s         0     C2M
CONNECT  AD701  O2  sp2       0     C2M 0     C2
CONNECT  AD701 1H7M s         0     C7M
CONNECT  AD701 2H7M s         0     C7M
CONNECT  AD701 3H7M s         0     C7M
CONNECT  AD701  C7M sp3       0     O7  0    1H7M 0    2H7M 0    3H7M

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
#opt ub3lyp/lanl2dz nosymm geom=connectivity  pop=chelpg scf(maxcycle=600)   Agnes 10/09
CHARGE   AD701  C1  -0.26
CHARGE   AD701  H1   0.15
CHARGE   AD701  C2   0.51
CHARGE   AD701  C2M  0.14
CHARGE   AD701  C3  -0.30
CHARGE   AD701  H3   0.16
CHARGE   AD701  C4  -0.03
CHARGE   AD701  H4   0.12
CHARGE   AD701  C5   0.04
CHARGE   AD701  H5   0.10
CHARGE   AD701  C6  -0.33
CHARGE   AD701  H6   0.17
CHARGE   AD701  C7   0.51
CHARGE   AD701  O7  -0.48
CHARGE   AD701  C8  -0.22
CHARGE   AD701  H8   0.13
CHARGE   AD701  C9   0.48
CHARGE   AD701  C10  0.47
CHARGE   AD701  O9  -0.49
CHARGE   AD701  O10 -0.52
CHARGE   AD701  C11 -0.05
CHARGE   AD701  C12 -0.12
CHARGE   AD701  C13 -0.00
CHARGE   AD701  C14 -0.12
CHARGE   AD701 1H2M  0.04
CHARGE   AD701 2H2M  0.04
CHARGE   AD701 3H2M  0.07
CHARGE   AD701  O2  -0.48
CHARGE   AD701 1H7M  0.05
CHARGE   AD701 2H7M  0.05
CHARGE   AD701 3H7M  0.07
CHARGE   AD701  C7M  0.10




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


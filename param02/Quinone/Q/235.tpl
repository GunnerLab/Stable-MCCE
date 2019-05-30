#2,3-dimethoxy, 5-methyl-Benzoquinone  by Agnes May 2008      
  
CONFLIST 235        235BK 23501 235DM

NATOM    235DM      0
NATOM    235BK      0
NATOM    23501      23

IATOM    23501  C6  0
IATOM    23501  C5  1
IATOM    23501  CM5 2
IATOM    23501  O3  3
IATOM    23501  C4  4
IATOM    23501  C1  5
IATOM    23501  C2  6
IATOM    23501  C3  7
IATOM    23501  CM3 8
IATOM    23501  O1  9
IATOM    23501  O4  10
IATOM    23501  CM2 11
IATOM    23501 2H1  12
IATOM    23501 2H2  13
IATOM    23501 2H3  14
IATOM    23501 3H1  15
IATOM    23501 3H2  16
IATOM    23501 3H3  17
IATOM    23501  O2  18
IATOM    23501  H6  19
IATOM    23501 5H1  20
IATOM    23501 5H2  21
IATOM    23501 5H3  22


ATOMNAME 23501    0  C6 
ATOMNAME 23501    1  C5 
ATOMNAME 23501    2  CM5
ATOMNAME 23501    3  O3 
ATOMNAME 23501    4  C4 
ATOMNAME 23501    5  C1 
ATOMNAME 23501    6  C2 
ATOMNAME 23501    7  C3 
ATOMNAME 23501    8  CM3 
ATOMNAME 23501    9  O1 
ATOMNAME 23501   10  O4 
ATOMNAME 23501   11  CM2 
ATOMNAME 23501   12 2H1
ATOMNAME 23501   13 2H2
ATOMNAME 23501   14 2H3
ATOMNAME 23501   15 3H1
ATOMNAME 23501   16 3H2
ATOMNAME 23501   17 3H3
ATOMNAME 23501   18  O2
ATOMNAME 23501   19  H6
ATOMNAME 23501   20 5H1
ATOMNAME 23501   21 5H2
ATOMNAME 23501   22 5H3


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   23501      0

PKA      23501      0.0

ELECTRON 23501      0

EM       23501      0.0

RXN      23501      -2.859


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  23501  C6  sp2        0    C5   0    C1   0    H6
CONNECT  23501  C5  sp2        0    C6   0    C4   0    CM5
CONNECT  23501  CM5 sp3        0    C5   0   5H1   0   5H2   0   5H3
CONNECT  23501  O3  sp2        0    C3   0    CM3
CONNECT  23501  C4  sp2        0    C5   0    C3   0    O4
CONNECT  23501  C1  sp2        0    C6   0    C2   0    O1
CONNECT  23501  C2  sp2        0    C1   0    C3   0    O2
CONNECT  23501  C3  sp2        0    C4   0    C2   0    O3
CONNECT  23501  CM3 sp3        0    O3   0   3H1   0   3H2   0   3H3  
CONNECT  23501  CM2 sp3        0    O2   0   2H1   0   2H2   0   2H3   
CONNECT  23501  O1  s          0    C1
CONNECT  23501  O4  s          0    C4
CONNECT  23501 2H1  s          0    CM2
CONNECT  23501 2H2  s          0    CM2
CONNECT  23501 2H3  s          0    CM2
CONNECT  23501 3H1  s          0    CM3
CONNECT  23501 3H2  s          0    CM3
CONNECT  23501 3H3  s          0    CM3
CONNECT  23501  O2  sp2        0    C2   0    CM2
CONNECT  23501  H6  s          0    C6   
CONNECT  23501 5H1  s          0    CM5  
CONNECT  23501 5H2  s          0    CM5   
CONNECT  23501 5H3  s          0    CM5   


#de-protonated---------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   235    C6  1.70
RADIUS   235    C5  1.70
RADIUS   235    CM5 1.70
RADIUS   235    O3  1.40
RADIUS   235    C4  1.70
RADIUS   235    C1  1.70
RADIUS   235    C2  1.70
RADIUS   235    C3  1.70
RADIUS   235    CM3 1.70
RADIUS   235    O1  1.40
RADIUS   235    O4  1.40
RADIUS   235    CM2 1.70
RADIUS   235   2H1  1.00
RADIUS   235   2H2  1.00
RADIUS   235   2H3  1.00
RADIUS   235   3H1  1.00
RADIUS   235   3H2  1.00
RADIUS   235   3H3  1.00
RADIUS   235    O2  1.40
RADIUS   235    H6  1.00
RADIUS   235   5H1  1.00
RADIUS   235   5H2  1.00
RADIUS   235   5H3  1.00



#NEUTRAL------
#23456789A123456789B123456789C
#p ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=800)   Agnes May, 2008
CHARGE   23501  O1  -0.49
CHARGE   23501  C1   0.58
CHARGE   23501  C2   0.07
CHARGE   23501  O2  -0.25
CHARGE   23501  CM2  0.05
CHARGE   23501  C3   0.04
CHARGE   23501  CM3  0.15
CHARGE   23501  C4   0.44
CHARGE   23501  O4  -0.43
CHARGE   23501  C5   0.02
CHARGE   23501  CM5  0.02
CHARGE   23501 2H1   0.07
CHARGE   23501 2H2   0.06
CHARGE   23501 2H3   0.05
CHARGE   23501 3H1   0.05 
CHARGE   23501 3H2   0.03 
CHARGE   23501 3H3   0.01 
CHARGE   23501  O3  -0.32
CHARGE   23501 5H1   0.02 
CHARGE   23501 5H2   0.03 
CHARGE   23501 5H3   0.02 
CHARGE   23501  C6  -0.40
CHARGE   23501  H6   0.18
 

#ParaNam|Res  |Atom|Param/toggle
TRANS    235          t

TORSION  23501  CM2  O2   C2   C3   f      0.0           2    180.00
TORSION  23501  CM3  O3   C3   C2   f      0.0           2    180.00
TORSION  23501  CM2  O2   C2   C1   f      0.215         6    180.00     -0.195        4    180.00    -0.427         3    180.00     4.208         2    180.00    -0.632         1    180.00
TORSION  23501  2H1  CM2  O2   C2   f     -0.180         3      0.00
TORSION  23501  3H1  CM3  O3   C3   f     -0.180         3      0.00
TORSION  23501  5H1  CM5  C5   C6   f     -0.180         3      0.00


#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  235   0     C2 - O2   CM2
ROTAMER  235   1     C3 - O3   CM3
#=========================================================================


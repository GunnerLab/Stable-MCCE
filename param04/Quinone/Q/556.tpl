#2,3-dimethoxy, 5,6-methyl-Benzoquinone  by Agnes May 2008      
  
CONFLIST 556        556BK 55601 556DM

NATOM    556DM      0
NATOM    556BK      0
NATOM    55601      26

IATOM    55601  C6  0
IATOM    55601  C5  1
IATOM    55601  CM5 2
IATOM    55601  O3  3
IATOM    55601  C4  4
IATOM    55601  C1  5
IATOM    55601  C2  6
IATOM    55601  C3  7
IATOM    55601  CM3 8
IATOM    55601  O1  9
IATOM    55601  O4  10
IATOM    55601  CM2 11
IATOM    55601 2H1  12
IATOM    55601 2H2  13
IATOM    55601 2H3  14
IATOM    55601 3H1  15
IATOM    55601 3H2  16
IATOM    55601 3H3  17
IATOM    55601  O2  18
IATOM    55601  CM6 19
IATOM    55601 5H1  20
IATOM    55601 5H2  21
IATOM    55601 5H3  22
IATOM    55601 6H1  23
IATOM    55601 6H2  24
IATOM    55601 6H3  25


ATOMNAME 55601    0  C6 
ATOMNAME 55601    1  C5 
ATOMNAME 55601    2  CM5
ATOMNAME 55601    3  O3 
ATOMNAME 55601    4  C4 
ATOMNAME 55601    5  C1 
ATOMNAME 55601    6  C2 
ATOMNAME 55601    7  C3 
ATOMNAME 55601    8  CM3 
ATOMNAME 55601    9  O1 
ATOMNAME 55601   10  O4 
ATOMNAME 55601   11  CM2 
ATOMNAME 55601   12 2H1
ATOMNAME 55601   13 2H2
ATOMNAME 55601   14 2H3
ATOMNAME 55601   15 3H1
ATOMNAME 55601   16 3H2
ATOMNAME 55601   17 3H3
ATOMNAME 55601   18  O2
ATOMNAME 55601   19  CM6
ATOMNAME 55601   20 5H1
ATOMNAME 55601   21 5H2
ATOMNAME 55601   22 5H3
ATOMNAME 55601   23 6H1
ATOMNAME 55601   24 6H2
ATOMNAME 55601   25 6H3


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   55601      0

PKA      55601      0.0

ELECTRON 55601      0

EM       55601      0.0

RXN      55601      -2.859


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  55601  C6  sp2        0    C5   0    C1   0    CM6
CONNECT  55601  C5  sp2        0    C6   0    C4   0    CM5
CONNECT  55601  CM5 sp3        0    C5   0   5H1   0   5H2   0   5H3
CONNECT  55601  O3  sp2        0    C3   0    CM3
CONNECT  55601  C4  sp2        0    C5   0    C3   0    O4
CONNECT  55601  C1  sp2        0    C6   0    C2   0    O1
CONNECT  55601  C2  sp2        0    C1   0    C3   0    O2
CONNECT  55601  C3  sp2        0    C4   0    C2   0    O3
CONNECT  55601  CM3 sp3        0    O3   0   3H1   0   3H2   0   3H3  
CONNECT  55601  CM2 sp3        0    O2   0   2H1   0   2H2   0   2H3   
CONNECT  55601  O1  s          0    C1
CONNECT  55601  O4  s          0    C4
CONNECT  55601 2H1  s          0    CM2
CONNECT  55601 2H2  s          0    CM2
CONNECT  55601 2H3  s          0    CM2
CONNECT  55601 3H1  s          0    CM3
CONNECT  55601 3H2  s          0    CM3
CONNECT  55601 3H3  s          0    CM3
CONNECT  55601  O2  sp2        0    C2   0    CM2
CONNECT  55601  CM6 sp3        0    C6   0   6H1   0   6H2   0   6H3
CONNECT  55601 5H1  s          0    CM5 
CONNECT  55601 5H2  s          0    CM5  
CONNECT  55601 5H3  s          0    CM5  
CONNECT  55601 6H1  s          0    CM6 
CONNECT  55601 6H2  s          0    CM6  
CONNECT  55601 6H3  s          0    CM6  


#de-protonated---------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   556    C6  1.70
RADIUS   556    C5  1.70
RADIUS   556    CM5 1.70
RADIUS   556    O3  1.40
RADIUS   556    C4  1.70
RADIUS   556    C1  1.70
RADIUS   556    C2  1.70
RADIUS   556    C3  1.70
RADIUS   556    CM3 1.70
RADIUS   556    O1  1.40
RADIUS   556    O4  1.40
RADIUS   556    CM2 1.70
RADIUS   556   2H1  1.00
RADIUS   556   2H2  1.00
RADIUS   556   2H3  1.00
RADIUS   556   3H1  1.00
RADIUS   556   3H2  1.00
RADIUS   556   3H3  1.00
RADIUS   556    O2  1.40
RADIUS   556    CM6 1.70
RADIUS   556   5H1  1.00
RADIUS   556   5H2  1.00
RADIUS   556   5H3  1.00
RADIUS   556   6H1  1.00
RADIUS   556   6H2  1.00
RADIUS   556   6H3  1.00



#NEUTRAL------
#23456789A123456789B123456789C
#p ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=800)   Agnes May, 2008
CHARGE   55601  O1  -0.41
CHARGE   55601  C1   0.40
CHARGE   55601  C2   0.08
CHARGE   55601  O2  -0.32
CHARGE   55601  CM2  0.17
CHARGE   55601  C3   0.09
CHARGE   55601  CM3  0.16
CHARGE   55601  C4   0.42
CHARGE   55601  O4  -0.41
CHARGE   55601  C5  -0.12
CHARGE   55601  CM5  0.08
CHARGE   55601 2H1   0.02
CHARGE   55601 2H2   0.02
CHARGE   55601 2H3   0.04
CHARGE   55601 3H1   0.04 
CHARGE   55601 3H2   0.02 
CHARGE   55601 3H3   0.02 
CHARGE   55601  O3  -0.32
CHARGE   55601 5H1   0.01 
CHARGE   55601 5H2   0.01 
CHARGE   55601 5H3   0.01 
CHARGE   55601  C6  -0.14
CHARGE   55601  CM6  0.10
CHARGE   55601 6H1   0.01
CHARGE   55601 6H2   0.01
CHARGE   55601 6H3   0.01
 
 

#ParaNam|Res  |Atom|Param/toggle
TRANS    556          t

#TORSION  55601  CM2  O2   C2   C1   f        0.0         2    0.00
#TORSION  55601  CM3  O3   C3   C4   f        0.0         2    0.00

TORSION  55601  CM2  O2   C2   C3   f      0.0           2    180.00
TORSION  55601  CM3  O3   C3   C2   f      0.0           2    180.00
TORSION  55601  CM2  O2   C2   C1   f      0.215         6    180.00     -0.195        4    180.00    -0.427         3    180.00     4.208         2    180.00    -0.632         1    180.00
TORSION  55601  CM3  O3   C3   C4   f      0.215         6    180.00     -0.195        4    180.00    -0.427         3    180.00     4.208         2    180.00    -0.632         1    180.00
TORSION  55601  2H1  CM2  O2   C2   f     -0.180         3      0.00
TORSION  55601  3H1  CM3  O3   C3   f     -0.180         3      0.00
TORSION  55601  6H1  CM6  C6   C5   f     -0.180         3      0.00
TORSION  55601  5H1  CM5  C5   C6   f     -0.180         3      0.00


#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  556   0     C3 - O3   CM3
ROTAMER  556   1     C2 - O2   CM2
#=========================================================================



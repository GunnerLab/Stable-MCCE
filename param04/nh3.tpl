CONFLIST NH3        NH3BK NH3+1 NH301 NH3-1 NH3DM

NATOM    NH3BK      0
NATOM    NH3+1      5
NATOM    NH301      4
NATOM    NH3-1      3
NATOM    NH3DM      0

IATOM    NH3+1  N   0
IATOM    NH3+1 1H   1
IATOM    NH3+1 2H   2
IATOM    NH3+1 3H   3
IATOM    NH3+1 4H   4

IATOM    NH301  N   0
IATOM    NH301 1H   1
IATOM    NH301 2H   2
IATOM    NH301 3H   3

IATOM    NH3-1  N   0
IATOM    NH3-1 1H   1
IATOM    NH3-1 2H   2


ATOMNAME NH3+1    0  N  
ATOMNAME NH3+1    1 1H  
ATOMNAME NH3+1    2 2H  
ATOMNAME NH3+1    3 3H  
ATOMNAME NH3+1    4 4H  

ATOMNAME NH301    0  N
ATOMNAME NH301    1 1H
ATOMNAME NH301    2 2H
ATOMNAME NH301    3 3H

ATOMNAME NH3-1    0  N
ATOMNAME NH3-1    1 1H
ATOMNAME NH3-1    2 2H

#Created using mcce2.5.1
#ammonia species
#Done by Witol Szejgis and Manoj Mandal

#-------|-----|----|----

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
RXN      NH3+1      -23.44
PROTON   NH3+1      1
ELECTRON NH3+1      0
RXN      NH301      -2.682
PROTON   NH301      0
ELECTRON NH301      0
RXN      NH3-1      -28.253
PROTON   NH3-1      -1
ELECTRON NH3-1      1
RXN      NH3DM      0
PROTON   NH3DM      0
ELECTRON NH3DM      0
PKA      NH3+1      9.25
PKA      NH301      0
PKA      NH3-1      36
PKA      NH3DM      0

#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  NH3+1  N   sp3       0    1H   0    2H   0    3H   0    4H           
CONNECT  NH3+1 1H   s         0     N        
CONNECT  NH3+1 2H   s         0     N        
CONNECT  NH3+1 3H   s         0     N        
CONNECT  NH3+1 4H   s         0     N         

CONNECT  NH301  N   sp3       0    1H   0    2H   0    3H   
CONNECT  NH301 1H   s         0     N
CONNECT  NH301 2H   s         0     N
CONNECT  NH301 3H   s         0     N

CONNECT  NH3-1  N   sp3       0    1H   0    2H   
CONNECT  NH3-1 1H   s         0     N
CONNECT  NH3-1 2H   s         0     N


CHARGE   NH3+1  N     -0.744
CHARGE   NH3+1  1H    0.436
CHARGE   NH3+1  2H    0.436
CHARGE   NH3+1  3H    0.436
CHARGE   NH3+1  4H    0.436

CHARGE   NH301  N     -1.125
CHARGE   NH301  1H    0.375
CHARGE   NH301  2H    0.375
CHARGE   NH301  3H    0.375

CHARGE   NH3-1  N     -1.50
CHARGE   NH3-1  1H    0.250
CHARGE   NH3-1  2H    0.250

#Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   NH3    N      1.50
RADIUS   NH3    1H     1.00
RADIUS   NH3    2H     1.00
RADIUS   NH3    3H     1.00
RADIUS   NH3    4H     1.00

#........#.....
TRANS    NH3         t

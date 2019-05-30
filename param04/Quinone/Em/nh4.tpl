CONFLIST NH4        NH4BK NH4+1 

NATOM    NH4BK      0
NATOM    NH4+1      5

IATOM    NH4+1  N   0
IATOM    NH4+1 1H   1
IATOM    NH4+1 2H   2
IATOM    NH4+1 3H   3
IATOM    NH4+1 4H   4

ATOMNAME NH4+1    0  N  
ATOMNAME NH4+1    1 1H  
ATOMNAME NH4+1    2 2H  
ATOMNAME NH4+1    3 3H  
ATOMNAME NH4+1    4 4H  
#Created using mcce2.4
#ammonium ion +1 charge
#Done by Gail 1/31/06

#-------|-----|----|----





#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
RXN      NH4+1      -23.579
PROTON   NH4+1      1

#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  NH4+1  N   sp3       0    1H   0    2H   0    3H   0    4H           
CONNECT  NH4+1 1H   s         0     N        
CONNECT  NH4+1 2H   s         0     N        
CONNECT  NH4+1 3H   s         0     N        
CONNECT  NH4+1 4H   s         0     N         

CHARGE   NH4+1  N     -0.744
CHARGE   NH4+1  1H    0.436
CHARGE   NH4+1  2H    0.436
CHARGE   NH4+1  3H    0.436
CHARGE   NH4+1  4H    0.436

#Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   NH4    N     1.50
RADIUS   NH4    1H     1.00
RADIUS   NH4    2H     1.00
RADIUS   NH4    3H     1.00
RADIUS   NH4    4H     1.00

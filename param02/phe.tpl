CONFLIST PHE        PHEBK PHE01

NATOM    PHEBK      6
NATOM    PHE01      14

IATOM    PHEBK  N   0
IATOM    PHEBK  H   1
IATOM    PHEBK  CA  2
IATOM    PHEBK  HA  3
IATOM    PHEBK  C   4
IATOM    PHEBK  O   5
IATOM    PHE01  CB  0
IATOM    PHE01 1HB  1
IATOM    PHE01 2HB  2
IATOM    PHE01  CG  3
IATOM    PHE01  CD1 4
IATOM    PHE01  HD1 5
IATOM    PHE01  CD2 6
IATOM    PHE01  HD2 7
IATOM    PHE01  CE1 8
IATOM    PHE01  HE1 9
IATOM    PHE01  CE2 10
IATOM    PHE01  HE2 11
IATOM    PHE01  CZ  12
IATOM    PHE01  HZ  13

ATOMNAME PHEBK    0  N  
ATOMNAME PHEBK    1  H  
ATOMNAME PHEBK    2  CA 
ATOMNAME PHEBK    3  HA 
ATOMNAME PHEBK    4  C  
ATOMNAME PHEBK    5  O  
ATOMNAME PHE01    0  CB 
ATOMNAME PHE01    1 1HB 
ATOMNAME PHE01    2 2HB 
ATOMNAME PHE01    3  CG 
ATOMNAME PHE01    4  CD1
ATOMNAME PHE01    5  HD1
ATOMNAME PHE01    6  CD2
ATOMNAME PHE01    7  HD2
ATOMNAME PHE01    8  CE1
ATOMNAME PHE01    9  HE1
ATOMNAME PHE01   10  CE2
ATOMNAME PHE01   11  HE2
ATOMNAME PHE01   12  CZ 
ATOMNAME PHE01   13  HZ






#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   PHE01      0
PKA      PHE01      0.0
ELECTRON PHE01      0
EM       PHE01      0.0
RXN      PHE01      0.00

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  PHEBK  N   sp2       -1    C   0     CA  0     H
CONNECT  PHEBK  H   s         0     N
CONNECT  PHEBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  PHEBK  HA  s         0     CA
CONNECT  PHEBK  C   sp2       0     CA  0     O   1     N
CONNECT  PHEBK  O   sp2       0     C
CONNECT  PHE01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  PHE01 1HB  s         0     CB
CONNECT  PHE01 2HB  s         0     CB
CONNECT  PHE01  CG  sp2       0     CB  0     CD1 0     CD2
CONNECT  PHE01  CD1 sp2       0     CG  0     CE1 0     HD1
CONNECT  PHE01  HD1 s         0     CD1
CONNECT  PHE01  CD2 sp2       0     CG  0     CE2 0     HD2
CONNECT  PHE01  HD2 s         0     CD2
CONNECT  PHE01  CE1 sp2       0     CD1 0     CZ  0     HE1
CONNECT  PHE01  HE1 s         0     CE1
CONNECT  PHE01  CE2 sp2       0     CD2 0     CZ  0     HE2
CONNECT  PHE01  HE2 s         0     CE2
CONNECT  PHE01  CZ  sp2       0     CE1 0     CE2 0     HZ
CONNECT  PHE01  HZ  s         0     CZ



#3.Atom Parameters: Partial Charges and Radii
CHARGE   PHEBK  N    -0.350
CHARGE   PHEBK  H     0.250
CHARGE   PHEBK  CA    0.100
CHARGE   PHEBK  C     0.550
CHARGE   PHEBK  O    -0.550

# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   PHE    N   1.50
RADIUS   PHE    H   1.00
RADIUS   PHE    CA  2.00
RADIUS   PHE    HA  0.00
RADIUS   PHE    C   1.70
RADIUS   PHE    O   1.40
RADIUS   PHE    CB  2.00
RADIUS   PHE   1HB  0.00
RADIUS   PHE   2HB  0.00
RADIUS   PHE    CG  1.70
RADIUS   PHE    CD1 1.70
RADIUS   PHE    HD1 1.00
RADIUS   PHE    CD2 1.70
RADIUS   PHE    HD2 1.00
RADIUS   PHE    CE1 1.70
RADIUS   PHE    HE1 1.00
RADIUS   PHE    CE2 1.70
RADIUS   PHE    HE2 1.00
RADIUS   PHE    CZ  1.70
RADIUS   PHE    HZ  1.00


#4.Rotomer
#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  PHE   0     CA - CB   CG   CD1  CD2  CE1  CE2  CZ
ROTAMER  PHE   1     CB - CG   CD1  CD2  CE1  CE2  CZ
#=========================================================================



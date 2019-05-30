#CYS making disulfide bond to another CYS
CONFLIST CYD        CYDBK CYD01

NATOM    CYDBK      6
NATOM    CYD01      4

IATOM    CYDBK  N   0
IATOM    CYDBK  H   1
IATOM    CYDBK  CA  2
IATOM    CYDBK  HA  3
IATOM    CYDBK  C   4
IATOM    CYDBK  O   5
IATOM    CYD01  CB  0
IATOM    CYD01 1HB  1
IATOM    CYD01 2HB  2
IATOM    CYD01  SG  3

ATOMNAME CYDBK    0  N  
ATOMNAME CYDBK    1  H  
ATOMNAME CYDBK    2  CA 
ATOMNAME CYDBK    3  HA 
ATOMNAME CYDBK    4  C  
ATOMNAME CYDBK    5  O  
ATOMNAME CYD01    0  CB 
ATOMNAME CYD01    1 1HB 
ATOMNAME CYD01    2 2HB 
ATOMNAME CYD01    3  SG 

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   CYD01      0
PKA      CYD01      0.0
ELECTRON CYD01      0
EM       CYD01      0.0
EM       CYD01      0.0
RXN      CYD01      0.0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  CYDBK  N   sp2       -1    C   0     CA  0     H
CONNECT  CYDBK  H   s         0     N
CONNECT  CYDBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  CYDBK  HA  s         0     CA
CONNECT  CYDBK  C   sp2       0     CA  0     O   1     N
CONNECT  CYDBK  O   sp2       0     C

CONNECT  CYD01  CB  sp3       0     CA  0     SG  0    1HB  0    2HB
CONNECT  CYD01 1HB  s         0     CB
CONNECT  CYD01 2HB  s         0     CB
CONNECT  CYD01  SG  sp2       0     CB  LIG   SG

#3.Atom Parameters: Partial Charges and Radii
CHARGE   CYDBK  N    -0.350
CHARGE   CYDBK  H     0.250
CHARGE   CYDBK  CA    0.100
CHARGE   CYDBK  C     0.550
CHARGE   CYDBK  O    -0.550

# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   CYD    N   1.55
RADIUS   CYD    H   1.20
RADIUS   CYD    CA  1.70
RADIUS   CYD    HA  1.20
RADIUS   CYD    C   1.70
RADIUS   CYD    O   1.52
RADIUS   CYD    CB  1.70
RADIUS   CYD   1HB  1.20
RADIUS   CYD   2HB  1.20
RADIUS   CYD    SG  1.80

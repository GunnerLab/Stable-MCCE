###Aspartate

CONFLIST ASP        ASPBK ASP01 ASP02 ASP-1

NATOM    ASPBK      6
NATOM    ASP01      7
NATOM    ASP02      7
NATOM    ASP-1      6

IATOM    ASPBK  N   0
IATOM    ASPBK  H   1
IATOM    ASPBK  CA  2
IATOM    ASPBK  HA  3
IATOM    ASPBK  C   4
IATOM    ASPBK  O   5
IATOM    ASP01  CB  0
IATOM    ASP01 1HB  1
IATOM    ASP01 2HB  2
IATOM    ASP01  CG  3
IATOM    ASP01  OD1 4
IATOM    ASP01  HD1 5
IATOM    ASP01  OD2 6
IATOM    ASP02  CB  0
IATOM    ASP02 1HB  1
IATOM    ASP02 2HB  2
IATOM    ASP02  CG  3
IATOM    ASP02  OD1 4
IATOM    ASP02  OD2 5
IATOM    ASP02  HD2 6
IATOM    ASP-1  CB  0
IATOM    ASP-1 1HB  1
IATOM    ASP-1 2HB  2
IATOM    ASP-1  CG  3
IATOM    ASP-1  OD1 4
IATOM    ASP-1  OD2 5

ATOMNAME ASPBK    0  N  
ATOMNAME ASPBK    1  H  
ATOMNAME ASPBK    2  CA 
ATOMNAME ASPBK    3  HA 
ATOMNAME ASPBK    4  C  
ATOMNAME ASPBK    5  O  
ATOMNAME ASP01    0  CB 
ATOMNAME ASP01    1 1HB 
ATOMNAME ASP01    2 2HB 
ATOMNAME ASP01    3  CG 
ATOMNAME ASP01    4  OD1
ATOMNAME ASP01    5  HD1
ATOMNAME ASP01    6  OD2
ATOMNAME ASP02    0  CB 
ATOMNAME ASP02    1 1HB 
ATOMNAME ASP02    2 2HB 
ATOMNAME ASP02    3  CG 
ATOMNAME ASP02    4  OD1
ATOMNAME ASP02    5  OD2
ATOMNAME ASP02    6  HD2
ATOMNAME ASP-1    0  CB 
ATOMNAME ASP-1    1 1HB 
ATOMNAME ASP-1    2 2HB 
ATOMNAME ASP-1    3  CG 
ATOMNAME ASP-1    4  OD1
ATOMNAME ASP-1    5  OD2

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   ASP01      0
PROTON   ASP02      0
PROTON   ASP-1      -1
PKA      ASP01      0.0
PKA      ASP02      0.0
PKA      ASP-1      4.75
ELECTRON ASP01      0
ELECTRON ASP02      0
ELECTRON ASP-1      0
EM       ASP01      0.0
EM       ASP02      0.0
EM       ASP-1      0.0
RXN      ASP01      -6.32
RXN      ASP02      -6.24
RXN      ASP-1      -41.68

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  ASPBK  N   sp2       -1    C   0     CA  0     H
CONNECT  ASPBK  H   s         0     N
CONNECT  ASPBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  ASPBK  HA  s         0     CA
CONNECT  ASPBK  C   sp2       0     CA  0     O   1     N
CONNECT  ASPBK  O   sp2       0     C
CONNECT  ASP01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  ASP01 1HB  s         0     CB
CONNECT  ASP01 2HB  s         0     CB
CONNECT  ASP01  CG  sp2       0     CB  0     OD1 0     OD2
CONNECT  ASP01  OD1 sp3       0     CG  0     HD1
CONNECT  ASP01  HD1 s         0     OD1
CONNECT  ASP01  OD2 sp2       0     CG
CONNECT  ASP02  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  ASP02 1HB  s         0     CB
CONNECT  ASP02 2HB  s         0     CB
CONNECT  ASP02  CG  sp2       0     CB  0     OD1 0     OD2
CONNECT  ASP02  OD1 sp2       0     CG
CONNECT  ASP02  OD2 sp3       0     CG  0     HD2
CONNECT  ASP02  HD2 s         0     OD2
CONNECT  ASP-1  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  ASP-1 1HB  s         0     CB
CONNECT  ASP-1 2HB  s         0     CB
CONNECT  ASP-1  CG  sp2       0     CB  0     OD1 0     OD2
CONNECT  ASP-1  OD1 sp2       0     CG
CONNECT  ASP-1  OD2 sp2       0     CG

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   ASP    N   1.50
RADIUS   ASP    H   1.00
RADIUS   ASP    CA  2.00
RADIUS   ASP    HA  0.00
RADIUS   ASP    C   1.70
RADIUS   ASP    O   1.40
RADIUS   ASP    CB  2.00
RADIUS   ASP   1HB  0.00
RADIUS   ASP   2HB  0.00
RADIUS   ASP    CG  1.70
RADIUS   ASP    OD1 1.40
RADIUS   ASP    HD1 1.00
RADIUS   ASP    OD2 1.40
RADIUS   ASP    HD2 1.00


CHARGE   ASPBK  N    -0.350
CHARGE   ASPBK  H     0.250
CHARGE   ASPBK  CA    0.100
CHARGE   ASPBK  C     0.550
CHARGE   ASPBK  O    -0.550
CHARGE   ASP-1  CG    0.100
CHARGE   ASP-1  OD1  -0.550
CHARGE   ASP-1  OD2  -0.550
CHARGE   ASP01  CG    0.550
CHARGE   ASP01  OD1  -0.495
CHARGE   ASP01  OD2  -0.490
CHARGE   ASP01  HD1   0.435
CHARGE   ASP02  CG    0.550
CHARGE   ASP02  OD1  -0.490
CHARGE   ASP02  OD2  -0.495
CHARGE   ASP02  HD2   0.435

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  ASP   0     CA - CB   CG   OD1  OD2
ROTAMER  ASP   1     CB - CG   OD1  OD2
#=========================================================================

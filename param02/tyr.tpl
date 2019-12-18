CONFLIST TYR        TYRBK TYR01 TYR-1

NATOM    TYRBK      6
NATOM    TYR01      15
NATOM    TYR-1      14

IATOM    TYRBK  N   0
IATOM    TYRBK  H   1
IATOM    TYRBK  CA  2
IATOM    TYRBK  HA  3
IATOM    TYRBK  C   4
IATOM    TYRBK  O   5
IATOM    TYR01  CB  0
IATOM    TYR01 1HB  1
IATOM    TYR01 2HB  2
IATOM    TYR01  CG  3
IATOM    TYR01  CD1 4
IATOM    TYR01  HD1 5
IATOM    TYR01  CD2 6
IATOM    TYR01  HD2 7
IATOM    TYR01  CE1 8
IATOM    TYR01  HE1 9
IATOM    TYR01  CE2 10
IATOM    TYR01  HE2 11
IATOM    TYR01  CZ  12
IATOM    TYR01  OH  13
IATOM    TYR01  HH  14
IATOM    TYR-1  CB  0
IATOM    TYR-1 1HB  1
IATOM    TYR-1 2HB  2
IATOM    TYR-1  CG  3
IATOM    TYR-1  CD1 4
IATOM    TYR-1  HD1 5
IATOM    TYR-1  CD2 6
IATOM    TYR-1  HD2 7
IATOM    TYR-1  CE1 8
IATOM    TYR-1  HE1 9
IATOM    TYR-1  CE2 10
IATOM    TYR-1  HE2 11
IATOM    TYR-1  CZ  12
IATOM    TYR-1  OH  13

ATOMNAME TYRBK    0  N  
ATOMNAME TYRBK    1  H  
ATOMNAME TYRBK    2  CA 
ATOMNAME TYRBK    3  HA
ATOMNAME TYRBK    4  C  
ATOMNAME TYRBK    5  O  
ATOMNAME TYR01    0  CB 
ATOMNAME TYR01    1 1HB 
ATOMNAME TYR01    2 2HB 
ATOMNAME TYR01    3  CG 
ATOMNAME TYR01    4  CD1
ATOMNAME TYR01    5  HD1
ATOMNAME TYR01    6  CD2
ATOMNAME TYR01    7  HD2
ATOMNAME TYR01    8  CE1
ATOMNAME TYR01    9  HE1
ATOMNAME TYR01   10  CE2
ATOMNAME TYR01   11  HE2
ATOMNAME TYR01   12  CZ 
ATOMNAME TYR01   13  OH 
ATOMNAME TYR01   14  HH 
ATOMNAME TYR-1    0  CB 
ATOMNAME TYR-1    1 1HB 
ATOMNAME TYR-1    2 2HB 
ATOMNAME TYR-1    3  CG 
ATOMNAME TYR-1    4  CD1
ATOMNAME TYR-1    5  HD1
ATOMNAME TYR-1    6  CD2
ATOMNAME TYR-1    7  HD2
ATOMNAME TYR-1    8  CE1
ATOMNAME TYR-1    9  HE1
ATOMNAME TYR-1   10  CE2
ATOMNAME TYR-1   11  HE2
ATOMNAME TYR-1   12  CZ 
ATOMNAME TYR-1   13  OH 






#1.Basic Conformer Information: name, pka, em, rxn.

#Tyr has a large sidechain self vdw energy (vdw0) due to the 1-4 interactions within the ring (eg. CG - CZ)
#this interaction is ionization/position independent
#In the AMBER program, bond length energy is used to keep the ring intact

#23456789A123456789B123456789C
PROTON   TYR01      0
PKA      TYR01      0.0
ELECTRON TYR01      0
EM       TYR01      0.0
#RXN      TYR01      -1.24
RXN      TYR01      -5.88

PROTON   TYR-1      -1
PKA      TYR-1      10.2
ELECTRON TYR-1      0
EM       TYR-1      0.0
#RXN      TYR-1      -14.53
#RXN      TYR-1      -19.7
RXN      TYR-1      -41.04

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  TYRBK  N   sp2       -1    C   0     CA  0     H
CONNECT  TYRBK  H   s         0     N
CONNECT  TYRBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  TYRBK  HA  s         0     CA
CONNECT  TYRBK  C   sp2       0     CA  0     O   1     N
CONNECT  TYRBK  O   sp2       0     C

CONNECT  TYR01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  TYR01 1HB  s         0     CB
CONNECT  TYR01 2HB  s         0     CB
CONNECT  TYR01  CG  sp2       0     CB  0     CD1 0     CD2
CONNECT  TYR01  CD1 sp2       0     CG  0     CE1 0     HD1
CONNECT  TYR01  HD1 s         0     CD1
CONNECT  TYR01  CD2 sp2       0     CG  0     CE2 0     HD2
CONNECT  TYR01  HD2 s         0     CD2
CONNECT  TYR01  CE1 sp2       0     CD1 0     CZ  0     HE1
CONNECT  TYR01  HE1 s         0     CE1
CONNECT  TYR01  CE2 sp2       0     CD2 0     CZ  0     HE2
CONNECT  TYR01  HE2 s         0     CE2
CONNECT  TYR01  CZ  sp2       0     CE1 0     CE2 0     OH
CONNECT  TYR01  OH  sp3       0     CZ  0     HH
CONNECT  TYR01  HH  s         0     OH

CONNECT  TYR-1  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  TYR-1 1HB  s         0     CB
CONNECT  TYR-1 2HB  s         0     CB
CONNECT  TYR-1  CG  sp2       0     CB  0     CD1 0     CD2
CONNECT  TYR-1  CD1 sp2       0     CG  0     CE1 0     HD1
CONNECT  TYR-1  HD1 s         0     CD1
CONNECT  TYR-1  CD2 sp2       0     CG  0     CE2 0     HD2
CONNECT  TYR-1  HD2 s         0     CD2
CONNECT  TYR-1  CE1 sp2       0     CD1 0     CZ  0     HE1
CONNECT  TYR-1  HE1 s         0     CE1
CONNECT  TYR-1  CE2 sp2       0     CD2 0     CZ  0     HE2
CONNECT  TYR-1  HE2 s         0     CE2
CONNECT  TYR-1  CZ  sp2       0     CE1 0     CE2 0     OH
CONNECT  TYR-1  OH  sp3       0     CZ 
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#-------|-----|----|----|----|----|----|---------|---------|---------|----
#        CONF  ATOM ATOM ATOM ATOM      phi0(min)  n_fold   Amplitude(barrier,kcal/mol)    

DONOR    TYR01  HH   OH 
ACCEPTOR TYR-1  OH   CZ 

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   TYR    N   1.50
RADIUS   TYR    H   1.00
RADIUS   TYR    CA  2.00
RADIUS   TYR    HA  0.00
RADIUS   TYR    C   1.70
RADIUS   TYR    O   1.40
RADIUS   TYR    CB  2.00
RADIUS   TYR   1HB  0.00
RADIUS   TYR   2HB  0.00
RADIUS   TYR    CG  1.70
RADIUS   TYR    CD1 1.70
RADIUS   TYR    HD1 1.00
RADIUS   TYR    CD2 1.70
RADIUS   TYR    HD2 1.00
RADIUS   TYR    CE1 1.70
RADIUS   TYR    HE1 1.00
RADIUS   TYR    CE2 1.70
RADIUS   TYR    HE2 1.00
RADIUS   TYR    CZ  1.70
RADIUS   TYR    OH  1.40
RADIUS   TYR    HH  1.00

CHARGE   TYRBK  N    -0.350
CHARGE   TYRBK  H     0.250
CHARGE   TYRBK  CA    0.100
CHARGE   TYRBK  C     0.550
CHARGE   TYRBK  O    -0.550
CHARGE   TYR01  CB    0.125
CHARGE   TYR01  CG   -0.125
CHARGE   TYR01  CD1  -0.125
CHARGE   TYR01  HD1   0.125
CHARGE   TYR01  CE1  -0.125
CHARGE   TYR01  HE1   0.125
CHARGE   TYR01  CZ    0.055
CHARGE   TYR01  OH   -0.490
CHARGE   TYR01  HH    0.435
CHARGE   TYR01  CE2  -0.125
CHARGE   TYR01  HE2   0.125
CHARGE   TYR01  CD2  -0.125
CHARGE   TYR01  HD2   0.125
CHARGE   TYR-1  CB    0.125
CHARGE   TYR-1  CG   -0.195
CHARGE   TYR-1  CD1  -0.195
CHARGE   TYR-1  HD1   0.125
CHARGE   TYR-1  CE1  -0.195
CHARGE   TYR-1  HE1   0.125
CHARGE   TYR-1  CZ   -0.070
CHARGE   TYR-1  OH   -0.580
CHARGE   TYR-1  CE2  -0.195
CHARGE   TYR-1  HE2   0.125
CHARGE   TYR-1  CD2  -0.195
CHARGE   TYR-1  HD2   0.125

#4.Rotomer
#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  TYR   0     CA - CB   CG   CD1  CD2  CE1  CE2  CZ   OH
ROTAMER  TYR   1     CB - CG   CD1  CD2  CE1  CE2  CZ   OH
#=========================================================================

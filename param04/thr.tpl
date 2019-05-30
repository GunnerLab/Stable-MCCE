CONFLIST THR        THRBK THR01

NATOM    THRBK      6
NATOM    THR01      8

IATOM    THRBK  N   0
IATOM    THRBK  H   1
IATOM    THRBK  CA  2
IATOM    THRBK  HA  3
IATOM    THRBK  C   4
IATOM    THRBK  O   5
IATOM    THR01  CB  0
IATOM    THR01  HB  1
IATOM    THR01  OG1 2
IATOM    THR01  HG1 3
IATOM    THR01  CG2 4
IATOM    THR01 1HG2 5
IATOM    THR01 2HG2 6
IATOM    THR01 3HG2 7

ATOMNAME THRBK    0  N  
ATOMNAME THRBK    1  H  
ATOMNAME THRBK    2  CA 
ATOMNAME THRBK    3  HA 
ATOMNAME THRBK    4  C  
ATOMNAME THRBK    5  O  
ATOMNAME THR01    0  CB 
ATOMNAME THR01    1  HB 
ATOMNAME THR01    2  OG1
ATOMNAME THR01    3  HG1
ATOMNAME THR01    4  CG2
ATOMNAME THR01    5 1HG2
ATOMNAME THR01    6 2HG2
ATOMNAME THR01    7 3HG2






#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   THR01      0
PKA      THR01      0.0
ELECTRON THR01      0
EM       THR01      0.0
RXN      THR01      -2.0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  THRBK  N   sp2       -1    C   0     CA  0     H
CONNECT  THRBK  H   s         0     N
CONNECT  THRBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  THRBK  HA  s         0     CA
CONNECT  THRBK  C   sp2       0     CA  0     O   1     N
CONNECT  THRBK  O   sp2       0     C
CONNECT  THR01  CB  sp3       0     CA  0     OG1 0     CG2 0     HB
CONNECT  THR01  HB  s         0     CB
CONNECT  THR01  OG1 sp3       0     CB  0     HG1
CONNECT  THR01  HG1 s         0     OG1
CONNECT  THR01  CG2 sp3       0     CB  0    1HG2 0    2HG2 0    3HG2
CONNECT  THR01 1HG2 s         0     CG2
CONNECT  THR01 2HG2 s         0     CG2
CONNECT  THR01 3HG2 s         0     CG2
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#-------|-----|----|----|----|----|----|---------|---------|---------|----
#        CONF  ATOM ATOM ATOM ATOM      phi0(min)  n_fold   Amplitude(barrier,kcal/mol)


DONOR    THR01  HG1  OG1

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   THR    N   1.50
RADIUS   THR    H   1.00
RADIUS   THR    CA  2.00
RADIUS   THR    HA  0.00
RADIUS   THR    C   1.70
RADIUS   THR    O   1.40
RADIUS   THR    CB  2.00
RADIUS   THR    HB  0.00
RADIUS   THR    OG1 1.40
RADIUS   THR    HG1 1.00
RADIUS   THR    CG2 2.00
RADIUS   THR   1HG2 0.00
RADIUS   THR   2HG2 0.00
RADIUS   THR   3HG2 0.00

CHARGE   THRBK  N    -0.350
CHARGE   THRBK  H     0.250
CHARGE   THRBK  CA    0.100
CHARGE   THRBK  C     0.550
CHARGE   THRBK  O    -0.550
CHARGE   THR01  CB    0.000
CHARGE   THR01  OG1  -0.490
CHARGE   THR01  HG1   0.490
CHARGE   THR01  CG2   0.000

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  THR   0     CA - CB   OG1  CG2
#-------|---|----|-|---------|---------|---------|---------|---------|---------|
#ROT_SWAP THR   0     OG1- CG2
#=========================================================================

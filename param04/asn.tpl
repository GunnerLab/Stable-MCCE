### Asparagine

CONFLIST ASN        ASNBK ASN01 

NATOM    ASNBK      6
NATOM    ASN01      8

IATOM    ASNBK  N   0
IATOM    ASNBK  H   1
IATOM    ASNBK  CA  2
IATOM    ASNBK  HA  3
IATOM    ASNBK  C   4
IATOM    ASNBK  O   5
IATOM    ASN01  CB  0
IATOM    ASN01 1HB  1
IATOM    ASN01 2HB  2
IATOM    ASN01  CG  3
IATOM    ASN01  OD1 4
IATOM    ASN01  ND2 5
IATOM    ASN01 1HD2 6
IATOM    ASN01 2HD2 7

ATOMNAME ASNBK    0  N  
ATOMNAME ASNBK    1  H  
ATOMNAME ASNBK    2  CA 
ATOMNAME ASNBK    3  HA 
ATOMNAME ASNBK    4  C  
ATOMNAME ASNBK    5  O  
ATOMNAME ASN01    0  CB 
ATOMNAME ASN01    1 1HB 
ATOMNAME ASN01    2 2HB 
ATOMNAME ASN01    3  CG 
ATOMNAME ASN01    4  OD1
ATOMNAME ASN01    5  ND2
ATOMNAME ASN01    6 1HD2
ATOMNAME ASN01    7 2HD2






#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   ASN01      0
PKA      ASN01      0.0
ELECTRON ASN01      0
EM       ASN01      0.0
RXN      ASN01      -3.50

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  ASNBK  N   sp2       -1    C   0     CA  0     H
CONNECT  ASNBK  H   s         0     N
CONNECT  ASNBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  ASNBK  HA  s         0     CA
CONNECT  ASNBK  C   sp2       0     CA  0     O   1     N
CONNECT  ASNBK  O   sp2       0     C
CONNECT  ASN01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  ASN01 1HB  s         0     CB
CONNECT  ASN01 2HB  s         0     CB
CONNECT  ASN01  CG  sp2       0     CB  0     OD1 0     ND2
CONNECT  ASN01  OD1 sp2       0     CG
CONNECT  ASN01  ND2 sp2       0     CG  0    1HD2 0    2HD2
CONNECT  ASN01 1HD2 s         0     ND2
CONNECT  ASN01 2HD2 s         0     ND2
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
ACCEPTOR ASN01  OD1  CG 
DONOR    ASN01 1HD2  ND2
DONOR    ASN01 2HD2  ND2

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   ASN    N   1.50
RADIUS   ASN    H   1.00
RADIUS   ASN    CA  2.00
RADIUS   ASN    HA  0.00
RADIUS   ASN    C   1.70
RADIUS   ASN    O   1.40
RADIUS   ASN    CB  2.00
RADIUS   ASN   1HB  0.00
RADIUS   ASN   2HB  0.00
RADIUS   ASN    CG  1.70
RADIUS   ASN    OD1 1.40
RADIUS   ASN    ND2 1.50
RADIUS   ASN   1HD2 1.00
RADIUS   ASN   2HD2 1.00

CHARGE   ASNBK  N    -0.350
CHARGE   ASNBK  H     0.250
CHARGE   ASNBK  CA    0.100
CHARGE   ASNBK  C     0.550
CHARGE   ASNBK  O    -0.550
CHARGE   ASN01  CG    0.550
CHARGE   ASN01  OD1  -0.550
CHARGE   ASN01  ND2  -0.780
CHARGE   ASN01 1HD2   0.390
CHARGE   ASN01 2HD2   0.390

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  ASN   0     CA - CB   CG   OD1  ND2
ROTAMER  ASN   1     CB - CG   OD1  ND2
#=========================================================================

#-------|---|----|-|---------|---------|---------|---------|---------|---------|
ROT_SWAP ASN   0     OD1- ND2


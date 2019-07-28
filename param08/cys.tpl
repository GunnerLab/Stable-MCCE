CONFLIST CYS        CYSBK CYS01 CYS-1

NATOM    CYSBK      6
NATOM    CYS01      5
NATOM    CYS-1      4

IATOM    CYSBK  N   0
IATOM    CYSBK  H   1
IATOM    CYSBK  CA  2
IATOM    CYSBK  HA  3
IATOM    CYSBK  C   4
IATOM    CYSBK  O   5
IATOM    CYS01  CB  0
IATOM    CYS01 1HB  1
IATOM    CYS01 2HB  2
IATOM    CYS01  SG  3
IATOM    CYS01  HG  4
IATOM    CYS-1  CB  0
IATOM    CYS-1 1HB  1
IATOM    CYS-1 2HB  2
IATOM    CYS-1  SG  3

ATOMNAME CYSBK    0  N
ATOMNAME CYSBK    1  H
ATOMNAME CYSBK    2  CA
ATOMNAME CYSBK    3  HA
ATOMNAME CYSBK    4  C
ATOMNAME CYSBK    5  O
ATOMNAME CYS01    0  CB
ATOMNAME CYS01    1 1HB
ATOMNAME CYS01    2 2HB
ATOMNAME CYS01    3  SG
ATOMNAME CYS01    4  HG
ATOMNAME CYS-1    0  CB
ATOMNAME CYS-1    1 1HB
ATOMNAME CYS-1    2 2HB
ATOMNAME CYS-1    3  SG

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   CYS01      0
PROTON   CYS-1      -1
PKA      CYS01      0.0
PKA      CYS-1      9.1
ELECTRON CYS01      0
ELECTRON CYS-1      0
EM       CYS01      0.0
EM       CYS-1      0.0
RXN      CYS01      -0.46 # +-0.01
RXN      CYS-1      -8.75 # +-0.05

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  CYSBK  N   sp2       -1    C   0     CA  0     H
CONNECT  CYSBK  H   s         0     N
CONNECT  CYSBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  CYSBK  HA  s         0     CA
CONNECT  CYSBK  C   sp2       0     CA  0     O   1     N
CONNECT  CYSBK  O   sp2       0     C
CONNECT  CYS01  CB  sp3       0     CA  0     SG  0    1HB  0    2HB
CONNECT  CYS01 1HB  s         0     CB
CONNECT  CYS01 2HB  s         0     CB
CONNECT  CYS01  SG  sp2       0     CB  0     HG
CONNECT  CYS01  HG  s         0     SG
CONNECT  CYS-1  CB  sp3       0     CA  0     SG  0    1HB  0    2HB
CONNECT  CYS-1 1HB  s         0     CB  
CONNECT  CYS-1 2HB  s         0     CB  
CONNECT  CYS-1  SG  sp2       0     CB  

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   CYS    N   1.50
RADIUS   CYS    H   1.00
RADIUS   CYS    CA  2.00
RADIUS   CYS    HA  0.00
RADIUS   CYS    C   1.70
RADIUS   CYS    O   1.40
RADIUS   CYS    CB  2.00
RADIUS   CYS   1HB  0.00
RADIUS   CYS   2HB  2.00
RADIUS   CYS    SG  1.85
RADIUS   CYS    HG  1.00

CHARGE   CYSBK  N    -0.350
CHARGE   CYSBK  H     0.250
CHARGE   CYSBK  CA    0.100
CHARGE   CYSBK  C     0.550
CHARGE   CYSBK  O    -0.550
CHARGE   CYS01  SG   -0.290
CHARGE   CYS01  HG    0.290
CHARGE   CYS-1  CB   -0.080
CHARGE   CYS-1  SG   -0.920

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  CYS   0     CA - CB   SG
#=========================================================================


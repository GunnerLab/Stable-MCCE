CONFLIST SER        SERBK SER01

NATOM    SERBK      6
NATOM    SER01      5

IATOM    SERBK  N   0
IATOM    SERBK  H   1
IATOM    SERBK  CA  2
IATOM    SERBK  HA  3
IATOM    SERBK  C   4
IATOM    SERBK  O   5
IATOM    SER01  CB  0
IATOM    SER01 1HB  1
IATOM    SER01 2HB  2
IATOM    SER01  OG  3
IATOM    SER01  HG  4

ATOMNAME SERBK    0  N  
ATOMNAME SERBK    1  H  
ATOMNAME SERBK    2  CA 
ATOMNAME SERBK    3  HA 
ATOMNAME SERBK    4  C  
ATOMNAME SERBK    5  O  
ATOMNAME SER01    0  CB 
ATOMNAME SER01    1 1HB 
ATOMNAME SER01    2 2HB 
ATOMNAME SER01    3  OG 
ATOMNAME SER01    4  HG 



#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   SER01      0
PKA      SER01      0.0
ELECTRON SER01      0
EM       SER01      0.0
RXN      SER01      -0.96 # +- 0.08

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  SERBK  N   sp2       -1    C   0     CA  0     H
CONNECT  SERBK  H   s         0     N
CONNECT  SERBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  SERBK  HA  s         0     CA
CONNECT  SERBK  C   sp2       0     CA  0     O   1     N
CONNECT  SERBK  O   sp2       0     C
CONNECT  SER01  CB  sp3       0     CA  0     OG  0    1HB  0    2HB
CONNECT  SER01 1HB  s         0     CB
CONNECT  SER01 2HB  s         0     CB
CONNECT  SER01  OG  sp3       0     CB  0     HG
CONNECT  SER01  HG  s         0     OG
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#-------|-----|----|----|----|----|----|---------|---------|---------|----
#        CONF  ATOM ATOM ATOM ATOM      phi0(min)  n_fold   Amplitude(barrier,kcal/mol)

#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
DONOR    SER01  HG   OG1

#3.Atom Parameters: Partial Charges and Radii
CHARGE   SERBK  N    -0.350
CHARGE   SERBK  H     0.250
CHARGE   SERBK  CA    0.100
CHARGE   SERBK  C     0.550
CHARGE   SERBK  O    -0.550
CHARGE   SER01  CB    0.000
CHARGE   SER01  OG   -0.490
CHARGE   SER01  HG    0.490

# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   SER    N   1.50
RADIUS   SER    H   1.00
RADIUS   SER    CA  2.00
RADIUS   SER    HA  0.00
RADIUS   SER    C   1.70
RADIUS   SER    O   1.40
RADIUS   SER    CB  2.00
RADIUS   SER   1HB  0.00
RADIUS   SER   2HB  0.00
RADIUS   SER    OG  1.40
RADIUS   SER    HG  1.00

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  SER   0     CA - CB   OG
#=========================================================================

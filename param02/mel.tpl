CONFLIST MEL        MELBK MEL01

NATOM    MELBK      6
NATOM    MEL01      11

IATOM    MELBK  N   0
IATOM    MELBK  H   1
IATOM    MELBK  CA  2
IATOM    MELBK  HA  3
IATOM    MELBK  C   4
IATOM    MELBK  O   5

IATOM    MEL01  CB  0
IATOM    MEL01 1HB  1
IATOM    MEL01 2HB  2
IATOM    MEL01  CG  3
IATOM    MEL01 1HG  4
IATOM    MEL01 2HG  5
IATOM    MEL01  SD  6
IATOM    MEL01  CE  7
IATOM    MEL01 1HE  8
IATOM    MEL01 2HE  9
IATOM    MEL01 3HE  10

ATOMNAME MELBK    0  N  
ATOMNAME MELBK    1  H  
ATOMNAME MELBK    2  CA 
ATOMNAME MELBK    3  HA 
ATOMNAME MELBK    4  C  
ATOMNAME MELBK    5  O  
ATOMNAME MEL01    0  CB 
ATOMNAME MEL01    1 1HB 
ATOMNAME MEL01    2 2HB 
ATOMNAME MEL01    3  CG 
ATOMNAME MEL01    4 1HG 
ATOMNAME MEL01    5 2HG 
ATOMNAME MEL01    6  SD 
ATOMNAME MEL01    7  CE 
ATOMNAME MEL01    8 1HE 
ATOMNAME MEL01    9 2HE 
ATOMNAME MEL01   10 3HE 


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   MEL01      0
PKA      MEL01      0.0
ELECTRON MEL01      0
EM       MEL01      0.0
RXN      MEL01      0.0

PROTON   MEL02      0
PKA      MEL02      0.0
ELECTRON MEL02      0
EM       MEL02      0.0
RXN      MEL02      0.0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  MELBK  N   sp2       -1    C   0     CA  0     H
CONNECT  MELBK  H   s         0     N
CONNECT  MELBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  MELBK  HA  s         0     CA
CONNECT  MELBK  C   sp2       0     CA  0     O   1     N
CONNECT  MELBK  O   sp2       0     C

CONNECT  MEL01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  MEL01 1HB  s         0     CB
CONNECT  MEL01 2HB  s         0     CB
CONNECT  MEL01  CG  sp3       0     CB  0     SD  0    1HG  0    2HG
CONNECT  MEL01 1HG  s         0     CG
CONNECT  MEL01 2HG  s         0     CG
CONNECT  MEL01  SD  sp3       0     CG  0     CE
CONNECT  MEL01  CE  sp3       0     SD  0    1HE  0    2HE  0    3HE
CONNECT  MEL01 1HE  s         0     CE
CONNECT  MEL01 2HE  s         0     CE
CONNECT  MEL01 3HE  s         0     CE


#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
H_DIHED  MELBK  1    H    N    CA   C   60
H_DIHED  MEL01  1   1HE   CE   SD   CG  60
H_DIHED  MEL02  1   1HE   CE   SD   CG  60

#3.Atom Parameters: Partial Charges and Radii
CHARGE   MEL01  CG   0.06
CHARGE   MEL01  SD  -0.12
CHARGE   MEL01  CE   0.06

# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   MEL    N   1.55
RADIUS   MEL    H   1.20
RADIUS   MEL    CA  1.70
RADIUS   MEL    HA  1.20
RADIUS   MEL    C   1.70
RADIUS   MEL    O   1.52
RADIUS   MEL    CB  1.70
RADIUS   MEL   1HB  1.20
RADIUS   MEL   2HB  1.20
RADIUS   MEL    CG  1.70
RADIUS   MEL   1HG  1.20
RADIUS   MEL   2HG  1.20
RADIUS   MEL    SD  1.80
RADIUS   MEL    CE  1.70
RADIUS   MEL   1HE  1.20
RADIUS   MEL   2HE  1.20
RADIUS   MEL   3HE  1.20

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
#=========================================================================

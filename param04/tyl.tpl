# Modified by Beicer C. Tapia
# 12/22/08
# Protein: Cytochrome b6f complex
# PDB: 1VF5.pdb

# NOTE - REMOVED EXTRA HYDROGEN IN THE BACKBONE TO BE LINKED IN THE Hef.tpl TO THE Fe

# NOTE - This needs to be added in name.txt as well as the the Tyrososine part
# 4. HEC, Heme f including one his ligand has FE(II) and FE(III) states
#*****HEM*C*304  *****HEF*C*304
#*N  *HIS*C  26  *****BKB******   backbone residue
#*CA *HIS*C  26  *****BKB******   backbone residue
#*C  *HIS*C  26  *****BKB******   backbone residue
#*O  *HIS*C  26  *****BKB******   backbone residue
#*****HIS*C  26  a****HEF*C*304

# NOTE - THIS PART IS ALSO NEEDED SINCE FORMS PART OF THE N-TER OF THE hef.tpl
# 5. Tyrosine C -1 - Named TYL
#*****TYR*C   1  *****TYL*C***1


CONFLIST TYL        TYLBK TYL01 TYL-1 

NATOM    TYLBK      5
NATOM    TYL01      15
NATOM    TYL-1      14

IATOM    TYLBK  N   0
IATOM    TYLBK  CA  1
IATOM    TYLBK  HA  2
IATOM    TYLBK  C   3
IATOM    TYLBK  O   4
IATOM    TYL01  CB  0
IATOM    TYL01 1HB  1
IATOM    TYL01 2HB  2
IATOM    TYL01  CG  3
IATOM    TYL01  CD1 4
IATOM    TYL01  HD1 5
IATOM    TYL01  CD2 6
IATOM    TYL01  HD2 7
IATOM    TYL01  CE1 8
IATOM    TYL01  HE1 9
IATOM    TYL01  CE2 10
IATOM    TYL01  HE2 11
IATOM    TYL01  CZ  12
IATOM    TYL01  OH  13
IATOM    TYL01  HH  14
IATOM    TYL-1  CB  0
IATOM    TYL-1 1HB  1
IATOM    TYL-1 2HB  2
IATOM    TYL-1  CG  3
IATOM    TYL-1  CD1 4
IATOM    TYL-1  HD1 5
IATOM    TYL-1  CD2 6
IATOM    TYL-1  HD2 7
IATOM    TYL-1  CE1 8
IATOM    TYL-1  HE1 9
IATOM    TYL-1  CE2 10
IATOM    TYL-1  HE2 11
IATOM    TYL-1  CZ  12
IATOM    TYL-1  OH  13

ATOMNAME TYLBK    0  N  
ATOMNAME TYLBK    1  CA 
ATOMNAME TYLBK    2  HA 
ATOMNAME TYLBK    3  C  
ATOMNAME TYLBK    4  O  
ATOMNAME TYL01    0  CB 
ATOMNAME TYL01    1 1HB 
ATOMNAME TYL01    2 2HB 
ATOMNAME TYL01    3  CG 
ATOMNAME TYL01    4  CD1
ATOMNAME TYL01    5  HD1
ATOMNAME TYL01    6  CD2
ATOMNAME TYL01    7  HD2
ATOMNAME TYL01    8  CE1
ATOMNAME TYL01    9  HE1
ATOMNAME TYL01   10  CE2
ATOMNAME TYL01   11  HE2
ATOMNAME TYL01   12  CZ 
ATOMNAME TYL01   13  OH 
ATOMNAME TYL01   14  HH 
ATOMNAME TYL-1    0  CB 
ATOMNAME TYL-1    1 1HB 
ATOMNAME TYL-1    2 2HB 
ATOMNAME TYL-1    3  CG 
ATOMNAME TYL-1    4  CD1
ATOMNAME TYL-1    5  HD1
ATOMNAME TYL-1    6  CD2
ATOMNAME TYL-1    7  HD2
ATOMNAME TYL-1    8  CE1
ATOMNAME TYL-1    9  HE1
ATOMNAME TYL-1   10  CE2
ATOMNAME TYL-1   11  HE2
ATOMNAME TYL-1   12  CZ 
ATOMNAME TYL-1   13  OH 

#1.Basic Conformer Information: name, pka, em, rxn.

#TYL has a large sidechain self vdw energy (vdw0) due to the 1-4 interactions within the ring (eg. CG - CZ)
#this interaction is ionization/position independent
#In the AMBER program, bond length energy is used to keep the ring intact

#23456789A123456789B123456789C
PROTON   TYL01      0
PKA      TYL01      0.0
ELECTRON TYL01      0
EM       TYL01      0.0
#RXN      TYL01      -1.24
RXN      TYL01      -2.7

PROTON   TYL-1      -1
PKA      TYL-1      10.2
ELECTRON TYL-1      0
EM       TYL-1      0.0
#RXN      TYL-1      -14.53
#RXN      TYL-1      -19.7
RXN      TYL-1      -19.0 # adjusted by jmao

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  TYLBK  N   sp2       -1    C   0     CA  LIG  FE
CONNECT  TYLBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  TYLBK  HA  s         0     CA
CONNECT  TYLBK  C   sp2       0     CA  0     O   1     N
CONNECT  TYLBK  O   sp2       0     C

CONNECT  TYL01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  TYL01 1HB  s         0     CB
CONNECT  TYL01 2HB  s         0     CB
CONNECT  TYL01  CG  sp2       0     CB  0     CD1 0     CD2
CONNECT  TYL01  CD1 sp2       0     CG  0     CE1 0     HD1
CONNECT  TYL01  HD1 s         0     CD1
CONNECT  TYL01  CD2 sp2       0     CG  0     CE2 0     HD2
CONNECT  TYL01  HD2 s         0     CD2
CONNECT  TYL01  CE1 sp2       0     CD1 0     CZ  0     HE1
CONNECT  TYL01  HE1 s         0     CE1
CONNECT  TYL01  CE2 sp2       0     CD2 0     CZ  0     HE2
CONNECT  TYL01  HE2 s         0     CE2
CONNECT  TYL01  CZ  sp2       0     CE1 0     CE2 0     OH
CONNECT  TYL01  OH  sp3       0     CZ  0     HH
CONNECT  TYL01  HH  s         0     OH

CONNECT  TYL-1  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  TYL-1 1HB  s         0     CB
CONNECT  TYL-1 2HB  s         0     CB
CONNECT  TYL-1  CG  sp2       0     CB  0     CD1 0     CD2
CONNECT  TYL-1  CD1 sp2       0     CG  0     CE1 0     HD1
CONNECT  TYL-1  HD1 s         0     CD1
CONNECT  TYL-1  CD2 sp2       0     CG  0     CE2 0     HD2
CONNECT  TYL-1  HD2 s         0     CD2
CONNECT  TYL-1  CE1 sp2       0     CD1 0     CZ  0     HE1
CONNECT  TYL-1  HE1 s         0     CE1
CONNECT  TYL-1  CE2 sp2       0     CD2 0     CZ  0     HE2
CONNECT  TYL-1  HE2 s         0     CE2
CONNECT  TYL-1  CZ  sp2       0     CE1 0     CE2 0     OH
CONNECT  TYL-1  OH  sp3       0     CZ 
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#-------|-----|----|----|----|----|----|---------|---------|---------|----
#        CONF  ATOM ATOM ATOM ATOM      phi0(min)  n_fold   Amplitude(barrier,kcal/mol)    

DONOR    TYL01  HH   OH 
ACCEPTOR TYL-1  OH   CZ 

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   TYL    N   1.50
RADIUS   TYL    H   1.00
RADIUS   TYL    CA  2.00
RADIUS   TYL    HA  0.00
RADIUS   TYL    C   1.70
RADIUS   TYL    O   1.40
RADIUS   TYL    CB  2.00
RADIUS   TYL   1HB  0.00
RADIUS   TYL   2HB  0.00
RADIUS   TYL    CG  1.70
RADIUS   TYL    CD1 1.70
RADIUS   TYL    HD1 1.00
RADIUS   TYL    CD2 1.70
RADIUS   TYL    HD2 1.00
RADIUS   TYL    CE1 1.70
RADIUS   TYL    HE1 1.00
RADIUS   TYL    CE2 1.70
RADIUS   TYL    HE2 1.00
RADIUS   TYL    CZ  1.70
RADIUS   TYL    OH  1.40
RADIUS   TYL    HH  1.00

CHARGE   TYLBK  N    -0.350
CHARGE   TYLBK  H     0.250
CHARGE   TYLBK  CA    0.100
CHARGE   TYLBK  C     0.550
CHARGE   TYLBK  O    -0.550
CHARGE   TYL01  CB    0.125
CHARGE   TYL01  CG   -0.125
CHARGE   TYL01  CD1  -0.125
CHARGE   TYL01  HD1   0.125
CHARGE   TYL01  CE1  -0.125
CHARGE   TYL01  HE1   0.125
CHARGE   TYL01  CZ    0.055
CHARGE   TYL01  OH   -0.490
CHARGE   TYL01  HH    0.435
CHARGE   TYL01  CE2  -0.125
CHARGE   TYL01  HE2   0.125
CHARGE   TYL01  CD2  -0.125
CHARGE   TYL01  HD2   0.125
CHARGE   TYL-1  CB    0.125
CHARGE   TYL-1  CG   -0.195
CHARGE   TYL-1  CD1  -0.195
CHARGE   TYL-1  HD1   0.125
CHARGE   TYL-1  CE1  -0.195
CHARGE   TYL-1  HE1   0.125
CHARGE   TYL-1  CZ   -0.070
CHARGE   TYL-1  OH   -0.580
CHARGE   TYL-1  CE2  -0.195
CHARGE   TYL-1  HE2   0.125
CHARGE   TYL-1  CD2  -0.195
CHARGE   TYL-1  HD2   0.125

#4.Rotomer
#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  TYL   0     CA - CB   CG   CD1  CD2  CE1  CE2  CZ   OH
ROTAMER  TYL   1     CB - CG   CD1  CD2  CE1  CE2  CZ   OH
#=========================================================================

TORSION  TYL01  HH   OH   CZ   CE1  f        0.0         3      0.00

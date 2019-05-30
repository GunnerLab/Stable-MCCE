CONFLIST OCS        OCSBK OCS01 OCS-1 

NATOM    OCSBK      6
NATOM    OCS01      8
NATOM    OCS-1      7

IATOM    OCSBK  N   0
IATOM    OCSBK  H   1
IATOM    OCSBK  CA  2
IATOM    OCSBK  HA  3
IATOM    OCSBK  C   4
IATOM    OCSBK  O   5
IATOM    OCS01  CB  0
IATOM    OCS01  SG  1
IATOM    OCS01  OD1 2
IATOM    OCS01  OD2 3
IATOM    OCS01  OD3 4
IATOM    OCS01  HB2 5
IATOM    OCS01  HB3 6
IATOM    OCS01  HD2 7
IATOM    OCS-1  CB  0
IATOM    OCS-1  SG  1
IATOM    OCS-1  OD1 2
IATOM    OCS-1  OD2 3
IATOM    OCS-1  OD3 4
IATOM    OCS-1  HB2 5
IATOM    OCS-1  HB3 6

ATOMNAME OCSBK    0  N  
ATOMNAME OCSBK    1  H  
ATOMNAME OCSBK    2  CA 
ATOMNAME OCSBK    3  HA 
ATOMNAME OCSBK    4  C  
ATOMNAME OCSBK    5  O  
ATOMNAME OCS01    0  CB 
ATOMNAME OCS01    1  SG 
ATOMNAME OCS01    2  OD1
ATOMNAME OCS01    3  OD2
ATOMNAME OCS01    4  OD3
ATOMNAME OCS01    5  HB2
ATOMNAME OCS01    6  HB3
ATOMNAME OCS01    7  HD2
ATOMNAME OCS-1    0  CB 
ATOMNAME OCS-1    1  SG 
ATOMNAME OCS-1    2  OD1
ATOMNAME OCS-1    3  OD2
ATOMNAME OCS-1    4  OD3
ATOMNAME OCS-1    5  HB2
ATOMNAME OCS-1    6  HB3

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
#pKa taken from taurine
PROTON   OCS01      0
PROTON   OCS-1     -1
PKA      OCS01      0.0
PKA      OCS-1      -1.8
ELECTRON OCS01      0
ELECTRON OCS-1      0
EM       OCS01     -3.716
EM       OCS-1     -13.365

#2.Structure Connectivity
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  OCSBK  N         sp2   -1  C      0  CA     0  H      0  HN2 
CONNECT  OCSBK  H           s    0  N   
CONNECT  OCSBK  CA        sp3    0  N      0  CB     0  C      0  HA  
CONNECT  OCSBK  HA          s    0  CA  
CONNECT  OCSBK  C         sp2    0  CA     0  O 
CONNECT  OCSBK  O         sp2    0  C   
CONNECT  OCS01  CB        sp3    0  CA     0  SG     0  HB2    0  HB3 
CONNECT  OCS01  SG        sp3    0  CB     0  OD1    0  OD2    0  OD3 
CONNECT  OCS01  OD1       sp2    0  SG  
CONNECT  OCS01  OD2       sp3    0  SG     0  HD2 
CONNECT  OCS01  OD3       sp2    0  SG  
CONNECT  OCS01  HB2         s    0  CB  
CONNECT  OCS01  HB3         s    0  CB 
CONNECT  OCS01  HD2         s    0  OD2 
CONNECT  OCS-1  CB        sp3    0  CA     0  SG     0  HB2    0  HB3 
CONNECT  OCS-1  SG        sp3    0  CB     0  OD1    0  OD2    0  OD3 
CONNECT  OCS-1  OD1       sp2    0  SG  
CONNECT  OCS-1  OD2       sp2    0  SG     0  HD2 
CONNECT  OCS-1  OD3       sp2    0  SG   
CONNECT  OCS-1  HB2         s    0  CB  
CONNECT  OCS-1  HB3         s    0  CB 

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   OCS    N   1.50
RADIUS   OCS    H   1.00
RADIUS   OCS    CA  2.00
RADIUS   OCS    HA  0.00
RADIUS   OCS    C   1.70
RADIUS   OCS    O   1.40
RADIUS   OCS    CB  2.00
RADIUS   OCS    SG  1.85
RADIUS   OCS    OD1 1.40
RADIUS   OCS    OD2 1.40
RADIUS   OCS    OD3 1.40
RADIUS   OCS    HB2 0.00
RADIUS   OCS    HB3 0.00
RADIUS   OCS    HD2 1.00

#charges taken from Gaussian
CHARGE   OCSBK  N    -0.329442
CHARGE   OCSBK  H     0.214649
CHARGE   OCSBK  CA   -0.173269
CHARGE   OCSBK  HA    0.201465
CHARGE   OCSBK  C     0.199299
CHARGE   OCSBK  O    -0.152021
CHARGE   OCS-1  CB   -0.385504
CHARGE   OCS-1  SG    1.016656
CHARGE   OCS-1  OD1  -0.598772
CHARGE   OCS-1  OD2  -0.520687
CHARGE   OCS-1  OD3  -0.592502
CHARGE   OCS-1  HB2   0.134582
CHARGE   OCS-1  HB3   0.114829
CHARGE   OCS01  CB   -0.371138
CHARGE   OCS01  SG    1.025132
CHARGE   OCS01  OD1  -0.495268
CHARGE   OCS01  OD2  -0.284510
CHARGE   OCS01  OD3  -0.516042
CHARGE   OCS01  HB2   0.201809
CHARGE   OCS01  HB3   0.172128
CHARGE   OCS01  HD2   0.307209

#-------|-----|----|----|----|----|----|---------|---------|---------|---------|---------|---------|---------|---------|---------|
#        CONF  ATOM ATOM ATOM ATOM relx Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma
TORSION  OCS    CA   CB   SG   OD1  f        0.0         1      0.00
TORSION  OCS    CA   CB   SG   OD2  f        0.0         1      0.00
TORSION  OCS    CA   CB   SG   OD3  f        0.0         1      0.00
TORSION  OCS01  CB   SG  OD2   HD2  f        0.0         1      0.00

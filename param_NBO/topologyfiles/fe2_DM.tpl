CONFLIST FE2        FE2BK FE2+2 FE2+3  

NATOM    FE2BK      0
NATOM    FE2+2      1
NATOM    FE2+3      1

IATOM    FE2+2 FE   0
IATOM    FE2+3 FE   0

ATOMNAME FE2+3    0 FE  
ATOMNAME FE2+2    0 FE


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   FE2+2      0
PKA      FE2+2      0.0
ELECTRON FE2+2      0
EM       FE2+2      0.0
RXN      FE2+2      -78.85

PROTON   FE2+3      0
PKA      FE2+3      0.0
ELECTRON FE2+3      -1
EM       FE2+3      770
RXN      FE2+3      -177.4

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  FE2+2 FE   ion
CONNECT  FE2+3 FE   ion

DONOR    FEX+2 1     FE

#3.Atom Parameters: Partial Charges and Radii  THIS NEEDS TO BE CHECKED
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   FE2   FE   2.00

CHARGE   FE2+2 FE   2.00
CHARGE   FE2+3 FE   3.00

#TRANS    FE2         t

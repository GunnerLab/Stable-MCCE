CONFLIST FE2        FE2BK FE201 FE2+2  

NATOM    FE2BK      0
NATOM    FE201      1
NATOM    FE2+2      1

IATOM    FE201 FE   0
IATOM    FE2+2 FE   0

ATOMNAME FE201    0 FE  
ATOMNAME FE2+2    0 FE


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   FEX01      0
PKA      FEX01      0.0
ELECTRON FEX01      0
EM       FEX01      0.0
RXN      FEX01      0

PROTON   FEX+2      0
PKA      FEX+2      0.0
ELECTRON FEX+2      2
EM       FEX+2      0.0
RXN      FEX+2      0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  FE201 FE   ion
CONNECT  FE2+2 FE   ion

DONOR    FEX+2 1     FE

#3.Atom Parameters: Partial Charges and Radii  THIS NEEDS TO BE CHECKED
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   FE2   FE2  2.23

CHARGE   FE201 FE   0.00
CHARGE   FE2+2 FE   2.00


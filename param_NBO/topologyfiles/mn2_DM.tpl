CONFLIST MN2        MN2BK MN2+2 MN2+3  

NATOM    MN2BK      0
NATOM    MN2+2      1
NATOM    MN2+3      1

IATOM    MN2+2 MN   0
IATOM    MN2+3 MN   0

ATOMNAME MN2+3    0 MN  
ATOMNAME MN2+2    0 MN


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   MN2+2      0
PKA      MN2+2      0.0
ELECTRON MN2+2      0
EM       MN2+2      0.0
RXN      MN2+2      -68.635
#RXN      MN2+2      -71.755
#RXN      MN2+2      -78.85

PROTON   MN2+3      0
PKA      MN2+3      0.0
ELECTRON MN2+3      -1
EM       MN2+3      1500
RXN      MN2+3      -154.428
#RXN      MN2+3      -161.448
#RXN      MN2+3      -177.4

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  MN2+2 MN   ion
CONNECT  MN2+3 MN   ion

#3.Atom Parameters: Partial Charges and Radii  THIS NEEDS TO BE CHECKED
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   MN2   MN   2.30

CHARGE   MN2+2 MN   2.00
CHARGE   MN2+3 MN   3.00

VDW_RAD  MN2+2 MN   0
VDW_EPS  MN2+2 MN   0
VDW_RAD  MN2+3 MN   0
VDW_EPS  MN2+3 MN   0


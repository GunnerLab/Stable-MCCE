CONFLIST CA2        CA2BK CA2+2 CA2DM

NATOM    CA2BK      0
NATOM    CA2+2      1
NATOM    CA2DM      0

IATOM    CA2+2 CA2  0

ATOMNAME CA2+2   0  CA2 


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   CA2+2      0
PKA      CA2+2      0.0
ELECTRON CA2+2      0
EM       CA2+2      0.0
RXN      CA2+2      -70.79

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  CA2+2 CA2  ion

#3.Atom Parameters: Partial Charges and Radii
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   CA2   CA2  2.23

CHARGE   CA2+2 CA2  2.00
#-------|-----|----|--------------
#ParaNam|Res |Atom|Param/toggle
TRANS    CA2         t


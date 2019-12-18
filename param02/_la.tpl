CONFLIST _LA        _LABK _LA+3 _LADM

NATOM    _LABK      0
NATOM    _LA+3      1
NATOM    _LADM      0

IATOM    _LA+3 LA   0

ATOMNAME _LA+3    0 LA  


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   _LA+3      0
PKA      _LA+3      0.0
ELECTRON _LA+3      0
EM       _LA+3      0.0
RXN      _LA+3      -182.146

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  _LA+3 LA   ion

#3.Atom Parameters: Partial Charges and Radii
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   _LA   LA   1.95

CHARGE   _LA+3 LA   3.00


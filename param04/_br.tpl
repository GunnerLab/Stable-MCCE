CONFLIST _BR        _BRBK _BR-1 _BRDM

NATOM    _BRBK      0
NATOM    _BR-1      1
NATOM    _BRDM      0

IATOM    _BR-1 BR   0

ATOMNAME _BR-1      0 BR  


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   _BR-1      0
PKA      _BR-1      0.0
ELECTRON _BR-1      0
EM       _BR-1      0.0
RXN      _BR-1      0.0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  _BR-1 BR   ion

#3.Atom Parameters: Partial Charges and Radii
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   _BR   BR   1.96

CHARGE   _BR-1 BR   -1.00

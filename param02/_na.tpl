#11/27/06 Gail Schneider
CONFLIST _NA        _NABK _NA+2 _NADM

NATOM    _NABK      0
NATOM    _NA+2      1
NATOM    _NADM      0

IATOM    _NA+2 NA   0

ATOMNAME _NA+2    0 NA  


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   _NA+2      0
PKA      _NA+2      0.0
ELECTRON _NA+2      0
EM       _NA+2      0.0
RXN      _NA+2      -70.79

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  _NA+2 NA   ion

#3.Atom Parameters: Partial Charges and Radii
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   _NA   NA   2.27

CHARGE   _NA+2 NA   2.00

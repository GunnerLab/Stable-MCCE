CONFLIST _CU        _CUBK _CU+1 _CU+2 _CUDM

NATOM    _CUBK      0
NATOM    _CU+1      1
NATOM    _CU+2      1
NATOM    _CUDM      0

IATOM    _CU+1 CU   0
IATOM    _CU+2 CU   0

ATOMNAME _CU+1    0 CU  
ATOMNAME _CU+2    0 CU  

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   _CU+2      0
PKA      _CU+2      0
ELECTRON _CU+2      -1
EM       _CU+2      0
RXN      _CU+2      -51.54  # with radii 1.45, edited by Jun  # 1/2 of existing value

PROTON   _CU+1      0
PKA      _CU+1      0
ELECTRON _CU+1      0
EM       _CU+1      0
RXN      _CU+1      -12.89  # with radii 1.45, edited by Jun # 1/2 of existing value

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  _CU+1 CU   ion
CONNECT  _CU+2 CU   ion

#3.Atom Parameters: Partial Charges and Radii
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   _CU   CU   1.45

CHARGE   _CU+1 CU   1.00
CHARGE   _CU+2 CU   2.00

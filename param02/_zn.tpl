# Copied and modified from _mg.tpl file. -Yifan 11/13/06
# dummy atom added 6/5/07 by marilyn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
#ONFLIST _CL        _CLBK _CL-1 _CLDM
CONFLIST _ZN        _ZNBK _ZN+2 _ZNDM

NATOM    _ZNBK      0
NATOM    _ZN+2      1
NATOM    _ZNDM      0

IATOM    _ZN+2 ZN   0

ATOMNAME _ZN+2    0 ZN  

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   _ZN+2      0
PKA      _ZN+2      0.0
ELECTRON _ZN+2      0
EM       _ZN+2      0.0
RXN      _ZN+2      -91.248


#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  _ZN+2 ZN   ion 

#3.Atom Parameters: Partial Charges and Radii
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   _ZN   ZN   1.73

CHARGE   _ZN+2 ZN   2.00

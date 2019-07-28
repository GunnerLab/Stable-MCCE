CONFLIST CTR        CTRBK CTR01 CTR02 CTR-1

NATOM    CTRBK      0
NATOM    CTR01      4
NATOM    CTR02      4
NATOM    CTR-1      3

IATOM    CTR01  C   0
IATOM    CTR01  O   1
IATOM    CTR01  HO  2
IATOM    CTR01  OXT 3
IATOM    CTR02  C   0
IATOM    CTR02  O   1
IATOM    CTR02  OXT 2
IATOM    CTR02  HXT 3
IATOM    CTR-1  C   0
IATOM    CTR-1  O   1
IATOM    CTR-1  OXT 2

ATOMNAME CTR01    0  C  
ATOMNAME CTR01    1  O  
ATOMNAME CTR01    2  HO 
ATOMNAME CTR01    3  OXT
ATOMNAME CTR02    0  C  
ATOMNAME CTR02    1  O  
ATOMNAME CTR02    2  OXT
ATOMNAME CTR02    3  HXT
ATOMNAME CTR-1    0  C  
ATOMNAME CTR-1    1  O  
ATOMNAME CTR-1    2  OXT







#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
#RXN      CTR01       -1.54
#RXN      CTR02       -1.54
#RXN      CTR-1       -17.72
RXN      CTR01       -1.11 # +-0.14
RXN      CTR02       -1.11 # +-0.14
RXN      CTR-1       -9.95
PKA      CTR01        0.0
PKA      CTR02        0.0
PKA      CTR-1        3.75
PROTON   CTR01        0
PROTON   CTR02        0
PROTON   CTR-1       -1
ELECTRON CTR01        0
ELECTRON CTR02        0
ELECTRON CTR-1        0
EM       CTR01        0
EM       CTR02        0
EM       CTR-1        0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H1$
CONNECT  CTR01  C   sp2       -1    CA  0     O   0     OXT
CONNECT  CTR01  O   sp3       0     C   0     HO 
CONNECT  CTR01  HO  s         0     O  
CONNECT  CTR01  OXT sp2       0     C  
CONNECT  CTR02  C   sp2       -1    CA  0     O   0     OXT
CONNECT  CTR02  O   sp2       0     C 
CONNECT  CTR02  OXT sp3       0     C   0     HXT
CONNECT  CTR02  HXT s         0     OXT
CONNECT  CTR-1  C   sp2       -1    CA  0     O   0     OXT
CONNECT  CTR-1  O   sp2       0     C 
CONNECT  CTR-1  OXT sp2       0     C 

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   CTR    C   1.70
RADIUS   CTR    O   1.40
RADIUS   CTR    HO  1.00
RADIUS   CTR    OXT 1.40
RADIUS   CTR    HXT 1.00

CHARGE   CTR01  C     0.000
CHARGE   CTR01  OXT  -0.200
CHARGE   CTR01  O    -0.200
CHARGE   CTR01  HO    0.400
CHARGE   CTR02  C     0.000
CHARGE   CTR02  OXT  -0.200
CHARGE   CTR02  O    -0.200
CHARGE   CTR02  HXT   0.400
CHARGE   CTR-1  C     0.140
CHARGE   CTR-1  OXT  -0.570
CHARGE   CTR-1  O    -0.570


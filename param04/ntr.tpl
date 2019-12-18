CONFLIST NTR        NTRBK NTR01 NTR+1

NATOM    NTRBK      0
NATOM    NTR01      5
NATOM    NTR+1      6

IATOM    NTR01  CA  0
IATOM    NTR01  HA  1
IATOM    NTR01  N   2
IATOM    NTR01 1H   3
IATOM    NTR01 2H   4
IATOM    NTR+1  CA  0
IATOM    NTR+1  HA  1
IATOM    NTR+1  N   2
IATOM    NTR+1 1H   3
IATOM    NTR+1 2H   4
IATOM    NTR+1 3H   5

ATOMNAME NTR01    0  CA
ATOMNAME NTR01    1  HA
ATOMNAME NTR01    2  N
ATOMNAME NTR01    3 1H
ATOMNAME NTR01    4 2H
ATOMNAME NTR+1    0  CA
ATOMNAME NTR+1    1  HA
ATOMNAME NTR+1    2  N
ATOMNAME NTR+1    3 1H
ATOMNAME NTR+1    4 2H
ATOMNAME NTR+1    5 3H


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   NTR01      0
PKA      NTR01      0.0
ELECTRON NTR01      0
EM       NTR01      0.0
#RXN      NTR01      -1.77
RXN      NTR01      0.0

PROTON   NTR+1      1
PKA      NTR+1      8.0
ELECTRON NTR+1      0
EM       NTR+1      0.0
#RXN      NTR+1      -23.4
RXN      NTR+1      -22.4 # adjusted by jmao


#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H1$
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  NTR01  CA  sp3       0     N   +1    C   +1    CB  0     HA
CONNECT  NTR01  HA  s         0     CA
CONNECT  NTR01  N   sp2       0     CA  0    1H   0    2H
CONNECT  NTR01 1H   s         0     N
CONNECT  NTR01 2H   s         0     N

CONNECT  NTR+1  CA  sp3       0     N   +1    C   +1    CB  0     HA
CONNECT  NTR+1  HA  s         0     CA
CONNECT  NTR+1  N   sp3       0     CA  0    1H   0    2H   0    3H
CONNECT  NTR+1 1H   s         0     N 
CONNECT  NTR+1 2H   s         0     N 
CONNECT  NTR+1 3H   s         0     N 

#3.Atom Parameters: Partial Charges and Radii
CHARGE   NTR01  N    -0.003
CHARGE   NTR01  CA    0.001
CHARGE   NTR01 1H     0.001
CHARGE   NTR01 2H     0.001

CHARGE   NTR+1  N    -0.100
CHARGE   NTR+1  CA    0.050
CHARGE   NTR+1 1H     0.350
CHARGE   NTR+1 2H     0.350
CHARGE   NTR+1 3H     0.350

# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   NTR    CA  2.00
RADIUS   NTR    HA  0.00
RADIUS   NTR    N   1.50
RADIUS   NTR   1H   1.00
RADIUS   NTR   2H   1.00
RADIUS   NTR   3H   1.00


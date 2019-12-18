CONFLIST NTG        NTGBK NTG01 NTG+1

NATOM    NTGBK      0
NATOM    NTG01      6
NATOM    NTG+1      7

IATOM    NTG01  CA  0
IATOM    NTG01 1HA  1
IATOM    NTG01 2HA  2
IATOM    NTG01  N   3
IATOM    NTG01 1H   4
IATOM    NTG01 2H   5
IATOM    NTG+1  CA  0
IATOM    NTG+1 1HA  1
IATOM    NTG+1 2HA  2
IATOM    NTG+1  N   3
IATOM    NTG+1 1H   4
IATOM    NTG+1 2H   5
IATOM    NTG+1 3H   6

ATOMNAME NTG01    0  CA
ATOMNAME NTG01    1 1HA
ATOMNAME NTG01    2 2HA
ATOMNAME NTG01    3  N
ATOMNAME NTG01    4 1H
ATOMNAME NTG01    5 2H
ATOMNAME NTG+1    0  CA
ATOMNAME NTG+1    1 1HA
ATOMNAME NTG+1    2 2HA
ATOMNAME NTG+1    3  N
ATOMNAME NTG+1    4 1H
ATOMNAME NTG+1    5 2H
ATOMNAME NTG+1    6 3H





# This file is for NTR of GLY, in which CA has 1HA and 2HA

#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   NTG01      0
PKA      NTG01      0.0
ELECTRON NTG01      0
EM       NTG01      0.0
RXN      NTG01      0.0

PROTON   NTG+1      1
PKA      NTG+1      8.0
ELECTRON NTG+1      0
EM       NTG+1      0.0
#RXN      NTG+1      -23.4 
RXN      NTG+1      -22.4  # adjusted by jmao

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H1$
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  NTG01  CA  sp3       0     N   +1    C   0    1HA  0    2HA
CONNECT  NTG01 1HA  s         0     CA
CONNECT  NTG01 2HA  s         0     CA
CONNECT  NTG01  N   sp2       0     CA  0    1H   0    2H
CONNECT  NTG01 1H   s         0     N
CONNECT  NTG01 2H   s         0     N

CONNECT  NTG+1  CA  sp3       0     N   +1    C   0    1HA  0    2HA
CONNECT  NTG+1 1HA  s         0     CA
CONNECT  NTG+1 2HA  s         0     CA
CONNECT  NTG+1  N   sp3       0     CA  0    1H   0    2H   0    3H
CONNECT  NTG+1 1H   s         0     N
CONNECT  NTG+1 2H   s         0     N
CONNECT  NTG+1 3H   s         0     N

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   NTG    CA  2.00
RADIUS   NTG   1HA  0.00
RADIUS   NTG   2HA  0.00
RADIUS   NTG    N   1.50
RADIUS   NTG   1H   1.00
RADIUS   NTG   2H   1.00
RADIUS   NTG   3H   1.00

CHARGE   NTG01  N    -0.003
CHARGE   NTG01  CA    0.001
CHARGE   NTG01 1H     0.001
CHARGE   NTG01 2H     0.001
CHARGE   NTG+1  N    -0.100
CHARGE   NTG+1  CA    0.050
CHARGE   NTG+1 1H     0.350
CHARGE   NTG+1 2H     0.350
CHARGE   NTG+1 3H     0.350



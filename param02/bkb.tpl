#Backbone part of HIS, used when HIS side chain is a ligand
CONFLIST BKB        BKBBK
NATOM    BKBBK      6
IATOM    BKBBK  N   0
IATOM    BKBBK  H   1
IATOM    BKBBK  CA  2
IATOM    BKBBK  HA  3
IATOM    BKBBK  C   4
IATOM    BKBBK  O   5
ATOMNAME BKBBK    0  N
ATOMNAME BKBBK    1  H
ATOMNAME BKBBK    2  CA
ATOMNAME BKBBK    3  HA
ATOMNAME BKBBK    4  C
ATOMNAME BKBBK    5  O

#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  BKBBK  N   sp2       -1    C   0     CA  0     H
CONNECT  BKBBK  H   s         0     N
CONNECT  BKBBK  CA  sp3       0     N   0     C   LIG   ?   0     HA
CONNECT  BKBBK  HA  s         0     CA
CONNECT  BKBBK  C   sp2       0     CA  0     O   +1    N
CONNECT  BKBBK  O   sp2       0     C

CHARGE   BKBBK  N    -0.350
CHARGE   BKBBK  H     0.250
CHARGE   BKBBK  CA    0.100
CHARGE   BKBBK  C     0.550
CHARGE   BKBBK  O    -0.550

#RADIUS   BKB    CA  2.00
#RADIUS   BKB    C   1.70
#RADIUS   BKB    O   1.40
RADIUS   BKB    CA    1.80            #repair mg 7/28
RADIUS   BKB    C     1.80            #repair mg 7/28
RADIUS   BKB    O     1.52            #repair mg 7/28
RADIUS   BKB    N   1.50
RADIUS   BKB    H   1.00
RADIUS   BKB    HA  0.00

#Backbone part of HIS, used when HIS side chain is a ligand
CONFLIST HIB        HIBBK
NATOM    HIBBK      6
IATOM    HIBBK  N   0
IATOM    HIBBK  H   1
IATOM    HIBBK  CA  2
IATOM    HIBBK  HA  3
IATOM    HIBBK  C   4
IATOM    HIBBK  O   5
ATOMNAME HIBBK    0  N  
ATOMNAME HIBBK    1  H  
ATOMNAME HIBBK    2  CA 
ATOMNAME HIBBK    3  HA 
ATOMNAME HIBBK    4  C  
ATOMNAME HIBBK    5  O  

CONNECT  HIBBK  N   sp2       -1    C   0     CA  0     H
CONNECT  HIBBK  H   s         0     N
CONNECT  HIBBK  CA  sp3       0     N   0     C   LIG   ?   0     HA
CONNECT  HIBBK  HA  s         0     CA
CONNECT  HIBBK  C   sp2       0     CA  0     O   1     N
CONNECT  HIBBK  O   sp2       0     C

CHARGE   HIBBK  N    -0.350
CHARGE   HIBBK  H     0.250
CHARGE   HIBBK  CA    0.100
CHARGE   HIBBK  C     0.550
CHARGE   HIBBK  O    -0.550

RADIUS   HIB    N   1.50
RADIUS   HIB    H   1.00
RADIUS   HIB    CA  2.00
RADIUS   HIB    HA  0.00
RADIUS   HIB    C   1.70
RADIUS   HIB    O   1.40

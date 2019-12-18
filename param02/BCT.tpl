CONFLIST BCT        BCTBK BCT01 BCT-1 

NATOM    BCTBK      2
NATOM    BCT01      2
NATOM    BCT-1      2

IATOM    BCTBK  C   0
IATOM    BCTBK  O1  1
IATOM    BCT01  O2  0
IATOM    BCT01  O3  1
IATOM    BCT-1  O2  0
IATOM    BCT-1  O3  1

ATOMNAME BCTBK    0  C  
ATOMNAME BCTBK    1  O1 
ATOMNAME BCT01    0  O2 
ATOMNAME BCT01    1  O3 
ATOMNAME BCT-1    0  O2 
ATOMNAME BCT-1    1  O3 





#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  BCTBK  C   ion        0    O1   0    O2   0    O3
CONNECT  BCTBK  O1  ion        0    C

CONNECT  BCT01  O2  ion        0    C
CONNECT  BCT01  O3  ion        0    C
CONNECT  BCT-1  O2  ion        0    C
CONNECT  BCT-1  O3  ion        0    C

CHARGE   BCT01  O2   0.00
CHARGE   BCT01  O3   0.00

CHARGE   BCT-1  O2  -0.50
CHARGE   BCT-1  O3  -0.50


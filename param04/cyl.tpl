#mg 8/3/06
# from cys
# when sg is a ligand; used for hemes
# charge removed from SG, no titration; no rotation
#
CONFLIST CYL        CYLBK CYL01 

NATOM    CYLBK      6
NATOM    CYL01      4

IATOM    CYLBK  N   0
IATOM    CYLBK  H   1
IATOM    CYLBK  CA  2
IATOM    CYLBK  HA  3
IATOM    CYLBK  C   4
IATOM    CYLBK  O   5
IATOM    CYL01  CB  0
IATOM    CYL01 1HB  1
IATOM    CYL01 2HB  2
IATOM    CYL01  SG  3

ATOMNAME CYLBK    0  N  
ATOMNAME CYLBK    1  H  
ATOMNAME CYLBK    2  CA 
ATOMNAME CYLBK    3  HA 
ATOMNAME CYLBK    4  C  
ATOMNAME CYLBK    5  O  
ATOMNAME CYL01    0  CB 
ATOMNAME CYL01    1 1HB 
ATOMNAME CYL01    2 2HB 
ATOMNAME CYL01    3  SG 

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  CYLBK  N   sp2       -1    C   0     CA  0     H
CONNECT  CYLBK  H   s         0     N
CONNECT  CYLBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  CYLBK  HA  s         0     CA
CONNECT  CYLBK  C   sp2       0     CA  0     O   1     N
CONNECT  CYLBK  O   sp2       0     C
CONNECT  CYL01  CB  sp3       0     CA  0     SG  0    1HB  0    2HB
CONNECT  CYL01 1HB  s         0     CB
CONNECT  CYL01 2HB  s         0     CB
CONNECT  CYL01  SG  sp2       0     CB  LIG   ?  

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   CYL    N   1.50
RADIUS   CYL    H   1.00
RADIUS   CYL    CA  2.00
RADIUS   CYL    HA  0.00
RADIUS   CYL    C   1.70
RADIUS   CYL    O   1.40
RADIUS   CYL    CB  2.00
RADIUS   CYL   1HB  0.00
RADIUS   CYL   2HB  2.00
RADIUS   CYL    SG  1.85

CHARGE   CYLBK  N    -0.350
CHARGE   CYLBK  H     0.250
CHARGE   CYLBK  CA    0.100
CHARGE   CYLBK  C     0.550
CHARGE   CYLBK  O    -0.550
CHARGE   CYL01  SG    0.000



### This is a temporary parameter file made for residue BKB ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST BKB        BKBBK 

NATOM    BKBBK      4

IATOM    BKBBK  N      0
IATOM    BKBBK  CA     1
IATOM    BKBBK  C      2
IATOM    BKBBK  O      3

ATOMNAME BKBBK    0  N  
ATOMNAME BKBBK    1  CA 
ATOMNAME BKBBK    2  C  
ATOMNAME BKBBK    3  O  

CONNECT  BKBBK  N   ion        0    CA 
CONNECT  BKBBK  CA  ion        0    N    0    C  
CONNECT  BKBBK  C   ion        0    CA   0    O  
CONNECT  BKBBK  O   ion        0    C  

### This is a temporary parameter file made for residue HEA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST HEA        HEABK 

NATOM    HEABK      44

IATOM    HEABK aCB     0
IATOM    HEABK aCG     1
IATOM    HEABK aND1    2
IATOM    HEABK aCD2    3
IATOM    HEABK aCE1    4
IATOM    HEABK aNE2    5
IATOM    HEABK bCB     6
IATOM    HEABK bCG     7
IATOM    HEABK bND1    8
IATOM    HEABK bCD2    9
IATOM    HEABK bCE1   10
IATOM    HEABK bNE2   11
IATOM    HEABK FE     12
IATOM    HEABK  CHA   13
IATOM    HEABK  CHB   14
IATOM    HEABK  CHC   15
IATOM    HEABK  CHD   16
IATOM    HEABK  N A   17
IATOM    HEABK  C1A   18
IATOM    HEABK  C2A   19
IATOM    HEABK  C3A   20
IATOM    HEABK  C4A   21
IATOM    HEABK  CMA   22
IATOM    HEABK  OMA   23
IATOM    HEABK  N B   24
IATOM    HEABK  C1B   25
IATOM    HEABK  C2B   26
IATOM    HEABK  C3B   27
IATOM    HEABK  C4B   28
IATOM    HEABK  CMB   29
IATOM    HEABK  N C   30
IATOM    HEABK  C1C   31
IATOM    HEABK  C2C   32
IATOM    HEABK  C3C   33
IATOM    HEABK  C4C   34
IATOM    HEABK  CMC   35
IATOM    HEABK  CAC   36
IATOM    HEABK  CBC   37
IATOM    HEABK  N D   38
IATOM    HEABK  C1D   39
IATOM    HEABK  C2D   40
IATOM    HEABK  C3D   41
IATOM    HEABK  C4D   42
IATOM    HEABK  CMD   43

ATOMNAME HEABK    0 aCB 
ATOMNAME HEABK    1 aCG 
ATOMNAME HEABK    2 aND1
ATOMNAME HEABK    3 aCD2
ATOMNAME HEABK    4 aCE1
ATOMNAME HEABK    5 aNE2
ATOMNAME HEABK    6 bCB 
ATOMNAME HEABK    7 bCG 
ATOMNAME HEABK    8 bND1
ATOMNAME HEABK    9 bCD2
ATOMNAME HEABK   10 bCE1
ATOMNAME HEABK   11 bNE2
ATOMNAME HEABK   12 FE  
ATOMNAME HEABK   13  CHA
ATOMNAME HEABK   14  CHB
ATOMNAME HEABK   15  CHC
ATOMNAME HEABK   16  CHD
ATOMNAME HEABK   17  N A
ATOMNAME HEABK   18  C1A
ATOMNAME HEABK   19  C2A
ATOMNAME HEABK   20  C3A
ATOMNAME HEABK   21  C4A
ATOMNAME HEABK   22  CMA
ATOMNAME HEABK   23  OMA
ATOMNAME HEABK   24  N B
ATOMNAME HEABK   25  C1B
ATOMNAME HEABK   26  C2B
ATOMNAME HEABK   27  C3B
ATOMNAME HEABK   28  C4B
ATOMNAME HEABK   29  CMB
ATOMNAME HEABK   30  N C
ATOMNAME HEABK   31  C1C
ATOMNAME HEABK   32  C2C
ATOMNAME HEABK   33  C3C
ATOMNAME HEABK   34  C4C
ATOMNAME HEABK   35  CMC
ATOMNAME HEABK   36  CAC
ATOMNAME HEABK   37  CBC
ATOMNAME HEABK   38  N D
ATOMNAME HEABK   39  C1D
ATOMNAME HEABK   40  C2D
ATOMNAME HEABK   41  C3D
ATOMNAME HEABK   42  C4D
ATOMNAME HEABK   43  CMD

CONNECT  HEABK aCB  ion        0   aCG 
CONNECT  HEABK aCG  ion        0   aCB   0   aND1  0   aCD2  0   aCE1  0   aNE2
CONNECT  HEABK aND1 ion        0   aCG   0   aCD2  0   aCE1  0   aNE2
CONNECT  HEABK aCD2 ion        0   aCG   0   aND1  0   aCE1  0   aNE2
CONNECT  HEABK aCE1 ion        0   aCG   0   aND1  0   aCD2  0   aNE2
CONNECT  HEABK aNE2 ion        0   aCG   0   aND1  0   aCD2  0   aCE1  0   FE  
CONNECT  HEABK bCB  ion        0   bCG 
CONNECT  HEABK bCG  ion        0   bCB   0   bND1  0   bCD2
CONNECT  HEABK bND1 ion        0   bCG   0   bCD2  0   bCE1  0   bNE2
CONNECT  HEABK bCD2 ion        0   bCG   0   bND1  0   bCE1  0   bNE2
CONNECT  HEABK bCE1 ion        0   bND1  0   bCD2  0   bNE2
CONNECT  HEABK bNE2 ion        0   bND1  0   bCD2  0   bCE1  0   FE  
CONNECT  HEABK FE   ion        0   aNE2  0   bNE2  0    N A  0    N B  0    N C  0    N D
CONNECT  HEABK  CHA ion        0    C1A  0    C4D
CONNECT  HEABK  CHB ion        0    C4A  0    C1B
CONNECT  HEABK  CHC ion        0    C4B  0    C1C
CONNECT  HEABK  CHD ion        0    C4C  0    C1D
CONNECT  HEABK  N A ion        0   FE    0    C1A  0    C4A
CONNECT  HEABK  C1A ion        0    CHA  0    N A  0    C2A  0    C4A
CONNECT  HEABK  C2A ion        0    C1A  0    C3A
CONNECT  HEABK  C3A ion        0    C2A  0    C4A  0    CMA
CONNECT  HEABK  C4A ion        0    CHB  0    N A  0    C1A  0    C3A
CONNECT  HEABK  CMA ion        0    C3A  0    OMA
CONNECT  HEABK  OMA ion        0    CMA
CONNECT  HEABK  N B ion        0   FE    0    C1B  0    C4B
CONNECT  HEABK  C1B ion        0    CHB  0    N B  0    C2B  0    C4B
CONNECT  HEABK  C2B ion        0    C1B  0    C3B  0    CMB
CONNECT  HEABK  C3B ion        0    C2B  0    C4B
CONNECT  HEABK  C4B ion        0    CHC  0    N B  0    C1B  0    C3B
CONNECT  HEABK  CMB ion        0    C2B
CONNECT  HEABK  N C ion        0   FE    0    C1C  0    C4C
CONNECT  HEABK  C1C ion        0    CHC  0    N C  0    C2C  0    C3C  0    C4C
CONNECT  HEABK  C2C ion        0    C1C  0    C3C  0    CMC
CONNECT  HEABK  C3C ion        0    C1C  0    C2C  0    C4C  0    CAC
CONNECT  HEABK  C4C ion        0    CHD  0    N C  0    C1C  0    C3C
CONNECT  HEABK  CMC ion        0    C2C
CONNECT  HEABK  CAC ion        0    C3C  0    CBC
CONNECT  HEABK  CBC ion        0    CAC
CONNECT  HEABK  N D ion        0   FE    0    C1D  0    C4D
CONNECT  HEABK  C1D ion        0    CHD  0    N D  0    C2D  0    C4D
CONNECT  HEABK  C2D ion        0    C1D  0    C3D  0    CMD
CONNECT  HEABK  C3D ion        0    C2D  0    C4D
CONNECT  HEABK  C4D ion        0    CHA  0    N D  0    C1D  0    C3D
CONNECT  HEABK  CMD ion        0    C2D

### This is a temporary parameter file made for residue CUB ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST CUB        CUBBK 

NATOM    CUBBK      19

IATOM    CUBBK aCB     0
IATOM    CUBBK aCG     1
IATOM    CUBBK aND1    2
IATOM    CUBBK aCD2    3
IATOM    CUBBK aCE1    4
IATOM    CUBBK aNE2    5
IATOM    CUBBK bCB     6
IATOM    CUBBK bCG     7
IATOM    CUBBK bND1    8
IATOM    CUBBK bCD2    9
IATOM    CUBBK bCE1   10
IATOM    CUBBK bNE2   11
IATOM    CUBBK cCB    12
IATOM    CUBBK cCG    13
IATOM    CUBBK cND1   14
IATOM    CUBBK cCD2   15
IATOM    CUBBK cCE1   16
IATOM    CUBBK cNE2   17
IATOM    CUBBK CU     18

ATOMNAME CUBBK    0 aCB 
ATOMNAME CUBBK    1 aCG 
ATOMNAME CUBBK    2 aND1
ATOMNAME CUBBK    3 aCD2
ATOMNAME CUBBK    4 aCE1
ATOMNAME CUBBK    5 aNE2
ATOMNAME CUBBK    6 bCB 
ATOMNAME CUBBK    7 bCG 
ATOMNAME CUBBK    8 bND1
ATOMNAME CUBBK    9 bCD2
ATOMNAME CUBBK   10 bCE1
ATOMNAME CUBBK   11 bNE2
ATOMNAME CUBBK   12 cCB 
ATOMNAME CUBBK   13 cCG 
ATOMNAME CUBBK   14 cND1
ATOMNAME CUBBK   15 cCD2
ATOMNAME CUBBK   16 cCE1
ATOMNAME CUBBK   17 cNE2
ATOMNAME CUBBK   18 CU  

CONNECT  CUBBK aCB  ion        0   aCG 
CONNECT  CUBBK aCG  ion        0   aCB   0   aND1  0   aCD2  0   aNE2
CONNECT  CUBBK aND1 ion        0   aCG   0   aCD2  0   aCE1  0   aNE2  0   CU  
CONNECT  CUBBK aCD2 ion        0   aCG   0   aND1  0   aNE2
CONNECT  CUBBK aCE1 ion        0   aND1  0   aNE2
CONNECT  CUBBK aNE2 ion        0   aCG   0   aND1  0   aCD2  0   aCE1
CONNECT  CUBBK bCB  ion        0   bCG 
CONNECT  CUBBK bCG  ion        0   bCB   0   bND1  0   bCD2  0   bNE2
CONNECT  CUBBK bND1 ion        0   bCG   0   bCD2  0   bCE1  0   bNE2
CONNECT  CUBBK bCD2 ion        0   bCG   0   bND1  0   bCE1  0   bNE2
CONNECT  CUBBK bCE1 ion        0   bND1  0   bCD2  0   bNE2
CONNECT  CUBBK bNE2 ion        0   bCG   0   bND1  0   bCD2  0   bCE1  0   CU  
CONNECT  CUBBK cCB  ion        0   cCG 
CONNECT  CUBBK cCG  ion        0   cCB   0   cND1  0   cCD2  0   cCE1  0   cNE2
CONNECT  CUBBK cND1 ion        0   cCG   0   cCD2  0   cCE1  0   cNE2
CONNECT  CUBBK cCD2 ion        0   cCG   0   cND1  0   cCE1  0   cNE2
CONNECT  CUBBK cCE1 ion        0   cCG   0   cND1  0   cCD2  0   cNE2
CONNECT  CUBBK cNE2 ion        0   cCG   0   cND1  0   cCD2  0   cCE1  0   CU  
CONNECT  CUBBK CU   ion        0   aND1  0   bNE2  0   cNE2

### This is a temporary parameter file made for residue TYF ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST TYF        TYFBK 

NATOM    TYFBK      12

IATOM    TYFBK  N      0
IATOM    TYFBK  CA     1
IATOM    TYFBK  C      2
IATOM    TYFBK  O      3
IATOM    TYFBK  CB     4
IATOM    TYFBK  CG     5
IATOM    TYFBK  CD1    6
IATOM    TYFBK  CD2    7
IATOM    TYFBK  CE1    8
IATOM    TYFBK  CE2    9
IATOM    TYFBK  CZ    10
IATOM    TYFBK  OH    11

ATOMNAME TYFBK    0  N  
ATOMNAME TYFBK    1  CA 
ATOMNAME TYFBK    2  C  
ATOMNAME TYFBK    3  O  
ATOMNAME TYFBK    4  CB 
ATOMNAME TYFBK    5  CG 
ATOMNAME TYFBK    6  CD1
ATOMNAME TYFBK    7  CD2
ATOMNAME TYFBK    8  CE1
ATOMNAME TYFBK    9  CE2
ATOMNAME TYFBK   10  CZ 
ATOMNAME TYFBK   11  OH 

CONNECT  TYFBK  N   ion        0    CA 
CONNECT  TYFBK  CA  ion        0    N    0    C    0    CB 
CONNECT  TYFBK  C   ion        0    CA   0    O  
CONNECT  TYFBK  O   ion        0    C  
CONNECT  TYFBK  CB  ion        0    CA   0    CG 
CONNECT  TYFBK  CG  ion        0    CB   0    CD1  0    CD2
CONNECT  TYFBK  CD1 ion        0    CG   0    CE1
CONNECT  TYFBK  CD2 ion        0    CG   0    CE2
CONNECT  TYFBK  CE1 ion        0    CD1  0    CZ 
CONNECT  TYFBK  CE2 ion        0    CD2  0    CZ 
CONNECT  TYFBK  CZ  ion        0    CE1  0    CE2  0    OH 
CONNECT  TYFBK  OH  ion        0    CZ 

### This is a temporary parameter file made for residue HE3 ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST HE3        HE3BK 

NATOM    HE3BK      38

IATOM    HE3BK  CB     0
IATOM    HE3BK  CG     1
IATOM    HE3BK  ND1    2
IATOM    HE3BK  CD2    3
IATOM    HE3BK  CE1    4
IATOM    HE3BK  NE2    5
IATOM    HE3BK FE      6
IATOM    HE3BK  CHA    7
IATOM    HE3BK  CHB    8
IATOM    HE3BK  CHC    9
IATOM    HE3BK  CHD   10
IATOM    HE3BK  N A   11
IATOM    HE3BK  C1A   12
IATOM    HE3BK  C2A   13
IATOM    HE3BK  C3A   14
IATOM    HE3BK  C4A   15
IATOM    HE3BK  CMA   16
IATOM    HE3BK  OMA   17
IATOM    HE3BK  N B   18
IATOM    HE3BK  C1B   19
IATOM    HE3BK  C2B   20
IATOM    HE3BK  C3B   21
IATOM    HE3BK  C4B   22
IATOM    HE3BK  CMB   23
IATOM    HE3BK  N C   24
IATOM    HE3BK  C1C   25
IATOM    HE3BK  C2C   26
IATOM    HE3BK  C3C   27
IATOM    HE3BK  C4C   28
IATOM    HE3BK  CMC   29
IATOM    HE3BK  CAC   30
IATOM    HE3BK  CBC   31
IATOM    HE3BK  N D   32
IATOM    HE3BK  C1D   33
IATOM    HE3BK  C2D   34
IATOM    HE3BK  C3D   35
IATOM    HE3BK  C4D   36
IATOM    HE3BK  CMD   37

ATOMNAME HE3BK    0  CB 
ATOMNAME HE3BK    1  CG 
ATOMNAME HE3BK    2  ND1
ATOMNAME HE3BK    3  CD2
ATOMNAME HE3BK    4  CE1
ATOMNAME HE3BK    5  NE2
ATOMNAME HE3BK    6 FE  
ATOMNAME HE3BK    7  CHA
ATOMNAME HE3BK    8  CHB
ATOMNAME HE3BK    9  CHC
ATOMNAME HE3BK   10  CHD
ATOMNAME HE3BK   11  N A
ATOMNAME HE3BK   12  C1A
ATOMNAME HE3BK   13  C2A
ATOMNAME HE3BK   14  C3A
ATOMNAME HE3BK   15  C4A
ATOMNAME HE3BK   16  CMA
ATOMNAME HE3BK   17  OMA
ATOMNAME HE3BK   18  N B
ATOMNAME HE3BK   19  C1B
ATOMNAME HE3BK   20  C2B
ATOMNAME HE3BK   21  C3B
ATOMNAME HE3BK   22  C4B
ATOMNAME HE3BK   23  CMB
ATOMNAME HE3BK   24  N C
ATOMNAME HE3BK   25  C1C
ATOMNAME HE3BK   26  C2C
ATOMNAME HE3BK   27  C3C
ATOMNAME HE3BK   28  C4C
ATOMNAME HE3BK   29  CMC
ATOMNAME HE3BK   30  CAC
ATOMNAME HE3BK   31  CBC
ATOMNAME HE3BK   32  N D
ATOMNAME HE3BK   33  C1D
ATOMNAME HE3BK   34  C2D
ATOMNAME HE3BK   35  C3D
ATOMNAME HE3BK   36  C4D
ATOMNAME HE3BK   37  CMD

CONNECT  HE3BK  CB  ion        0    CG 
CONNECT  HE3BK  CG  ion        0    CB   0    ND1  0    CD2
CONNECT  HE3BK  ND1 ion        0    CG   0    CD2  0    CE1  0    NE2
CONNECT  HE3BK  CD2 ion        0    CG   0    ND1  0    CE1  0    NE2
CONNECT  HE3BK  CE1 ion        0    ND1  0    CD2  0    NE2
CONNECT  HE3BK  NE2 ion        0    ND1  0    CD2  0    CE1
CONNECT  HE3BK FE   ion        0    N A  0    N B  0    N C  0    N D
CONNECT  HE3BK  CHA ion        0    C1A  0    C4D
CONNECT  HE3BK  CHB ion        0    C4A  0    C1B
CONNECT  HE3BK  CHC ion        0    C4B  0    C1C
CONNECT  HE3BK  CHD ion        0    C4C  0    C1D
CONNECT  HE3BK  N A ion        0   FE    0    C1A  0    C4A
CONNECT  HE3BK  C1A ion        0    CHA  0    N A  0    C2A  0    C4A
CONNECT  HE3BK  C2A ion        0    C1A  0    C3A
CONNECT  HE3BK  C3A ion        0    C2A  0    C4A  0    CMA
CONNECT  HE3BK  C4A ion        0    CHB  0    N A  0    C1A  0    C3A
CONNECT  HE3BK  CMA ion        0    C3A  0    OMA
CONNECT  HE3BK  OMA ion        0    CMA
CONNECT  HE3BK  N B ion        0   FE    0    C1B  0    C4B
CONNECT  HE3BK  C1B ion        0    CHB  0    N B  0    C2B  0    C4B
CONNECT  HE3BK  C2B ion        0    C1B  0    C3B  0    CMB
CONNECT  HE3BK  C3B ion        0    C2B  0    C4B
CONNECT  HE3BK  C4B ion        0    CHC  0    N B  0    C1B  0    C3B
CONNECT  HE3BK  CMB ion        0    C2B
CONNECT  HE3BK  N C ion        0   FE    0    C1C  0    C4C
CONNECT  HE3BK  C1C ion        0    CHC  0    N C  0    C2C  0    C4C
CONNECT  HE3BK  C2C ion        0    C1C  0    C3C  0    CMC
CONNECT  HE3BK  C3C ion        0    C2C  0    C4C  0    CAC
CONNECT  HE3BK  C4C ion        0    CHD  0    N C  0    C1C  0    C3C
CONNECT  HE3BK  CMC ion        0    C2C
CONNECT  HE3BK  CAC ion        0    C3C  0    CBC
CONNECT  HE3BK  CBC ion        0    CAC
CONNECT  HE3BK  N D ion        0   FE    0    C1D  0    C4D
CONNECT  HE3BK  C1D ion        0    CHD  0    N D  0    C2D  0    C4D
CONNECT  HE3BK  C2D ion        0    C1D  0    C3D  0    CMD
CONNECT  HE3BK  C3D ion        0    C2D  0    C4D
CONNECT  HE3BK  C4D ion        0    CHA  0    N D  0    C1D  0    C3D
CONNECT  HE3BK  CMD ion        0    C2D

### This is a temporary parameter file made for residue CUA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST CUA        CUABK 

NATOM    CUABK      18

IATOM    CUABK bCB     0
IATOM    CUABK bCG     1
IATOM    CUABK bND1    2
IATOM    CUABK bCD2    3
IATOM    CUABK bCE1    4
IATOM    CUABK bNE2    5
IATOM    CUABK cCB     6
IATOM    CUABK cSG     7
IATOM    CUABK dCB     8
IATOM    CUABK dSG     9
IATOM    CUABK aCB    10
IATOM    CUABK aCG    11
IATOM    CUABK aND1   12
IATOM    CUABK aCD2   13
IATOM    CUABK aCE1   14
IATOM    CUABK aNE2   15
IATOM    CUABK CU1    16
IATOM    CUABK CU2    17

ATOMNAME CUABK    0 bCB 
ATOMNAME CUABK    1 bCG 
ATOMNAME CUABK    2 bND1
ATOMNAME CUABK    3 bCD2
ATOMNAME CUABK    4 bCE1
ATOMNAME CUABK    5 bNE2
ATOMNAME CUABK    6 cCB 
ATOMNAME CUABK    7 cSG 
ATOMNAME CUABK    8 dCB 
ATOMNAME CUABK    9 dSG 
ATOMNAME CUABK   10 aCB 
ATOMNAME CUABK   11 aCG 
ATOMNAME CUABK   12 aND1
ATOMNAME CUABK   13 aCD2
ATOMNAME CUABK   14 aCE1
ATOMNAME CUABK   15 aNE2
ATOMNAME CUABK   16 CU1 
ATOMNAME CUABK   17 CU2 

CONNECT  CUABK bCB  ion        0   bCG 
CONNECT  CUABK bCG  ion        0   bCB   0   bND1  0   bCD2  0   bNE2
CONNECT  CUABK bND1 ion        0   bCG   0   bCD2  0   bCE1  0   bNE2  0   CU2 
CONNECT  CUABK bCD2 ion        0   bCG   0   bND1  0   bNE2
CONNECT  CUABK bCE1 ion        0   bND1  0   bNE2
CONNECT  CUABK bNE2 ion        0   bCG   0   bND1  0   bCD2  0   bCE1
CONNECT  CUABK cCB  ion        0   cSG 
CONNECT  CUABK cSG  ion        0   cCB   0   CU1 
CONNECT  CUABK dCB  ion        0   dSG 
CONNECT  CUABK dSG  ion        0   dCB 
CONNECT  CUABK aCB  ion        0   aCG 
CONNECT  CUABK aCG  ion        0   aCB   0   aND1  0   aCD2  0   aNE2
CONNECT  CUABK aND1 ion        0   aCG   0   aCE1  0   aNE2  0   CU1 
CONNECT  CUABK aCD2 ion        0   aCG   0   aNE2
CONNECT  CUABK aCE1 ion        0   aND1  0   aNE2
CONNECT  CUABK aNE2 ion        0   aCG   0   aND1  0   aCD2  0   aCE1
CONNECT  CUABK CU1  ion        0   cSG   0   aND1
CONNECT  CUABK CU2  ion        0   bND1

### This is a temporary parameter file made for residue _MG ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _MG        _MGBK 

NATOM    _MGBK      1

IATOM    _MGBK MG      0

ATOMNAME _MGBK    0 MG  

CONNECT  _MGBK MG   ion      

### This is a temporary parameter file made for residue _CA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST _CA        _CABK 

NATOM    _CABK      1

IATOM    _CABK CA      0

ATOMNAME _CABK    0 CA  

CONNECT  _CABK CA   ion      

### This is a temporary parameter file made for residue PAA ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST PAA        PAABK 

NATOM    PAABK      5

IATOM    PAABK  CAA    0
IATOM    PAABK  CBA    1
IATOM    PAABK  CGA    2
IATOM    PAABK  O1A    3
IATOM    PAABK  O2A    4

ATOMNAME PAABK    0  CAA
ATOMNAME PAABK    1  CBA
ATOMNAME PAABK    2  CGA
ATOMNAME PAABK    3  O1A
ATOMNAME PAABK    4  O2A

CONNECT  PAABK  CAA ion        0    CBA
CONNECT  PAABK  CBA ion        0    CAA  0    CGA
CONNECT  PAABK  CGA ion        0    CBA  0    O1A  0    O2A
CONNECT  PAABK  O1A ion        0    CGA  0    O2A
CONNECT  PAABK  O2A ion        0    CGA  0    O1A

### This is a temporary parameter file made for residue PDD ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST PDD        PDDBK 

NATOM    PDDBK      5

IATOM    PDDBK  CAD    0
IATOM    PDDBK  CBD    1
IATOM    PDDBK  CGD    2
IATOM    PDDBK  O1D    3
IATOM    PDDBK  O2D    4

ATOMNAME PDDBK    0  CAD
ATOMNAME PDDBK    1  CBD
ATOMNAME PDDBK    2  CGD
ATOMNAME PDDBK    3  O1D
ATOMNAME PDDBK    4  O2D

CONNECT  PDDBK  CAD ion        0    CBD
CONNECT  PDDBK  CBD ion        0    CAD  0    CGD
CONNECT  PDDBK  CGD ion        0    CBD  0    O1D  0    O2D
CONNECT  PDDBK  O1D ion        0    CGD  0    O2D
CONNECT  PDDBK  O2D ion        0    CGD  0    O1D

### This is a temporary parameter file made for residue FAL ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST FAL        FALBK 

NATOM    FALBK      18

IATOM    FALBK  C11    0
IATOM    FALBK  O11    1
IATOM    FALBK  C12    2
IATOM    FALBK  C13    3
IATOM    FALBK  C14    4
IATOM    FALBK  C15    5
IATOM    FALBK  C16    6
IATOM    FALBK  C17    7
IATOM    FALBK  C18    8
IATOM    FALBK  C19    9
IATOM    FALBK  C20   10
IATOM    FALBK  C21   11
IATOM    FALBK  C22   12
IATOM    FALBK  C23   13
IATOM    FALBK  C24   14
IATOM    FALBK  C25   15
IATOM    FALBK  C26   16
IATOM    FALBK  C27   17

ATOMNAME FALBK    0  C11
ATOMNAME FALBK    1  O11
ATOMNAME FALBK    2  C12
ATOMNAME FALBK    3  C13
ATOMNAME FALBK    4  C14
ATOMNAME FALBK    5  C15
ATOMNAME FALBK    6  C16
ATOMNAME FALBK    7  C17
ATOMNAME FALBK    8  C18
ATOMNAME FALBK    9  C19
ATOMNAME FALBK   10  C20
ATOMNAME FALBK   11  C21
ATOMNAME FALBK   12  C22
ATOMNAME FALBK   13  C23
ATOMNAME FALBK   14  C24
ATOMNAME FALBK   15  C25
ATOMNAME FALBK   16  C26
ATOMNAME FALBK   17  C27

CONNECT  FALBK  C11 ion        0    O11  0    C12
CONNECT  FALBK  O11 ion        0    C11
CONNECT  FALBK  C12 ion        0    C11  0    C13
CONNECT  FALBK  C13 ion        0    C12  0    C14
CONNECT  FALBK  C14 ion        0    C13  0    C15
CONNECT  FALBK  C15 ion        0    C14  0    C16  0    C26
CONNECT  FALBK  C16 ion        0    C15  0    C17
CONNECT  FALBK  C17 ion        0    C16  0    C18
CONNECT  FALBK  C18 ion        0    C17  0    C19
CONNECT  FALBK  C19 ion        0    C18  0    C20  0    C27
CONNECT  FALBK  C20 ion        0    C19  0    C21
CONNECT  FALBK  C21 ion        0    C20  0    C22
CONNECT  FALBK  C22 ion        0    C21  0    C23
CONNECT  FALBK  C23 ion        0    C22  0    C24  0    C25
CONNECT  FALBK  C24 ion        0    C23
CONNECT  FALBK  C25 ion        0    C23
CONNECT  FALBK  C26 ion        0    C15
CONNECT  FALBK  C27 ion        0    C19

### This is a temporary parameter file made for residue PEH ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST PEH        PEHBK 

NATOM    PEHBK      51

IATOM    PEHBK  P      0
IATOM    PEHBK  N      1
IATOM    PEHBK  O11    2
IATOM    PEHBK  O12    3
IATOM    PEHBK  O13    4
IATOM    PEHBK  O14    5
IATOM    PEHBK  C11    6
IATOM    PEHBK  C12    7
IATOM    PEHBK  C1     8
IATOM    PEHBK  C2     9
IATOM    PEHBK  C3    10
IATOM    PEHBK  O31   11
IATOM    PEHBK  O32   12
IATOM    PEHBK  C31   13
IATOM    PEHBK  C32   14
IATOM    PEHBK  C33   15
IATOM    PEHBK  C34   16
IATOM    PEHBK  C35   17
IATOM    PEHBK  C36   18
IATOM    PEHBK  C37   19
IATOM    PEHBK  C38   20
IATOM    PEHBK  C39   21
IATOM    PEHBK  C3A   22
IATOM    PEHBK  C3B   23
IATOM    PEHBK  C3C   24
IATOM    PEHBK  C3D   25
IATOM    PEHBK  C3E   26
IATOM    PEHBK  C3F   27
IATOM    PEHBK  C3G   28
IATOM    PEHBK  C3H   29
IATOM    PEHBK  C3I   30
IATOM    PEHBK  O21   31
IATOM    PEHBK  O22   32
IATOM    PEHBK  C21   33
IATOM    PEHBK  C22   34
IATOM    PEHBK  C23   35
IATOM    PEHBK  C24   36
IATOM    PEHBK  C25   37
IATOM    PEHBK  C26   38
IATOM    PEHBK  C27   39
IATOM    PEHBK  C28   40
IATOM    PEHBK  C29   41
IATOM    PEHBK  C2A   42
IATOM    PEHBK  C2B   43
IATOM    PEHBK  C2C   44
IATOM    PEHBK  C2D   45
IATOM    PEHBK  C2E   46
IATOM    PEHBK  C2F   47
IATOM    PEHBK  C2G   48
IATOM    PEHBK  C2H   49
IATOM    PEHBK  C2I   50

ATOMNAME PEHBK    0  P  
ATOMNAME PEHBK    1  N  
ATOMNAME PEHBK    2  O11
ATOMNAME PEHBK    3  O12
ATOMNAME PEHBK    4  O13
ATOMNAME PEHBK    5  O14
ATOMNAME PEHBK    6  C11
ATOMNAME PEHBK    7  C12
ATOMNAME PEHBK    8  C1 
ATOMNAME PEHBK    9  C2 
ATOMNAME PEHBK   10  C3 
ATOMNAME PEHBK   11  O31
ATOMNAME PEHBK   12  O32
ATOMNAME PEHBK   13  C31
ATOMNAME PEHBK   14  C32
ATOMNAME PEHBK   15  C33
ATOMNAME PEHBK   16  C34
ATOMNAME PEHBK   17  C35
ATOMNAME PEHBK   18  C36
ATOMNAME PEHBK   19  C37
ATOMNAME PEHBK   20  C38
ATOMNAME PEHBK   21  C39
ATOMNAME PEHBK   22  C3A
ATOMNAME PEHBK   23  C3B
ATOMNAME PEHBK   24  C3C
ATOMNAME PEHBK   25  C3D
ATOMNAME PEHBK   26  C3E
ATOMNAME PEHBK   27  C3F
ATOMNAME PEHBK   28  C3G
ATOMNAME PEHBK   29  C3H
ATOMNAME PEHBK   30  C3I
ATOMNAME PEHBK   31  O21
ATOMNAME PEHBK   32  O22
ATOMNAME PEHBK   33  C21
ATOMNAME PEHBK   34  C22
ATOMNAME PEHBK   35  C23
ATOMNAME PEHBK   36  C24
ATOMNAME PEHBK   37  C25
ATOMNAME PEHBK   38  C26
ATOMNAME PEHBK   39  C27
ATOMNAME PEHBK   40  C28
ATOMNAME PEHBK   41  C29
ATOMNAME PEHBK   42  C2A
ATOMNAME PEHBK   43  C2B
ATOMNAME PEHBK   44  C2C
ATOMNAME PEHBK   45  C2D
ATOMNAME PEHBK   46  C2E
ATOMNAME PEHBK   47  C2F
ATOMNAME PEHBK   48  C2G
ATOMNAME PEHBK   49  C2H
ATOMNAME PEHBK   50  C2I

CONNECT  PEHBK  P   ion        0    O11  0    O12  0    O13  0    O14
CONNECT  PEHBK  N   ion        0    C12
CONNECT  PEHBK  O11 ion        0    P    0    C1 
CONNECT  PEHBK  O12 ion        0    P    0    C11
CONNECT  PEHBK  O13 ion        0    P  
CONNECT  PEHBK  O14 ion        0    P  
CONNECT  PEHBK  C11 ion        0    O12  0    C12
CONNECT  PEHBK  C12 ion        0    N    0    C11
CONNECT  PEHBK  C1  ion        0    O11  0    C2 
CONNECT  PEHBK  C2  ion        0    C1   0    C3   0    O21
CONNECT  PEHBK  C3  ion        0    C2   0    O31
CONNECT  PEHBK  O31 ion        0    C3   0    C31
CONNECT  PEHBK  O32 ion        0    C31
CONNECT  PEHBK  C31 ion        0    O31  0    O32  0    C32
CONNECT  PEHBK  C32 ion        0    C31  0    C33
CONNECT  PEHBK  C33 ion        0    C32  0    C34
CONNECT  PEHBK  C34 ion        0    C33  0    C35
CONNECT  PEHBK  C35 ion        0    C34  0    C36
CONNECT  PEHBK  C36 ion        0    C35  0    C37
CONNECT  PEHBK  C37 ion        0    C36  0    C38
CONNECT  PEHBK  C38 ion        0    C37  0    C39
CONNECT  PEHBK  C39 ion        0    C38  0    C3A
CONNECT  PEHBK  C3A ion        0    C39  0    C3B
CONNECT  PEHBK  C3B ion        0    C3A  0    C3C
CONNECT  PEHBK  C3C ion        0    C3B  0    C3D
CONNECT  PEHBK  C3D ion        0    C3C  0    C3E
CONNECT  PEHBK  C3E ion        0    C3D  0    C3F
CONNECT  PEHBK  C3F ion        0    C3E  0    C3G
CONNECT  PEHBK  C3G ion        0    C3F  0    C3H
CONNECT  PEHBK  C3H ion        0    C3G  0    C3I
CONNECT  PEHBK  C3I ion        0    C3H
CONNECT  PEHBK  O21 ion        0    C2   0    C21
CONNECT  PEHBK  O22 ion        0    C21
CONNECT  PEHBK  C21 ion        0    O21  0    O22  0    C22
CONNECT  PEHBK  C22 ion        0    C21  0    C23
CONNECT  PEHBK  C23 ion        0    C22  0    C24
CONNECT  PEHBK  C24 ion        0    C23  0    C25
CONNECT  PEHBK  C25 ion        0    C24  0    C26
CONNECT  PEHBK  C26 ion        0    C25  0    C27
CONNECT  PEHBK  C27 ion        0    C26  0    C28
CONNECT  PEHBK  C28 ion        0    C27  0    C29
CONNECT  PEHBK  C29 ion        0    C28  0    C2A
CONNECT  PEHBK  C2A ion        0    C29  0    C2B
CONNECT  PEHBK  C2B ion        0    C2A  0    C2C
CONNECT  PEHBK  C2C ion        0    C2B  0    C2D
CONNECT  PEHBK  C2D ion        0    C2C  0    C2E
CONNECT  PEHBK  C2E ion        0    C2D  0    C2F
CONNECT  PEHBK  C2F ion        0    C2E  0    C2G
CONNECT  PEHBK  C2G ion        0    C2F  0    C2H
CONNECT  PEHBK  C2H ion        0    C2G  0    C2I
CONNECT  PEHBK  C2I ion        0    C2H


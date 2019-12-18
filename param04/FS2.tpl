###------#
#
#  FE2-S2 Cluster - BC1 complex
#
# Charges -2/-1
#                          by Jun & cf 2006
###------#

CONFLIST FS2        FS2BK FS201 FS202

DEL_HYDR FS2          f

NATOM    FS2BK      0
NATOM    FS201      32
NATOM    FS202      34

IATOM    FS201 aCB     0
IATOM    FS201 aHB1    1
IATOM    FS201 aHB2    2
IATOM    FS201 aCG     3
IATOM    FS201 aND1    4
IATOM    FS201 aCD2    5
IATOM    FS201 aCE1    6
IATOM    FS201 aNE2    7
IATOM    FS201 aHD2    8
IATOM    FS201 aHE1    9
IATOM    FS201 bCB    10
IATOM    FS201 bHB1   11
IATOM    FS201 bHB2   12
IATOM    FS201 bCG    13
IATOM    FS201 bND1   14
IATOM    FS201 bCD2   15
IATOM    FS201 bCE1   16
IATOM    FS201 bNE2   17
IATOM    FS201 bHD2   18
IATOM    FS201 bHE1   19
IATOM    FS201 cCB    20
IATOM    FS201 cHB1   21
IATOM    FS201 cHB2   22
IATOM    FS201 cSG    23
IATOM    FS201 dCB    24
IATOM    FS201 dHB1   25
IATOM    FS201 dHB2   26
IATOM    FS201 dSG    27
IATOM    FS201  S1    28
IATOM    FS201  FE1   29
IATOM    FS201  S2    30
IATOM    FS201  FE2   31

ATOMNAME FS201    0 aCB  
ATOMNAME FS201    1 aHB1 
ATOMNAME FS201    2 aHB2 
ATOMNAME FS201    3 aCG  
ATOMNAME FS201    4 aND1 
ATOMNAME FS201    5 aCD2 
ATOMNAME FS201    6 aCE1 
ATOMNAME FS201    7 aNE2 
ATOMNAME FS201    8 aHD2 
ATOMNAME FS201    9 aHE1 
ATOMNAME FS201   10 bCB  
ATOMNAME FS201   11 bHB1 
ATOMNAME FS201   12 bHB2 
ATOMNAME FS201   13 bCG  
ATOMNAME FS201   14 bND1 
ATOMNAME FS201   15 bCD2 
ATOMNAME FS201   16 bCE1 
ATOMNAME FS201   17 bNE2 
ATOMNAME FS201   18 bHD2 
ATOMNAME FS201   19 bHE1                 
ATOMNAME FS201   20 cCB  
ATOMNAME FS201   21 cHB1 
ATOMNAME FS201   22 cHB2 
ATOMNAME FS201   23 cSG  
ATOMNAME FS201   24 dCB  
ATOMNAME FS201   25 dHB1 
ATOMNAME FS201   26 dHB2 
ATOMNAME FS201   27 dSG  
ATOMNAME FS201   28  S1  
ATOMNAME FS201   29  FE1 
ATOMNAME FS201   30  S2  
ATOMNAME FS201   31  FE2 

#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS201 dCB  sp3        0   dSG   0   dHB1  0   dHB2  LIG  CA
CONNECT  FS201 dHB1 s          0   dCB
CONNECT  FS201 dHB2 s          0   dCB
CONNECT  FS201 dSG  ion        0   dCB   0    FE1
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS201 cCB  sp3        0   cSG   0   cHB1  0   cHB2  LIG  CA
CONNECT  FS201 cHB1 s          0   cCB
CONNECT  FS201 cHB2 s          0   cCB
CONNECT  FS201 cSG  ion        0   cCB   0    FE1
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS201 aCB  sp3        0   aCG   0   aHB1  0   aHB2  LIG  CA
CONNECT  FS201 aHB1 s          0   aCB
CONNECT  FS201 aHB2 s          0   aCB
CONNECT  FS201 aCG  sp2        0   aCB   0   aND1  0   aCD2
CONNECT  FS201 aND1 sp3        0   aCG   0    FE2  0   aCE1  
CONNECT  FS201 aCD2 sp2        0   aCG   0   aHD2  0   aNE2
CONNECT  FS201 aCE1 sp2        0   aND1  0   aHE1  0   aNE2
CONNECT  FS201 aNE2 sp2        0   aCD2  0   aCE1
CONNECT  FS201 aHD2 s          0   aCD2
CONNECT  FS201 aHE1 s          0   aCE1
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS201 bCB  sp3        0   bCG   0   bHB1  0   bHB2  LIG  CA
CONNECT  FS201 bHB1 s          0   bCB
CONNECT  FS201 bHB2 s          0   bCB
CONNECT  FS201 bCG  sp2        0   bCB   0   bND1  0   bCD2
CONNECT  FS201 bND1 sp3        0   bCG   0    FE2  0   bCE1  
CONNECT  FS201 bCD2 sp2        0   bCG   0   bHD2  0   bNE2
CONNECT  FS201 bCE1 sp2        0   bND1  0   bHE1  0   bNE2
CONNECT  FS201 bNE2 sp2        0   bCD2  0   bCE1 
CONNECT  FS201 bHD2 s          0   bCD2
CONNECT  FS201 bHE1 s          0   bCE1
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS201  FE1 sp3        0   cSG   0   dSG   0    S1   0    S2
CONNECT  FS201  FE2 sp3        0   aND1  0   bND1  0    S1   0    S2
CONNECT  FS201  S1  ion        0    FE1  0    FE2
CONNECT  FS201  S2  ion        0    FE1  0    FE2

IATOM    FS202 aCB     0
IATOM    FS202 aHB1    1
IATOM    FS202 aHB2    2
IATOM    FS202 aCG     3
IATOM    FS202 aND1    4
IATOM    FS202 aCD2    5
IATOM    FS202 aCE1    6
IATOM    FS202 aNE2    7
IATOM    FS202 aHE2    8
IATOM    FS202 aHD2    9
IATOM    FS202 aHE1   10
IATOM    FS202 bCB    11
IATOM    FS202 bHB1   12
IATOM    FS202 bHB2   13
IATOM    FS202 bCG    14
IATOM    FS202 bND1   15
IATOM    FS202 bCD2   16
IATOM    FS202 bCE1   17
IATOM    FS202 bNE2   18
IATOM    FS202 bHE2   19
IATOM    FS202 bHD2   20
IATOM    FS202 bHE1   21
IATOM    FS202 cCB    22
IATOM    FS202 cHB1   23
IATOM    FS202 cHB2   24
IATOM    FS202 cSG    25
IATOM    FS202 dCB    26
IATOM    FS202 dHB1   27
IATOM    FS202 dHB2   28
IATOM    FS202 dSG    29
IATOM    FS202  S1    30
IATOM    FS202  FE1   31
IATOM    FS202  S2    32
IATOM    FS202  FE2   33

                      
ATOMNAME FS202    0 aCB 
ATOMNAME FS202    1 aHB1
ATOMNAME FS202    2 aHB2
ATOMNAME FS202    3 aCG 
ATOMNAME FS202    4 aND1
ATOMNAME FS202    5 aCD2
ATOMNAME FS202    6 aCE1
ATOMNAME FS202    7 aNE2
ATOMNAME FS202    8 aHE2
ATOMNAME FS202    9 aHD2
ATOMNAME FS202   10 aHE1
                        
ATOMNAME FS202   11 bCB 
ATOMNAME FS202   12 bHB1
ATOMNAME FS202   13 bHB2
ATOMNAME FS202   14 bCG 
ATOMNAME FS202   15 bND1
ATOMNAME FS202   16 bCD2
ATOMNAME FS202   17 bCE1
ATOMNAME FS202   18 bNE2
ATOMNAME FS202   19 bHE2
ATOMNAME FS202   20 bHD2
ATOMNAME FS202   21 bHE1
                        
ATOMNAME FS202   22 cCB 
ATOMNAME FS202   23 cHB1
ATOMNAME FS202   24 cHB2
ATOMNAME FS202   25 cSG 
                        
ATOMNAME FS202   26 dCB 
ATOMNAME FS202   27 dHB1
ATOMNAME FS202   28 dHB2
ATOMNAME FS202   29 dSG 

ATOMNAME FS202   30  S1
ATOMNAME FS202   31  FE1
ATOMNAME FS202   32  S2
ATOMNAME FS202   33  FE2

#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS202 dCB  sp3        0   dSG   0   dHB1  0   dHB2  LIG  CA
CONNECT  FS202 dHB1 s          0   dCB
CONNECT  FS202 dHB2 s          0   dCB
CONNECT  FS202 dSG  ion        0   dCB   0    FE1
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS202 cCB  sp3        0   cSG   0   cHB1  0   cHB2  LIG  CA
CONNECT  FS202 cHB1 s          0   cCB
CONNECT  FS202 cHB2 s          0   cCB
CONNECT  FS202 cSG  ion        0   cCB   0    FE1
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS202 aCB  sp3        0   aCG   0   aHB1  0   aHB2  LIG  CA
CONNECT  FS202 aHB1 s          0   aCB
CONNECT  FS202 aHB2 s          0   aCB
CONNECT  FS202 aCG  sp2        0   aCB   0   aND1  0   aCD2
CONNECT  FS202 aND1 sp3        0   aCG   0    FE2  0   aCE1  
CONNECT  FS202 aCD2 sp2        0   aCG   0   aHD2  0   aNE2
CONNECT  FS202 aCE1 sp2        0   aND1  0   aHE1  0   aNE2
CONNECT  FS202 aNE2 sp2        0   aCD2  0   aCE1  0   aHE2
CONNECT  FS202 aHE2 s          0   aNE2
CONNECT  FS202 aHD2 s          0   aCD2
CONNECT  FS202 aHE1 s          0   aCE1
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS202 bCB  sp3        0   bCG   0   bHB1  0   bHB2  LIG  CA
CONNECT  FS202 bHB1 s          0   bCB
CONNECT  FS202 bHB2 s          0   bCB
CONNECT  FS202 bCG  sp2        0   bCB   0   bND1  0   bCD2
CONNECT  FS202 bND1 sp3        0   bCG   0    FE2  0   bCE1  
CONNECT  FS202 bCD2 sp2        0   bCG   0   bHD2  0   bNE2
CONNECT  FS202 bCE1 sp2        0   bND1  0   bHE1  0   bNE2
CONNECT  FS202 bNE2 sp2        0   bCD2  0   bCE1  0   bHE2 
CONNECT  FS202 bHE2 s          0   bNE2
CONNECT  FS202 bHD2 s          0   bCD2
CONNECT  FS202 bHE1 s          0   bCE1
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  FS202  FE1 sp3        0   cSG   0   dSG   0    S1   0    S2
CONNECT  FS202  FE2 sp3        0   aND1  0   bND1  0    S1   0    S2
CONNECT  FS202  S1  ion        0    FE1  0    FE2
CONNECT  FS202  S2  ion        0    FE1  0    FE2

PROTON   FS201      0
PROTON   FS202      2

PKA      FS201      0.0
PKA      FS202      12.3

ELECTRON FS201      0
ELECTRON FS202      1

EM       FS201      0.0
EM       FS202      0.0

RXN      FS201      0.0
RXN      FS202      0.0

#
# Taken from Ullmann et al 2002 J. Biol Inorg. Chemem 7; 632-639
# And modified ... (to many atoms around their sphere ;))
# so increase_charge.py on the rest of the atoms
#
CHARGE   FS201 aCB   0.0946
CHARGE   FS201 aCG  -0.2594
CHARGE   FS201 aND1 -0.2634
CHARGE   FS201 aCD2  0.0016
CHARGE   FS201 aCE1  0.1186
CHARGE   FS201 aNE2 -0.5594
CHARGE   FS201 aHB1 -0.1084
CHARGE   FS201 aHB2 -0.0124
CHARGE   FS201 aHD2  0.0996
CHARGE   FS201 aHE1  0.0766
CHARGE   FS201 bCB  -0.1534
CHARGE   FS201 bCG   0.0566
CHARGE   FS201 bND1 -0.4004
CHARGE   FS201 bCD2 -0.0124
CHARGE   FS201 bCE1  0.0906
CHARGE   FS201 bNE2 -0.5494
CHARGE   FS201 bHB1  0.0096
CHARGE   FS201 bHB2  0.0156
CHARGE   FS201 bHD2  0.0946
CHARGE   FS201 bHE1  0.0876
CHARGE   FS201 cCB   0.0536
CHARGE   FS201 cSG  -0.2874
CHARGE   FS201 cHB1 -0.0514
CHARGE   FS201 cHB2 -0.0444
CHARGE   FS201 dCB   0.0696
CHARGE   FS201 dSG  -0.5174
CHARGE   FS201 dHB1 -0.0344
CHARGE   FS201 dHB2 -0.0174
CHARGE   FS201  FE1  0.5454
CHARGE   FS201  FE2  0.6644
CHARGE   FS201  S1  -0.3577
CHARGE   FS201  S2  -0.4497

CHARGE   FS202 aCB   0.0315
CHARGE   FS202 aCG  -0.0475
CHARGE   FS202 aND1 -0.1525
CHARGE   FS202 aCD2 -0.2815
CHARGE   FS202 aCE1 -0.0605
CHARGE   FS202 aNE2 -0.2735
CHARGE   FS202 aHB1 -0.1115
CHARGE   FS202 aHB2  0.0305
CHARGE   FS202 aHD2  0.2075
CHARGE   FS202 aHE1  0.0975
CHARGE   FS202 aHE2  0.3485
CHARGE   FS202 bCB  -0.1855
CHARGE   FS202 bCG   0.2735
CHARGE   FS202 bND1 -0.1975
CHARGE   FS202 bCD2 -0.4165
CHARGE   FS202 bCE1 -0.2055
CHARGE   FS202 bNE2 -0.1545
CHARGE   FS202 bHB1  0.0295
CHARGE   FS202 bHB2  0.0175
CHARGE   FS202 bHD2  0.2335
CHARGE   FS202 bHE1  0.1875
CHARGE   FS202 bHE2  0.3165
CHARGE   FS202 cCB   0.0725
CHARGE   FS202 cSG  -0.3075
CHARGE   FS202 cHB1 -0.0505
CHARGE   FS202 cHB2 -0.0335
CHARGE   FS202 dCB   0.0845
CHARGE   FS202 dSG  -0.5235
CHARGE   FS202 dHB1 -0.0285
CHARGE   FS202 dHB2 -0.0045
CHARGE   FS202  FE1  0.5795
CHARGE   FS202  FE2  0.3965
CHARGE   FS202  S1  -0.3785
CHARGE   FS202  S2  -0.4935
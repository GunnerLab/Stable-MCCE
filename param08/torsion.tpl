# torsion energy = Vn/2 * [1 + cos(n_fold * torsion_angle - gamma)]
#123456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#                               relax hydroxyl
#-------|-----|----|----|----|----|----|---------|---------|---------|---------|---------|---------|
#        CONF  ATOM ATOM ATOM ATOM relx Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma  
TORSION  NTR   1H    N    CA   C    f      1.400         3      0.00
TORSION  CTR01  HO   O    C    CA   t      4.600         2    180.00
TORSION  CTR02  HXT  OXT  C    CA   t      4.600         2    180.00

TORSION  ASP    CG   CB   CA   N    f      1.400         3      0.00
TORSION  ASP    OD1  CG   CB   CA   f      0.000         2      0.00 
TORSION  ASP01  HD1  OD1  CG   CB   t      4.600         2    180.00
TORSION  ASP02  HD2  OD2  CG   CB   t      4.600         2    180.00 

TORSION  ASN    CG   CB   CA   N    f      1.400         3      0.00
TORSION  ASN    ND2  CG   CB   CA   f      0.100         4      0.00      0.070        2      0.00
TORSION  ASN   1HD2  ND2  CG   CB   f     10.000         2    180.00

TORSION  GLU    CG   CB   CA   N    f      1.400         3      0.00
TORSION  GLU    CD   CG   CB   CA   f      1.400         3      0.00 
TORSION  GLU    OE1  CD   CG   CB   f      0.000         2      0.00 
TORSION  GLU    OE2  CD   CG   CB   f      0.000         2      0.00
TORSION  GLU01  HE1  OE1  CD   CG   t      4.600         2    180.00
TORSION  GLU02  HE2  OE2  CD   CG   t      4.600         2    180.00 

TORSION  GLN    CG   CB   CA   N    f      1.400         3      0.00
TORSION  GLN    CD   CG   CB   CA   f      1.400         3      0.00
TORSION  GLN    NE2  CD   CG   CB   f      0.100         4      0.00      0.070        2      0.00
TORSION  GLN   1HE2  NE2  CD   CG   f     10.000         2    180.00

#-------|-----|----|----|----|----|----|---------|---------|---------|---------|---------|---------|---------|---------|---------|
#        CONF  ATOM ATOM ATOM ATOM relx Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma
TORSION  LYS    CG   CB   CA   N    f      1.400         3      0.00 
TORSION  LYS    CD   CG   CB   CA   f      0.180         3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  LYS    CE   CD   CG   CB   f      0.180         3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  LYS    NZ   CE   CD   CG   f      1.400         3      0.00
TORSION  LYS   1HZ   NZ   CE   CD   f      1.400         3      0.00 

TORSION  ARG    CG   CB   CA   N    f      1.400         3      0.00 
TORSION  ARG    CD   CG   CB   CA   f      0.180         3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  ARG    NE   CD   CG   CB   f      1.400         3      0.00 
TORSION  ARG    CZ   NE   CD   CG   f      0.000         3      0.00 
TORSION  ARG    NH1  CZ   NE   CD   f      9.600         2    180.00 
TORSION  ARG   1HH1  NH1  CZ   NE   f      9.600         2    180.00
TORSION  ARG   1HH2  NH2  CZ   NE   f      9.600         2    180.00 

TORSION  ALA   1HB   CB   CA   N    f      1.400         3      0.00 

TORSION  VAL    CG1  CB   CA   N    f      1.400         3      0.00 
TORSION  VAL   1HG1  CG1  CB   CA   f      0.160         3      0.00
TORSION  VAL   1HG2  CG2  CB   CA   f      0.160         3      0.00

TORSION  LEU    CG   CB   CA   N    f      1.400         3      0.00
TORSION  LEU    CD1  CG   CB   CA   f      0.180         3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  LEU   1HD1  CD1  CG   CB   f      0.160         3      0.00 
TORSION  LEU   1HD2  CD2  CG   CB   f      0.160         3      0.00 

TORSION  ILE    CG1  CB   CA   N    f      1.400         3      0.00
TORSION  ILE    CD1  CG1  CB   CA   f      0.180         3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  ILE   1HD1  CD1  CG1  CB   f      0.160         3      0.00
TORSION  ILE   1HG2  CG2  CB   CA   f      0.160         3      0.00      

TORSION  SER    OG   CB   CA   N    f      1.400         3      0.00 
TORSION  SER    HG   OG   CB   CA   t      0.160         3      0.00      0.250        1      0.00

TORSION  THR    OG1  CB   CA   N    f      1.400         3      0.00
TORSION  THR    HG1  OG1  CB   CA   t      0.160         3      0.00      0.250        1      0.00
TORSION  THR   1HG2  CG2  CB   CA   f      0.160         3      0.00

TORSION  CYS    SG   CB   CA   N    f      1.400         3      0.00
TORSION  CYS01  HG   SG   CB   CA   t      0.750         3      0.00

TORSION  MET    CG   CB   CA   N    f      1.400         3      0.00 
TORSION  MET    SD   CG   CB   CA   f      1.400         3      0.00
TORSION  MET    CE   SD   CG   CB   f      1.000         3      0.00
TORSION  MET   1HE   CE   SD   CG   f      1.000         3      0.00

TORSION  HIS    CG   CB   CA   N    f      1.400         3      0.00 
TORSION  HIS    CE1  ND1  CG   CB   f      1.400         3      0.00 

TORSION  TYR    HH   OH   CZ   CE1  t      1.800         2    180.00
TORSION  TYF    HH   OH   CZ   CE1  t      1.800         2    180.00

TORSION  RSB   1H16  C16  C1   C2   f      0.160         3      0.00
TORSION  RSB   1H17  C17  C1   C2   f      0.160         3      0.00
TORSION  RSB   1H18  C18  C5   C4   f      0.160         3      0.00
TORSION  RSB   1H19  C19  C9   C8   f      0.160         3      0.00
TORSION  RSB   1H20  C20  C13  C12  f      0.160         3      0.00

TORSION  UbQ    HO2  O2   C2   C1   t      1.800         2    180.00
TORSION  UbQ    HO5  O5   C5   C4   t      1.800         2    180.00

TORSION  HEA   1HMB  CMB  C2B  C1B  f      1.400         3      0.00   # like 1HB in VAL
TORSION  HEA   1HMC  CMC  C2C  C1C  f      1.400         3      0.00
TORSION  HEA   1HMD  CMD  C2D  C1D  f      1.400         3      0.00
TORSION  HEA   1HBC  CBC  CAC  C3C  f      1.400         3      0.00
TORSION  HA3   1HMB  CMB  C2B  C1B  f      1.400         3      0.00
TORSION  HA3   1HMC  CMC  C2C  C1C  f      1.400         3      0.00
TORSION  HA3   1HMD  CMD  C2D  C1D  f      1.400         3      0.00
TORSION  HA3   1HBC  CBC  CAC  C3C  f      1.400         3      0.00

TORSION  PAA01  H1A  O1A  CGA  CBA  t      4.600         2    180.00   # like asp hydroxyl
TORSION  PAA02  H2A  O2A  CGA  CBA  t      4.600         2    180.00
TORSION  PDD01  H1D  O1D  CGD  CBD  t      4.600         2    180.00
TORSION  PDD02  H2D  O2D  CGD  CBD  t      4.600         2    180.00

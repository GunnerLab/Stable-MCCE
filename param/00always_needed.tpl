### This file has AMBER L-J interaction parameter, covalent radii parameter,
### torsion parameter, and toggles for residues to be moved in relaxation subroutine -Yifan

# VDW C12 and C6 parameters (AMBER AutoDock V3.0, Self-consistent Lennard-Jones 12-6 parameters)
# VDW = C12/r^12 - C12/r^6
#
# Comments by Junjun:
# VDW reaches minimum when r = reqm
# VDW reaches the second 0 when r = (R1+R2), where R1 and R2 are atomic radii.
# reqm = 1.122*(R1 + R2)

#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H1$
VDWAMBER C12   C-C  2516582.400     # reqm = 4.00, eij = -0.150 kcal/mol
VDWAMBER C6    C-C  1228.800000
VDWAMBER C12   C-N  1198066.249     # reqm = 3.75, eij = -0.155 kcal/mol
VDWAMBER C6    C-N   861.634784
VDWAMBER C12   C-O   820711.722     # reqm = 3.60, eij = -0.173 kcal/mol
VDWAMBER C6    C-O   754.059521
VDWAMBER C12   C-S  2905899.052     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    C-S  1418.896022
VDWAMBER C12   C-H    29108.222     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    C-H    79.857949

VDWAMBER C12   N-C  1198066.249     # reqm = 3.75, eij = -0.155 kcal/mol
VDWAMBER C6    N-C   861.634784
VDWAMBER C12   N-N   540675.281     # reqm = 3.50, eij = -0.160 kcal/mol
VDWAMBER C6    N-N   588.245000
VDWAMBER C12   N-O   357365.541     # reqm = 3.35, eij = -0.179 kcal/mol
VDWAMBER C6    N-O   505.677729
VDWAMBER C12   N-S  1383407.742     # reqm = 3.75, eij = -0.179 kcal/mol
VDWAMBER C6    N-S   994.930149
VDWAMBER C12   N-H    10581.989     # reqm = 2.75, eij = -0.057 kcal/mol
VDWAMBER C6    N-H    48.932922

VDWAMBER C12   O-C   820711.722     # reqm = 3.60, eij = -0.173 kcal/mol
VDWAMBER C6    O-C   754.059521
VDWAMBER C12   O-N   357365.541     # reqm = 3.35, eij = -0.179 kcal/mol
VDWAMBER C6    O-N   505.677729
VDWAMBER C12   O-O   230584.301     # reqm = 3.20, eij = -0.200 kcal/mol
VDWAMBER C6    O-O   429.496730
VDWAMBER C12   O-S   947676.268     # reqm = 3.60, eij = -0.200 kcal/mol
VDWAMBER C6    O-S   870.712934
VDWAMBER C12   O-H     6035.457     # reqm = 2.60, eij = -0.063 kcal/mol
VDWAMBER C6    O-H    39.075098

VDWAMBER C12   S-C  2905899.052     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    S-C  1418.896022
VDWAMBER C12   S-N  1383407.742     # reqm = 3.75, eij = -0.179 kcal/mol
VDWAMBER C6    S-N   994.930149
VDWAMBER C12   S-O   947676.268     # reqm = 3.60, eij = -0.200 kcal/mol
VDWAMBER C6    S-O   870.712934
VDWAMBER C12   S-S  3355443.200     # reqm = 4.00, eij = -0.200 kcal/mol
VDWAMBER C6    S-S  1638.400000
VDWAMBER C12   S-H    33611.280     # reqm = 4.00, eij = -0.200 kcal/mol
VDWAMBER C6    S-H    92.212017

VDWAMBER C12   H-C    29108.222     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    H-C    79.857949
VDWAMBER C12   H-N    10581.989     # reqm = 2.75, eij = -0.057 kcal/mol
VDWAMBER C6    H-N    48.932922
VDWAMBER C12   H-O     6035.457     # reqm = 2.60, eij = -0.063 kcal/mol
VDWAMBER C6    H-O    39.075098
VDWAMBER C12   H-S    33611.280     # reqm = 4.00, eij = -0.200 kcal/mol
VDWAMBER C6    H-S    92.212017
VDWAMBER C12   H-H       81.920     # reqm = 2.00, eij = -0.020 kcal/mol
VDWAMBER C6    H-H     2.560000

# default is the same as S-C
VDWAMBER C12   X-X  2905899.052     # reqm = 4.00, eij = -0.173 kcal/mol
VDWAMBER C6    X-X  1418.896022

#Covalant radii, used for getting bond length information, currently all bonds are assumed
#to be single bond - Yifan, 05/09/2006
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
RADCOVAL SINGL  H   0.30

RADCOVAL SINGL  C   0.772
RADCOVAL DOUBL  C   0.667
RADCOVAL TRIPL  C   0.603

RADCOVAL SINGL  N   0.70

RADCOVAL SINGL  O   0.66

RADCOVAL SINGL  F   0.64

RADCOVAL SINGL  S   1.04

RADCOVAL SINGL  E   1.29   #Fe-O radius 1.95, Yifan
RADCOVAL SINGL  U   1.29   #Cu-O radius 1.95, Yifan

#Bond angles
#123456789012345678901234567890123456789012345678901234567890
#                  |one block of def   |
#-------|-----|----|----|----|---------|----|----|---------|
BOND_ANG HIS    CG   CB   ND1  126.0     CB   CD2  126.0   ND1  CD2  108.0
BOND_ANG HIS    ND1  CG   CE1  108.0
BOND_ANG HIS    CD2  CG   NE2  108.0
BOND_ANG HIS    CE1  ND1  NE2  108.0
BOND_ANG HIS    NE2  CD2  CE1  108.0

BOND_ANG TRP    CG   CB   CD1  126.0     CB   CD2  126.0   CD1  CD2  108.0
BOND_ANG TRP    CD1  CG   NE1  108.0
BOND_ANG TRP    CD2  CG   CE2  108.0     CE2  CE3  120.0
BOND_ANG TRP    NE1  CD1  CE2  108.0
BOND_ANG TRP    CE2  CD2  NE1  108.0     CD2  CZ2  120.0

### list of residues that can be moved in the heavy atom relaxation subroutine
RELAX    ALA        t
RELAX    ARG        t
RELAX    ASN        t
RELAX    ASP        t
RELAX    CYS        t
RELAX    GLN        t
RELAX    GLU        t
RELAX    HIS        t
RELAX    ILE        t
RELAX    LEU        t
RELAX    LYS        t
RELAX    MET        t
RELAX    PAA        t
RELAX    PDD        t
RELAX    PHE        t
RELAX    SER        t
RELAX    THR        t
RELAX    TRP        t
RELAX    TYR        t
RELAX    VAL        t
RELAX    HOH        t
RELAX    _CL        t

### parameters for calculate torsion energy
# torsion energy = Vn/2 * [1 + cos(n_fold * torsion_angle - gamma)]
#123456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#                               relax hydroxyl
#-------|-----|----|----|----|----|----|---------|---------|---------|---------|---------|---------|---------|---------|---------
#        CONF  ATOM ATOM ATOM ATOM relx Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma
TORSION  NTR    H2   N    CA   C    f      1.400         3      0.00
TORSION  CTR01  HO   O    C    CA   t      1.800         2    180.00     1.360         1      0.00
TORSION  CTR02  HXT  OXT  C    CA   t      1.800         2    180.00     1.360         1      0.00

#-------|-----|----|----|----|----|----|---------|---------|---------|---------|---------|---------|---------|---------|---------|
#        CONF  ATOM ATOM ATOM ATOM relx Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma
TORSION  ASP    CG   CB   CA   N    f      1.400         3      0.00
TORSION  ASP    OD1  CG   CB   CA   f      0.000        12    180.00      0.500        2    180.00
TORSION  ASP01  HD1  OD1  CG   CB   t      1.800         2    180.00     1.360         1      0.00
TORSION  ASP02  HD2  OD2  CG   CB   t      1.800         2    180.00     1.360         1      0.00


TORSION  ASN    CG   CB   CA   N    f      1.400         3      0.00
TORSION  ASN    OD1  CG   CB   CA   f      0.000        12    180.00      0.100        4      0.00      0.500        2    180.00
TORSION  ASN   HD21  ND2  CG   CB   f     10.000         2    180.00

TORSION  GLU    CG   CB   CA   N    f      1.400         3    360.00
TORSION  GLU    CD   CG   CB   CA   f      1.400         3    360.00
TORSION  GLU    OE1  CD   CG   CB   f      0.100         4    360.00      0.500        2    180.00
TORSION  GLU01  HE1  OE1  CD   CG   t      1.800         2    180.00     1.360         1      0.00
TORSION  GLU02  HE2  OE2  CD   CG   t      1.800         2    180.00     1.360         1      0.00

TORSION  GLN    CG   CB   CA   N    f      1.400         3      0.00
TORSION  GLN    CD   CG   CB   CA   f      0.000        12    180.00      1.400        3      0.00
TORSION  GLN    OE1  CD   CG   CB   f      0.000        12    180.00      0.100        4      0.00      0.500        2    180.00
TORSION  GLN   HE21  NE2  CD   CG   f     10.000         2    180.00

#-------|-----|----|----|----|----|----|---------|---------|---------|---------|---------|---------|---------|---------|---------|
#        CONF  ATOM ATOM ATOM ATOM relx Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma    Vn/2(kcal)  n_fold   gamma
TORSION  LYS    CG   CB   CA   N    f      1.400         3      0.00
TORSION  LYS    CD   CG   CB   CA   f      0.180         3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  LYS    CE   CD   CG   CB   f      0.180         3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  LYS    NZ   CE   CD   CG   f      1.400         3      0.00
TORSION  LYS    HZ1  NZ   CE   CD   f      1.400         3      0.00

TORSION  ARG    CG   CB   CA   N    f      0.000        36    180.00      1.400        3      0.00
TORSION  ARG    CD   CG   CB   CA   f      0.000        12    180.00      0.180        3      0.00     0.250         2    180.00     0.200         1    180.00
TORSION  ARG    NE   CD   CG   CB   f      1.400         3      0.00
TORSION  ARG    CZ   NE   CD   CG   f      0.000        12    180.00
TORSION  ARG    NH1  CZ   NE   CD   f      9.600         2    180.00
TORSION  ARG   HH11  NH1  CZ   NE   f      9.600         2    180.00
TORSION  ARG   HH21  NH2  CZ   NE   f      9.600         2    180.00

TORSION  ALA    HB1  CB   CA   N    f      1.400         3      0.00

TORSION  VAL    CG1  CB   CA   N    f      1.400         3      0.00
TORSION  VAL   HG11  CG1  CB   CA   f      0.160         3      0.00
TORSION  VAL   HG21  CG2  CB   CA   f      0.160         3      0.00

TORSION  LEU    CG   CB   CA   N    f      0.000        36    180.00      1.400        3      0.00
TORSION  LEU    CD1  CG   CB   CA   f      0.000        12    180.00      0.180        3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  LEU   1HD1  CD1  CG   CB   f      0.160         3      0.00
TORSION  LEU   1HD2  CD2  CG   CB   f      0.160         3      0.00

TORSION  ILE    CG1  CB   CA   N    f      1.400         3      0.00
TORSION  ILE    CD1  CG1  CB   CA   f      0.000        12    180.00      0.180        3      0.00      0.250        2    180.00     0.200         1    180.00
TORSION  ILE   HD11  CD1  CG1  CB   f      0.160         3      0.00
TORSION  ILE   HG21  CG2  CB   CA   f      0.160         3      0.00

TORSION  SER    OG   CB   CA   N    f      1.400         3      0.00
TORSION  SER    HG   OG   CB   CA   t      0.160         3      0.00      0.250        1      0.00

TORSION  THR    OG1  CB   CA   N    f      1.400         3      0.00
TORSION  THR    HG1  OG1  CB   CA   t      0.160         3      0.00      0.250        1      0.00
TORSION  THR   HG21  CG2  CB   CA   f      0.160         3      0.00

TORSION  CYS    SG   CB   CA   N    f      1.400         3      0.00
TORSION  CYS01  HG   SG   CB   CA   t      0.750         3      0.00

TORSION  MET    CG   CB   CA   N    f      1.400         3      0.00
TORSION  MET    SD   CG   CB   CA   f      1.400         3      0.00
TORSION  MET    CE   SD   CG   CB   f      1.000         3      0.00
TORSION  MET    HE1  CE   SD   CG   f      1.000         3      0.00

#-------|-----|----|----|----|----|----|---------|---------|---------|---------|---------|---------|---------|---------|---------|
TORSION  HIS    CG   CB   CA   N    f      1.400         3      0.00
TORSION  HIS    ND1  CG   CB   CA   f        0.0        12    180.00   # X-CC-CT-X from amber94
TORSION  HIS    CE1  ND1  CG   CB   f        0.0         1      0.00
TORSION  HIS    NE2  CE1  ND1  CG   f        0.0         1    180.00

TORSION  PHE    CG   CB   CA   N    f      1.400         3      0.00
TORSION  PHE    CD1  CG   CB   CA   f        0.0         6    180.00
TORSION  PHE    CE1  CD1  CG   CB   f        0.0         1      0.00
TORSION  PHE    CE2  CD2  CG   CB   f        0.0         1      0.00
TORSION  PHE    CZ   CE1  CD1  CG   f        0.0         1    180.00

TORSION  TRP    CG   CB   CA   N    f        0.0         6    180.00     1.400         3      0.00
TORSION  TRP    CD1  CG   CB   CA   f        0.0         6    180.00
TORSION  TRP    NE1  CD1  CG   CB   f        0.0         1      0.00
TORSION  TRP    CE2  CD2  CG   CB   f        0.0         1      0.00
TORSION  TRP    CZ3  CE3  CD2  CG   f        0.0         1      0.00
TORSION  TRP    CH2  CZ2  CE2  CD2  f        0.0         1    180.00

TORSION  TYR    CG   CB   CA   N    f        0.0        12    180.00     1.400         3      0.00
TORSION  TYR    CD1  CG   CB   CA   f        0.0        12    180.00
TORSION  TYR    CE1  CD1  CG   CB   f        0.0         1      0.00
TORSION  TYR    CE2  CD2  CG   CB   f        0.0         1      0.00
TORSION  TYR    CZ   CE1  CD1  CG   f        0.0         1    180.00
TORSION  TYR    HH   OH   CZ   CE1  t      1.800         2    180.00
TORSION  TYF    HH   OH   CZ   CE1  t      1.800         2    180.00

TORSION  RSB   1H16  C16  C1   C2   f      0.160         3      0.00
TORSION  RSB   1H17  C17  C1   C2   f      0.160         3      0.00
TORSION  RSB   1H18  C18  C5   C4   f      0.160         3      0.00
TORSION  RSB   1H19  C19  C9   C8   f      0.160         3      0.00
TORSION  RSB   1H20  C20  C13  C12  f      0.160         3      0.00

TORSION  UbQ    HO2  O2   C2   C1   t      1.800         2    180.00
TORSION  UbQ    HO5  O5   C5   C4   t      1.800         2    180.00

TORSION  HEA    HMB  CMB  C2B  C1B  f      1.400         3      0.00   # like 1HB in VAL
TORSION  HEA    HMC  CMC  C2C  C1C  f      1.400         3      0.00
TORSION  HEA    HMD  CMD  C2D  C1D  f      1.400         3      0.00
TORSION  HEA    HBC  CBC  CAC  C3C  f      1.400         3      0.00
TORSION  HA3    HMB  CMB  C2B  C1B  f      1.400         3      0.00
TORSION  HA3    HMC  CMC  C2C  C1C  f      1.400         3      0.00
TORSION  HA3    HMD  CMD  C2D  C1D  f      1.400         3      0.00
TORSION  HA3    HBC  CBC  CAC  C3C  f      1.400         3      0.00

TORSION  PAA01  H1A  O1A  CGA  CBA  t      1.800         2    180.00     1.360         1      0.00
TORSION  PAA02  H2A  O2A  CGA  CBA  t      1.800         2    180.00     1.360         1      0.00
TORSION  PDD01  H1D  O1D  CGD  CBD  t      1.800         2    180.00     1.360         1      0.00
TORSION  PDD02  H2D  O2D  CGD  CBD  t      1.800         2    180.00     1.360         1      0.00

### IPECE information, list of atoms used to calculate scores in ipece subroutine
IPECE_SC ASP    OD1   0.5
IPECE_SC ASP    OD2   0.5
IPECE_SC GLU    OE1   0.5
IPECE_SC GLU    OE2   0.5
IPECE_SC ARG    NH1   0.5
IPECE_SC ARG    NH2   0.5
IPECE_SC LYS    NZ    1.0

#This file defines "floating" cofactors that can be deleted if the solvent
#accessible is over the number defined in the following line (<=256 char)
#0.05     cut off water if % SAS exceeds this number         (H2O_SASCUTOFF)
#-------------------Floating cofactor names start from this column, space separated
FLOATING COFACTOR   HOH SO4 NO3 ZN CA NA CL DOD VO4


#-------|-----|----|----|----|----|----|----|----|
#EBUILD  RES         A list of atoms to copy for rebuilding
REBUILD  ALA         CB
REBUILD  ARG         CB
REBUILD  ASN         CB
REBUILD  ASP         CB
REBUILD  CYS         CB
REBUILD  GLU         CB
REBUILD  GLN         CB
REBUILD  HIS         CB
REBUILD  ILE         CB
REBUILD  LEU         CB
REBUILD  LYS         CB
REBUILD  MET         CB
REBUILD  PHE         CB
REBUILD  SER         CB
REBUILD  THR         CB
REBUILD  TRP         CB
REBUILD  TYR         CB
REBUILD  VAL         CB
REBUILD  PAA         CAA
REBUILD  PDD         CAD

#-------|-----|----|------------------------------
DEF_CONF ASP        ASP-1
DEF_CONF GLU        GLU-1
DEF_CONF PAA        PAA-1
DEF_CONF PDD        PDD-1
DEF_CONF ARG        ARG+1
DEF_CONF LYS        LYS+1
DEF_CONF TYR        TYR01
DEF_CONF HIS        HIS+1

#-------|-----|----|----|----|----|----|----|----|
BOND_LEN TYR    CG   CD1 1.40  CD2 1.40
BOND_LEN TYR    CD1  CE1 1.40
BOND_LEN TYR    CD2  CE2 1.40
BOND_LEN TYR    CE1  CZ  1.40
BOND_LEN TYR    CE2  CZ  1.40

BOND_LEN TRP    CG   CD1 1.34  CD2 1.40
BOND_LEN TRP    CD1  NE1 1.43
BOND_LEN TRP    CD2  CE2 1.33  CE3 1.40
BOND_LEN TRP    NE1  CE2 1.31
BOND_LEN TRP    CE2  CZ2 1.40
BOND_LEN TRP    CE3  CZ3 1.40
BOND_LEN TRP    CZ2  CH2 1.40
BOND_LEN TRP    CZ3  CH2 1.35

BOND_LEN PHE    CG   CD1 1.40  CD2 1.40
BOND_LEN TRP    CD1  CE1 1.40
BOND_LEN TRP    CD2  CE2 1.40
BOND_LEN TRP    CE1  CZ  1.40
BOND_LEN TRP    CE2  CZ  1.40


DONOR    ARG01 HH11  NH1
DONOR    ARG01 HH12  NH1
DONOR    ARG01 HH21  NH2
DONOR    ARG01 HH22  NH2
DONOR    ARG02 HH11  NH1
DONOR    ARG02 HH21  NH2
DONOR    ARG02 HH22  NH2
DONOR    ARG02  HE   NE
DONOR    ARG03 HH11  NH1
DONOR    ARG03 HH12  NH1
DONOR    ARG03 HH21  NH2
DONOR    ARG03  HE   NE
DONOR    ARG+1 HH11  NH1
DONOR    ARG+1 HH12  NH1
DONOR    ARG+1 HH21  NH2
DONOR    ARG+1 HH22  NH2
DONOR    ARG+1  HE   NE
DONOR    ASN01 HD22  ND2
DONOR    ASN01 HD22  ND2
DONOR    FEX+2 1     FE
DONOR    FS4+2 1     FE1
DONOR    GLN01 HE21  NE2
DONOR    GLN01 HE22  NE2
DONOR    GLU01  HE1  OE1
DONOR    GLU02  HE2  OE2
DONOR    HOH01  H1   O
DONOR    HOH01  H2   O
DONOR    PTR01 HO2P  O2P
DONOR    PTR01 HO3P  O3P
DONOR    PTR-1 HO2P  O2P
DONOR    SER01  HG   OG1
DONOR    THR01  HG1  OG1
DONOR    TPO01 HO2P  O2P
DONOR    TPO01 HO3P  O3P
DONOR    TPO-1 HO2P  O2P
DONOR    TRP01  HE1  NE1
DONOR    TYL01  HH   OH
DONOR    TYR01  HH   OH
DONOR    UbQP1  HO2  O2
DONOR    UbQP2  HO5  O5
DONOR    UbQS1  HO2  O2
DONOR    UbQS2  HO5  O5
ACCEPTOR ASN01  OD1  CG
ACCEPTOR GLN01  OE1  CD
ACCEPTOR GLU01  OE2  CD
ACCEPTOR GLU02  OE1  CD
ACCEPTOR GLU-1  OE1  CD
ACCEPTOR GLU-1  OE2  CD
ACCEPTOR HOH01  O    H1
ACCEPTOR PTR-1  O2P  P
ACCEPTOR PTR-2  O2P  P
ACCEPTOR PTR-2  O3P  P
ACCEPTOR TPO-1  O2P  P
ACCEPTOR TPO-2  O2P  P
ACCEPTOR TPO-2  O3P  P
ACCEPTOR TYL-1  OH   CZ
ACCEPTOR TYR-1  OH   CZ

#23456789A123456789B123456789C123456789D123456789
#Edit by Cai-20160701-used in final data calculation.
HDONOR   ASP01       HD1
HDONOR   ASP02       HD2

HDONOR   ARG01      1HH1 2HH1 1HH2 2HH2
HDONOR   ARG02      1HH1 1HH2 2HH2  HE
HDONOR   ARG03      1HH1 2HH1 1HH2  HE
HDONOR   ARG+1      1HH1 2HH1 1HH2 2HH2  HE

HDONOR   GLU01       HE1
HDONOR   GLU02       HE2

HDONOR   HIS01       HE2
HDONOR   HIS02       HD1
HDONOR   HIS+1       HD1  HE2

HDONOR   SER01       HG 
#add acceptor atom for neutral ser
HACCEPT  SER01       OG

HDONOR   THR01       HG1
#add acceptor atom for neutral thr
HACCEPT  THR01       OG1

##add ASN
HDONOR   ASN01      1HD2 2HD2
HACCEPT  ASN01       OD1

##add GLN
HDONOR   GLN01      1HE2 2HE2
HACCEPT  GLN01       OE1

##add CYS
HDONOR   CYS01       HG
HACCEPT  CYS01       SG 
HACCEPT  CYS-1       SG

##add MET
HACCEPT  MET01       SD

HDONOR   TYR01       HH

##add TRP
HDONOR   TRP01       HE1

HDONOR   LYS01      1HZ  2HZ
HDONOR   LYS+1      1HZ  2HZ  3HZ
## add lys neutral donor######
HACCEPT  LYS01       NZ
############################

HDONOR   PAA01       H1A
HDONOR   PAA02       H2A

HDONOR   PDD01       H1D
HDONOR   PDD02       H2D
################################################
# hydrogen bond acceptor #######################
################################################
HACCEPT  ARG01       NE
HACCEPT  ARG02       NH1
HACCEPT  ARG03       NH2

HACCEPT  ASP01       OD1  OD2
HACCEPT  ASP02       OD1  OD2
HACCEPT  ASP-1       OD1  OD2

HACCEPT  GLU01       OE1  OE2
HACCEPT  GLU02       OE1  OE2
HACCEPT  GLU-1       OE1  OE2

HACCEPT  HIS01       ND1
HACCEPT  HIS02       NE2 

#add acceptor atom for neutral TYR 
HACCEPT  TYR01       OH
HACCEPT  TYR-1       OH

HACCEPT  PAA01       O1A  O2A
HACCEPT  PAA02       O1A  O2A
HACCEPT  PAA-1       O1A  O2A

HACCEPT  PDD01       O1D  O2D
HACCEPT  PDD02       O1D  O2D
HACCEPT  PDD-1       O1D  O2D

#for water
HDONOR   HOH01      1H   2H		# water molecule
HDONOR   HOH-1       H               	# water molecule
HDONOR   HOH+1      1H   2H   3H        # water molecule

HACCEPT  HOH01       O			# water molecule
HACCEPT  HOH-1       O			# water molecule
HACCEPT  HOH+1       O			# water molecule

# for HE3
#HDONOR   HE3+O       HD1 
HACCEPT  HE3+O       O    #OMA
HDONOR   HE30W      1H   2H
HACCEPT  HE30W       O    #OMA
HDONOR   HE3+W      1H   2H
HACCEPT  HE3+W       O    #OMA
HDONOR   HE30H       H  
HACCEPT  HE30H       O    #OMA
HDONOR   HE3+H       H
HACCEPT  HE3+H       O    #OMA


# for CUB
HDONOR   CUB1W      1H   2H   #cHD1
HACCEPT  CUB1W       O
HDONOR   CUB2W      1H   2H   #cHD1
HACCEPT  CUB2W       O
HDONOR   CUB2I      1H   2H
HACCEPT  CUB2I       O   #cND1
HDONOR   CUB1H       H   #cHD1
HACCEPT  CUB1H       O    
HDONOR   CUB2H       H   #cHD1
HACCEPT  CUB2H       O

# for TYF
HDONOR   TYF01       HH
HACCEPT  TYF01       OH
HACCEPT  TYF-1       OH
HDONOR   TYF+1       HH
HACCEPT  TYF+1       OH
HACCEPT  TYF02       OH

# for HLI
HDONOR   HLI01       HD1
HACCEPT  HLI-1       ND1

# for FA3
HDONOR   FA301       HO1
HACCEPT  FA3-1       O11
HACCEPT  FA301       O11

# for HEA
HACCEPT  HEA01       OMA
HACCEPT  HEA+1       OMA

# for FAL
HDONOR   FAL01       HO1
HACCEPT  FAL01       O11
HACCEPT  FAL-1       O11

#23456789A123456789B123456789C123456789D123456789
#Edit by Cai-20160701-used in final data calculation.
HDONOR   ASP01       HD1
HDONOR   ASP02       HD2

HDONOR   ARG01      1HH1 2HH1 1HH2 2HH2
HDONOR   ARG02      1HH1 1HH2 2HH2  HE
HDONOR   ARG03      1HH1 2HH1 1HH2  HE
HDONOR   ARG+1      1HH1 2HH1 1HH2 2HH2  HE

HDONOR   GLU01       HE1
HDONOR   GLU02       HE2

HDONOR   HIS01       HE2
HDONOR   HIS02       HD1
HDONOR   HIS+1       HD1  HE2

HDONOR   SER01       HG 
#add acceptor atom for neutral ser
HACCEPT  SER01       OG

HDONOR   THR01       HG1
#add acceptor atom for neutral thr
HACCEPT  THR01       OG1

##add ASN
HDONOR   ASN01      1HD2 2HD2
HACCEPT  ASN01       OD1

##add GLN
HDONOR   GLN01      1HE2 2HE2
HACCEPT  GLN01       OE1

##add CYS
HDONOR   CYS01       HG
HACCEPT  CYS01       SG 
HACCEPT  CYS-1       SG

##add MET
HACCEPT  MET01       SD

HDONOR   TYR01       HH

##add TRP
HDONOR   TRP01       HE1

HDONOR   LYS01      1HZ  2HZ
HDONOR   LYS+1      1HZ  2HZ  3HZ
## add lys neutral donor######
HACCEPT  LYS01       NZ
############################

HDONOR   PAA01       H1A
HDONOR   PAA02       H2A

HDONOR   PDD01       H1D
HDONOR   PDD02       H2D
################################################
# hydrogen bond acceptor #######################
################################################
HACCEPT  ARG01       NE
HACCEPT  ARG02       NH1
HACCEPT  ARG03       NH2

HACCEPT  ASP01       OD1  OD2
HACCEPT  ASP02       OD1  OD2
HACCEPT  ASP-1       OD1  OD2

HACCEPT  GLU01       OE1  OE2
HACCEPT  GLU02       OE1  OE2
HACCEPT  GLU-1       OE1  OE2

HACCEPT  HIS01       ND1
HACCEPT  HIS02       NE2 

#add acceptor atom for neutral TYR 
HACCEPT  TYR01       OH
HACCEPT  TYR-1       OH

HACCEPT  PAA01       O1A  O2A
HACCEPT  PAA02       O1A  O2A
HACCEPT  PAA-1       O1A  O2A

HACCEPT  PDD01       O1D  O2D
HACCEPT  PDD02       O1D  O2D
HACCEPT  PDD-1       O1D  O2D

#for water
HDONOR   HOH01      1H   2H		# water molecule
HDONOR   HOH-1       H               	# water molecule
HDONOR   HOH+1      1H   2H   3H        # water molecule

HACCEPT  HOH01       O			# water molecule
HACCEPT  HOH-1       O			# water molecule
HACCEPT  HOH+1       O			# water molecule

# for HE3
#HDONOR   HE3+O       HD1 
HACCEPT  HE3+O       O    #OMA
HDONOR   HE30W      1H   2H
HACCEPT  HE30W       O    #OMA
HDONOR   HE3+W      1H   2H
HACCEPT  HE3+W       O    #OMA
HDONOR   HE30H       H  
HACCEPT  HE30H       O    #OMA
HDONOR   HE3+H       H
HACCEPT  HE3+H       O    #OMA


# for CUB
HDONOR   CUB1W      1H   2H   #cHD1
HACCEPT  CUB1W       O
HDONOR   CUB2W      1H   2H   #cHD1
HACCEPT  CUB2W       O
HDONOR   CUB2I      1H   2H
HACCEPT  CUB2I       O   #cND1
HDONOR   CUB1H       H   #cHD1
HACCEPT  CUB1H       O    
HDONOR   CUB2H       H   #cHD1
HACCEPT  CUB2H       O

# for TYF
HDONOR   TYF01       HH
HACCEPT  TYF01       OH
HACCEPT  TYF-1       OH
HDONOR   TYF+1       HH
HACCEPT  TYF+1       OH
HACCEPT  TYF02       OH

# for HLI
HDONOR   HLI01       HD1
HACCEPT  HLI-1       ND1

# for FA3
HDONOR   FA301       HO1
HACCEPT  FA3-1       O11
HACCEPT  FA301       O11

# for HEA
HACCEPT  HEA01       OMA
HACCEPT  HEA+1       OMA

# for FAL
HDONOR   FAL01       HO1
HACCEPT  FAL01       O11
HACCEPT  FAL-1       O11

param04/00always_needed.tpl:VDW_RAD  NTR01  CA  1.908
param04/00always_needed.tpl:VDW_RAD  NTR01  HA  1.1
param04/00always_needed.tpl:VDW_RAD  NTR01  N   1.824
param04/00always_needed.tpl:VDW_RAD  NTR01 1H   0.6
param04/00always_needed.tpl:VDW_RAD  NTR01 2H   0.6
param04/00always_needed.tpl:VDW_RAD  NTR+1  CA  1.908
param04/00always_needed.tpl:VDW_RAD  NTR+1  HA  1.1
param04/00always_needed.tpl:VDW_RAD  NTR+1  N   1.824
param04/00always_needed.tpl:VDW_RAD  NTR+1 1H   0.6
param04/00always_needed.tpl:VDW_RAD  NTR+1 2H   0.6
param04/00always_needed.tpl:VDW_RAD  NTR+1 3H   0.6
param04/00always_needed.tpl:VDW_RAD  CTR01  C   1.908
param04/00always_needed.tpl:VDW_RAD  CTR01  O   1.6612
param04/00always_needed.tpl:VDW_RAD  CTR01  OXT 1.6612
param04/00always_needed.tpl:VDW_RAD  CTR02  C   1.908
param04/00always_needed.tpl:VDW_RAD  CTR02  O   1.6612
param04/00always_needed.tpl:VDW_RAD  CTR02  OXT 1.6612
param04/00always_needed.tpl:VDW_RAD  CTR-1  C   1.908
param04/00always_needed.tpl:VDW_RAD  CTR-1  O   1.6612
param04/00always_needed.tpl:VDW_RAD  CTR-1  OXT 1.6612
param04/00always_needed.tpl:VDW_RAD  ALABK  N   1.824
param04/00always_needed.tpl:VDW_RAD  ALABK  H   0.6
param04/00always_needed.tpl:VDW_RAD  ALABK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  ALABK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  ALABK  C   1.908
param04/00always_needed.tpl:VDW_RAD  ALABK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  ALA01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ALA01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ALA01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ALA01 3HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARGBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  ARGBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  ARGBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  ARGBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  ARGBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  ARGBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  ARG01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ARG01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARG01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARG01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  ARG01 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  ARG01 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  ARG01  CD  1.908
param04/00always_needed.tpl:VDW_RAD  ARG01 1HD  1.387
param04/00always_needed.tpl:VDW_RAD  ARG01 2HD  1.387
param04/00always_needed.tpl:VDW_RAD  ARG01  NE  1.824
param04/00always_needed.tpl:VDW_RAD  ARG01  CZ  1.908
param04/00always_needed.tpl:VDW_RAD  ARG01  NH1 1.824
param04/00always_needed.tpl:VDW_RAD  ARG01 1HH1 0.6
param04/00always_needed.tpl:VDW_RAD  ARG01 2HH1 0.6
param04/00always_needed.tpl:VDW_RAD  ARG01  NH2 1.824
param04/00always_needed.tpl:VDW_RAD  ARG01 1HH2 0.6
param04/00always_needed.tpl:VDW_RAD  ARG01 2HH2 0.6
param04/00always_needed.tpl:VDW_RAD  ARG02  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ARG02 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARG02 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARG02  CG  1.908
param04/00always_needed.tpl:VDW_RAD  ARG02 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  ARG02 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  ARG02  CD  1.908
param04/00always_needed.tpl:VDW_RAD  ARG02 1HD  1.387
param04/00always_needed.tpl:VDW_RAD  ARG02 2HD  1.387
param04/00always_needed.tpl:VDW_RAD  ARG02  NE  1.824
param04/00always_needed.tpl:VDW_RAD  ARG02  HE  0.6
param04/00always_needed.tpl:VDW_RAD  ARG02  CZ  1.908
param04/00always_needed.tpl:VDW_RAD  ARG02  NH1 1.824
param04/00always_needed.tpl:VDW_RAD  ARG02 1HH1 0.6
param04/00always_needed.tpl:VDW_RAD  ARG02  NH2 1.824
param04/00always_needed.tpl:VDW_RAD  ARG02 1HH2 0.6
param04/00always_needed.tpl:VDW_RAD  ARG02 2HH2 0.6
param04/00always_needed.tpl:VDW_RAD  ARG03  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ARG03 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARG03 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARG03  CG  1.908
param04/00always_needed.tpl:VDW_RAD  ARG03 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  ARG03 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  ARG03  CD  1.908
param04/00always_needed.tpl:VDW_RAD  ARG03 1HD  1.387
param04/00always_needed.tpl:VDW_RAD  ARG03 2HD  1.387
param04/00always_needed.tpl:VDW_RAD  ARG03  NE  1.824
param04/00always_needed.tpl:VDW_RAD  ARG03  HE  0.6
param04/00always_needed.tpl:VDW_RAD  ARG03  CZ  1.908
param04/00always_needed.tpl:VDW_RAD  ARG03  NH1 1.824
param04/00always_needed.tpl:VDW_RAD  ARG03 1HH1 0.6
param04/00always_needed.tpl:VDW_RAD  ARG03 2HH1 0.6
param04/00always_needed.tpl:VDW_RAD  ARG03  NH2 1.824
param04/00always_needed.tpl:VDW_RAD  ARG03 1HH2 0.6
param04/00always_needed.tpl:VDW_RAD  ARG+1  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ARG+1 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARG+1 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ARG+1  CG  1.908
param04/00always_needed.tpl:VDW_RAD  ARG+1 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  ARG+1 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  ARG+1  CD  1.908
param04/00always_needed.tpl:VDW_RAD  ARG+1 1HD  1.387
param04/00always_needed.tpl:VDW_RAD  ARG+1 2HD  1.387
param04/00always_needed.tpl:VDW_RAD  ARG+1  NE  1.824
param04/00always_needed.tpl:VDW_RAD  ARG+1  HE  0.6
param04/00always_needed.tpl:VDW_RAD  ARG+1  CZ  1.908
param04/00always_needed.tpl:VDW_RAD  ARG+1  NH1 1.824
param04/00always_needed.tpl:VDW_RAD  ARG+1 1HH1 0.6
param04/00always_needed.tpl:VDW_RAD  ARG+1 2HH1 0.6
param04/00always_needed.tpl:VDW_RAD  ARG+1  NH2 1.824
param04/00always_needed.tpl:VDW_RAD  ARG+1 1HH2 0.6
param04/00always_needed.tpl:VDW_RAD  ARG+1 2HH2 0.6
param04/00always_needed.tpl:VDW_RAD  ASNBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  ASNBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  ASNBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  ASNBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  ASNBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  ASNBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  ASN01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ASN01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ASN01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ASN01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  ASN01  OD1 1.6612
param04/00always_needed.tpl:VDW_RAD  ASN01  ND2 1.824
param04/00always_needed.tpl:VDW_RAD  ASN01 1HD2 0.6
param04/00always_needed.tpl:VDW_RAD  ASN01 2HD2 0.6
param04/00always_needed.tpl:VDW_RAD  ASPBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  ASPBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  ASPBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  ASPBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  ASPBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  ASPBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  ASP01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ASP01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ASP01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ASP01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  ASP01  OD1 1.721
param04/00always_needed.tpl:VDW_RAD  ASP01  HD1 0
param04/00always_needed.tpl:VDW_RAD  ASP01  OD2 1.6612
param04/00always_needed.tpl:VDW_RAD  ASP02  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ASP02 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ASP02 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ASP02  CG  1.908
param04/00always_needed.tpl:VDW_RAD  ASP02  OD1 1.6612
param04/00always_needed.tpl:VDW_RAD  ASP02  OD2 1.721
param04/00always_needed.tpl:VDW_RAD  ASP02  HD2 0
param04/00always_needed.tpl:VDW_RAD  ASP-1  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ASP-1 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  ASP-1 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  ASP-1  CG  1.908
param04/00always_needed.tpl:VDW_RAD  ASP-1  OD1 1.6612
param04/00always_needed.tpl:VDW_RAD  ASP-1  OD2 1.6612
param04/00always_needed.tpl:VDW_RAD  CYSBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  CYSBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  CYSBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  CYSBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  CYSBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  CYSBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  CYS01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  CYS01 1HB  1.387
param04/00always_needed.tpl:VDW_RAD  CYS01 2HB  1.387
param04/00always_needed.tpl:VDW_RAD  CYS01  SG  2
param04/00always_needed.tpl:VDW_RAD  CYS01  HG  0.6
param04/00always_needed.tpl:VDW_RAD  CYS-1  CB  1.908
param04/00always_needed.tpl:VDW_RAD  CYS-1 1HB  1.387
param04/00always_needed.tpl:VDW_RAD  CYS-1 2HB  1.387
param04/00always_needed.tpl:VDW_RAD  CYS-1  SG  2
param04/00always_needed.tpl:VDW_RAD  CYDBK  N   1.824                      #copied from CYS for di-sulfide bridge; marilyn 6/5/07
param04/00always_needed.tpl:VDW_RAD  CYDBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  CYDBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  CYDBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  CYDBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  CYDBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  CYD01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  CYD01 1HB  1.387
param04/00always_needed.tpl:VDW_RAD  CYD01 2HB  1.387
param04/00always_needed.tpl:VDW_RAD  CYD01  SG  2
param04/00always_needed.tpl:VDW_RAD  CYD01  HG  0.6
param04/00always_needed.tpl:VDW_RAD  GLNBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  GLNBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  GLNBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  GLNBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  GLNBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  GLNBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  GLN01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  GLN01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  GLN01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  GLN01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  GLN01 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  GLN01 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  GLN01  CD  1.908
param04/00always_needed.tpl:VDW_RAD  GLN01  OE1 1.6612
param04/00always_needed.tpl:VDW_RAD  GLN01  NE2 1.824
param04/00always_needed.tpl:VDW_RAD  GLN01 1HE2 0.6
param04/00always_needed.tpl:VDW_RAD  GLN01 2HE2 0.6
param04/00always_needed.tpl:VDW_RAD  GLUBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  GLUBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  GLUBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  GLUBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  GLUBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  GLUBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  GLU01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  GLU01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  GLU01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  GLU01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  GLU01 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  GLU01 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  GLU01  CD  1.908
param04/00always_needed.tpl:VDW_RAD  GLU01  OE1 1.721
param04/00always_needed.tpl:VDW_RAD  GLU01  HE1 0
param04/00always_needed.tpl:VDW_RAD  GLU01  OE2 1.6612
param04/00always_needed.tpl:VDW_RAD  GLU02  CB  1.908
param04/00always_needed.tpl:VDW_RAD  GLU02 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  GLU02 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  GLU02  CG  1.908
param04/00always_needed.tpl:VDW_RAD  GLU02 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  GLU02 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  GLU02  CD  1.908
param04/00always_needed.tpl:VDW_RAD  GLU02  OE1 1.6612
param04/00always_needed.tpl:VDW_RAD  GLU02  OE2 1.721
param04/00always_needed.tpl:VDW_RAD  GLU02  HE2 0
param04/00always_needed.tpl:VDW_RAD  GLU-1  CB  1.908
param04/00always_needed.tpl:VDW_RAD  GLU-1 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  GLU-1 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  GLU-1  CG  1.908
param04/00always_needed.tpl:VDW_RAD  GLU-1 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  GLU-1 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  GLU-1  CD  1.908
param04/00always_needed.tpl:VDW_RAD  GLU-1  OE1 1.6612
param04/00always_needed.tpl:VDW_RAD  GLU-1  OE2 1.6612
param04/00always_needed.tpl:VDW_RAD  GLYBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  GLYBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  GLYBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  GLYBK 1HA  1.387
param04/00always_needed.tpl:VDW_RAD  GLYBK 2HA  1.387
param04/00always_needed.tpl:VDW_RAD  GLYBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  GLYBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  HISBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  HISBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  HISBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  HISBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  HISBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  HISBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  HIS01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  HIS01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIS01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIS01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  HIS01  ND1 1.824
param04/00always_needed.tpl:VDW_RAD  HIS01  CE1 1.908
param04/00always_needed.tpl:VDW_RAD  HIS01  HE1 1.359
param04/00always_needed.tpl:VDW_RAD  HIS01  NE2 1.824
param04/00always_needed.tpl:VDW_RAD  HIS01  HE2 0.6
param04/00always_needed.tpl:VDW_RAD  HIS01  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  HIS01  HD2 1.409
param04/00always_needed.tpl:VDW_RAD  HIS02  CB  1.908
param04/00always_needed.tpl:VDW_RAD  HIS02 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIS02 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIS02  CG  1.908
param04/00always_needed.tpl:VDW_RAD  HIS02  ND1 1.824
param04/00always_needed.tpl:VDW_RAD  HIS02  HD1 0.6
param04/00always_needed.tpl:VDW_RAD  HIS02  CE1 1.908
param04/00always_needed.tpl:VDW_RAD  HIS02  HE1 1.359
param04/00always_needed.tpl:VDW_RAD  HIS02  NE2 1.824
param04/00always_needed.tpl:VDW_RAD  HIS02  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  HIS02  HD2 1.409
param04/00always_needed.tpl:VDW_RAD  HIS+1  CB  1.908
param04/00always_needed.tpl:VDW_RAD  HIS+1 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIS+1 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIS+1  CG  1.908
param04/00always_needed.tpl:VDW_RAD  HIS+1  ND1 1.824
param04/00always_needed.tpl:VDW_RAD  HIS+1  HD1 0.6
param04/00always_needed.tpl:VDW_RAD  HIS+1  CE1 1.908
param04/00always_needed.tpl:VDW_RAD  HIS+1  HE1 1.359
param04/00always_needed.tpl:VDW_RAD  HIS+1  NE2 1.824
param04/00always_needed.tpl:VDW_RAD  HIS+1  HE2 0.6
param04/00always_needed.tpl:VDW_RAD  HIS+1  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  HIS+1  HD2 1.409
param04/00always_needed.tpl:VDW_RAD  ILEBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  ILEBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  ILEBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  ILEBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  ILEBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  ILEBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  ILE01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  ILE01  HB  1.487
param04/00always_needed.tpl:VDW_RAD  ILE01  CG1 1.908
param04/00always_needed.tpl:VDW_RAD  ILE01 1HG1 1.487
param04/00always_needed.tpl:VDW_RAD  ILE01 2HG1 1.487
param04/00always_needed.tpl:VDW_RAD  ILE01  CG2 1.908
param04/00always_needed.tpl:VDW_RAD  ILE01 1HG2 1.487
param04/00always_needed.tpl:VDW_RAD  ILE01 2HG2 1.487
param04/00always_needed.tpl:VDW_RAD  ILE01 3HG2 1.487
param04/00always_needed.tpl:VDW_RAD  ILE01  CD1 1.908
param04/00always_needed.tpl:VDW_RAD  ILE01 1HD1 1.487
param04/00always_needed.tpl:VDW_RAD  ILE01 2HD1 1.487
param04/00always_needed.tpl:VDW_RAD  ILE01 3HD1 1.487
param04/00always_needed.tpl:VDW_RAD  LEUBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  LEUBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  LEUBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  LEUBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  LEUBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  LEUBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  LEU01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  LEU01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  LEU01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  LEU01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  LEU01  HG  1.487
param04/00always_needed.tpl:VDW_RAD  LEU01  CD1 1.908
param04/00always_needed.tpl:VDW_RAD  LEU01 1HD1 1.487
param04/00always_needed.tpl:VDW_RAD  LEU01 2HD1 1.487
param04/00always_needed.tpl:VDW_RAD  LEU01 3HD1 1.487
param04/00always_needed.tpl:VDW_RAD  LEU01  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  LEU01 1HD2 1.487
param04/00always_needed.tpl:VDW_RAD  LEU01 2HD2 1.487
param04/00always_needed.tpl:VDW_RAD  LEU01 3HD2 1.487
param04/00always_needed.tpl:VDW_RAD  LYSBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  LYSBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  LYSBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  LYSBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  LYSBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  LYSBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  LYS01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  LYS01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  LYS01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  LYS01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  LYS01 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  LYS01 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  LYS01  CD  1.908
param04/00always_needed.tpl:VDW_RAD  LYS01 1HD  1.487
param04/00always_needed.tpl:VDW_RAD  LYS01 2HD  1.487
param04/00always_needed.tpl:VDW_RAD  LYS01  CE  1.908
param04/00always_needed.tpl:VDW_RAD  LYS01 1HE  1.1
param04/00always_needed.tpl:VDW_RAD  LYS01 2HE  1.1
param04/00always_needed.tpl:VDW_RAD  LYS01  NZ  1.824
param04/00always_needed.tpl:VDW_RAD  LYS01 1HZ  0.6
param04/00always_needed.tpl:VDW_RAD  LYS01 2HZ  0.6
param04/00always_needed.tpl:VDW_RAD  LYS+1  CB  1.908
param04/00always_needed.tpl:VDW_RAD  LYS+1 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  LYS+1 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  LYS+1  CG  1.908
param04/00always_needed.tpl:VDW_RAD  LYS+1 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  LYS+1 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  LYS+1  CD  1.908
param04/00always_needed.tpl:VDW_RAD  LYS+1 1HD  1.487
param04/00always_needed.tpl:VDW_RAD  LYS+1 2HD  1.487
param04/00always_needed.tpl:VDW_RAD  LYS+1  CE  1.908
param04/00always_needed.tpl:VDW_RAD  LYS+1 1HE  1.1
param04/00always_needed.tpl:VDW_RAD  LYS+1 2HE  1.1
param04/00always_needed.tpl:VDW_RAD  LYS+1  NZ  1.824
param04/00always_needed.tpl:VDW_RAD  LYS+1 1HZ  0.6
param04/00always_needed.tpl:VDW_RAD  LYS+1 2HZ  0.6
param04/00always_needed.tpl:VDW_RAD  LYS+1 3HZ  0.6
param04/00always_needed.tpl:VDW_RAD  METBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  METBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  METBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  METBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  METBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  METBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  MET01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  MET01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  MET01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  MET01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  MET01 1HG  1.387
param04/00always_needed.tpl:VDW_RAD  MET01 2HG  1.387
param04/00always_needed.tpl:VDW_RAD  MET01  SD  2
param04/00always_needed.tpl:VDW_RAD  MET01  CE  1.908
param04/00always_needed.tpl:VDW_RAD  MET01 1HE  1.387
param04/00always_needed.tpl:VDW_RAD  MET01 2HE  1.387
param04/00always_needed.tpl:VDW_RAD  MET01 3HE  1.387
param04/00always_needed.tpl:VDW_RAD  PHEBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  PHEBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  PHEBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  PHEBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  PHEBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  PHEBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  PHE01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  PHE01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  PHE01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  PHE01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  PHE01  CD1 1.908
param04/00always_needed.tpl:VDW_RAD  PHE01  HD1 1.459
param04/00always_needed.tpl:VDW_RAD  PHE01  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  PHE01  HD2 1.459
param04/00always_needed.tpl:VDW_RAD  PHE01  CE1 1.908
param04/00always_needed.tpl:VDW_RAD  PHE01  HE1 1.459
param04/00always_needed.tpl:VDW_RAD  PHE01  CE2 1.908
param04/00always_needed.tpl:VDW_RAD  PHE01  HE2 1.459
param04/00always_needed.tpl:VDW_RAD  PHE01  CZ  1.908
param04/00always_needed.tpl:VDW_RAD  PHE01  HZ  1.459
param04/00always_needed.tpl:VDW_RAD  PROBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  PROBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  PROBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  PROBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  PROBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  PRO01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  PRO01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  PRO01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  PRO01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  PRO01 1HG  1.487
param04/00always_needed.tpl:VDW_RAD  PRO01 2HG  1.487
param04/00always_needed.tpl:VDW_RAD  PRO01  CD  1.908
param04/00always_needed.tpl:VDW_RAD  PRO01 1HD  1.387
param04/00always_needed.tpl:VDW_RAD  PRO01 2HD  1.387
param04/00always_needed.tpl:VDW_RAD  SERBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  SERBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  SERBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  SERBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  SERBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  SERBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  SER01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  SER01 1HB  1.387
param04/00always_needed.tpl:VDW_RAD  SER01 2HB  1.387
param04/00always_needed.tpl:VDW_RAD  SER01  OG  1.721
param04/00always_needed.tpl:VDW_RAD  SER01  HG  0
param04/00always_needed.tpl:VDW_RAD  THRBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  THRBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  THRBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  THRBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  THRBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  THRBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  THR01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  THR01  HB  1.387
param04/00always_needed.tpl:VDW_RAD  THR01  OG1 1.721
param04/00always_needed.tpl:VDW_RAD  THR01  HG1 0
param04/00always_needed.tpl:VDW_RAD  THR01  CG2 1.908
param04/00always_needed.tpl:VDW_RAD  THR01 1HG2 1.487
param04/00always_needed.tpl:VDW_RAD  THR01 2HG2 1.487
param04/00always_needed.tpl:VDW_RAD  THR01 3HG2 1.487
param04/00always_needed.tpl:VDW_RAD  TRPBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  TRPBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  TRPBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  TRPBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  TRPBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  TRPBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  TRP01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  TRP01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  TRP01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  TRP01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  TRP01  CD1 1.908
param04/00always_needed.tpl:VDW_RAD  TRP01  HD1 1.409
param04/00always_needed.tpl:VDW_RAD  TRP01  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  TRP01  NE1 1.824
param04/00always_needed.tpl:VDW_RAD  TRP01  HE1 0.6
param04/00always_needed.tpl:VDW_RAD  TRP01  CE2 1.908
param04/00always_needed.tpl:VDW_RAD  TRP01  CE3 1.908
param04/00always_needed.tpl:VDW_RAD  TRP01  HE3 1.459
param04/00always_needed.tpl:VDW_RAD  TRP01  CZ2 1.908
param04/00always_needed.tpl:VDW_RAD  TRP01  HZ2 1.459
param04/00always_needed.tpl:VDW_RAD  TRP01  CZ3 1.908
param04/00always_needed.tpl:VDW_RAD  TRP01  HZ3 1.459
param04/00always_needed.tpl:VDW_RAD  TRP01  CH2 1.908
param04/00always_needed.tpl:VDW_RAD  TRP01  HH2 1.459
param04/00always_needed.tpl:VDW_RAD  TYRBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  TYRBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  TYRBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  TYRBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  TYRBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  TYRBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  TYR01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  TYR01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  TYR01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  TYR01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  TYR01  CD1 1.908
param04/00always_needed.tpl:VDW_RAD  TYR01  HD1 1.459
param04/00always_needed.tpl:VDW_RAD  TYR01  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  TYR01  HD2 1.459
param04/00always_needed.tpl:VDW_RAD  TYR01  CE1 1.908
param04/00always_needed.tpl:VDW_RAD  TYR01  HE1 1.459
param04/00always_needed.tpl:VDW_RAD  TYR01  CE2 1.908
param04/00always_needed.tpl:VDW_RAD  TYR01  HE2 1.459
param04/00always_needed.tpl:VDW_RAD  TYR01  CZ  1.908
param04/00always_needed.tpl:VDW_RAD  TYR01  OH  1.721
param04/00always_needed.tpl:VDW_RAD  TYR01  HH  0
param04/00always_needed.tpl:VDW_RAD  TYR-1  CB  1.908
param04/00always_needed.tpl:VDW_RAD  TYR-1 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  TYR-1 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  TYR-1  CG  1.908
param04/00always_needed.tpl:VDW_RAD  TYR-1  CD1 1.908
param04/00always_needed.tpl:VDW_RAD  TYR-1  HD1 1.459
param04/00always_needed.tpl:VDW_RAD  TYR-1  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  TYR-1  HD2 1.459
param04/00always_needed.tpl:VDW_RAD  TYR-1  CE1 1.908
param04/00always_needed.tpl:VDW_RAD  TYR-1  HE1 1.459
param04/00always_needed.tpl:VDW_RAD  TYR-1  CE2 1.908
param04/00always_needed.tpl:VDW_RAD  TYR-1  HE2 1.459
param04/00always_needed.tpl:VDW_RAD  TYR-1  CZ  1.908
param04/00always_needed.tpl:VDW_RAD  TYR-1  OH  1.721
param04/00always_needed.tpl:VDW_RAD  VALBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  VALBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  VALBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  VALBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  VALBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  VALBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  VAL01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  VAL01  HB  1.487
param04/00always_needed.tpl:VDW_RAD  VAL01  CG1 1.908
param04/00always_needed.tpl:VDW_RAD  VAL01 1HG1 1.487
param04/00always_needed.tpl:VDW_RAD  VAL01 2HG1 1.487
param04/00always_needed.tpl:VDW_RAD  VAL01 3HG1 1.487
param04/00always_needed.tpl:VDW_RAD  VAL01  CG2 1.908
param04/00always_needed.tpl:VDW_RAD  VAL01 1HG2 1.487
param04/00always_needed.tpl:VDW_RAD  VAL01 2HG2 1.487
param04/00always_needed.tpl:VDW_RAD  VAL01 3HG2 1.487
param04/00always_needed.tpl:VDW_RAD  HOH01  O   1.7682
param04/00always_needed.tpl:VDW_RAD  HOH01 1H   0
param04/00always_needed.tpl:VDW_RAD  HOH01 2H   0
param04/00always_needed.tpl:VDW_RAD  HOH-1  O   1.7682
param04/00always_needed.tpl:VDW_RAD  HOH-1 1H   0
param04/00always_needed.tpl:VDW_RAD  HOH+1  O   1.7682
param04/00always_needed.tpl:VDW_RAD  HOH+1 1H   0
param04/00always_needed.tpl:VDW_RAD  HOH+1 2H   0
param04/00always_needed.tpl:VDW_RAD  HOH+1 3H   0
param04/00always_needed.tpl:VDW_RAD  RSBBK  N     1.824
param04/00always_needed.tpl:VDW_RAD  RSBBK  CA    1.908
param04/00always_needed.tpl:VDW_RAD  RSBBK  C     1.908
param04/00always_needed.tpl:VDW_RAD  RSBBK  O     1.70
param04/00always_needed.tpl:VDW_RAD  RSB01  CB    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  CG    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  CD    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  CE    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  NZ    1.824
param04/00always_needed.tpl:VDW_RAD  RSB01  C1    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C2    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C3    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C4    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C5    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C6    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C7    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C8    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C9    1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C10   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C11   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C12   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C13   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C14   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C15   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C16   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C17   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C18   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C19   1.908
param04/00always_needed.tpl:VDW_RAD  RSB01  C20   1.908
param04/00always_needed.tpl:VDW_RAD  CTR01  HO    1.1
param04/00always_needed.tpl:VDW_RAD  CTR02  HXT   1.1
param04/00always_needed.tpl:VDW_RAD  RSBBK  H     1.1
param04/00always_needed.tpl:VDW_RAD  RSBBK  HA    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1HB    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2HB    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1HG    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2HG    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1HD    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2HD    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1HE    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2HE    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1H2    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2H2    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1H3    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2H3    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1H4    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2H4    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01  H7    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01  H8    1.1
param04/00always_needed.tpl:VDW_RAD  RSB01  H10   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01  H11   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01  H12   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01  H14   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01  H15   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1H16   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2H16   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 3H16   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1H17   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2H17   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 3H17   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1H18   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2H18   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 3H18   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1H19   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2H19   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 3H19   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 1H20   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 2H20   1.1
param04/00always_needed.tpl:VDW_RAD  RSB01 3H20   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  CB    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1HB    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2HB    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  CG    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1HG    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2HG    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  CD    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1HD    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2HD    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  CE    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1HE    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2HE    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  NZ    1.824
param04/00always_needed.tpl:VDW_RAD  RSB+1  HZ    0.6
param04/00always_needed.tpl:VDW_RAD  RSB+1  C1    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  C2    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1H2    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2H2    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C3    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1H3    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2H3    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C4    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1H4    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2H4    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C5    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  C6    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  C7    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  H7    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C8    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  H8    1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C9    1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  C10   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  H10   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C11   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  H11   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C12   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  H12   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C13   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  C14   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  H14   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C15   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1  H15   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C16   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1H16   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2H16   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 3H16   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C17   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1H17   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2H17   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 3H17   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C18   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1H18   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2H18   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 3H18   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C19   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1H19   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2H19   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 3H19   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1  C20   1.908
param04/00always_needed.tpl:VDW_RAD  RSB+1 1H20   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 2H20   1.1
param04/00always_needed.tpl:VDW_RAD  RSB+1 3H20   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 FE     1.40
param04/00always_needed.tpl:VDW_RAD  HEM01  CHA   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CHB   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CHC   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CHD   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  N A   1.824
param04/00always_needed.tpl:VDW_RAD  HEM01  C1A   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C2A   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C3A   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C4A   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CMA   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  N B   1.824
param04/00always_needed.tpl:VDW_RAD  HEM01  C1B   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C2B   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C3B   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C4B   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CMB   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CAB   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CBB   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  N C   1.824
param04/00always_needed.tpl:VDW_RAD  HEM01  C1C   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C2C   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C3C   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C4C   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CMC   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CAC   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CBC   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  N D   1.824
param04/00always_needed.tpl:VDW_RAD  HEM01  C1D   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C2D   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C3D   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  C4D   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  CMD   1.908
param04/00always_needed.tpl:VDW_RAD  HEM01  HHA   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01  HHB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01  HHC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01  HHD   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 1HMA   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 2HMA   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 3HMA   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 1HMB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 2HMB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 3HMB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01  HAB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 1HBB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 2HBB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 3HBB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 1HMC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 2HMC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 3HMC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01  HAC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 1HBC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 2HBC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 3HBC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 1HMD   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 2HMD   1.1
param04/00always_needed.tpl:VDW_RAD  HEM01 3HMD   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 FE     1.40
param04/00always_needed.tpl:VDW_RAD  HEM+1  CHA   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CHB   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CHC   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CHD   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  N A   1.824
param04/00always_needed.tpl:VDW_RAD  HEM+1  C1A   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C2A   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C3A   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C4A   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CMA   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  N B   1.824
param04/00always_needed.tpl:VDW_RAD  HEM+1  C1B   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C2B   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C3B   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C4B   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CMB   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CAB   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CBB   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  N C   1.824
param04/00always_needed.tpl:VDW_RAD  HEM+1  C1C   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C2C   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C3C   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C4C   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CMC   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CAC   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CBC   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  N D   1.824
param04/00always_needed.tpl:VDW_RAD  HEM+1  C1D   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C2D   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C3D   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  C4D   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  CMD   1.908
param04/00always_needed.tpl:VDW_RAD  HEM+1  HHA   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1  HHB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1  HHC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1  HHD   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 1HMA   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 2HMA   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 3HMA   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 1HMB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 2HMB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 3HMB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1  HAB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 1HBB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 2HBB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 3HBB   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 1HMC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 2HMC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 3HMC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1  HAC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 1HBC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 2HBC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 3HBC   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 1HMD   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 2HMD   1.1
param04/00always_needed.tpl:VDW_RAD  HEM+1 3HMD   1.1
param04/00always_needed.tpl:VDW_RAD  PAA01  CAA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA01 1HAA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA01 2HAA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA01  CBA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA01 1HBA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA01 2HBA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA01  CGA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA01  O1A 1.721
param04/00always_needed.tpl:VDW_RAD  PAA01  H1A 0
param04/00always_needed.tpl:VDW_RAD  PAA01  O2A 1.6612
param04/00always_needed.tpl:VDW_RAD  PAA02  CAA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA02 1HAA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA02 2HAA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA02  CBA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA02 1HBA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA02 2HBA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA02  CGA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA02  O1A 1.6612
param04/00always_needed.tpl:VDW_RAD  PAA02  O2A 1.721
param04/00always_needed.tpl:VDW_RAD  PAA02  H2A 0
param04/00always_needed.tpl:VDW_RAD  PAA-1  CAA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA-1 1HAA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA-1 2HAA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA-1  CBA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA-1 1HBA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA-1 2HBA 1.487
param04/00always_needed.tpl:VDW_RAD  PAA-1  CGA 1.908
param04/00always_needed.tpl:VDW_RAD  PAA-1  O1A 1.6612
param04/00always_needed.tpl:VDW_RAD  PAA-1  O2A 1.6612
param04/00always_needed.tpl:VDW_RAD  PDD01  CAD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD01 1HAD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD01 2HAD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD01  CBD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD01 1HBD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD01 2HBD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD01  CGD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD01  O1D 1.721
param04/00always_needed.tpl:VDW_RAD  PDD01  H1D 0
param04/00always_needed.tpl:VDW_RAD  PDD01  O2D 1.6612
param04/00always_needed.tpl:VDW_RAD  PDD02  CAD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD02 1HAD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD02 2HAD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD02  CBD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD02 1HBD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD02 2HBD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD02  CGD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD02  O1D 1.6612
param04/00always_needed.tpl:VDW_RAD  PDD02  O2D 1.721
param04/00always_needed.tpl:VDW_RAD  PDD02  H2D 0
param04/00always_needed.tpl:VDW_RAD  PDD-1  CAD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD-1 1HAD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD-1 2HAD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD-1  CBD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD-1 1HBD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD-1 2HBD 1.487
param04/00always_needed.tpl:VDW_RAD  PDD-1  CGD 1.908
param04/00always_needed.tpl:VDW_RAD  PDD-1  O1D 1.6612
param04/00always_needed.tpl:VDW_RAD  PDD-1  O2D 1.6612
param04/00always_needed.tpl:VDW_RAD  NTG01  CA  1.908
param04/00always_needed.tpl:VDW_RAD  NTG01 1HA  1.1
param04/00always_needed.tpl:VDW_RAD  NTG01 2HA  1.1
param04/00always_needed.tpl:VDW_RAD  NTG01  N   1.824
param04/00always_needed.tpl:VDW_RAD  NTG01 1H   0.6
param04/00always_needed.tpl:VDW_RAD  NTG01 2H   0.6
param04/00always_needed.tpl:VDW_RAD  NTG+1  CA  1.908
param04/00always_needed.tpl:VDW_RAD  NTG+1 1HA  1.1
param04/00always_needed.tpl:VDW_RAD  NTG+1 2HA  1.1
param04/00always_needed.tpl:VDW_RAD  NTG+1  N   1.824
param04/00always_needed.tpl:VDW_RAD  NTG+1 1H   0.6
param04/00always_needed.tpl:VDW_RAD  NTG+1 2H   0.6
param04/00always_needed.tpl:VDW_RAD  NTG+1 3H   0.6
param04/00always_needed.tpl:VDW_RAD  HILBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  HILBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  HILBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  HILBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  HILBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  HILBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  HIL01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  HIL01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIL01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIL01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  HIL01  ND1 1.824
param04/00always_needed.tpl:VDW_RAD  HIL01  CE1 1.908
param04/00always_needed.tpl:VDW_RAD  HIL01  HE1 1.359
param04/00always_needed.tpl:VDW_RAD  HIL01  NE2 1.824
param04/00always_needed.tpl:VDW_RAD  HIL01  HE2 0.6
param04/00always_needed.tpl:VDW_RAD  HIL01  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  HIL01  HD2 1.409
param04/00always_needed.tpl:VDW_RAD  HIL02  CB  1.908
param04/00always_needed.tpl:VDW_RAD  HIL02 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIL02 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  HIL02  CG  1.908
param04/00always_needed.tpl:VDW_RAD  HIL02  ND1 1.824
param04/00always_needed.tpl:VDW_RAD  HIL02  HD1 0.6
param04/00always_needed.tpl:VDW_RAD  HIL02  CE1 1.908
param04/00always_needed.tpl:VDW_RAD  HIL02  HE1 1.359
param04/00always_needed.tpl:VDW_RAD  HIL02  NE2 1.824
param04/00always_needed.tpl:VDW_RAD  HIL02  CD2 1.908
param04/00always_needed.tpl:VDW_RAD  HIL02  HD2 1.409
param04/00always_needed.tpl:VDW_RAD  MELBK  N   1.824
param04/00always_needed.tpl:VDW_RAD  MELBK  H   0.6
param04/00always_needed.tpl:VDW_RAD  MELBK  CA  1.908
param04/00always_needed.tpl:VDW_RAD  MELBK  HA  1.387
param04/00always_needed.tpl:VDW_RAD  MELBK  C   1.908
param04/00always_needed.tpl:VDW_RAD  MELBK  O   1.6612
param04/00always_needed.tpl:VDW_RAD  MEL01  CB  1.908
param04/00always_needed.tpl:VDW_RAD  MEL01 1HB  1.487
param04/00always_needed.tpl:VDW_RAD  MEL01 2HB  1.487
param04/00always_needed.tpl:VDW_RAD  MEL01  CG  1.908
param04/00always_needed.tpl:VDW_RAD  MEL01 1HG  1.387
param04/00always_needed.tpl:VDW_RAD  MEL01 2HG  1.387
param04/00always_needed.tpl:VDW_RAD  MEL01  SD  2
param04/00always_needed.tpl:VDW_RAD  MEL01  CE  1.908
param04/00always_needed.tpl:VDW_RAD  MEL01 1HE  1.387
param04/00always_needed.tpl:VDW_RAD  MEL01 2HE  1.387
param04/00always_needed.tpl:VDW_RAD  MEL01 3HE  1.387
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 0C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 1C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 2C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 3C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 4C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 5C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 6C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 7C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 8C99   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C00   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C01   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C02   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C03   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C04   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C05   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C06   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C07   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C08   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C09   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C10   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C11   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C12   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C13   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C14   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C15   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C16   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C17   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C18   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C19   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C20   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C21   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C22   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C23   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C24   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C25   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C26   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C27   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C28   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C29   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C30   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C31   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C32   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C33   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C34   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C35   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C36   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C37   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C38   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C39   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C40   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C41   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C42   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C43   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C44   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C45   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C46   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C47   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C48   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C49   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C50   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C51   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C52   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C53   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C54   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C55   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C56   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C57   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C58   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C59   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C60   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C61   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C62   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C63   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C64   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C65   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C66   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C67   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C68   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C69   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C70   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C71   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C72   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C73   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C74   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C75   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C76   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C77   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C78   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C79   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C80   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C81   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C82   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C83   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C84   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C85   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C86   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C87   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C88   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C89   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C90   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C91   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C92   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C93   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C94   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C95   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C96   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C97   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C98   1.91
param04/00always_needed.tpl:VDW_RAD  MEMBK 9C99   1.91
param04/00always_needed.tpl:#VDW_RAD  _CL-1 CL   2.21
param04/00always_needed.tpl:VDW_RAD  _CL-1 CL   2.25 

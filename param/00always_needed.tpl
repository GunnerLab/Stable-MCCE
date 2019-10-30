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

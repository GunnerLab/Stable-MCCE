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

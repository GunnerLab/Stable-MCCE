#include <stdio.h>
#include <math.h>
#include "mcce.h"

/* This fuction returns a floating point number to account for the electrostatic 
 * interaction of a hydrogen bond, without actual hydrogen.
 * This function is used to overcome the lack of electrostatic term in repacking
 * The optimal distance for O-O, or O-N is tuned to be 3.0 with vdw.
 */

/*
d   vdwO-N vdwO-O Ehbd=10/r minO-N
2.5   3.92   2.11  -4.00   -0.08 	-1.89
2.6   2.11   1.03  -3.85   -1.74 	-2.82
2.7   1.08   0.43  -3.70   -2.63 	-3.28
2.8   0.49   0.10  -3.57   -3.08 	-3.47
2.9   0.16  -0.07  -3.45   -3.29 	-3.52*
3.0  -0.02  -0.16  -3.33   -3.35* 	-3.49
3.1  -0.12  -0.19  -3.23   -3.34 	-3.42
3.2  -0.16  -0.20  -3.13   -3.29 	-3.33
3.3  -0.18  -0.19  -3.03   -3.21 	-3.22
3.4  -0.18  -0.18  -2.94   -3.12 	-3.12
3.5  -0.17  -0.17  -2.86   -3.03 	-3.02
3.6  -0.16  -0.15  -2.78   -2.93 	-2.93
3.7  -0.14  -0.13  -2.70   -2.85 	-2.84
3.8  -0.13  -0.12  -2.63   -2.76 	-2.75
3.9  -0.11  -0.10  -2.56   -2.68 	-2.67
4.0  -0.10  -0.09  -2.50   -2.60 	-2.59
*/

float hbond_extra(CONF a, CONF b)
{  int iatom, jatom;
   float E = 0.0;
   float d;

   for (iatom=0; iatom<a.n_atom; iatom++) {
      if (!a.atom[iatom].on) continue;
      if (a.atom[iatom].name[1] != 'N' && a.atom[iatom].name[1] != 'O') continue;
      for (jatom=0; jatom<b.n_atom; jatom++) {
         if (!b.atom[jatom].on) continue;
         if (b.atom[jatom].name[1] != 'N' && b.atom[jatom].name[1] != 'O') continue;
         d = dvv(a.atom[iatom].xyz, b.atom[jatom].xyz);
         if (d > 0.001) E += -10.0/d;
      }
   }

   return E;
}


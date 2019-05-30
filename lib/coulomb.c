#include <stdio.h>
#include <math.h>
#include "mcce.h"

#define CUTOFF_NEAR  1

float coulomb(ATOM atom1, ATOM atom2)
{  float e = 0.;
   float d;
   VECTOR v1, v2;

   if (!atom1.on) return e;
   if (!atom2.on) return e;
   if (atom1.crg < 1e-6 && atom1.crg > -1e-6) return e;
   if (atom2.crg < 1e-6 && atom2.crg > -1e-6) return e;

   v1 = atom1.xyz;
   v2 = atom2.xyz;

   d = dvv(v1, v2);

   if (d < 0.8) d = 0.8;
   e = 331.5*atom1.crg*atom2.crg/(env.epsilon_coulomb * d);
   return e;
}

VECTOR coulomb_frc(VECTOR v1, VECTOR v2, float crg1, float crg2) {
    VECTOR frc;
    float d, d2;
    float cutoff_near2 = CUTOFF_NEAR * CUTOFF_NEAR;
    
    if (crg1 < 1e-6 && crg1 > -1e-6) {frc.x=0;frc.y=0;frc.z=0;return frc;}
    if (crg2 < 1e-6 && crg2 > -1e-6) {frc.x=0;frc.y=0;frc.z=0;return frc;}

    d2 = ddvv(v1, v2);
    if (d2 < cutoff_near2) d2 = cutoff_near2;
    d = sqrt(d2);
    frc = vector_rescale(vector_vminusv(v2, v1), -331.5*crg1*crg2/(env.epsilon_coulomb * d * d2));
    return frc;
}


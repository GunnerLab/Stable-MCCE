#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

void geom_roll(GEOM *op, float phi, LINE axis)
{  /*--- "Rotate about an arbitrary axis"  CRC Math handbook 4th edition */
   VECTOR v;
   float t, SIN, COS;
   int i, j;
   float m1[4][4], m2[4][4], m3[4][4];

   /* extreme contidions */
   t = fabs(axis.t.x) + fabs(axis.t.y) + fabs(axis.t.z);
   if ((fabs(phi) < 1.0E-8) || t <1.0E-8) return;


   /* translate to (0,0,0) */
   m1[0][0] = 1.0; m1[0][1] = 0.0; m1[0][2] = 0.0; m1[0][3] = -axis.p0.x;
   m1[1][0] = 0.0; m1[1][1] = 1.0; m1[1][2] = 0.0; m1[1][3] = -axis.p0.y;
   m1[2][0] = 0.0; m1[2][1] = 0.0; m1[2][2] = 1.0; m1[2][3] = -axis.p0.z;
   m1[3][0] = 0.0; m1[3][1] = 0.0; m1[3][2] = 0.0; m1[3][3] = 1.0;
   mxm4(m1, op->M, m2);

   /* rotate */
   v = axis.t;

   SIN = sin(phi); COS = cos(phi);

   m1[0][0] = v.x*v.x*(1.0-COS) + COS;
   m1[0][1] = v.x*v.y*(1.0-COS) - v.z*SIN;
   m1[0][2] = v.x*v.z*(1.0-COS) + v.y*SIN;

   m1[1][0] = v.x*v.y*(1.0-COS) + v.z*SIN;
   m1[1][1] = v.y*v.y*(1.0-COS) + COS;
   m1[1][2] = v.y*v.z*(1.0-COS) - v.x*SIN;

   m1[2][0] = v.x*v.z*(1.0-COS) - v.y*SIN;
   m1[2][1] = v.y*v.z*(1.0-COS) + v.x*SIN;
   m1[2][2] = v.z*v.z*(1.0-COS) + COS;

   m1[0][3] = 0.0; m1[1][3] = 0.0; m1[2][3] = 0.0; m1[3][3] = 1.0;
   m1[3][0] = 0.0; m1[3][1] = 0.0; m1[3][2] = 0.0;

   mxm4(m1, m2, m3);

   /* translate back */
   m1[0][0] = 1.0; m1[0][1] = 0.0; m1[0][2] = 0.0; m1[0][3] = axis.p0.x;
   m1[1][0] = 0.0; m1[1][1] = 1.0; m1[1][2] = 0.0; m1[1][3] = axis.p0.y;
   m1[2][0] = 0.0; m1[2][1] = 0.0; m1[2][2] = 1.0; m1[2][3] = axis.p0.z;
   m1[3][0] = 0.0; m1[3][1] = 0.0; m1[3][2] = 0.0; m1[3][3] = 1.0;
   mxm4(m1, m3, m2);

   /* load the result in m2 */
   for (i=0; i<4; i++)
      for (j=0; j<4; j++)
         op->M[i][j] = m2[i][j];

   return;
}

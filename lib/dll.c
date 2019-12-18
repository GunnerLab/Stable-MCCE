/*******************************************************************************
 * NAME
 *        dll - distance between two lines
 *
 * SYNOPSIS
 *        #include <math.h>
 *        #include <mcce.h>
 *
 *        float dll(LINE line1, LINE line2);
 *
 * DESCRIPTION
 *        The dll() function returns the distance between two lines.
 *
 * RETURN VALUE
 *        The dll() function returns a floating point number.
 *
 * SEE ALSO
 *        dpp
 *
 * DEPENDENCY
 *        det3
 *
 * EXAMPLE
 *        #include <math.h>
 *        #include <mcce.h>
 *        ...
 *        LINE line1, line2;
 *        float d;
 *        ...
 *        d=dll(line1, line2);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/

#include <math.h>
#include "mcce.h"

float dll(LINE line1, LINE line2)
/* distance between two lines */
{  float M[3][3];
   float t1, t2, t3, t4;

   M[0][0] = line2.p0.x - line1.p0.x;
   M[0][1] = line2.p0.y - line1.p0.y;
   M[0][2] = line2.p0.z - line1.p0.z;
   M[1][0] = line1.t.x;
   M[1][1] = line1.t.y;
   M[1][2] = line1.t.z;
   M[2][0] = line2.t.x;
   M[2][1] = line2.t.y;
   M[2][2] = line2.t.z;

   t1=det3(M);
   t2=line1.t.y * line2.t.z - line2.t.y * line1.t.z;
   t3=line1.t.z * line2.t.x - line2.t.z * line1.t.x;
   t4=line1.t.x * line2.t.y - line2.t.x * line1.t.y;

   return fabs(t1/sqrt(t2*t2+t3*t3+t4*t4));
}

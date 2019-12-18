/*******************************************************************************
 * NAME
 *        avv - angle between two vectors
 *
 * SYNOPSIS
 *        #include <math.h>
 *        #include <mcce.h>
 *
 *        double avv(VECTOR v1, VECTOR v2);
 *
 * DESCRIPTION
 *        The avv() function returns  the angle  in radians  between two vectors.
 *
 * RETURN VALUE
 *        The returned floating point number is between 0 and pi.
 *
 * SEE ALSO
 *        vdotv
 *
 * DEPENDENCY
 *        vector_normalize
 *
 * EXAMPLE
 *        #include <math.h>
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2;
 *        float phi;
 *        ...
 *        phi=avv(v1, v2);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/

#include <math.h>
#include "mcce.h"

double avv(VECTOR v1, VECTOR v2)
{  double t;
   v1 = vector_normalize(v1); v2 = vector_normalize(v2);
   t = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
   if (t>1.0) t = 1.0;
   else if (t < -1.0) t = -1.0;
   return acos(t);
}


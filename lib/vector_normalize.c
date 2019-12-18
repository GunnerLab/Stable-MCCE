/*******************************************************************************
 * NAME
 *        vector_normalize - normalize a vector
 *
 * SYNOPSIS
 *        #include <math.h>
 *        #include <mcce.h>
 *
 *        VECTOR vector_normalize(VECTOR v);
 *
 * DESCRIPTION
 *        The vector_normalize() function returns the normalized vector of the
 *        passed in vector, or returns (0,0,0) when error occurs.
 *
 * RETURN VALUE
 *        The returned value is in VECTOR type.
 *
 * EXAMPLE
 *        #include <math.h>
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2;
 *        ...
 *        v2=vector_normalize(v1);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/
#include <math.h>
#include "mcce.h"

VECTOR vector_normalize(VECTOR v)
{  VECTOR vn;
   float d;
   d = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
   if (d < 1.0e-20) {
      vn.x = 0.0;
      vn.y = 0.0;
      vn.z = 0.0;
   }
   else {
      vn.x = v.x/d;
      vn.y = v.y/d;
      vn.z = v.z/d;
   }
   return vn;
}

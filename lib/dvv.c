/*******************************************************************************
 * NAME
 *        dvv - distance between two vectors (points)
 *
 * SYNOPSIS
 *        #include <math.h>
 *        #include <mcce.h>
 *
 *        double dvv(VECTOR v1, VECTOR v2);
 *
 * DESCRIPTION
 *        The dvv() function returns the distance between two vectors (points).
 *        A vector is a structure as:
 *
 *                typedef struct {double x, y, z;} VECTOR;
 *
 * RETURN VALUE
 *        Returned floating point number is the distance of two points.
 *
 * SEE ALSO
 *        ddvv()
 *
 * EXAMPLE
 *        #include <math.h>
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2;
 *        float d;
 *        ...
 *        d=dvv(v1, v2);
 *
 * AUTHOR
 *        Junjun Mao, 05/27/2003
 *******************************************************************************/

#include <math.h>
#include "mcce.h"

double dvv(VECTOR v1, VECTOR v2)
{  double dx, dy, dz;
   dx = v2.x- v1.x;
   dy = v2.y- v1.y;
   dz = v2.z- v1.z;
   return sqrt(dx*dx+dy*dy+dz*dz);
}

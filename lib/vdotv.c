/*******************************************************************************
 * NAME
 *        vdotv - dot product of two vectors
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        double vdotv(VECTOR v1, VECTOR v2);
 *
 * DESCRIPTION
 *        The vdotv() returns the dot product of two vectors.  If the two vectors
 *        are  normalized  vectors,  the  dot product  is the cosine of the angle
 *        between two vectors.
 *
 * RETURN VALUE
 *        The returned floating point number is the dot product of two vectors.
 *
 * SEE ALSO
 *        vector_vxv
 *
 * EXAMPLE
 *        #include <stdio.h>
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2;
 *        float phi_cos;
 *        ...
 *        v1 = vector_normalize(v1);
 *        v2 = vector_normalize(v2);
 *        phi_cos = vdotv(v1, v2);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/

#include "mcce.h"

double vdotv(VECTOR v1, VECTOR v2)
{  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

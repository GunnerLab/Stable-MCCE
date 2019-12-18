/*******************************************************************************
 * NAME
 *        line_3v - make a plane from three points
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        PLANE plane_3v(VECTOR v1, VECTOR v2, VECTOR v3);
 *
 * DESCRIPTION
 *        The plane_3v() function returns a PLANE object passing three points.
 *        The normal direction is defined by the cross product of direction v1
 *        to v2 and the direction of v2 to v3.
 *
 * RETURN VALUE
 *        The plane_3v() function returns a structure in PLANE type.
 *
 * SEE ALSO
 *        line_2v
 *
 * DEPENDENCY
 *       vector_normalize, vector_vxv, vector_vminusv
 *
 * EXAMPLE
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2, v3;
 *        PLANE p;
 *        ...
 *        p=plane_3p(v1, v2, v3);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/

#include "mcce.h"

PLANE plane_3v(VECTOR v1, VECTOR v2, VECTOR v3)
{  PLANE plane;
   plane.p0 = v1;
   plane.t = vector_normalize(vector_vxv(vector_vminusv(v2, v1),vector_vminusv(v3, v2)));
   return plane;
}

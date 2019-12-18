/*******************************************************************************
 * NAME
 *        vector_vxv - cross product of two vectors
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        VECTOR vector_vxv(VECTOR v1, VECTOR v2);
 *
 * DESCRIPTION
 *        The vector_vxv() function returns the resulted cross product of vector
 *        of v1 and vector v2.
 *
 * RETURN VALUE
 *        The returned value is in VECTOR type.
 *
 * SEE ALSO
 *        vector_vplusv, vector_vminusv
 *
 * EXAMPLE
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2, v3;
 *        ...
 *        v3=vector_vxv(v1, v2);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/
#include "mcce.h"

VECTOR vector_vxv(VECTOR v1, VECTOR v2)
{  VECTOR z;
   z.x = v1.y*v2.z - v1.z*v2.y;
   z.y = v1.z*v2.x - v1.x*v2.z;
   z.z = v1.x*v2.y - v1.y*v2.x;
   return z;
}

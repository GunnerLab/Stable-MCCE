/*******************************************************************************
 * NAME
 *        vector_vplusv - vector plus vector
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        VECTOR vector_vplusv(VECTOR v1, VECTOR v2);
 *
 * DESCRIPTION
 *        The vector_vplusv() function returns the resulted vector of v1 plus v2.
 *
 * RETURN VALUE
 *        The returned value is in VECTOR type.
 *
 * SEE ALSO
 *        vector_vminusv, vector_vxv
 *
 * EXAMPLE
 *        #include <math.h>
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2, v3;
 *        ...
 *        v3=vector_vplusv(v1, v2);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/
#include "mcce.h"

VECTOR vector_vplusv(VECTOR v1, VECTOR v2)
{  VECTOR z;
   z.x = v1.x + v2.x;
   z.y = v1.y + v2.y;
   z.z = v1.z + v2.z;
   return z;
}

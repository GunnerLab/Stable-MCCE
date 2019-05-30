/*******************************************************************************
 * NAME
 *        avvvv - dihedral angle calculated from 4 points
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        double avvvv(VECTOR v1, VECTOR v2, VECTOR v3, VECTOR v4);
 *
 * DESCRIPTION
 *        The avvvv() function returns the dihedral angle in radians between two
 *        planes composed by  v1, v2, v3  and  v2, v3, v4.   The order of points
 *        the normal direction of the plane by the righ-hand rule.
 *
 * RETURN VALUE
 *        The returned floating point number is between 0 and pi.
 *
 * SEE ALSO
 *        avv, all
 *
 * DEPENDENCY
 *        app, plane_3v
 *
 * EXAMPLE
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2, v3, v4;
 *        float phi;
 *        ...
 *        phi=avvvv(v1, v2, v3, v4);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/

#include "mcce.h"

double avvvv(VECTOR v1, VECTOR v2, VECTOR v3, VECTOR v4)
{ return app(plane_3v(v1, v2, v3), plane_3v(v2, v3, v4));
}

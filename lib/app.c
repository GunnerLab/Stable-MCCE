/*******************************************************************************
 * NAME
 *        app - angle between two planes
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        float app(PLANE p1, PLANE p2);
 *
 * DESCRIPTION
 *        The app() function returns the angle in radians between two planes.
 *
 * RETURN VALUE
 *        The returned floating point number is between 0 and pi.
 *
 * SEE ALSO
 *        avv, all
 *
 * DEPENDENCY
 *        avv
 *
 * EXAMPLE
 *        #include <mcce.h>
 *        ...
 *        PLANE p1, p2;
 *        float phi;
 *        ...
 *        phi=app(p1, p2);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/
#include "mcce.h"

float app(PLANE plane1, PLANE plane2)
{  return avv(plane1.t, plane2.t);
}

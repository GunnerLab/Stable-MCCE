/*******************************************************************************
 * NAME
 *        all - angle between two lines
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        float all(LINE line1, LINE line2);
 *
 * DESCRIPTION
 *        The all() function returns the angle in radians between two lines.
 *
 * RETURN VALUE
 *        The returned floating point number is between 0 and pi.
 *
 * SEE ALSO
 *        avv, app
 *
 * DEPENDENCY
 *        avv
 *
 * EXAMPLE
 *        #include <mcce.h>
 *        ...
 *        LINE line1, line2;
 *        float phi;
 *        ...
 *        phi=all(line1, line2);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/
#include "mcce.h"

float all(LINE line1, LINE line2)
{  return avv(line1.t, line2.t);
}

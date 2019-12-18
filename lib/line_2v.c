/*******************************************************************************
 * NAME
 *        line_2v - make a line from two vectors
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        LINE line_2v(VECTOR v1, VECTOR v2);
 *
 * DESCRIPTION
 *        The line_2v() function returns a LINE object passing two points, from
 *        v1 to v2.
 *
 * RETURN VALUE
 *        The line_2v() function returns a structure in LINE type.
 *
 * SEE ALSO
 *        plane_3v
 *
 * DEPENDENCY
 *       vector_normalize
 *
 * EXAMPLE
 *        #include <mcce.h>
 *        ...
 *        VECTOR v1, v2;
 *        LINE line;
 *        ...
 *        line=line_2v(v1, v2);
 *
 * AUTHOR
 *        Junjun Mao, 05/28/2003
 *******************************************************************************/
#include "mcce.h"

LINE line_2v(VECTOR v1, VECTOR v2)
{  LINE line;
   line.p0 = v1;
   line.t = vector_normalize(vector_vminusv(v2, v1));
   return line;
}

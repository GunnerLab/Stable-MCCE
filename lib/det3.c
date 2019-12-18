/*******************************************************************************
 * NAME
 *        det3 - determinant of a 3 by 3 matrix
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        float det3(float m[][3]);
 *
 * DESCRIPTION
 *        The det3() function returns the determinant of a 3 by 3 matrix.
 *
 * SEE ALSO
 *        det4
 *
 * EXAMPLE
 *        float det;
 *        float M[3][3];
 *        ...
 *        det = det3(M);
 *
 * AUTHOR
 *        Junjun Mao, 05/26/2003
 *******************************************************************************/

float det3(float m[][3])
{  return m[0][0]*m[1][1]*m[2][2]
        - m[0][0]*m[1][2]*m[2][1]
        + m[0][1]*m[1][2]*m[2][0]
        - m[0][1]*m[1][0]*m[2][2]
        + m[0][2]*m[1][0]*m[2][1]
        - m[0][2]*m[1][1]*m[2][0];
}

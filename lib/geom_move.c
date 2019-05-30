/*******************************************************************************
 * NAME
 *        geom_move - geometry translation function
 *
 * SYNOPSIS
 *        #include <stdio.h>
 *        #include <string.h>
 *        #include <mcce.h>
 *
 *        void geom_move(GEOM *op, VECTOR v);
 *
 * DESCRIPTION
 *        The geom_move()  function records a translation  defined by vector v to
 *        the geometry operation recorder op so the translation can be applied to
 *        multiple points later by geom_apply.
 *
 * SEE ALSO
 *        geom_reset, geom_dump, geom_roll, geom_inverse, geom_3v_onto_3v,
 *        geom_apply
 *
 * EXAMPLE
 *        #include <stdio.h>
 *        #include <string.h>
 *        #include <mcce.h>
 *
 *        int main()
 *        {  GEOM recorder;
 *           VECTOR v;
 *
 *           v.x = 1.0; v.y = 2.0; v.z = 3.0;
 *
 *           geom_reset(&recorder);
 *           geom_move(&recorder, v);
 *           geom_dump(recorder);
 *
 *           return 0;
 *        }
 *
 * DEPENDENCY
 *        mxm4
 *
 * AUTHOR
 *        Junjun Mao, 05/26/2003
 *******************************************************************************/

#include <stdio.h>
#include <string.h>
#include "mcce.h"

void geom_move(GEOM *op, VECTOR v)
{  int i, j;
   float m1[4][4], m2[4][4];

   /* translation */
   m1[0][0] = 1.0; m1[0][1] = 0.0; m1[0][2] = 0.0; m1[0][3] = v.x;
   m1[1][0] = 0.0; m1[1][1] = 1.0; m1[1][2] = 0.0; m1[1][3] = v.y;
   m1[2][0] = 0.0; m1[2][1] = 0.0; m1[2][2] = 1.0; m1[2][3] = v.z;
   m1[3][0] = 0.0; m1[3][1] = 0.0; m1[3][2] = 0.0; m1[3][3] = 1.0;

   /* apply translation */
   mxm4(m1, op->M, m2);

   /* load the geom operation */
   for (i=0; i<4; i++)
      for (j=0; j<4; j++)
         op->M[i][j] = m2[i][j];

   return;
}

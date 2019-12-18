/*******************************************************************************
 * NAME
 *        geom_apply - apply geometry transformation to a vector
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        void geom_apply(GEOM op, VECTOR *v)
 *
 * DESCRIPTION
 *        The geom_apply() applies the transformation recorded in transformation
 *        recorder op to a vector. To apply the transformation to atoms, use
 *        geom_apply_atoms()
 *
 * SEE ALSO
 *        geom_reset, geom_dump, geom_move, geom_roll, geom_3v_onto_3v,
 *        geom_apply_atoms()
 *
 * EXAMPLE
 *        #include <stdio.h>
 *        #include <mcce.h>
 *
 *        int main()
 *        {  GEOM recorder;
 *           VECTOR v1, v2;
 *
 *           v.x = 1.0; v.y = 1.0; v.z = 1.0;
 *           geom_reset(&recorder);
 *           geom_move(&recorder, v);
 *           geom_apply(&v);
 *
 *           return 0;
 *        }
 *
 * AUTHOR
 *        Junjun Mao, 05/27/2003
 *******************************************************************************/

#include "mcce.h"

void geom_apply(GEOM op, VECTOR *v)
{  float vh[4];

   vh[0] = v->x;
   vh[1] = v->y;
   vh[2] = v->z;
   vh[3] = 1.0;

   v->x = op.M[0][0]*vh[0] + op.M[0][1]*vh[1] + op.M[0][2]*vh[2] + op.M[0][3]*vh[3];
   v->y = op.M[1][0]*vh[0] + op.M[1][1]*vh[1] + op.M[1][2]*vh[2] + op.M[1][3]*vh[3];
   v->z = op.M[2][0]*vh[0] + op.M[2][1]*vh[1] + op.M[2][2]*vh[2] + op.M[2][3]*vh[3];

   return;
}

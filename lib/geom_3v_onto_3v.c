/*******************************************************************************
 * NAME
 *        geom_3v_onto_3v - geometry superimpose function
 *
 * SYNOPSIS
 *        #include <stdio.h>
 *        #include <string.h>
 *        #include <mcce.h>
 *
 *        GEOM geom_3v_onto_3v(VECTOR v1, VECTOR v2, VECTOR v3, VECTOR t1,
 *                             VECTOR t2, VECTOR t3);
 *
 * DESCRIPTION
 *        The geom_3v_onto_3v() function superimposes point v1 to t1, then line
 *        v1-v2 to t1-t2,  then plane  v1-v2-v3  to  t1-t2-t3,  and returns the
 *        transformation in a GEOM object.
 *
 * RETUERN VALUE
 *        The geom_3v_onto_3v() function returns a GEOM object.
 *
 * SEE ALSO
 *        geom_reset, geom_dump, geom_move, geom_roll, geom_inverse, geom_apply,
 *        all, avv, line_3v, plane_3v
 *
 * EXAMPLE
 *        #include <stdio.h>
 *        #include <string.h>
 *        #include <mcce.h>
 *
 *        int main()
 *        {  GEOM recorder;
 *           VECTOR v1, v2, v3, t1, t2, t3;
 *           v1.x = 0.0; v1.y = 0.0; v1.z = 0.0;
 *           v2.x = 1.0; v2.y = 0.0; v2.z = 0.0;
 *           v3.x = 0.0; v3.y = 1.0; v3.z = 0.0;
 * 
 *           t1.x = 0.0; t1.y = 0.0; t1.z = 0.0;
 *           t2.x = 0.0; t2.y = 0.0; t2.z = 1.0;
 *           t3.x = 0.0; t3.y = 1.0; t3.z = 0.0;
 *
 *           recorder = geom_3v_onto_3v(v1, v2, v3, t1, t2, t3);
 *           geom_dump(recorder);
 *
 *           return 0;
 *        }
 *        
 * DEPENDENCY
 *        geom_reset, geom_move, geom_roll, geom_apply, all, avv, line_2v,
 *        plane_3v
 *        
 * AUTHOR
 *        Junjun Mao, 05/27/2003
 *******************************************************************************/
 

#include <stdio.h>
#include <string.h>
#include "mcce.h"

GEOM geom_3v_onto_3v(VECTOR v1, VECTOR v2, VECTOR v3, VECTOR t1, VECTOR t2, VECTOR t3)
{  GEOM op;
   float angle;
   LINE axis;
   VECTOR v123;
   PLANE plane_v123, plane_t123;

   /* step 1, superimpose v1 to t1 */
   geom_reset(&op);
   geom_move(&op, vector_vminusv(t1, v1));

   /* step 2, align v1-v2 to t1-t2 */
   angle = all(line_2v(v1, v2), line_2v(t1, t2));
   axis.p0 = t1;
   axis.t  = vector_normalize(vector_vxv(vector_vminusv(v2, v1), vector_vminusv(t2, t1)));
   geom_roll(&op, angle, axis);

   /* normal direction of plane v123 should be updated */
   geom_apply(op, &v1);
   geom_apply(op, &v2);
   geom_apply(op, &v3);
   plane_v123 = plane_3v(v1, v2, v3);
   plane_t123 = plane_3v(t1, t2, t3);
   v123 = plane_v123.t;

   /*step 3, align v1-v2-v3 to t1-t2-t3 */
   angle = avv(v123, plane_t123.t);
   axis.p0 = t1;
   axis.t  = vector_normalize(vector_vxv(v123, plane_t123.t));
   geom_roll(&op, angle, axis);

   return op;
}


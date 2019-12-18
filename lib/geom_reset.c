/*******************************************************************************
 * NAME
 *        geom_reset - reset geometry operation recorder
 *
 * SYNOPSIS
 *        #include <string.h>
 *        #include <mcce.h>
 *
 *        void geom_reset(GEOM *op)
 *
 * DESCRIPTION
 *        The geom_reset() function resets the transformation matrix and history
 *        records of the geometry  operation recorder op.  It is  recommended to
 *        call this function  for newly defined GEOM object  because c  does not
 *        initialize a non-static varaible.
 *
 * SEE ALSO
 *        geom_dump, geom_move, geom_roll, geom_inverse, geom_3v_onto_3v, 
 *        geom_apply
 *
 * EXAMPLE
 *        #include <string.h>
 *        #include <mcce.h>
 *
 *        int main()
 *        {  GEOM recorder;
 *           geom_reset(&recorder);
 *           return 0;
 *        }
 *
 * AUTHOR
 *        Junjun Mao, 05/26/2003
 *******************************************************************************/

#include "mcce.h"

void geom_reset(GEOM *op)
{  op->M[0][0] = 1.0; op->M[0][1] = 0.0; op->M[0][2] = 0.0; op->M[0][3] = 0.0;
   op->M[1][0] = 0.0; op->M[1][1] = 1.0; op->M[1][2] = 0.0; op->M[1][3] = 0.0;
   op->M[2][0] = 0.0; op->M[2][1] = 0.0; op->M[2][2] = 1.0; op->M[2][3] = 0.0;
   op->M[3][0] = 0.0; op->M[3][1] = 0.0; op->M[3][2] = 0.0; op->M[3][3] = 1.0;
   return;
}

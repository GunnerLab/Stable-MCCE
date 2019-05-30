#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int cpy_res(RES *tgt, RES *src)
{  int i_conf;

   /* faithful copy of the residue */
   memcpy(tgt, src, sizeof(RES));
   
   /* prepare the lower level array */
   tgt->n_conf = 0;
   tgt->conf   = NULL;
   if (src->n_conf) {
       for (i_conf = 0; i_conf<src->n_conf; i_conf++) {
           ins_conf(tgt, i_conf, src->conf[i_conf].n_atom);
           cpy_conf(&tgt->conf[i_conf], &src->conf[i_conf]);
       }
   }
   
   return 0;
}

#include <string.h>
#include <stdlib.h>
#include "mcce.h"

int cpy_prot(PROT *tgt, PROT *src)
{  int i;

   if (tgt->n_res) printf("cpy_prot(): can not copy to a non-empty target .\n");

   /* faithful copy of the residue */
   memcpy(tgt, src, sizeof(PROT));

   /* copy the lower level array */
   if (src->n_res) {
       tgt->res = (RES *) malloc(src->n_res * sizeof(RES));
       for (i = 0; i<src->n_res; i++) cpy_res(&tgt->res[i], &src->res[i]);
   }
   
   return 0;
}

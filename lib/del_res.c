#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int del_res(PROT *prot, int pos)
{  int i;

   if (pos >= prot->n_res) {
      printf("del_res(): off range deletion.\n");
      return USERERR;
   }

   /* free all the lower level data under this group */

   if (prot->res[pos].n_conf) {
      for (i = 0; i < prot->res[pos].n_conf; i++)
         if (prot->res[pos].conf[i].n_atom) free(prot->res[pos].conf[i].atom);
      free(prot->res[pos].conf);
   }

   /* move contents after position "pos" backward by 1 unit */
   memmove(prot->res+pos, prot->res+pos+1, (prot->n_res - pos - 1) * sizeof(RES));

   /* shrink the memory, update number of atoms */
   prot->n_res--;
   prot->res = (RES *) realloc(prot->res, prot->n_res * sizeof(RES));

   return pos;
}

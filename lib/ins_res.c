#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int ins_res(PROT *prot, int ins)
{  if (ins > prot->n_res) {
      printf("ins_res(): off range insertion.\n");
      return USERERR;
   }

   /*--- increase the residue list size by 1 */
   if (!(prot->res = (RES *) realloc(prot->res, (prot->n_res + 1) * sizeof(RES)))) {
      printf("ins_res(): Fails resizing memory.\n");
      return USERERR;
   }
   prot->n_res++;

   /* move contents after position "ins" forward by 1 */
   memmove(prot->res+ins+1, prot->res+ins, (prot->n_res - ins - 1) * sizeof(RES));

   /* initialize the inserted residue */
   memset(&prot->res[ins], 0, sizeof(RES));

   return ins;
}

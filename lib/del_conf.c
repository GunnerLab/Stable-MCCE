#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int del_conf(RES *res_p, int pos)
{  if (pos >= res_p->n_conf) {
      printf("del_conf(): offrange deletion. residue %3s %c%04d%c, n_conf=%d, delete pos=%d\n",res_p->resName,res_p->chainID,res_p->resSeq,res_p->iCode,res_p->n_conf,pos);
      return USERERR;
   }

   /* free all the lower level data under this conformer */
   if (res_p->conf[pos].atom) {
      free(res_p->conf[pos].atom);
   }

   /* move contents after position "pos" backward by 1 unit */
   memmove(res_p->conf+pos, res_p->conf+pos+1, (res_p->n_conf - pos - 1) * sizeof(CONF));

   /* update number of conformers, shrink the memory */
   res_p->n_conf--;
   res_p->conf = (CONF *) realloc(res_p->conf, res_p->n_conf * sizeof(CONF));

   return pos;
}

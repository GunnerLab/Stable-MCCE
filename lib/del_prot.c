#include <string.h>
#include <stdlib.h>
#include "mcce.h"

int del_prot(PROT *prot)
{  int i, j;

   /* free all the lower level data under this prot */
   if (prot->n_res) {
      for (i=0; i<prot->n_res; i++) {
         if (prot->res[i].n_conf) {
            for (j=0; j<prot->res[i].n_conf; j++) {
               if (prot->res[i].conf[j].n_atom>0) free(prot->res[i].conf[j].atom);/* free atoms */
            }
            free(prot->res[i].conf);	/* free conformers */
         }
      }
      free(prot->res);			/* free residues */
   }

   memset(prot, 0, sizeof(PROT));

   return 0;
}

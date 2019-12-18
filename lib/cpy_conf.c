#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int cpy_conf(CONF *tgt, CONF *src)
{  int i;
   ATOM *swap;

   if (src->n_atom != tgt->n_atom) {
      printf("cpy_conf(): conformers do not have same number of atom slots.\n");
      return USERERR;
   }

   /* atom array was already declared, so preserve atoms (number of atoms are the same) */
   swap=tgt->atom;

   /* faithful copy of the conformer */
   memcpy(tgt, src, sizeof(CONF));

   /* copy sublevel contents */
   tgt->atom = swap;
   if (tgt->n_atom && tgt->atom) {
       for (i = 0; i<src->n_atom; i++) tgt->atom[i] = src->atom[i];
   }
   return 0;
}

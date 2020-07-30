#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int atom_element(ATOM *atom)
{   if (atom->name[0] == 'H' && atom->name[3] != ' ') { // "HG  " for mercury, "HG12" for gamma H
        strncpy(atom->element, " H", 2);
        atom->element[2] = '\0';
    }
    else {
        strncpy(atom->element, atom->name, 2);
        atom->element[2] = '\0';
    }
    //printf("name = %2s\n", atom->element);
    return 0;
}

int prot_atom_element(PROT prot)
{  int i, j, k;
   for (i=0; i<prot.n_res; i++) {
      for (j=0; j<prot.res[i].n_conf; j++) {
         for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
            //if (!(prot.res[i].conf[j].atom[k].on)) continue;
            atom_element(&(prot.res[i].conf[j].atom[k]));
            //printf("element = %s\n", prot.res[i].conf[j].atom[k].name);
         }
      }
   }
   return 0;
}
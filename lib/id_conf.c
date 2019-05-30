#include <stdio.h>
#include "mcce.h"

void id_conf(PROT prot)
{  int kr, kc;
   char ins;
   int C;

   for (kr=0; kr<prot.n_res; kr++) {
      C=0;
      if (prot.res[kr].iCode == '\0' || prot.res[kr].iCode == ' ')
         ins = '_';
      else ins = prot.res[kr].iCode;
      for (kc=0; kc<prot.res[kr].n_conf; kc++) {
         sprintf(prot.res[kr].conf[kc].uniqID, "%3s%c%c%c%04d%c%03d", prot.res[kr].resName,
                                                                     prot.res[kr].conf[kc].history[0],
                                                                     prot.res[kr].conf[kc].history[1],
                                                                     prot.res[kr].chainID,
                                                                     prot.res[kr].resSeq,
                                                                     ins, C);
         prot.res[kr].conf[kc].uniqID[sizeof(prot.res[kr].conf[kc].uniqID)-1] = '\0';
         C++;
      }
   }

   return;
}

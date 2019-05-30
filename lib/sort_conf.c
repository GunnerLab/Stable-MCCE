#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int sort_conf(PROT prot)
{   int i_res;
    int i_conf;
    int j_conf;
    CONF swap_conf;
    STRINGS conflist;
    int  i_conflist;

    for (i_res=0; i_res<prot.n_res; i_res++) {
        if (param_get("CONFLIST", prot.res[i_res].resName, "", &conflist)) {
            printf("sort_conf(): Can't find parameter \"CONFLIST\" of residue \"%s\"\n", prot.res[i_res].resName);
            continue;
        }

        i_conflist = 0;
        for (i_conf=0; i_conf < prot.res[i_res].n_conf; i_conf++) {
            if (i_conflist >= conflist.n) {
               printf("   WARNING: sort_conf(): some conformers in residue \"%s %03d %c\" can not be sorted.\n",
                      prot.res[i_res].resName,
                      prot.res[i_res].resSeq,
                      prot.res[i_res].chainID);
                      break;
            }
            if (!strcmp(prot.res[i_res].conf[i_conf].confName, conflist.strings[i_conflist])) continue;
            /* find the first confomer that can be switched and switch */
            for (j_conf=i_conf+1; j_conf<prot.res[i_res].n_conf; j_conf++) {
               if (!strcmp(prot.res[i_res].conf[j_conf].confName, conflist.strings[i_conflist])) {
                  memcpy(&swap_conf, &prot.res[i_res].conf[j_conf], sizeof(CONF));
                  memmove(&prot.res[i_res].conf[i_conf+1], &prot.res[i_res].conf[i_conf], (j_conf-i_conf)*sizeof(CONF)); /* move instead of switch. Yifan 03/04/01 */
                  memcpy(&prot.res[i_res].conf[i_conf], &swap_conf, sizeof(CONF));
                  break;
               }
            }
            if (j_conf>=prot.res[i_res].n_conf) {
               i_conflist ++;
               i_conf--;
            }
        }
    }
    return 0;
}


int cmp_Eself(const void *a, const void *b)
{  CONF *A, *B;

   A = (CONF *) a;
   B = (CONF *) b;
   if (A->E_self<B->E_self) return -1;
   else if (A->E_self==B->E_self) return 0;
   else return 1;
}

/*
int cmp_Eself(CONF *A, CONF *B)
{  if (A->E_self<B->E_self) return -1;
   else if (A->E_self==B->E_self) return 0;
   else return 1;
}
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int sort_res(PROT prot) {
    int i_res,j_res;
    RES swap_res;

    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (j_res=i_res+1; j_res<prot.n_res; j_res++) {
            if (prot.res[i_res].chainID > prot.res[j_res].chainID) {
                  memcpy(&swap_res, &prot.res[j_res], sizeof(RES));
                  memmove(&prot.res[i_res+1], &prot.res[i_res], (j_res-i_res)*sizeof(RES)); /* move instead of switch. Yifan 03/04/01 */
                  memcpy(&prot.res[i_res], &swap_res, sizeof(RES));
                  j_res = i_res;
            }
        }
    }

    /*
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (j_res=i_res+1; j_res<prot.n_res; j_res++) {
            if (prot.res[i_res].chainID != prot.res[j_res].chainID) continue;
            if (prot.res[i_res].groupID > prot.res[j_res].groupID) {
                  memcpy(&swap_res, &prot.res[j_res], sizeof(RES));
                  memmove(&prot.res[i_res+1], &prot.res[i_res], (j_res-i_res)*sizeof(RES));  //move instead of switch. Yifan 03/04/01 
                  memcpy(&prot.res[i_res], &swap_res, sizeof(RES));
                  j_res = i_res;
            }
        }
    }
    */

    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (j_res=i_res+1; j_res<prot.n_res; j_res++) {
            if ((prot.res[i_res].chainID != prot.res[j_res].chainID)) continue;
            if (prot.res[i_res].resSeq > prot.res[j_res].resSeq) {
                  memcpy(&swap_res, &prot.res[j_res], sizeof(RES));
                  memmove(&prot.res[i_res+1], &prot.res[i_res], (j_res-i_res)*sizeof(RES)); /* move instead of switch. Yifan 03/04/01 */
                  memcpy(&prot.res[i_res], &swap_res, sizeof(RES));
                  j_res = i_res;
            }
        }
    }
    return 0;
}

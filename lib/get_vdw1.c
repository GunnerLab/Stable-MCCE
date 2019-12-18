#include <stdio.h>
#include "mcce.h"

int     n_elem;
//float   C6_matrix[N_ELEM_MAX][N_ELEM_MAX];
//float   C12_matrix[N_ELEM_MAX][N_ELEM_MAX];

void get_vdw1(PROT prot)
{
    float e;
    int i, j, k;
    
    assign_rad(prot);
    assign_crg(prot);
    get_connect12(prot);
    //setup_C6_C12(prot);
    setup_vdw_fast(prot);
    
    for (i=0; i<prot.n_res; i++) {
        setup_connect_res(prot, i);
        prot.res[i].conf[0].E_vdw1 = 0.0; /* conformer 0 is defined to have 0 */
        
        for (j=1; j<prot.res[i].n_conf; j++) {
            e = 0.0;
            for (k=0; k<prot.n_res; k++) {
                e += vdw_conf_fast(i, j, k, 0, prot, 0);
                //printf("VDW interaction %s-%s %8.3f\n",prot.res[i].conf[j].uniqID, prot.res[k].conf[0].uniqID, vdw_conf_fast(i, j, k, 0, prot, 0));
            }
            prot.res[i].conf[j].E_vdw1 = e;
        }
        free_connect_res(prot, i);
    }
    return;
}

void get_vdw1_pr_bkb(PROT prot)
{
    float e;
    int i, j, k;
    
    assign_rad(prot);
    assign_crg(prot);
    get_connect12(prot);
    //setup_C6_C12(prot);
    setup_vdw_fast(prot);
    
    for (i=0; i<prot.n_res; i++) {
        setup_connect_res(prot, i);
        prot.res[i].conf[0].E_vdw1 = 0.0; /* conformer 0 is defined to have 0 */
        
        for (j=1; j<prot.res[i].n_conf; j++) {
            e = 0.0;
            for (k=0; k<prot.n_res; k++) {
                e += vdw_conf_fast(i, j, k, 0, prot, 0);
                printf("VDW interaction %s-%s %8.3f\n",
                prot.res[i].conf[j].uniqID, prot.res[k].conf[0].uniqID, vdw_conf_fast(i, j, k, 0, prot, 0));
            }
            prot.res[i].conf[j].E_vdw1 = e;
        }
        free_connect_res(prot, i);
    }
    return;
}

void get_vdw1_res(int i, PROT prot)
{
    float e;
    int j, k;
    
    setup_vdw_fast_res(i, prot);
    setup_connect_res(prot, i);
    prot.res[i].conf[0].E_vdw1 = 0.0; /* conformer 0 is backbone */
    
    for (j=1; j<prot.res[i].n_conf; j++) {
        e = 0.0;
        for (k=0; k<prot.n_res; k++) {
            e += vdw_conf_fast(i, j, k, 0, prot, 0);
        }
        prot.res[i].conf[j].E_vdw1 = e;
    }
    free_connect_res(prot, i);
}

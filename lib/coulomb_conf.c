#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

float coulomb_conf(int ires, int iconf, int jres, int jconf, PROT prot)
{
    float  e = 0.0;
    int    iatom, jatom;
    ATOM   *connect13[MAX_CONNECTED2],    *connect14[MAX_CONNECTED3];
    ATOM   *iatom_p, *jatom_p;
    int    connect13_res[MAX_CONNECTED2], connect14_res[MAX_CONNECTED3];
    int    iconnect, jconnect, kconnect, n_connect13, n_connect14;
    
    if (ires<0 || ires>prot.n_res-1) printf("   Error! coulomb_conf(): residue index out of range in protein\n"),exit(-1);
    if (jres<0 || jres>prot.n_res-1) printf("   Error! coulomb_conf(): residue index out of range in protein\n"),exit(-1);
    if (iconf<0 || iconf>prot.res[ires].n_conf-1) printf("   Error! coulomb_conf(): conformer index out of range in protein\n"),exit(-1);
    if (jconf<0 || jconf>prot.res[jres].n_conf-1) printf("   Error! coulomb_conf(): conformer index out of range in protein\n"),exit(-1);
    for (iatom=0; iatom<prot.res[ires].conf[iconf].n_atom; iatom++) {
        iatom_p = &prot.res[ires].conf[iconf].atom[iatom];
        if (!iatom_p->on) continue;
        
        n_connect13 = 0;
        n_connect14 = 0;
        memset(connect13,    0,MAX_CONNECTED2*sizeof(void *));
        memset(connect13_res,0,MAX_CONNECTED2*sizeof(int));
        memset(connect14,    0,MAX_CONNECTED3*sizeof(void *));
        memset(connect14_res,0,MAX_CONNECTED3*sizeof(int));
        
        for (iconnect = 0; iconnect < MAX_CONNECTED; iconnect++) {
            if (!iatom_p->connect12[iconnect]) break;
            n_connect13++;
            connect13[n_connect13-1] = iatom_p->connect12[iconnect];
            connect13_res[n_connect13-1] = iatom_p->connect12_res[iconnect];
            
            for (jconnect = 0; jconnect < MAX_CONNECTED; jconnect++) {
                if (!iatom_p->connect12[iconnect]->connect12[jconnect]) break;
                n_connect13++;
                connect13[n_connect13-1] = iatom_p->connect12[iconnect]->connect12[jconnect];
                connect13_res[n_connect13-1] = iatom_p->connect12[iconnect]->connect12_res[jconnect];

                for (kconnect = 0; kconnect < MAX_CONNECTED; kconnect++) {
                    if (!iatom_p->connect12[iconnect]->connect12[jconnect]->connect12[kconnect]) break;
                    n_connect14++;
                    connect14[n_connect14-1] = iatom_p->connect12[iconnect]->connect12[jconnect]->connect12[kconnect];
                    connect14_res[n_connect14-1] = iatom_p->connect12[iconnect]->connect12[jconnect]->connect12_res[kconnect];
                }
            }
        }
        
        for (jatom=0; jatom<prot.res[jres].conf[jconf].n_atom; jatom++) {
            jatom_p = &prot.res[jres].conf[jconf].atom[jatom];
            if (!jatom_p->on) continue;
            
            for (iconnect = 0; iconnect < n_connect13; iconnect++) {
                if (jres == connect13_res[iconnect] && !strcmp(jatom_p->name,connect13[iconnect]->name)) break;
            }
            if (iconnect < n_connect13) continue;

            for (iconnect = 0; iconnect < n_connect14; iconnect++) {
                if (jres == connect14_res[iconnect] && !strcmp(jatom_p->name,connect14[iconnect]->name)) break;
            }
            if (iconnect < n_connect13) {
                e += coulomb(*iatom_p, *jatom_p) * env.factor_14lj;
                continue;
            }
            
                e += coulomb(*iatom_p, *jatom_p);
                //printf("iatom %d,jatom %d\n",iatom,jatom);
        }
    }
    
    return e;
}


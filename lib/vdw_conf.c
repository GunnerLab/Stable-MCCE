#include <stdio.h>
#include <string.h>
#include "mcce.h"

float vdw_conf(int i_res, int i_conf, int j_res, int j_conf, PROT prot)
{
    float  e = 0.0;
    int    iatom, jatom;
    ATOM   *connect123[MAX_CONNECTED2],    *connect14[MAX_CONNECTED3];
    ATOM   *iatom_p, *jatom_p;
    int    connect123_res[MAX_CONNECTED2], connect14_res[MAX_CONNECTED3];
    int    iconnect, jconnect, kconnect, n_connect123, n_connect14, cal_vdw;

    if (!prot.res[i_res].cal_vdw) {
        if ( !param_get("CAL_VDW",prot.res[i_res].resName, "", &cal_vdw) ) {
            prot.res[i_res].cal_vdw = cal_vdw;
        }
        else prot.res[i_res].cal_vdw = 1;
    }
    if (!prot.res[j_res].cal_vdw) {
        if ( !param_get("CAL_VDW",prot.res[j_res].resName, "", &cal_vdw) ) {
            prot.res[j_res].cal_vdw = cal_vdw;
        }
        else prot.res[j_res].cal_vdw = 1;
    }
    if (!prot.res[i_res].cal_vdw || !prot.res[j_res].cal_vdw) return e;
    
    for (iatom=0; iatom<prot.res[i_res].conf[i_conf].n_atom; iatom++) {
        iatom_p = &prot.res[i_res].conf[i_conf].atom[iatom];
        if (!iatom_p->on) continue;
        
        /* skip atom if CAL_VDW is 0 */
        if (!prot.res[i_res].conf[i_conf].atom[iatom].cal_vdw) {
            if ( !param_get("CAL_VDW",prot.res[i_res].conf[i_conf].confName, prot.res[i_res].conf[i_conf].atom[iatom].name, &cal_vdw) ) {
                prot.res[i_res].conf[i_conf].atom[iatom].cal_vdw = cal_vdw;
            }
            else prot.res[i_res].conf[i_conf].atom[iatom].cal_vdw = 1;
        }
        if (!prot.res[i_res].conf[i_conf].atom[iatom].cal_vdw) continue;

        n_connect123 = 0;
        n_connect14 = 0;
        memset(connect123,    0,MAX_CONNECTED2*sizeof(ATOM *));
        memset(connect123_res,0,MAX_CONNECTED2*sizeof(int));
        memset(connect14,    0,MAX_CONNECTED3*sizeof(ATOM *));
        memset(connect14_res,0,MAX_CONNECTED3*sizeof(int));

        for (iconnect = 0; iconnect < MAX_CONNECTED; iconnect++) {
            if (!iatom_p->on) break;
            if (!iatom_p->connect12[iconnect]) break;
            n_connect123++;
            connect123[n_connect123-1] = iatom_p->connect12[iconnect];
            connect123_res[n_connect123-1] = iatom_p->connect12_res[iconnect];

            for (jconnect = 0; jconnect < MAX_CONNECTED; jconnect++) {
                if (!iatom_p->connect12[iconnect]->on) break;
                if (!iatom_p->connect12[iconnect]->connect12[jconnect]) break;
                n_connect123++;
                connect123[n_connect123-1] = iatom_p->connect12[iconnect]->connect12[jconnect];
                connect123_res[n_connect123-1] = iatom_p->connect12[iconnect]->connect12_res[jconnect];

                for (kconnect = 0; kconnect < MAX_CONNECTED; kconnect++) {
                    if (!iatom_p->connect12[iconnect]->connect12[jconnect]->on) break;
                    if (!iatom_p->connect12[iconnect]->connect12[jconnect]->connect12[kconnect]) break;
                    n_connect14++;
                    connect14[n_connect14-1] = iatom_p->connect12[iconnect]->connect12[jconnect]->connect12[kconnect];
                    connect14_res[n_connect14-1] = iatom_p->connect12[iconnect]->connect12[jconnect]->connect12_res[kconnect];
                }
            }
        }

        for (jatom=0; jatom<prot.res[j_res].conf[j_conf].n_atom; jatom++) {
            jatom_p = &prot.res[j_res].conf[j_conf].atom[jatom];
            if (iatom_p == jatom_p) continue;
            if (!jatom_p->on) continue;
            
            /* skip atom if CAL_VDW is 0 */
            if (!prot.res[j_res].conf[j_conf].atom[jatom].cal_vdw) {
                if ( !param_get("CAL_VDW",prot.res[j_res].conf[j_conf].confName, prot.res[j_res].conf[j_conf].atom[jatom].name, &cal_vdw) ) {
                    prot.res[j_res].conf[j_conf].atom[jatom].cal_vdw = cal_vdw;
                }
                else prot.res[j_res].conf[j_conf].atom[jatom].cal_vdw = 1;
            }
            if (!prot.res[j_res].conf[j_conf].atom[jatom].cal_vdw) continue;

            for (iconnect = 0; iconnect < n_connect123; iconnect++) {
                if (j_res == connect123_res[iconnect] && !strcmp(jatom_p->name,connect123[iconnect]->name)) break;
            }
            if (iconnect < n_connect123) continue;

            for (iconnect = 0; iconnect < n_connect14; iconnect++) {
                if (j_res == connect14_res[iconnect] && !strcmp(jatom_p->name,connect14[iconnect]->name)) break;
            }
            if (iconnect < n_connect14) {
                e += vdw(*iatom_p, *jatom_p) * env.factor_14lj;

                /* print large atom-atom vdw
                if (vdw(*iatom_p, *jatom_p)* env.factor_14lj > 0.5) {
                    if (!i_conf || !j_conf || (i_res==j_res && i_conf==j_conf)) {
                        printf("Large vdw btw \"%s %s %c%04d%c%03d\" and \"%s %s %c%04d%c%03d\": %8.3f, dist=%8.4f\n",
                            iatom_p->name, prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].chainID,i_conf,
                            jatom_p->name, prot.res[j_res].resName,prot.res[j_res].chainID,prot.res[j_res].resSeq,prot.res[j_res].chainID,j_conf,
                            vdw(*iatom_p, *jatom_p)* env.factor_14lj, sqrt(ddvv(iatom_p->xyz, jatom_p->xyz)));
                    }
                }
                 print end */
                
                continue;
            }
                e += vdw(*iatom_p, *jatom_p);
                
                /* print large atom-atom vdw 
                if (vdw(*iatom_p, *jatom_p) > 0.5) {
                    if (!i_conf || !j_conf || (i_res==j_res && i_conf==j_conf)) {
                        printf("Large vdw btw \"%s %s %c%04d%c%03d\" and \"%s %s %c%04d%c%03d\": %8.3f, dist=%8.4f\n",
                            iatom_p->name, prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].chainID,i_conf,
                            jatom_p->name, prot.res[j_res].resName,prot.res[j_res].chainID,prot.res[j_res].resSeq,prot.res[j_res].chainID,j_conf,
                            vdw(*iatom_p, *jatom_p), sqrt(ddvv(iatom_p->xyz, jatom_p->xyz)));
                    }
                }
                print end */
        }
    }

    if (i_res == j_res && i_conf == j_conf) return 0.5*e;
    else return e;
}

float vdw_conf_hv(int i_res, int i_conf, int j_res, int j_conf, PROT prot)
{
    float  e = 0.0;
    int    iatom, jatom;
    ATOM   *connect123[MAX_CONNECTED2],    *connect14[MAX_CONNECTED3];
    ATOM   *iatom_p, *jatom_p;
    int    connect123_res[MAX_CONNECTED2], connect14_res[MAX_CONNECTED3];
    int    iconnect, jconnect, kconnect, n_connect123, n_connect14, cal_vdw;

    if (!prot.res[i_res].conf[i_conf].n_atom) return e;
    if (!prot.res[j_res].conf[j_conf].n_atom) return e;
    
    if (!prot.res[i_res].cal_vdw) {
        if ( !param_get("CAL_VDW",prot.res[i_res].resName, "", &cal_vdw) ) {
            prot.res[i_res].cal_vdw = cal_vdw;
        }
        else prot.res[i_res].cal_vdw = 1;
    }
    if (!prot.res[j_res].cal_vdw) {
        if ( !param_get("CAL_VDW",prot.res[j_res].resName, "", &cal_vdw) ) {
            prot.res[j_res].cal_vdw = cal_vdw;
        }
        else prot.res[j_res].cal_vdw = 1;
    }
    if (!prot.res[i_res].cal_vdw || !prot.res[j_res].cal_vdw) return e;

    for (iatom=0; iatom<prot.res[i_res].conf[i_conf].n_atom; iatom++) {
        iatom_p = &prot.res[i_res].conf[i_conf].atom[iatom];
        if (!iatom_p->on) continue;
        if (iatom_p->name[1] == 'H') continue;

        n_connect123 = 0;
        n_connect14 = 0;
        memset(connect123,    0,MAX_CONNECTED2*sizeof(ATOM *));
        memset(connect123_res,0,MAX_CONNECTED2*sizeof(int));
        memset(connect14,    0,MAX_CONNECTED3*sizeof(ATOM *));
        memset(connect14_res,0,MAX_CONNECTED3*sizeof(int));

        for (iconnect = 0; iconnect < MAX_CONNECTED; iconnect++) {
            if (!iatom_p->on) break;
            if (!iatom_p->connect12[iconnect]) break;
            n_connect123++;
            connect123[n_connect123-1] = iatom_p->connect12[iconnect];
            connect123_res[n_connect123-1] = iatom_p->connect12_res[iconnect];

            for (jconnect = 0; jconnect < MAX_CONNECTED; jconnect++) {
                if (!iatom_p->connect12[iconnect]->on) break;
                if (!iatom_p->connect12[iconnect]->connect12[jconnect]) break;
                n_connect123++;
                connect123[n_connect123-1] = iatom_p->connect12[iconnect]->connect12[jconnect];
                connect123_res[n_connect123-1] = iatom_p->connect12[iconnect]->connect12_res[jconnect];

                for (kconnect = 0; kconnect < MAX_CONNECTED; kconnect++) {
                    if (!iatom_p->connect12[iconnect]->connect12[jconnect]->on) break;
                    if (!iatom_p->connect12[iconnect]->connect12[jconnect]->connect12[kconnect]) break;
                    n_connect14++;
                    connect14[n_connect14-1] = iatom_p->connect12[iconnect]->connect12[jconnect]->connect12[kconnect];
                    connect14_res[n_connect14-1] = iatom_p->connect12[iconnect]->connect12[jconnect]->connect12_res[kconnect];
                }
            }
        }

        for (jatom=0; jatom<prot.res[j_res].conf[j_conf].n_atom; jatom++) {
            jatom_p = &prot.res[j_res].conf[j_conf].atom[jatom];
            if (!jatom_p->on) continue;
            if (jatom_p->name[1] == 'H') continue;
            if (iatom_p == jatom_p) continue;

            for (iconnect = 0; iconnect < n_connect123; iconnect++) {
                if (j_res == connect123_res[iconnect] && !strcmp(jatom_p->name,connect123[iconnect]->name)) break;
            }
            if (iconnect < n_connect123) continue;

            for (iconnect = 0; iconnect < n_connect14; iconnect++) {
                if (j_res == connect14_res[iconnect] && !strcmp(jatom_p->name,connect14[iconnect]->name)) break;
            }
            if (iconnect < n_connect14) {
                e += vdw(*iatom_p, *jatom_p) * env.factor_14lj;
                continue;
            }
                e += vdw(*iatom_p, *jatom_p);
        }
    }

    if (i_res == j_res && i_conf == j_conf) return 0.5*e;
    else return e;
}


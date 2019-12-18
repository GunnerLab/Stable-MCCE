#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

#define PRINT_THR       4
//extern int   n_elem;
//extern float C6_matrix[N_ELEM_MAX][N_ELEM_MAX];
//extern float C12_matrix[N_ELEM_MAX][N_ELEM_MAX];
extern void collect_all_connect(int i_res, int i_conf, int i_atom, PROT prot, int *n_connect12, ATOM ***connect12, int *n_connect13, ATOM ***connect13, int *n_connect14, ATOM ***connect14);

/*
int setup_C6_C12(PROT prot) {
    int   i_elem, j_elem;
    int   i_res, i_conf, i_atom;
    char  elem[N_ELEM_MAX];

    // Setup C6,C12 matrix 
    n_elem = 0;
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                for (i_elem=0;i_elem<n_elem;i_elem++)
                    if (elem[i_elem] == prot.res[i_res].conf[i_conf].atom[i_atom].name[1]) break;
                if (i_elem == n_elem) {
                    n_elem++;
                    elem[i_elem] = prot.res[i_res].conf[i_conf].atom[i_atom].name[1];
                }
                prot.res[i_res].conf[i_conf].atom[i_atom].i_elem = i_elem;
            }
        }
    }
    for (i_elem=0;i_elem<n_elem;i_elem++) {
        for (j_elem=0;j_elem<n_elem;j_elem++) {
            char pair[4];
            pair[0] = elem[i_elem];
            pair[1] = '-';
            pair[2] = elem[j_elem];
            pair[3] = '\0';
            
            if(param_get("VDWAMBER", "C6",  pair, &C6_matrix[i_elem][j_elem])) {
                param_get("VDWAMBER", "C6",  "X-X", &C6_matrix[i_elem][j_elem]);
            }
            if (param_get("VDWAMBER", "C12", pair, &C12_matrix[i_elem][j_elem])) {
                param_get("VDWAMBER", "C12", "X-X", &C12_matrix[i_elem][j_elem]);
            }
        }
    }
    
    return n_elem;
}
*/

void setup_vdw_fast(PROT prot) {
    int i_res;
    for (i_res=0;i_res<prot.n_res;i_res++) {
        setup_vdw_fast_res(i_res,prot);
    }
}

void setup_vdw_fast_res(int i_res, PROT prot) {
    /* This subroutine sets up back index, r_min, r_max for all residues and conformers */
    int i_conf,i_atom,cal_vdw;
    ATOM *atom_p;
    int res_initialized = 0;
    
    prot.res[i_res].i_res_prot = i_res;
    
    if ( !param_get("CAL_VDW",prot.res[i_res].resName, "", &cal_vdw) ) {
        prot.res[i_res].cal_vdw = cal_vdw;
    }
    else prot.res[i_res].cal_vdw = 1;
    
    /* get r_min, r_max for each res and conf */
    for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
        prot.res[i_res].conf[i_conf].i_res_prot = i_res;
        prot.res[i_res].conf[i_conf].i_conf_res = i_conf;
        
        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
            if (!atom_p->on) continue;
            prot.res[i_res].conf[i_conf].r_min = atom_p->xyz;
            prot.res[i_res].conf[i_conf].r_max = atom_p->xyz;
            
            if (!res_initialized) {
                prot.res[i_res].r_min = atom_p->xyz;
                prot.res[i_res].r_max = atom_p->xyz;
                res_initialized = 1;
            }
            break;
        }
    }
    
    for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
            atom_p->i_res_prot  = i_res;
            atom_p->i_conf_res  = i_conf;
            atom_p->i_atom_conf = i_atom;
            if (!atom_p->on) continue;
            
            if (atom_p->xyz.x < prot.res[i_res].r_min.x)
                prot.res[i_res].r_min.x = atom_p->xyz.x;
            if (atom_p->xyz.y < prot.res[i_res].r_min.y)
                prot.res[i_res].r_min.y = atom_p->xyz.y;
            if (atom_p->xyz.z < prot.res[i_res].r_min.z)
                prot.res[i_res].r_min.z = atom_p->xyz.z;
            if (atom_p->xyz.x > prot.res[i_res].r_max.x)
                prot.res[i_res].r_max.x = atom_p->xyz.x;
            if (atom_p->xyz.y > prot.res[i_res].r_max.y)
                prot.res[i_res].r_max.y = atom_p->xyz.y;
            if (atom_p->xyz.z > prot.res[i_res].r_max.z)
                prot.res[i_res].r_max.z = atom_p->xyz.z;
            
            if (atom_p->xyz.x < prot.res[i_res].conf[i_conf].r_min.x)
                prot.res[i_res].conf[i_conf].r_min.x = atom_p->xyz.x;
            if (atom_p->xyz.y < prot.res[i_res].conf[i_conf].r_min.y)
                prot.res[i_res].conf[i_conf].r_min.y = atom_p->xyz.y;
            if (atom_p->xyz.z < prot.res[i_res].conf[i_conf].r_min.z)
                prot.res[i_res].conf[i_conf].r_min.z = atom_p->xyz.z;
            if (atom_p->xyz.x > prot.res[i_res].conf[i_conf].r_max.x)
                prot.res[i_res].conf[i_conf].r_max.x = atom_p->xyz.x;
            if (atom_p->xyz.y > prot.res[i_res].conf[i_conf].r_max.y)
                prot.res[i_res].conf[i_conf].r_max.y = atom_p->xyz.y;
            if (atom_p->xyz.z > prot.res[i_res].conf[i_conf].r_max.z)
                prot.res[i_res].conf[i_conf].r_max.z = atom_p->xyz.z;
        }
    }
    
}

void setup_connect_res(PROT prot, int i_res) {
    float d2;
    int i_conf,i_atom,i_connect;
    ATOM *atom_p;
    FILE *debug_fp;

    /* setup connectivity for for all conformers in one residue */
    prot.res[i_res].n_connect12 = calloc(prot.res[i_res].n_conf, sizeof(int *));
    prot.res[i_res].connect12   = calloc(prot.res[i_res].n_conf, sizeof(int *));
    prot.res[i_res].n_connect13 = calloc(prot.res[i_res].n_conf, sizeof(int *));
    prot.res[i_res].connect13   = calloc(prot.res[i_res].n_conf, sizeof(int *));
    prot.res[i_res].n_connect14 = calloc(prot.res[i_res].n_conf, sizeof(int *));
    prot.res[i_res].connect14   = calloc(prot.res[i_res].n_conf, sizeof(int *));

    for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
        if (!prot.res[i_res].conf[i_conf].n_atom) continue;
        prot.res[i_res].n_connect12[i_conf] = calloc(prot.res[i_res].conf[i_conf].n_atom, sizeof(int));
        prot.res[i_res].connect12[i_conf]   = calloc(prot.res[i_res].conf[i_conf].n_atom, sizeof(ATOM **));
        prot.res[i_res].n_connect13[i_conf] = calloc(prot.res[i_res].conf[i_conf].n_atom, sizeof(int));
        prot.res[i_res].connect13[i_conf]   = calloc(prot.res[i_res].conf[i_conf].n_atom, sizeof(ATOM **));
        prot.res[i_res].n_connect14[i_conf] = calloc(prot.res[i_res].conf[i_conf].n_atom, sizeof(int));
        prot.res[i_res].connect14[i_conf]   = calloc(prot.res[i_res].conf[i_conf].n_atom, sizeof(ATOM **));
        
        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
            if (!atom_p->on) continue;
            collect_all_connect(i_res, i_conf, i_atom, prot,
            &prot.res[i_res].n_connect12[i_conf][i_atom], &prot.res[i_res].connect12[i_conf][i_atom],
            &prot.res[i_res].n_connect13[i_conf][i_atom], &prot.res[i_res].connect13[i_conf][i_atom],
            &prot.res[i_res].n_connect14[i_conf][i_atom], &prot.res[i_res].connect14[i_conf][i_atom]);
        }
        
        /* find max distance for each connectivity list */
        prot.res[i_res].r12sq_max = 0.;
        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
            if (!atom_p->on) continue;
            for (i_connect=0; i_connect<prot.res[i_res].n_connect12[i_conf][i_atom]; i_connect++) {
                if (!prot.res[i_res].connect12[i_conf][i_atom][i_connect]->on) continue;
                d2 = ddvv(atom_p->xyz, prot.res[i_res].connect12[i_conf][i_atom][i_connect]->xyz);
                if (d2 > prot.res[i_res].r12sq_max) prot.res[i_res].r12sq_max = d2;
                if (d2 > PRINT_THR) {
                    debug_fp = fopen(env.debug_log,"a");
                    if (!debug_fp) debug_fp = stdout;
                    fprintf(debug_fp, "Warning, bond length btw %s-%s in res %s %c%04d is %8.3f\n",
                    atom_p->name, prot.res[i_res].connect12[i_conf][i_atom][i_connect]->name,prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,sqrt(d2));
                    fclose(debug_fp);
                }
            }
        }
        prot.res[i_res].r13sq_max = 0.;
        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
            for (i_connect=0; i_connect<prot.res[i_res].n_connect13[i_conf][i_atom]; i_connect++) {
                d2 = ddvv(atom_p->xyz, prot.res[i_res].connect13[i_conf][i_atom][i_connect]->xyz);
                if (d2 > prot.res[i_res].r13sq_max) prot.res[i_res].r13sq_max = d2;
            }
        }
        prot.res[i_res].r14sq_max = 0.;
        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
            for (i_connect=0; i_connect<prot.res[i_res].n_connect14[i_conf][i_atom]; i_connect++) {
                d2 = ddvv(atom_p->xyz, prot.res[i_res].connect14[i_conf][i_atom][i_connect]->xyz);
                if (d2 > prot.res[i_res].r14sq_max) prot.res[i_res].r14sq_max = d2;
            }
        }
        prot.res[i_res].r12sq_max += prot.res[i_res].r12sq_max*env.hv_relax_shake_tol*2.;
        prot.res[i_res].r13sq_max += prot.res[i_res].r13sq_max*env.hv_relax_shake_tol*2.;
        prot.res[i_res].r14sq_max += prot.res[i_res].r14sq_max*env.hv_relax_shake_tol*2.;
    }
}

void free_connect_res(PROT prot, int i_res) {
    int i_conf,i_atom;
    ATOM *atom_p;

    for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
        if (!prot.res[i_res].conf[i_conf].n_atom) continue;
        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
            if (!atom_p->on) continue;
            free(prot.res[i_res].connect12[i_conf][i_atom]);
            free(prot.res[i_res].connect13[i_conf][i_atom]);
            free(prot.res[i_res].connect14[i_conf][i_atom]);
        }
        free(prot.res[i_res].n_connect12[i_conf]);
        free(prot.res[i_res].n_connect13[i_conf]);
        free(prot.res[i_res].n_connect14[i_conf]);
        free(prot.res[i_res].connect12[i_conf]);
        free(prot.res[i_res].connect13[i_conf]);
        free(prot.res[i_res].connect14[i_conf]);
    }
    free(prot.res[i_res].n_connect12);
    free(prot.res[i_res].n_connect13);
    free(prot.res[i_res].n_connect14);
    free(prot.res[i_res].connect12);
    free(prot.res[i_res].connect13);
    free(prot.res[i_res].connect14);
}    

float vdw_conf_fast(int i_res, int i_conf, int j_res, int j_conf, PROT prot, int handle_hv) {
    /* This is a fast version of vdw_conf, pre-setup is need to use this function and to make calculation fast,
    to setup, call the setup functions before get into the vdw loop. See example:
    handle_hv = 0: full vdw.
    handle_hv = 1: heavy atom vdw.
    */
    float pair_vdw = 0., d2, d6, d12, e, C6, C12;
    VECTOR v1,v2;
    int i_atom, j_atom, i_connect, done;
    ATOM *iatom_p,*jatom_p;
    float cutoff_near2 = VDW_CUTOFF_NEAR * VDW_CUTOFF_NEAR;
    float cutoff_far2  = VDW_CUTOFF_FAR  * VDW_CUTOFF_FAR;
    /* FILE *debug_fp; */

    if (!prot.res[i_res].conf[i_conf].n_atom) return pair_vdw;
    if (!prot.res[j_res].conf[j_conf].n_atom) return pair_vdw;
    if (!prot.res[i_res].cal_vdw) return pair_vdw;
    if (!prot.res[j_res].cal_vdw) return pair_vdw;
    
    for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
        iatom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
        if (!iatom_p->on) continue;
        if (handle_hv) {
            if (iatom_p->name[1] == 'H') continue;
        }
        v1 = iatom_p->xyz;
        for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
            jatom_p = &prot.res[j_res].conf[j_conf].atom[j_atom];
            if (!jatom_p->on) continue;
            if (handle_hv) {
                if (jatom_p->name[1] == 'H') continue;
            }
            v2 = jatom_p->xyz;
            
            if (iatom_p == jatom_p) continue;
            d2 = ddvv(v1, v2);
            if (d2 < cutoff_far2) {
                done = 0;
                
                if (d2 <= prot.res[i_res].r12sq_max + 1.0) {
                    for (i_connect=0;i_connect<prot.res[i_res].n_connect12[i_conf][i_atom];i_connect++) {
                        if (jatom_p == prot.res[i_res].connect12[i_conf][i_atom][i_connect]) {
                            done = 1;
                        }
                    }
                }
                if (done) continue;
                
                if (d2 <= prot.res[i_res].r13sq_max + 1.0) {
                    for (i_connect=0;i_connect<prot.res[i_res].n_connect13[i_conf][i_atom];i_connect++) {
                        if (jatom_p == prot.res[i_res].connect13[i_conf][i_atom][i_connect]) {
                            done = 1;
                        }
                    }
                }
                if (done) continue;
                
                /* calculate vdw */
                if (d2 < cutoff_near2)
                    e = VDW_ELIMIT_NEAR;            /* Cutoff */
                else {
                    float sig_min = iatom_p->vdw_rad + jatom_p->vdw_rad;
                    float eps = sqrt(iatom_p->vdw_eps*jatom_p->vdw_eps);
                    C12 = eps*pow(sig_min,12);
                    C6 = 2.*eps*pow(sig_min,6);

                    //C6 = C6_matrix[iatom_p->i_elem][jatom_p->i_elem];
                    //C12 = C12_matrix[iatom_p->i_elem][jatom_p->i_elem];
                    d6 = d2*d2*d2;
                    d12 = d6*d6;
                    e = C12/d12 - C6/d6;
                }
                
                if (d2 <= prot.res[i_res].r14sq_max + 1.0) {
                    for (i_connect=0;i_connect<prot.res[i_res].n_connect14[i_conf][i_atom];i_connect++) {
                        if (jatom_p == prot.res[i_res].connect14[i_conf][i_atom][i_connect]) {
                            
                            pair_vdw += e * env.factor_14lj;
                            
                            /* print large atom-atom vdw 
                            if (e * env.factor_14lj > 2.0) {
                                if (!i_conf || !j_conf || (i_res==j_res && i_conf==j_conf)) {
                                    //debug_fp = fopen(env.debug_log,"a");
                                    //if (!debug_fp) debug_fp = stdout;
                                    //fprintf(debug_fp, "Warning, interaction btw atom %s of res %s %c%04d and atom %s of res %s %c%04d is %8.3f, d=%11.7f, r12max=%11.7f\n",
                                        printf("Large vdw btw \"%s %s %c%04d\" and \"%s %s %c%04d\": %8.3f, dist=%8.4f\n",
                                            iatom_p->name, prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,
                                            jatom_p->name, prot.res[j_res].resName,prot.res[j_res].chainID,prot.res[j_res].resSeq, e * env.factor_14lj, sqrt(d2));
                                        //fclose(debug_fp);
                                }
                            }
                             print end */
                            
                            done = 1;
                            break;
                        }
                    }
                }
                if (done) continue;
                
                /* print large atom-atom vdw 
                if (e > 2.0) {
                    if (!i_conf || !j_conf || (i_res==j_res && i_conf==j_conf)) {
                        //debug_fp = fopen(env.debug_log,"a");
                        //if (!debug_fp) debug_fp = stdout;
                        //fprintf(debug_fp, "Warning, interaction btw atom %s of res %s %c%04d and atom %s of res %s %c%04d is %8.3f, d=%11.7f, r12max=%11.7f\n",
                        printf("Large vdw btw \"%s %s %c%04d\" and \"%s %s %c%04d\": %8.3f, dist=%8.4f\n",
                        iatom_p->name, prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,
                        jatom_p->name, prot.res[j_res].resName,prot.res[j_res].chainID,prot.res[j_res].resSeq, e, sqrt(d2));
                        //fclose(debug_fp);
                    }
                }
                 print end */
                
                pair_vdw += e;
                //printf("%s-%s: %8.3f\n",iatom_p->name,jatom_p->name,e);
            }
        }
    }
    
    if (i_res == j_res && i_conf == j_conf) return pair_vdw/2.;
    else return pair_vdw;
}

float vdw_conf_fast_print(int i_res, int i_conf, int j_res, int j_conf, PROT prot) {
    /* This is a fast version of vdw_conf, pre-setup is need to use this function and to make calculation fast,
    to setup, call the setup functions before get into the vdw loop. See example:
    */
    float pair_vdw = 0., d2, d6, d12, e, C6, C12;
    VECTOR v1,v2;
    int i_atom, j_atom, i_connect, done;
    ATOM *iatom_p,*jatom_p;
    float cutoff_near2 = VDW_CUTOFF_NEAR * VDW_CUTOFF_NEAR;
    float cutoff_far2  = VDW_CUTOFF_FAR  * VDW_CUTOFF_FAR;
    /* FILE *debug_fp; */

    if (!prot.res[i_res].conf[i_conf].n_atom) return pair_vdw;
    if (!prot.res[j_res].conf[j_conf].n_atom) return pair_vdw;
    if (!prot.res[i_res].cal_vdw) {
        printf("CAL_VDW parameter of residue %s is off\n",prot.res[i_res].resName);
        return pair_vdw;
    }
    if (!prot.res[j_res].cal_vdw) {
        printf("CAL_VDW parameter of residue %s is off\n",prot.res[j_res].resName);
        return pair_vdw;
    }

    
    for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
        iatom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
        if (!iatom_p->on) continue;

        v1 = iatom_p->xyz;
        for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
            jatom_p = &prot.res[j_res].conf[j_conf].atom[j_atom];
            if (!jatom_p->on) continue;

            v2 = jatom_p->xyz;
            
            if (iatom_p == jatom_p) continue;
            d2 = ddvv(v1, v2);
            if (d2 < cutoff_far2) {
                done = 0;
                
                if (d2 <= prot.res[i_res].r12sq_max + 1.0) {
                    for (i_connect=0;i_connect<prot.res[i_res].n_connect12[i_conf][i_atom];i_connect++) {
                        if (jatom_p == prot.res[i_res].connect12[i_conf][i_atom][i_connect]) {
                            done = 1;
                        }
                    }
                }
                if (done) continue;
                
                if (d2 <= prot.res[i_res].r13sq_max + 1.0) {
                    for (i_connect=0;i_connect<prot.res[i_res].n_connect13[i_conf][i_atom];i_connect++) {
                        if (jatom_p == prot.res[i_res].connect13[i_conf][i_atom][i_connect]) {
                            done = 1;
                        }
                    }
                }
                if (done) continue;
                
                /* calculate vdw */
                if (d2 < cutoff_near2)
                    e = VDW_ELIMIT_NEAR;            /* Cutoff */
                else {
                    float sig_min = iatom_p->vdw_rad + jatom_p->vdw_rad;
                    float eps = sqrt(iatom_p->vdw_eps*jatom_p->vdw_eps);
                    C12 = eps*pow(sig_min,12);
                    C6 = 2.*eps*pow(sig_min,6);

                    d6 = d2*d2*d2;
                    d12 = d6*d6;
                    e = C12/d12 - C6/d6;
                }
                
                /* 1,4 connectivity */
                if (d2 <= prot.res[i_res].r14sq_max + 1.0) {
                    for (i_connect=0;i_connect<prot.res[i_res].n_connect14[i_conf][i_atom];i_connect++) {
                        if (jatom_p == prot.res[i_res].connect14[i_conf][i_atom][i_connect]) {
                            
                            pair_vdw += e * env.factor_14lj;
                            
                            /* print large atom-atom vdw */
                            printf("VDW (1,4) btw \"%s %s\" and \"%s %s\": %8.3f, dist=%8.4f\n",
                                iatom_p->name, prot.res[i_res].conf[i_conf].uniqID,
                                jatom_p->name, prot.res[j_res].conf[j_conf].uniqID,
                                e * env.factor_14lj, sqrt(d2));
                            /* print end */
                            
                            done = 1;
                            break;
                        }
                    }
                }
                if (done) continue;
                
                /* print large atom-atom vdw */
                printf("VDW       btw \"%s %s\" and \"%s %s\": %8.3f, dist=%8.4f\n",
                    iatom_p->name, prot.res[i_res].conf[i_conf].uniqID,
                    jatom_p->name, prot.res[j_res].conf[j_conf].uniqID,
                    e, sqrt(d2));
                /* print end */
                
                pair_vdw += e;
            }
        }
    }
    
    if (i_res == j_res && i_conf == j_conf) return pair_vdw/2.;
    else return pair_vdw;
}

float coulomb_conf_fast(int i_res, int i_conf, int j_res, int j_conf, PROT prot) {
    /* This is a fast version of vdw_conf, pre-setup is need to use this function and to make calculation fast,
    to setup, call the setup functions before get into the vdw loop. See example:
    */
    float pair_coulomb = 0., d, d2, e;
    VECTOR v1,v2;
    int i_atom, j_atom, i_connect, done;
    ATOM *iatom_p,*jatom_p;
    
    if (!prot.res[i_res].conf[i_conf].n_atom) return pair_coulomb;
    if (!prot.res[j_res].conf[j_conf].n_atom) return pair_coulomb;
    
    for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
        iatom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
        if (!iatom_p->on) continue;
        v1 = iatom_p->xyz;
        for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
            jatom_p = &prot.res[j_res].conf[j_conf].atom[j_atom];
            if (!jatom_p->on) continue;
            v2 = jatom_p->xyz;
            
            if (iatom_p == jatom_p) continue;
            d2 = ddvv(v1, v2);
            done = 0;
            
            if (d2 <= prot.res[i_res].r12sq_max) {
                for (i_connect=0;i_connect<prot.res[i_res].n_connect12[i_conf][i_atom];i_connect++) {
                    if (jatom_p == prot.res[i_res].connect12[i_conf][i_atom][i_connect]) {
                        done = 1;
                    }
                }
            }
            if (done) continue;
            
            if (d2 <= prot.res[i_res].r13sq_max) {
                for (i_connect=0;i_connect<prot.res[i_res].n_connect13[i_conf][i_atom];i_connect++) {
                    if (jatom_p == prot.res[i_res].connect13[i_conf][i_atom][i_connect]) {
                        done = 1;
                    }
                }
            }
            if (done) continue;
            
            /* calculate coulumb */
            d = sqrt(d2);
            if (d < 0.8) d = 0.8;
            e = 331.5*iatom_p->crg*jatom_p->crg/(env.epsilon_coulomb * d);
            
            if (d2 <= prot.res[i_res].r14sq_max) {
                for (i_connect=0;i_connect<prot.res[i_res].n_connect14[i_conf][i_atom];i_connect++) {
                    if (jatom_p == prot.res[i_res].connect14[i_conf][i_atom][i_connect]) {
                        
                        pair_coulomb += e * env.factor_14lj;
                        done = 1;
                        break;
                    }
                }
            }
            if (done) continue;
            
            pair_coulomb += e;
        }
    }
    
    if (i_res == j_res && i_conf == j_conf) return pair_coulomb/2.;
    else return pair_coulomb;
}

int out_of_range(VECTOR i_min, VECTOR i_max, VECTOR j_min, VECTOR j_max, float range2) {
    VECTOR bigger_min, smaller_max, dist;
    
    /*
    .i:       |------------|
    .        r_min        r_max
    .j:                               |---------|
    .                                r_min     r_max
    .                              bigger_min
    .                smaller_max
    .                       --------->
    .                        distance
    */
    bigger_min = i_min;
    if (j_min.x > bigger_min.x) bigger_min.x = j_min.x;
    if (j_min.y > bigger_min.y) bigger_min.y = j_min.y;
    if (j_min.z > bigger_min.z) bigger_min.z = j_min.z;
    smaller_max = i_max;
    if (j_max.x < smaller_max.x) smaller_max.x = j_max.x;
    if (j_max.y < smaller_max.y) smaller_max.y = j_max.y;
    if (j_max.z < smaller_max.z) smaller_max.z = j_max.z;
    
    dist = vector_vminusv(bigger_min, smaller_max);
    if (dist.x < 0) dist.x = 0;
    if (dist.y < 0) dist.y = 0;
    if (dist.z < 0) dist.z = 0;
    
    if (vdotv(dist, dist) > range2) return 1;
    else return 0;
}


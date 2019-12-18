#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

void collect_all_connect12(int i_res, int i_conf, int i_atom, PROT prot, int (*n_connect12), ATOM **(*connect12)) {
    ATOM *atom_p;
    int j_res, j_conf, j_atom, i_connect;
    
    atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
    (*n_connect12) = 0;
    (*connect12) = NULL;
    
    /* Collect all 12 connnectivity */
    for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
        if (!atom_p->connect12[i_connect]) break;
        if (!atom_p->connect12[i_connect]->on) continue;
        (*n_connect12)++;
        (*connect12) = realloc((*connect12), (*n_connect12)*sizeof(ATOM *));
        ((*connect12))[(*n_connect12)-1] = atom_p->connect12[i_connect];
        
        /* search for connectivity duplicated in more than one conformer */
        if (i_res == atom_p->connect12[i_connect]->i_res_prot
            && i_conf == atom_p->connect12[i_connect]->i_conf_res) continue; /* within same conformer */
        if (atom_p->connect12[i_connect]->i_conf_res == 0) continue; /* connect to backbone */
        
        j_res = atom_p->connect12[i_connect]->i_res_prot;
        for (j_conf=1;j_conf<prot.res[j_res].n_conf;j_conf++) {
            if (j_conf == atom_p->connect12[i_connect]->i_conf_res) continue;
            if (!strcmp(prot.res[j_res].conf[j_conf].confName,
            prot.res[j_res].conf[atom_p->connect12[i_connect]->i_conf_res].confName)) {
                j_atom = atom_p->connect12[i_connect]->i_atom_conf;
                
                (*n_connect12)++;
                ((*connect12)) = realloc((*connect12), (*n_connect12)*sizeof(ATOM *));
                ((*connect12))[(*n_connect12)-1] = &prot.res[j_res].conf[j_conf].atom[j_atom];
            }
            else {
                j_atom = iatom(prot.res[j_res].conf[j_conf].confName, atom_p->connect12[i_connect]->name);
                
                if (j_atom != -1) {
                    (*n_connect12)++;
                    (*connect12) = realloc((*connect12), (*n_connect12)*sizeof(ATOM *));
                    ((*connect12))[(*n_connect12)-1] = &prot.res[j_res].conf[j_conf].atom[j_atom];
                }
            }
        }
    }
}

void collect_all_connect(int i_res, int i_conf, int i_atom, PROT prot, int (*n_connect12), ATOM ***connect12, int (*n_connect13), ATOM ***connect13, int (*n_connect14), ATOM ***connect14) {
    ATOM *atom_p;
    int j_res, j_conf, j_atom, i_connect,j_connect;
    int n_connect12_j;
    ATOM **connect12_j;

    atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
    (*n_connect12) = 0;
    (*connect12) = NULL;
    (*n_connect13) = 0;
    (*connect13) = NULL;
    (*n_connect14) = 0;
    (*connect14) = NULL;

    collect_all_connect12(i_res, i_conf, i_atom, prot, n_connect12, connect12);
    
    /* Collect all 13 connnectivity */
    for (i_connect=0; i_connect<(*n_connect12); i_connect++) {
        j_res  = ((*connect12))[i_connect]->i_res_prot;
        j_conf = ((*connect12))[i_connect]->i_conf_res;
        j_atom = ((*connect12))[i_connect]->i_atom_conf;

        collect_all_connect12(j_res, j_conf, j_atom, prot, &n_connect12_j, &connect12_j);
        
        (*n_connect13) += n_connect12_j;
        (*connect13) = realloc((*connect13), (*n_connect13)*sizeof(ATOM *));
        memcpy( (*connect13)+(*n_connect13)-(n_connect12_j), connect12_j, (n_connect12_j)*sizeof(ATOM *));
        free(connect12_j);
    }
    for (i_connect=(*n_connect13)-1; i_connect>=0; i_connect--) {
        if ((*connect13)[i_connect] == atom_p) {
            (*n_connect13)--;
            memmove((*connect13)+i_connect,(*connect13)+i_connect+1,((*n_connect13)-i_connect)*sizeof(ATOM *));
            continue;
        }
        
        if ( (*connect13)[i_connect]->i_res_prot == atom_p->i_res_prot) {
            if ((*connect13)[i_connect]->i_conf_res != atom_p->i_conf_res) {
                if ((*connect13)[i_connect]->i_conf_res*atom_p->i_conf_res) {
                    (*n_connect13)--;
                    memmove((*connect13)+i_connect,(*connect13)+i_connect+1,((*n_connect13)-i_connect)*sizeof(ATOM *));
                    continue;
                }
            }
        }
        
        for (j_connect=0; j_connect<(*n_connect12); j_connect++) {
            if ((*connect13)[i_connect] == ((*connect12))[j_connect]) break;
        }
        if (j_connect<(*n_connect12)) {
            (*n_connect13)--;
            memmove((*connect13)+i_connect,(*connect13)+i_connect+1,((*n_connect13)-i_connect)*sizeof(ATOM *));
            continue;
        }
    }
    (*connect13) = realloc((*connect13), (*n_connect13)*sizeof(ATOM *));
    
    /* Collect all 14 connnectivity */
    for (i_connect=0; i_connect<(*n_connect13); i_connect++) {
        j_res  = (*connect13)[i_connect]->i_res_prot;
        j_conf = (*connect13)[i_connect]->i_conf_res;
        j_atom = (*connect13)[i_connect]->i_atom_conf;
        
        collect_all_connect12(j_res, j_conf, j_atom, prot, &n_connect12_j, &connect12_j);
        
        (*n_connect14) += n_connect12_j;
        (*connect14) = realloc((*connect14), (*n_connect14)*sizeof(ATOM *));
        memcpy( ((*connect14))+(*n_connect14)-n_connect12_j, connect12_j, (n_connect12_j)*sizeof(ATOM *));
        free(connect12_j);
    }
    for (i_connect=(*n_connect14)-1; i_connect>=0; i_connect--) {
        if (((*connect14))[i_connect] == atom_p) {
            (*n_connect14)--;
            memmove((*connect14)+i_connect,(*connect14)+i_connect+1,((*n_connect14)-i_connect)*sizeof(ATOM *));
            continue;
        }
        
        if (((*connect14))[i_connect]->i_res_prot == atom_p->i_res_prot) {
            if (((*connect14))[i_connect]->i_conf_res != atom_p->i_conf_res) {
                if (((*connect14))[i_connect]->i_conf_res * atom_p->i_conf_res) {
                    (*n_connect14)--;
                    memmove((*connect14)+i_connect,(*connect14)+i_connect+1,((*n_connect14)-i_connect)*sizeof(ATOM *));
                    continue;
                }
            }
        }
        
        for (j_connect=0; j_connect<(*n_connect12); j_connect++) {
            if (((*connect14))[i_connect] == ((*connect12))[j_connect]) break;
        }
        if (j_connect<(*n_connect12)) {
            (*n_connect14)--;
            memmove((*connect14)+i_connect,(*connect14)+i_connect+1,((*n_connect14)-i_connect)*sizeof(ATOM *));
            continue;
        }
        
        for (j_connect=0; j_connect<(*n_connect13); j_connect++) {
            if (((*connect14))[i_connect] == (*connect13)[j_connect]) break;
        }
        if (j_connect<(*n_connect13)) {
            (*n_connect14)--;
            memmove((*connect14)+i_connect,(*connect14)+i_connect+1,((*n_connect14)-i_connect)*sizeof(ATOM *));
            continue;
        }
    }
    (*connect14) = realloc((*connect14), (*n_connect14)*sizeof(ATOM *));
}

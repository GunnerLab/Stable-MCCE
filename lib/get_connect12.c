#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

#define BOND_THR 2.4
#define MAX_LIGS 1000

int get_connect12(PROT prot) {
    int         i_res, i_conf;
    int         ret_val, err = 0;
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            ret_val = get_connect12_conf(i_res, i_conf, prot);
            if (ret_val == -1) {
                printf("   Error! get_connect12(): Error in CONNECT parameter, refer to debug log file: %s for detail\n", env.debug_log);
                return -1;
            }
            else err += ret_val;
        }
    }
    if (err)
        printf("   Error in connectivity, grep \"get_connect12\" from %s to find details!\n",env.debug_log);
    
    return 0;
}

int get_connect12_conf(int i_res, int i_conf, PROT prot) {
    CONNECT     connect;
    FILE        *debug_fp;
    int         i_atom;
    int         j_connect, k_connect, n_connect;
    int         j_res, j_conf, j_atom;
    int         k_lig, l_connect;
    int         n_ligs, n_conf;
    int         connect_found, lig_treated;
    RES         *res_p;
    CONF        *conf_p;
    ATOM        *atom_p,*jatom_p, *ligs[MAX_LIGS];
    int         ligs_res[MAX_LIGS];
    int         err = 0;
    
    res_p = &prot.res[i_res];
    conf_p = &prot.res[i_res].conf[i_conf];
    /* Looping over all atoms */
    for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
        atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
        if (!atom_p->on) continue;
        for (j_connect=0; j_connect<MAX_CONNECTED; j_connect++) {
            atom_p->connect12[j_connect] = NULL;
        }
        
        /* Error checking for connectivity parameter */
        memset(atom_p->connect12, 0, MAX_CONNECTED*sizeof(void *));
        if( param_get("CONNECT", conf_p->confName, atom_p->name, &connect) ) {
            debug_fp = fopen(env.debug_log, "a");
            fprintf(debug_fp, "   Error! get_connect12(): Can't find CONNECT parameter of conformer \"%s\" atom \"%s\"\n", conf_p->confName, atom_p->name);
            fclose(debug_fp);
            err++;
            continue;
        }
        if( connect.n > MAX_CONNECTED ) {
            debug_fp = fopen(env.debug_log, "a");
            fprintf(debug_fp, "   Error! get_connect12(): Error in CONNECT parameter for conformer \"%s\" atom \"%s\", number of connected atoms is bigger than array size\nCheck parameter file or fix mcce.h with a bigger MAX_CONNECTED \n", conf_p->confName, atom_p->name);
            fclose(debug_fp);
            return USERERR;
        }
        
        lig_treated = 0;        /* for each atom, ligand connectivy is checked at most once, no matter how many connectivities this atom has */
        /* Looping over each connectivity member */
        for (j_connect=0; j_connect<connect.n; j_connect++) {
            connect_found = 0;
            
            /* Ligand type connectivity or not */
            if (!connect.atom[j_connect].ligand) {      /* NOT ligand type */
                if (!connect.atom[j_connect].res_offset) { /* If within the same residue (off_set is 0) */
                    j_res = i_res;
                    if ( !param_get("IATOM", conf_p->confName, connect.atom[j_connect].name, &j_atom) ) {
                        /* First search atom in the same conformer */
                        j_conf = i_conf;
                        atom_p->connect12[j_connect] = &prot.res[j_res].conf[j_conf].atom[j_atom];
                        atom_p->connect12_res[j_connect] = j_res;
                        connect_found = 1;
                    }
                    else {
                        /* If not in same conformer, there are two situations: 
                        1. atom is in side chain conformer, connecting to backbone;
                        2. atom is in backbone conformer, connecting to side chain; */
                        for (j_conf = 0; j_conf < prot.res[j_res].n_conf; j_conf++) {
                            if (j_conf == i_conf) continue;
                            if ( !param_get("IATOM", prot.res[j_res].conf[j_conf].confName, connect.atom[j_connect].name, &j_atom) ) {
                                atom_p->connect12[j_connect] = &prot.res[j_res].conf[j_conf].atom[j_atom];
                                atom_p->connect12_res[j_connect] = j_res;
                                connect_found = 1;
                                break;
                            }
                        }
                        if (!connect_found) {
                            char err_msg1[MAXCHAR_LINE];
                            char err_msg2[MAXCHAR_LINE];
                            sprintf(err_msg1,"   Error! get_connect12(): connectivity of atom \"%s %s %c%04d\" is not complete",
                            atom_p->name, conf_p->confName, res_p->chainID, res_p->resSeq);
                            sprintf(err_msg2,"          get_connect12(): atom %s in the same residue is not found", connect.atom[j_connect].name);
                            if (param_exist(err_msg1, "", "")) {
                                if (param_exist(err_msg2, "", "")) {
                                    continue;
                                }
                            }

                            param_sav(err_msg1, "", "", "", 0);
                            param_sav(err_msg2, "", "", "", 0);
                                
                            debug_fp = fopen(env.debug_log, "a");
                            fprintf(debug_fp,"%s\n",err_msg1);
                            fprintf(debug_fp,"%s\n",err_msg2);
                            fclose(debug_fp);
                            err++;
                        }
                    }
                }
                else {        /* Not in the same reside. */
                    j_res = i_res + connect.atom[j_connect].res_offset;
                    if (j_res < 0 || j_res >= prot.n_res) {     /* j_res is out of residue list */
                        char err_msg1[MAXCHAR_LINE];
                        char err_msg2[MAXCHAR_LINE];
                        sprintf(err_msg1, "   Error! get_connect12(): connectivity of atom \"%s %s %c%04d\" is not complete:\n",
                        atom_p->name, res_p->resName, res_p->chainID, res_p->resSeq);
                        sprintf(err_msg2, "          get_connect12(): atom %s of residue i%+d is not found \n", connect.atom[j_connect].name,connect.atom[j_connect].res_offset);
                        if (param_exist(err_msg1, "", "")) {
                            if (param_exist(err_msg2, "", "")) {
                                continue;
                            }
                        }
                        
                        param_sav(err_msg1, "", "", "", 0);
                        param_sav(err_msg2, "", "", "", 0);
                        
                        debug_fp = fopen(env.debug_log, "a");
                        fprintf(debug_fp,"%s\n",err_msg1);
                        fprintf(debug_fp,"%s\n",err_msg2);
                        fclose(debug_fp);

                        err++;
                        continue;
                    }
                    
                    for (j_conf = 0; j_conf < prot.res[j_res].n_conf; j_conf++) {
                        if ( !param_get("IATOM", prot.res[j_res].conf[j_conf].confName, connect.atom[j_connect].name, &j_atom) ) {
                            if (!prot.res[j_res].conf[j_conf].atom[j_atom].on) {
                                atom_p->connect12[j_connect] = &prot.res[j_res].conf[j_conf].atom[j_atom];
                                atom_p->connect12_res[j_connect] = j_res;
                                connect_found = 1;

                                debug_fp = fopen(env.debug_log, "a");
                                fprintf(debug_fp, "   WARNING! get_connect12(): atom %s of residue i%+d is in connectivity list of atom \"%s %s %c%04d\", but it's not on\n",
                                connect.atom[j_connect].name,connect.atom[j_connect].res_offset,atom_p->name, res_p->resName, res_p->chainID, res_p->resSeq);
                                fclose(debug_fp);
                                break;
                            }
                            if (BOND_THR > dvv(atom_p->xyz, prot.res[j_res].conf[j_conf].atom[j_atom].xyz)) {
                                atom_p->connect12[j_connect] = &prot.res[j_res].conf[j_conf].atom[j_atom];
                                atom_p->connect12_res[j_connect] = j_res;
                                connect_found = 1;
                                break;
                            }
                        }
                    }
                    if (!connect_found) {
                        char err_msg1[MAXCHAR_LINE];
                        char err_msg2[MAXCHAR_LINE];
                        sprintf(err_msg1, "   Error! get_connect12(): connectivity of atom \"%s %s %c%04d\" is not complete:\n",
                        atom_p->name, res_p->resName, res_p->chainID, res_p->resSeq);
                        sprintf(err_msg2, "          get_connect12(): atom %s of residue i%+d is not found \n", connect.atom[j_connect].name,connect.atom[j_connect].res_offset);

                        if (param_exist(err_msg1, "", "")) {
                            if (param_exist(err_msg2, "", "")) {
                                continue;
                            }
                        }
                        
                        param_sav(err_msg1, "", "", "", 0);
                        param_sav(err_msg2, "", "", "", 0);
                        
                        debug_fp = fopen(env.debug_log, "a");
                        fprintf(debug_fp,"%s\n",err_msg1);
                        fprintf(debug_fp,"%s\n",err_msg2);
                        fclose(debug_fp);
                        
                        err++;
                    }
                }
            }
            else {      /* Ligand type connectivity */
                if (lig_treated) continue;
                
                /* Loop over all atoms and find all atom within bond threshold */
                n_ligs = 0;
                for (j_res=0; j_res < prot.n_res; j_res++) {
                    if (j_res == i_res) continue;
                    if ( prot.res[j_res].n_conf>2 ) n_conf=2;
                    else n_conf = prot.res[j_res].n_conf;
                    for (j_conf=0; j_conf < n_conf; j_conf++) {             /* searching is localized in backbone and first side chain */
                        for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
                            jatom_p = &prot.res[j_res].conf[j_conf].atom[j_atom];
                            if (!jatom_p->on) continue;
                            
                            if (BOND_THR > dvv(atom_p->xyz, jatom_p->xyz)) {
                                n_ligs++;
                                ligs[n_ligs-1] = jatom_p;
                                ligs_res[n_ligs-1] = j_res;
                            }
                        }
                    }
                }
                
                /* Go over connectivity list and put connected atoms in for those atom name is defined */
                for (k_connect=j_connect; k_connect<connect.n; k_connect++) {
                    if (!connect.atom[k_connect].ligand) continue;
                    if (atom_p->connect12[k_connect]) continue;
                    if (strchr(connect.atom[k_connect].name, '?')) continue; /* If atom name in parameter is a '?', then go to next one */
                    
                    for (k_lig=0; k_lig<n_ligs; k_lig++) {
                        if (strcmp(ligs[k_lig]->name, connect.atom[k_connect].name)) continue;
                        
                        /* Check if ligand atom is already in the connectivity list, this is for the case one atom connect to over one atom with same name */
                        for (l_connect=j_connect; l_connect<connect.n; l_connect++) {
                            if (ligs[k_lig] == atom_p->connect12[l_connect]) break;
                        }
                        if (l_connect < connect.n) continue;
                        
                        /* Adding ligs[k_lig] into list */
                        atom_p->connect12[k_connect] = ligs[k_lig];
                        atom_p->connect12_res[k_connect] = ligs_res[k_lig];
                        break;
                    }
                }
                
                /* Then go over connectivity list again and put atoms in for those atom name is not defined */
                for (k_connect=j_connect; k_connect<connect.n; k_connect++) {
                    int k_lig_add;
                    float lig_dist;
                    if (!connect.atom[k_connect].ligand) continue;
                    if (atom_p->connect12[k_connect]) continue;
                    if (!strchr(connect.atom[k_connect].name, '?')) continue;   /* If atom name in parameter is not a '?', then go to next one */
                    
                    k_lig_add = -1;
                    lig_dist = BOND_THR;
                    for (k_lig=0; k_lig<n_ligs; k_lig++) {
                        if (ligs[k_lig]->name[1] == 'H') continue; /* If it's a proton, then not considered */
                        
                        /* Check if ligand atom is already in the connectivity list, this is for the case one atom connect to over one atom with same name */
                        for (l_connect=j_connect; l_connect<connect.n; l_connect++) {
                            if (ligs[k_lig] == atom_p->connect12[l_connect]) break;
                        }
                        if (l_connect < connect.n) continue;
                        
                        if (dvv(atom_p->xyz,ligs[k_lig]->xyz) < lig_dist) {
                            lig_dist = dvv(atom_p->xyz,ligs[k_lig]->xyz);
                            k_lig_add = k_lig;
                        }
                    }
                    
                    /* Adding ligs[k_lig] into list */
                    if (k_lig_add >= 0) {
                        atom_p->connect12[k_connect] = ligs[k_lig_add];
                        atom_p->connect12_res[k_connect] = ligs_res[k_lig_add];
                    }
                }
                
                /* Go over connectivity list to check if there is empty slot */
                for (k_connect=j_connect; k_connect<connect.n; k_connect++) {
                    if (!connect.atom[k_connect].ligand) continue;
                    if (atom_p->connect12[k_connect]) continue;

                    char err_msg1[MAXCHAR_LINE];
                    sprintf(err_msg1, "   Warning! get_connect12(): An empty ligand connectivity slot found for atom %s in residue %s %d to atom %s\n", atom_p->name, res_p->resName, res_p->resSeq, connect.atom[k_connect].name);
                    
                    if (param_exist(err_msg1, "", "")) {
                        continue;
                    }
                    
                    param_sav(err_msg1, "", "", "", 0);
                    
                    debug_fp = fopen(env.debug_log, "a");
                    fprintf(debug_fp,"%s\n",err_msg1);
                    fclose(debug_fp);
                }
                
                /* Check if all ligand atoms are in the connectivity list */
                for (k_lig=0; k_lig<n_ligs; k_lig++) {
                    if (ligs[k_lig]->name[1] == 'H') continue; /* If it's a proton, then not considered */
                    
                    /* Check if ligand atom is already in the connectivity list, this is for the case one atom connect to over one atom with same name */
                    for (l_connect=j_connect; l_connect<connect.n; l_connect++) {
                        if (ligs[k_lig] == atom_p->connect12[l_connect]) break;
                    }
                    if (l_connect < connect.n) continue;

                    char err_msg1[MAXCHAR_LINE];
                    sprintf(err_msg1, "   Warning! get_connect12(): An atom (%s in residue %s %d) within bond threshold to atom %s in residue %s %d is not put in the connectivity list \n",ligs[k_lig]->name, ligs[k_lig]->resName, ligs[k_lig]->resSeq, atom_p->name, res_p->resName, res_p->resSeq);
                    
                    if (param_exist(err_msg1, "", "")) {
                        continue;
                    }
                    
                    param_sav(err_msg1, "", "", "", 0);
                    
                    debug_fp = fopen(env.debug_log, "a");
                    fprintf(debug_fp,"%s\n",err_msg1);
                    fclose(debug_fp);
                }
                
                lig_treated = 1;
                /* END of treating ligand type */
            }
            
            /* If connecting to a dummy atom, it could be the case connecting to NTR or CTR */
            if (!atom_p->connect12[j_connect]) continue;
            if ( atom_p->connect12[j_connect]->on) continue;
            if (strcmp(connect.atom[j_connect].name, " CA ") && strcmp(connect.atom[j_connect].name, " C  ")) continue;
            //printf("   Debugging! residue %s%4d,conformer %s, on=%d\n", res_p->resName,res_p->resSeq,conf_p->confName,atom_p->connect12[j_connect]->on);
            connect_found = 0;
            for (j_res=0; j_res<prot.n_res; j_res++) {
                if (strcmp(prot.res[j_res].resName, "NTR"))
                    if (strcmp(prot.res[j_res].resName, "NTG"))
                        if (strcmp(prot.res[j_res].resName, "CTR")) continue;
                for (j_conf=0; j_conf<prot.res[j_res].n_conf; j_conf++) {
                    for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
                        jatom_p = &prot.res[j_res].conf[j_conf].atom[j_atom];
                        if (!jatom_p->on) continue;
                        
                        if (strcmp(connect.atom[j_connect].name, jatom_p->name)) continue;
                        if (BOND_THR > dvv(atom_p->xyz, jatom_p->xyz)) {
                            atom_p->connect12[j_connect] = jatom_p;
                            atom_p->connect12_res[j_connect] = j_res;
                            connect_found = 1;
                            break;
                        }
                    }
                    if (connect_found) break;
                }
                if (connect_found) break;
            }
        }
        
        /* Move NULL pointer to the end of the array */
        n_connect = connect.n;
        for (j_connect=0; j_connect < n_connect-1; j_connect++) {
            if (!atom_p->connect12[j_connect]) {
                for (k_connect=j_connect; k_connect<n_connect-1; k_connect++) {
                    atom_p->connect12[k_connect] = atom_p->connect12[k_connect+1];
                    atom_p->connect12_res[k_connect] = atom_p->connect12_res[k_connect+1];
                }
                atom_p->connect12[n_connect-1] = NULL;
                atom_p->connect12_res[n_connect-1] = 0;
                
                j_connect--;
                n_connect--;
            }
        }
    }
    
    return 0;
}

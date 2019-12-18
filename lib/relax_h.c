/* This subroutine is a little primitive. It was developed before a lot of other subroutines.
There will be a lot of improvement need to be made. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"
#define  NGH_THR        5
#define  NGH_THR2       NGH_THR*NGH_THR
#define  N_ELEM_MAX     60
#define  NSTATE_MAX     1000000

typedef struct RELAX_ATOM_STRUCT{
    ATOM   *atom_p;
    ATOM   *torsion_atom[3];
    TORS   tors;
    int    i_res_prot;
    int    i_conf_res;
    int    i_atom_conf;
    int    move;  /* 0 is fixed, 1 rotate around an axis, 2 rotate around an origin, 3 moved by other atoms */
    int    *counter;
    VECTOR xyz_back;
    VECTOR xyz_new;
    
    VECTOR frc;
    LINE   axis;
    VECTOR origin;
    VECTOR torque;
    
    int    n_slaved;
    struct RELAX_ATOM_STRUCT *slaved[8];
    
    int    k_relax;
} RELAX_ATOM;

int counter;

//int  n_elem;
//char elem[N_ELEM_MAX];
//float **C6_matrix;
//float **C12_matrix;
float **factor_matrix;
float **pair_vdw;

int relax(int n_relax, RELAX_ATOM **relax_atoms, PROT prot);
int water_orient(VECTOR v0, VECTOR v1, VECTOR v2, float *theta, float *phi, float *psi);
extern int rm_comment(char *target, char *str);

int relax_h(PROT prot)
{
    int         i_res, i_subres, i_conf, i_atom;
    int         j_res, j_subres, j_conf, j_atom, n_conf;
    int         ic, jc;
    ATOM        *atom_p0, *atom_p1, *atom_p2, *atom_p3;
    RELAX_ATOM   *all_atoms=NULL;
    int         ka,na;
    int         loop_states;
    FILE        *progress_fp;
    //int         i_elem,j_elem;
    char        sbuff[MAXCHAR_LINE];
    TORS        tors;

    counter = 0;
    
    assign_crg(prot);
    get_connect12(prot);
    //n_elem = 0;
    //memset(elem,0,N_ELEM_MAX);
    
    prot.nc = 0;
    for (i_res=0; i_res<prot.n_res; i_res++) {
        int cal_vdw;
        prot.res[i_res].i_res_prot = i_res;
        if (!param_get("CAL_VDW",prot.res[i_res].resName, "", &cal_vdw)) {
            prot.res[i_res].cal_vdw = cal_vdw;
            if (!cal_vdw) continue;
        }
        else prot.res[i_res].cal_vdw = 1;
        
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            prot.res[i_res].conf[i_conf].i_conf_res = i_conf;
            prot.res[i_res].conf[i_conf].counter = 0;
            if (i_conf) {
                prot.nc++;
                prot.res[i_res].conf[i_conf].i_conf_prot = prot.nc-1;
            }
            /*
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
            */
        }
    }
    
    /*
    C6_matrix  = malloc(n_elem*sizeof(void *));
    C12_matrix = malloc(n_elem*sizeof(void *));
    for (i_elem=0;i_elem<n_elem;i_elem++) {
        C6_matrix[i_elem]  = malloc(n_elem*sizeof(float));
        C12_matrix[i_elem] = malloc(n_elem*sizeof(float));
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
            //if (j_elem<i_elem) printf("%f,%f\n",C6[i_elem][j_elem]-C6[j_elem][i_elem],C12[i_elem][j_elem]-C12[j_elem][i_elem]);
        }
    }
    */
    
    /* Decide which residues need optimize hydrogen */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        if (!strcmp(prot.res[i_res].resName, "HOH")) prot.res[i_res].opt_hyd = 1;
        else {
            prot.res[i_res].opt_hyd = 0;
            for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
                CONF *conf_p = &prot.res[i_res].conf[i_conf];
                for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                    int i_connect, n_connect;
                    ATOM *atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                    if (!atom_p->on) continue;
                    if ( torsion_atoms(conf_p, i_atom, &atom_p0, &atom_p1, &atom_p2, &atom_p3, &tors, 1) ) continue;
                    if (!tors.opt_hyd) continue;
                    
                    /* check if atom_p1 is connected to more than one atom besides atom_p0 */
                    // do optimization only if there are no more than 2 atoms connected to atom_p1, usually this means it's OH
                    n_connect = 0;
                    for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
                        if (!atom_p1->connect12[i_connect]) break;
                        if (!atom_p1->connect12[i_connect]->on) continue;
                        if (strcmp(atom_p1->connect12[i_connect]->name, atom_p0->name)) {
                            n_connect++;
                        }
                    }
                    if (n_connect > 1) continue;
                    
                    prot.res[i_res].opt_hyd = 1;
                    break;
                }
                if (prot.res[i_res].opt_hyd) break;
            }
        }
    }
    
    /* Group conformers of a residues into sub residues */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        RES *res_p;
        res_p = &prot.res[i_res];
        res_p->n_subres = 0;
        
        for (i_conf=1; i_conf<res_p->n_conf; i_conf++) {
            for (i_subres = res_p->n_subres-1; i_subres>=0; i_subres--) {
                if (strcmp(res_p->conf[i_conf].confName, res_p->subres[i_subres].conf[0]->confName)) continue;
                if (!cmp_conf_hv(res_p->conf[i_conf], *res_p->subres[i_subres].conf[0], 0.05)) break;
            }
            
            if (i_subres == -1) {
                res_p->n_subres++;
                res_p->subres = realloc(res_p->subres,res_p->n_subres*sizeof(SUBRES));
                i_subres = res_p->n_subres-1;
                memset(&res_p->subres[i_subres],0,sizeof(SUBRES));
            }
            
            res_p->subres[i_subres].n_conf++;
            n_conf = res_p->subres[i_subres].n_conf;
            res_p->subres[i_subres].conf = realloc(res_p->subres[i_subres].conf,n_conf*sizeof(void *));
            res_p->subres[i_subres].conf[n_conf-1] = &res_p->conf[i_conf];
            res_p->conf[i_conf].i_subres_res = i_subres;
            res_p->conf[i_conf].i_conf_subres = n_conf-1;
        }
    }
    
    /* heavy atom vdw */
    /** calculate a prot.nc * prot.nc matrix for the vdw interaction
     * only the interactions between conformers of the first conformer of sub residues are calculated
     */
    pair_vdw = malloc(prot.nc * sizeof(float *));
    memset(pair_vdw, 0, prot.nc * sizeof(float *));
    for (i_res=0; i_res<prot.n_res; i_res++) {
        if (!prot.res[i_res].cal_vdw) continue;
        for (i_subres=0; i_subres<prot.res[i_res].n_subres; i_subres++) {
            ic = prot.res[i_res].subres[i_subres].conf[0]->i_conf_prot;
            i_conf = prot.res[i_res].subres[i_subres].conf[0]->i_conf_res;

            pair_vdw[ic] = malloc(prot.nc * sizeof(float));
            for (j_res=0; j_res<prot.n_res; j_res++) {
                if (!prot.res[j_res].cal_vdw) continue;
                for (j_subres=0; j_subres<prot.res[j_res].n_subres; j_subres++) {
                    jc = prot.res[j_res].subres[j_subres].conf[0]->i_conf_prot;
                    j_conf = prot.res[j_res].subres[j_subres].conf[0]->i_conf_res;
                    pair_vdw[ic][jc] = vdw_conf_hv(i_res,i_conf,j_res,j_conf,prot);
                }
            }
        }
    }
    


    /* Find residue neighbor list for each residue */
    //printf("   Making residue-residue neighbor list...\n");
    for (i_res=0; i_res<prot.n_res; i_res++) {
        RES *ires_p;
        if (!prot.res[i_res].opt_hyd) continue;
        
        ires_p = &prot.res[i_res];

        // first neighbor is the residue itself.
        ires_p->n_ngh = 1;
        ires_p->ngh = malloc(sizeof(RES *));
        ires_p->ngh[0] = &prot.res[i_res];
        
        /* Get residue neighbor list */
        for (j_res=0; j_res<prot.n_res; j_res++) {
            RES *jres_p;
            int ngh_found, clash, all_clash;
            float d2;
            
            if (j_res == i_res) continue;
            if (!prot.res[j_res].cal_vdw) continue;
            jres_p = &prot.res[j_res];
            
            ngh_found = 0;
            all_clash = 1;
            for (i_subres=0; i_subres<ires_p->n_subres; i_subres++) {
                for (j_subres=0; j_subres<jres_p->n_subres; j_subres++) {
                    ic = prot.res[i_res].subres[i_subres].conf[0]->i_conf_prot;
                    jc = prot.res[j_res].subres[j_subres].conf[0]->i_conf_prot;
                    i_conf = prot.res[i_res].subres[i_subres].conf[0]->i_conf_res;
                    j_conf = prot.res[j_res].subres[j_subres].conf[0]->i_conf_res;
                    
                    if (pair_vdw[ic][jc] > env.relax_clash_thr) clash = 1;
                    else clash = 0;


                    if (!clash && prot.res[i_res].subres[i_subres].conf[0]->n_atom && prot.res[j_res].subres[j_subres].conf[0]->n_atom) all_clash = 0;
                    
                    if (ngh_found) continue;
                    for (i_atom=0; i_atom<ires_p->conf[i_conf].n_atom; i_atom++) {
                        ATOM *iatom_p = &ires_p->conf[i_conf].atom[i_atom];
                        if (!iatom_p->on) continue;
                        for (j_atom=0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) {
                            ATOM *jatom_p = &jres_p->conf[j_conf].atom[j_atom];
                            if (!jatom_p->on) continue;
                            
                            d2 = ddvv(iatom_p->xyz,jatom_p->xyz);
                            if (d2 < NGH_THR2) {
                                ngh_found = 1;
                                break;
                            }
                        }
                        if (ngh_found) break;
                    }
                }
            }
            
            // if two conformers always clash, they are not in each other's neighbor list.
            if (all_clash && (ires_p->n_conf-1) && (jres_p->n_conf-1)) {
                FILE *debug_fp = fopen(env.debug_log,"a");
                fprintf(debug_fp, " Residue %s %c%4d always clash with Residue %s %c%4d\n", 
                prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,
                prot.res[j_res].resName,prot.res[j_res].chainID,prot.res[j_res].resSeq);
                fclose(debug_fp);
                continue;
            }
            else if (ngh_found) {
                ires_p->n_ngh++;
                ires_p->ngh = realloc(ires_p->ngh, ires_p->n_ngh*sizeof(RES *));
                ires_p->ngh[ires_p->n_ngh-1] = jres_p;
            }
        }
        
        /*
        printf("%s %c%04d", ires_p->resName, ires_p->chainID, ires_p->resSeq);
        for (i_ngh=0;i_ngh<ires_p->n_ngh;i_ngh++) {
            printf("|%c%4d", ires_p->ngh[i_ngh]->chainID, ires_p->ngh[i_ngh]->resSeq);
        }
        printf("\n");
        */
    }
    
    /* Define configurations for relaxation */
    progress_fp = fopen(env.progress_log,"a");
    fprintf(progress_fp, "   Configure atoms for hydrogen relaxation...\n");
    fclose(progress_fp);


    // get all the atoms invovled in the relaxation, exclude the residues that don't calculating the vdw.
    na = 0;
    for (i_res=0; i_res<prot.n_res; i_res++) {
        int cal_vdw;
        if ( !param_get("CAL_VDW",prot.res[i_res].resName, "", &cal_vdw) ) {
            if (!cal_vdw) continue;
        }
        // all the atoms of all the conformers are chosen to relax.
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                na++;
                ka = na-1;
                all_atoms = realloc(all_atoms, na*sizeof(RELAX_ATOM));
                memset(&all_atoms[ka], 0, sizeof(RELAX_ATOM));
                prot.res[i_res].conf[i_conf].atom[i_atom].i_atom_prot = ka;
                
                all_atoms[ka].atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                all_atoms[ka].xyz_back = prot.res[i_res].conf[i_conf].atom[i_atom].xyz;
                all_atoms[ka].i_res_prot  = i_res;
                all_atoms[ka].i_conf_res = i_conf;
                all_atoms[ka].i_atom_conf = i_atom;
                all_atoms[ka].move   = 0;
            }
        }
    }
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        if(!prot.res[i_res].opt_hyd) continue;
        for (i_subres=0; i_subres<prot.res[i_res].n_subres; i_subres++) {
            for (i_conf=0; i_conf<prot.res[i_res].subres[i_subres].n_conf; i_conf++) {
                CONF *conf_p = prot.res[i_res].subres[i_subres].conf[i_conf];
                for (i_atom=0; i_atom<conf_p->n_atom; i_atom++) {
                    ATOM *atom_p;
                    if (!conf_p->atom[i_atom].on) continue;
                    atom_p = &conf_p->atom[i_atom];
                    ka = atom_p->i_atom_prot;

                    // deal with waters
                    if (!strcmp(prot.res[i_res].resName,"HOH")) {
                        if (atom_p->name[1]=='O') {
                            all_atoms[ka].move = 2;
                            //if (!i_conf) continue;
                            if (prot.res[i_res].n_hyd_pos) {
                                //printf("water prot.res[i_res].n_hyd_pos %d\n",prot.res[i_res].n_hyd_pos);
                                all_atoms[ka].counter = malloc(prot.res[i_res].n_hyd_pos * prot.res[i_res].n_hyd_pos * prot.res[i_res].n_hyd_pos * sizeof(int));
                                memset(all_atoms[ka].counter,0,prot.res[i_res].n_hyd_pos * prot.res[i_res].n_hyd_pos * prot.res[i_res].n_hyd_pos * sizeof(int));
                            }

                            all_atoms[ka].n_slaved = 0;
                            for (j_atom=0; j_atom <conf_p->n_atom; j_atom++) {
                                if (conf_p->atom[j_atom].name[1]=='O') continue;
                                all_atoms[ka].n_slaved++;
                                all_atoms[conf_p->atom[j_atom].i_atom_prot].origin = atom_p->xyz;
                                //if (i_atom == j_atom) continue;
                                all_atoms[ka].slaved[all_atoms[ka].n_slaved-1] = &all_atoms[conf_p->atom[j_atom].i_atom_prot];
                            }
                        }
                        else {
                            all_atoms[ka].move=3;   // H atoms..
                        }
                    }
                    else {        // amino acids
                        int i_connect;
                        if (torsion_atoms(conf_p, all_atoms[ka].i_atom_conf, &atom_p0, &atom_p1, &atom_p2, &atom_p3, &tors, 1)) continue;
                        //printf("%s opt=%d\n",atom_p0->name,tors.opt_hyd);
                        if (!tors.opt_hyd) continue;
                        int n_connect;

                        // atom_p0 is H, atom_p1 is O, in the 4 atoms involved in torsion
                        for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
                            if (!atom_p1->connect12[i_connect]) break;
                            if (!atom_p1->connect12[i_connect]->on) continue;
                            if (strcmp(atom_p1->connect12[i_connect]->name, atom_p0->name)) {
                                n_connect ++;
                            }
                        }
                        if (n_connect > 1) continue;
                        //if (!atom_p1->connect12[i_connect]) continue;
                        
                        all_atoms[ka].move = 1;        // hydroxyl
                        all_atoms[ka].tors = tors;
                        all_atoms[ka].axis = line_2v(atom_p2->xyz, atom_p1->xyz);
                        all_atoms[ka].torsion_atom[0] = atom_p1;
                        all_atoms[ka].torsion_atom[1] = atom_p2;
                        all_atoms[ka].torsion_atom[2] = atom_p3;
                        
                        if (!i_conf) {
                        	if (prot.res[i_res].n_hyd_pos) {
                        		//printf("prot.res[i_res].n_hyd_pos %d\n",prot.res[i_res].n_hyd_pos);
                        		all_atoms[ka].counter = malloc(prot.res[i_res].n_hyd_pos*sizeof(int));
                        		memset(all_atoms[ka].counter,0,prot.res[i_res].n_hyd_pos*sizeof(int));
                        	}
                        }
                    }
                }
                //printf("residue %s %c%4d atom %s ka=%d move %d\n", prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,all_atoms[ka].atom_p->name,ka,all_atoms[ka].move);
            }
        }
    }
    printf("   Total number of atoms %d\n", na);


    /* Relaxation */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        int i_ngh;
        long n_state;
        RES *ires_p = &prot.res[i_res];
        int n_trial;
        int j_ngh, n_atom;
        
        if (!prot.res[i_res].opt_hyd) continue;
        if (!prot.res[i_res].n_hyd_pos) continue;
        
        /* Calculate number of sub-microstates */
        // number of microstates depends on the number of sub residues of its neighbors
        n_state = 1;
        for (i_ngh=0;i_ngh<ires_p->n_ngh;i_ngh++) {
            if (ires_p->ngh[i_ngh]->n_subres<2) continue;
            n_state = n_state*(ires_p->ngh[i_ngh]->n_subres);
            
            if (n_state >= NSTATE_MAX) {
                n_state = NSTATE_MAX;
                break;
            }
        }

        progress_fp = fopen(env.progress_log,"a");
        if (n_state == NSTATE_MAX) {
            fprintf(progress_fp, "   %s %c%04d (n_subres=%2d,n_conf=%3d) has %3d neighbors and over %10ld states... ", prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].n_subres,prot.res[i_res].n_conf,prot.res[i_res].n_ngh,n_state);
        }
        else {
            fprintf(progress_fp, "   %s %c%04d (n_subres=%2d,n_conf=%3d) has %3d neighbors and      %10ld states... ", prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].n_subres,prot.res[i_res].n_conf,prot.res[i_res].n_ngh,n_state);
        }

        if (n_state < env.relax_nstates) {
            loop_states = 1;
            fprintf(progress_fp,"Loop over %5ld states...\n", n_state);
        }
        else {
            loop_states = 0;
            fprintf(progress_fp,"Randomly pick %5d states...\n", env.relax_nstates);
        }
        fclose(progress_fp);
        
        n_atom = 0;
        for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
            for (i_conf=0; i_conf<prot.res[i_res].ngh[i_ngh]->n_conf; i_conf++) {
                for (i_atom = 0; i_atom < prot.res[i_res].ngh[i_ngh]->conf[i_conf].n_atom; i_atom++) {
                    ATOM *iatom_p;
                    iatom_p = &prot.res[i_res].ngh[i_ngh]->conf[i_conf].atom[i_atom];
                    if (!iatom_p->on) continue;
                    n_atom++;
                    // all the atoms of neighbor residues are indexed for relaxation of this residue.
                    all_atoms[iatom_p->i_atom_prot].k_relax = n_atom-1;
                }
            }
        }

        // Get the factor matrix between every two atoms of the neighbor residues, which are to be relaxed.
        factor_matrix = malloc(n_atom * sizeof(void*));
        memset(factor_matrix, 0, n_atom * sizeof(void*));
        for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
            for (i_conf=0; i_conf<prot.res[i_res].ngh[i_ngh]->n_conf; i_conf++) {
                for (i_atom = 0; i_atom < prot.res[i_res].ngh[i_ngh]->conf[i_conf].n_atom; i_atom++) {
                    ATOM *iatom_p;
                    int   n_connect13;
                    int   n_connect14;
                    ATOM   *connect13[MAX_CONNECTED2],    *connect14[MAX_CONNECTED3];
                    int    connect13_res[MAX_CONNECTED2], connect14_res[MAX_CONNECTED3];
                    int    iconnect, jconnect, kconnect;
                    
                    iatom_p = &prot.res[i_res].ngh[i_ngh]->conf[i_conf].atom[i_atom];
                    if (!iatom_p->on) continue;
                    if (!all_atoms[iatom_p->i_atom_prot].move) continue;
                    factor_matrix[all_atoms[iatom_p->i_atom_prot].k_relax] = malloc(n_atom * sizeof(float));
                    
                    n_connect13 = 0;
                    n_connect14 = 0;
                    memset(connect13, 0, MAX_CONNECTED2*sizeof(void *));
                    memset(connect13_res,0,MAX_CONNECTED2*sizeof(int));
                    memset(connect14, 0, MAX_CONNECTED3*sizeof(void *));
                    memset(connect14_res,0,MAX_CONNECTED3*sizeof(int));

                    /** connect13 includes the atoms directly connecting to this atom
                     * and the atoms indirectly connecting to this atom by one atom.
                     * connect14 only includes the atoms that can connect to this atom through two atoms.
                     */
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
                    
                    for (j_ngh=0; j_ngh<prot.res[i_res].n_ngh; j_ngh++) {
                        for (j_conf=0; j_conf<prot.res[i_res].ngh[j_ngh]->n_conf; j_conf++) {
                            for (j_atom = 0; j_atom < prot.res[i_res].ngh[j_ngh]->conf[j_conf].n_atom; j_atom++) {
                                ATOM *jatom_p;
                                jatom_p = &prot.res[i_res].ngh[j_ngh]->conf[j_conf].atom[j_atom];
                                if (!jatom_p->on) continue;
                                
                                if (iatom_p  == jatom_p) {
                                	factor_matrix[all_atoms[iatom_p->i_atom_prot].k_relax][all_atoms[jatom_p->i_atom_prot].k_relax]=0;
                                	continue;
                                }
                                factor_matrix[all_atoms[iatom_p->i_atom_prot].k_relax][all_atoms[jatom_p->i_atom_prot].k_relax] = 1.;
                                //printf("HERE1,14=%d,13=%d,%d,%d,%f\n",n_connect14,n_connect13,all_atoms[iatom_p->i_atom_prot].k_relax,all_atoms[jatom_p->i_atom_prot].k_relax,factor_matrix[all_atoms[iatom_p->i_atom_prot].k_relax][all_atoms[jatom_p->i_atom_prot].k_relax]);
                                
                                for (iconnect = 0; iconnect < n_connect14; iconnect++) {
                                    if (all_atoms[jatom_p->i_atom_prot].i_res_prot == connect14_res[iconnect] && !strcmp(all_atoms[jatom_p->i_atom_prot].atom_p->name,connect14[iconnect]->name)) break;
                                }
                                if (iconnect < n_connect14) {
                                    factor_matrix[all_atoms[iatom_p->i_atom_prot].k_relax][all_atoms[jatom_p->i_atom_prot].k_relax] = env.factor_14lj;
                                }
                                //printf("iconnect= %d,14=%d\n",iconnect,n_connect14);
                                
                                for (iconnect = 0; iconnect < n_connect13; iconnect++) {
                                    if (all_atoms[jatom_p->i_atom_prot].i_res_prot == connect13_res[iconnect] && !strcmp(all_atoms[jatom_p->i_atom_prot].atom_p->name,connect13[iconnect]->name)) break;
                                }
                                if (iconnect < n_connect13) {
                                    factor_matrix[all_atoms[iatom_p->i_atom_prot].k_relax][all_atoms[jatom_p->i_atom_prot].k_relax] = 0.;
                                }
                                //printf("iconnect= %d,13=%d\n",iconnect,n_connect13);
                                //printf("HERE2,%d,%d,%f\n",all_atoms[iatom_p->i_atom_prot].k_relax,all_atoms[jatom_p->i_atom_prot].k_relax,factor_matrix[all_atoms[jatom_p->i_atom_prot].k_relax][all_atoms[jatom_p->i_atom_prot].k_relax]);
                            }
                        }
                    }
                }
            }
        }


        /* Initialize sub-microstates */
        n_trial=0;
        prot.res[i_res].i_subres_on = -1;
        for (i_ngh=1; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
            prot.res[i_res].ngh[i_ngh]->i_subres_on = 0;
        }
        
        /* Pick local microstates to relax h */
        while (1) {
            int n_relax,convg;
            float diff;
            int   i_bin,j_bin,k_bin,hyd_pos;
            RELAX_ATOM **relax_atoms = NULL;
            float  vdw_start, coulomb_start,vdw_end, coulomb_end;
            n_trial++;
            
            //printf("n_trial %d\n",n_trial);
        	//printf("loop states %d\n", loop_states);

            /* Get a sub-residue microstate */
            if (!loop_states) {
                int clash,clash_counter;
                if (n_trial>env.relax_nstates) break;
                clash = 1;
                clash_counter = 0;
                while (clash) {
                    clash_counter++;
                    for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
                        int k_rand;
                        if (!prot.res[i_res].ngh[i_ngh]->n_subres) continue;
                        k_rand = rand();
                        //printf("j_conf_on %d,%d\n",res_extra[j_res].n_subres,k_rand);
                        //prot.res[i_res].ngh[i_ngh]->i_subres_on = prot.res[i_res].ngh[i_ngh]->n_subres * ((float)k_rand/(RAND_MAX+1.0));
                        prot.res[i_res].ngh[i_ngh]->i_subres_on = k_rand % prot.res[i_res].ngh[i_ngh]->n_subres;
                        
                        //printf("%d/%d\n",prot.res[i_res].ngh[i_ngh]->i_subres_on, prot.res[i_res].ngh[i_ngh]->n_subres);
                    }
                    clash = 0;
                    for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
                        if (!prot.res[i_res].ngh[i_ngh]->n_subres) continue;
                        for (j_ngh=i_ngh+1; j_ngh<prot.res[i_res].n_ngh; j_ngh++) {
                            if (!prot.res[i_res].ngh[j_ngh]->n_subres) continue;
                            ic = prot.res[i_res].ngh[i_ngh]->subres[prot.res[i_res].ngh[i_ngh]->i_subres_on].conf[0]->i_conf_prot;
                            jc = prot.res[i_res].ngh[j_ngh]->subres[prot.res[i_res].ngh[j_ngh]->i_subres_on].conf[0]->i_conf_prot;
                            
                            /*
                            //debug_fp = fopen(env.debug_log,"a");
                            //fprintf(debug_fp,"i_ngh=%d,j_ngh=%d,n_trial=%d, e=%.3f\n",i_ngh,j_ngh,n_trial,vdw_conf_hv(ires,iconf,jres,jconf,prot));
                            //fclose(debug_fp);
                            */
                            if (pair_vdw[ic][jc] > env.relax_clash_thr) {
                                int flip, k_rand;
                                clash = 1;
                                if ( rand()/(RAND_MAX+1.0) < 0.5 ) flip = i_ngh;
                                else flip = j_ngh;
                                
                                k_rand = rand();
                                //printf("j_conf_on %d,%d\n",res_extra[j_res].n_subres,k_rand);
                                prot.res[i_res].ngh[flip]->i_subres_on = prot.res[i_res].ngh[flip]->n_subres * ((float)k_rand/(RAND_MAX+1.0));
                                break;
                            }
                        }
                        if (clash) break;
                    }
                    if (clash_counter>=1000) break;
                }
                //printf("clash = %d\n",clash_counter);
            }
            else {
                int clash;
                i_ngh = 0;
                prot.res[i_res].ngh[i_ngh]->i_subres_on++;
                while (prot.res[i_res].ngh[i_ngh]->i_subres_on >= prot.res[i_res].ngh[i_ngh]->n_subres) {
                    prot.res[i_res].ngh[i_ngh]->i_subres_on = 0;
                    i_ngh++;
                    if (i_ngh == prot.res[i_res].n_ngh) break;
                    if (!prot.res[i_res].ngh[i_ngh]->n_subres) continue;
                    prot.res[i_res].ngh[i_ngh]->i_subres_on++;
                    //printf("n_trial%d,j_ngh:%d,%d,%d,subres%d,%d\n",n_trial,j_ngh,res_extra[i_res].n_nghs,i_res,res_extra[j_res].j_subres_on,res_extra[j_res].n_subres);
                }
                if (i_ngh == prot.res[i_res].n_ngh) break;
                
                clash = 0;
                for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
                    if (!prot.res[i_res].ngh[i_ngh]->n_subres) continue;
                    for (j_ngh=i_ngh+1; j_ngh<prot.res[i_res].n_ngh; j_ngh++) {
                        if (!prot.res[i_res].ngh[j_ngh]->n_subres) continue;
                        ic = prot.res[i_res].ngh[i_ngh]->subres[prot.res[i_res].ngh[i_ngh]->i_subres_on].conf[0]->i_conf_prot;
                        jc = prot.res[i_res].ngh[j_ngh]->subres[prot.res[i_res].ngh[j_ngh]->i_subres_on].conf[0]->i_conf_prot;
                        if (pair_vdw[ic][jc] > env.relax_clash_thr) {
                            clash = 1;
                            break;
                        }
                    }
                    if (clash) break;
                }
                if (clash) continue;
            }
            
            /* Get a list of conformers for relaxation */
            for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
                int k_rand;
                if (prot.res[i_res].ngh[i_ngh]->n_conf < 2) continue;
                i_subres = prot.res[i_res].ngh[i_ngh]->i_subres_on;
                k_rand = rand();
                i_conf = k_rand / (RAND_MAX/prot.res[i_res].ngh[i_ngh]->subres[i_subres].n_conf+1.0);
                i_conf = prot.res[i_res].ngh[i_ngh]->subres[i_subres].n_conf * ((float)k_rand / (RAND_MAX+1.0));
                prot.res[i_res].ngh[i_ngh]->conf_w = prot.res[i_res].ngh[i_ngh]->subres[i_subres].conf[i_conf];
                //printf("i_conf %d, n_ngh %d\n",i_conf, prot.res[i_res].n_ngh);
            }
            
            /* Calculate starting energy */
            vdw_start = 0.;
            coulomb_start = 0;
            for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
                if (!prot.res[i_res].ngh[i_ngh]->n_subres) continue;
                for (j_ngh=0; j_ngh<prot.res[i_res].n_ngh; j_ngh++) {
                    if (!prot.res[i_res].ngh[j_ngh]->n_subres) continue;
                    vdw_start += vdw_conf(prot.res[i_res].ngh[i_ngh]->i_res_prot,prot.res[i_res].ngh[i_ngh]->conf_w->i_conf_res,
                    prot.res[i_res].ngh[j_ngh]->i_res_prot,prot.res[i_res].ngh[j_ngh]->conf_w->i_conf_res,prot);
                    coulomb_start += coulomb_conf(prot.res[i_res].ngh[i_ngh]->i_res_prot,prot.res[i_res].ngh[i_ngh]->conf_w->i_conf_res,
                    prot.res[i_res].ngh[j_ngh]->i_res_prot,prot.res[i_res].ngh[j_ngh]->conf_w->i_conf_res,prot);
                }
            }
            //printf("Start %8.3f,%8.3f\n",vdw_start,coulomb_start);
            
            /* Get a list of atoms for relaxation */
            // conf_w is the actual conformer in the relaxation
            n_relax=0;
            for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
                /* Add conformer 0 into list */
                for (i_atom =0; i_atom<prot.res[i_res].ngh[i_ngh]->conf[0].n_atom; i_atom++) {
                    ATOM *atom_p;
                    atom_p = &prot.res[i_res].ngh[i_ngh]->conf[0].atom[i_atom];
                    if (!atom_p->on) continue;
                    n_relax++;
                    relax_atoms = realloc(relax_atoms, n_relax*sizeof(RELAX_ATOM *));
                    relax_atoms[n_relax-1] = &all_atoms[atom_p->i_atom_prot];
                }
                
                /* Add sidechain into list */
                if (prot.res[i_res].ngh[i_ngh]->n_conf<2) continue;
                for (i_atom =0; i_atom<prot.res[i_res].ngh[i_ngh]->conf_w->n_atom; i_atom++) {
                    ATOM *atom_p;
                    atom_p = &prot.res[i_res].ngh[i_ngh]->conf_w->atom[i_atom];
                    if (!atom_p->on) continue;
                    n_relax++;
                    relax_atoms = realloc(relax_atoms, n_relax*sizeof(RELAX_ATOM *));
                    relax_atoms[n_relax-1] = &all_atoms[atom_p->i_atom_prot];
                }
            }
            //printf("n_relax=%d\n",n_relax);
            /*
            int  n_relax=0, j_ngh;
            RELAX_ATOM *relax_atoms=NULL;
            for (j_ngh=0; j_ngh<res_extra[i_res].n_nghs; j_ngh++) {
                j_res = res_extra[i_res].ngh_res[j_ngh].j_res;
                RES *jres_p = &prot.res[j_res];
                if(res_extra[j_res].hyd_opt) {
                    int n_relax_old = n_relax;
                    n_relax += jres_p->conf[0].n_atom;
                    relax_atoms = realloc(relax_atoms, n_relax*sizeof(RELAX_ATOM));
                    
                    for (j_atom =0; j_atom<jres_p->conf[0].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].atom_p = &jres_p->conf[0].atom[j_atom];
                    for (j_atom =0; j_atom<jres_p->conf[0].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].i_res_prot = j_res;
                    for (j_atom =0; j_atom<jres_p->conf[0].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].move=0;
                    
                    int k_rand=rand();
                    int j_conf = (int) ( (res_extra[j_res].n_conf_o-1.0) * k_rand/(RAND_MAX+1.0) + 1.0);
                    
                    ins_conf(jres_p, jres_p->n_conf, jres_p->conf[j_conf].n_atom);
                    if (cpy_conf(&jres_p->conf[jres_p->n_conf-1], &jres_p->conf[j_conf])) {printf("   Error! opt_h(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",jres_p->conf[j_conf].confName,jres_p->resName, jres_p->resSeq, jres_p->n_conf-1); return USERERR;}
                    j_conf = jres_p->n_conf-1;
                    get_connect12_conf(j_res,j_conf,prot);
                    res_extra[j_res].j_conf_on = j_conf;
                    
                    n_relax_old = n_relax;
                    n_relax += jres_p->conf[j_conf].n_atom;
                    relax_atoms = realloc(relax_atoms, n_relax*sizeof(RELAX_ATOM));
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].atom_p = &jres_p->conf[j_conf].atom[j_atom];
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].i_res_prot = j_res;
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) {
                        ATOM *jatom_p = &jres_p->conf[j_conf].atom[j_atom];
                        TORSION     torsion;
                        if ( param_get("TORSION",jres_p->conf[j_conf].confName, jatom_p->name, &torsion) ) {relax_atoms[n_relax_old+j_atom].move=0;continue;}
                        if ( !torsion.hyd_opt ) {relax_atoms[n_relax_old+j_atom].move=0;continue;}
                        relax_atoms[n_relax_old+j_atom].move=1;
                        int i_connect, j_connect;
                        for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
                            if (!jatom_p->connect12[i_connect]) break;
                            if (!strcmp(jatom_p->connect12[i_connect]->name, torsion.atom1)) break;
                        }
                        if (!jatom_p->connect12[i_connect]) {relax_atoms[n_relax_old+j_atom].move=0;continue;}
                        
                        for (j_connect=0; j_connect<MAX_CONNECTED; j_connect++) {
                            if (!jatom_p->connect12[i_connect]->connect12[j_connect]) break;
                            if (!strcmp(jatom_p->connect12[i_connect]->connect12[j_connect]->name, torsion.atom2)) break;
                        }
                        if (!jatom_p->connect12[i_connect]->connect12[j_connect]) {relax_atoms[n_relax_old+j_atom].move=0;continue;}
                        
                        relax_atoms[n_relax_old+j_atom].axis = line_2v(jatom_p->connect12[i_connect]->connect12[j_connect]->xyz,
                        jatom_p->connect12[i_connect]->xyz);
                    }
                }
                else if (res_extra[j_res].hoh_opt) {
                    int k_rand=rand();
                    int j_conf = (int) ( (res_extra[j_res].n_conf_o-1.0) * k_rand/(RAND_MAX+1.0) + 1.0);
                    
                    ins_conf(jres_p, jres_p->n_conf, jres_p->conf[j_conf].n_atom);
                    if (cpy_conf(&jres_p->conf[jres_p->n_conf-1], &jres_p->conf[j_conf])) {printf("   Error! opt_h(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",jres_p->conf[j_conf].confName,jres_p->resName, jres_p->resSeq, jres_p->n_conf-1); return USERERR;}
                    j_conf = jres_p->n_conf-1;
                    get_connect12_conf(j_res,j_conf,prot);
                    res_extra[j_res].j_conf_on = j_conf;
                    
                    int n_relax_old = n_relax;
                    n_relax += jres_p->conf[j_conf].n_atom;
                    relax_atoms = realloc(relax_atoms, n_relax*sizeof(RELAX_ATOM));
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].atom_p = &jres_p->conf[j_conf].atom[j_atom];
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].i_res_prot = j_res;
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) {
                        if (j_atom==0) relax_atoms[n_relax_old+j_atom].move=2;
                        else {
                            relax_atoms[n_relax_old+j_atom].move=3;
                            relax_atoms[n_relax_old+j_atom].origin = jres_p->conf[j_conf].atom[0].xyz;
                        }
                    }
                }
                else {
                    int n_relax_old = n_relax;
                    n_relax += jres_p->conf[0].n_atom;
                    relax_atoms = realloc(relax_atoms, n_relax*sizeof(RELAX_ATOM));
                    for (j_atom =0; j_atom<jres_p->conf[0].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].atom_p = &jres_p->conf[0].atom[j_atom];
                    for (j_atom =0; j_atom<jres_p->conf[0].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].i_res_prot = j_res;
                    for (j_atom =0; j_atom<jres_p->conf[0].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].move=0;
                    
                    if (jres_p->n_conf==1) continue;
                    
                    int k_rand=rand();
                    int j_conf = (int) ( (res_extra[j_res].n_conf_o-1.0) * k_rand/(RAND_MAX+1.0) + 1.0);
                    res_extra[j_res].j_conf_on = j_conf;
                    
                    n_relax_old = n_relax;
                    n_relax += jres_p->conf[j_conf].n_atom;
                    relax_atoms = realloc(relax_atoms, n_relax*sizeof(RELAX_ATOM));
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].atom_p = &jres_p->conf[j_conf].atom[j_atom];
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].i_res_prot = j_res;
                    for (j_atom =0; j_atom<jres_p->conf[j_conf].n_atom; j_atom++) relax_atoms[n_relax_old+j_atom].move=0;
                }
            }
            */
            
            /* Print out interactions 
            float check_vdw = 0, check_coulomb=0;
            for (j_ngh=0; j_ngh<res_extra[i_res].n_nghs; j_ngh++) {
                j_res = res_extra[i_res].ngh_res[j_ngh].j_res;
                j_conf = res_extra[j_res].j_conf_on;
                int k_res, k_conf, k_ngh;
                for (k_ngh=0; k_ngh<res_extra[i_res].n_nghs; k_ngh++) {
                    k_res = res_extra[i_res].ngh_res[k_ngh].j_res;
                    k_conf = res_extra[k_res].j_conf_on;
                    check_vdw += vdw_conf(j_res,j_conf,k_res,k_conf,prot);
                    check_coulomb += coulomb_conf(j_res,j_conf,k_res,k_conf,prot);
                }
            }
            printf("before relaxation vdw %8.3f, coulomb %8.3f; ", check_vdw/2.,check_coulomb/2.);
            */
            
            /* Print out each state to a PDB 
            FILE *temp_fp;
            char  sbuffer[MAXCHAR_LINE];
            sprintf(sbuffer,"%06d",counter);counter++;
            strcat(sbuffer,".pdb");
            temp_fp = fopen(sbuffer,"w");
            for (i_atom=0; i_atom< n_relax; i_atom++) {
                if (!relax_atoms[i_atom]->atom_p->on) continue;
                fprintf(temp_fp, "ATOM %d %3d  %4s %3s %c%4d%c   %8.3f%8.3f%8.3f %7.3f       %6.2f       \n",
                relax_atoms[i_atom]->move,
                relax_atoms[i_atom]->k_relax,
                relax_atoms[i_atom]->atom_p->name,
                prot.res[relax_atoms[i_atom]->i_res_prot].resName,
                prot.res[relax_atoms[i_atom]->i_res_prot].chainID,
                prot.res[relax_atoms[i_atom]->i_res_prot].resSeq,
                prot.res[relax_atoms[i_atom]->i_res_prot].iCode,
                relax_atoms[i_atom]->atom_p->xyz.x,
                relax_atoms[i_atom]->atom_p->xyz.y,
                relax_atoms[i_atom]->atom_p->xyz.z,
                relax_atoms[i_atom]->atom_p->crg,
                relax_atoms[i_atom]->atom_p->rad
                );
            }
            fclose(temp_fp);
            */
            //printf("before relax\n");
            convg = relax(n_relax, relax_atoms, prot);
            //printf("after relax, %d\n", convg);
            
            /* Calculate energy after relaxation */
            vdw_end = 0.;
            coulomb_end = 0;
            for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
                if (!prot.res[i_res].ngh[i_ngh]->n_subres) continue;
                for (j_ngh=0; j_ngh<prot.res[i_res].n_ngh; j_ngh++) {
                    if (!prot.res[i_res].ngh[j_ngh]->n_subres) continue;
                    vdw_end += vdw_conf(prot.res[i_res].ngh[i_ngh]->i_res_prot,prot.res[i_res].ngh[i_ngh]->conf_w->i_conf_res,
                    prot.res[i_res].ngh[j_ngh]->i_res_prot,prot.res[i_res].ngh[j_ngh]->conf_w->i_conf_res,prot);
                    coulomb_end += coulomb_conf(prot.res[i_res].ngh[i_ngh]->i_res_prot,prot.res[i_res].ngh[i_ngh]->conf_w->i_conf_res,
                    prot.res[i_res].ngh[j_ngh]->i_res_prot,prot.res[i_res].ngh[j_ngh]->conf_w->i_conf_res,prot);
                }
            }

            /*
            printf("ntrial=%d, niter = %4d/%d",n_trial,convg,env.relax_niter);
            printf("Start %8.3f,%8.3f",vdw_start,coulomb_start);
            printf("End   %8.3f,%8.3f\n",vdw_end,coulomb_end);
            */
            //printf("DEBUG=========================================1\n");
            /* If difference of energy lower than threshold, increase counter */
            diff = coulomb_end - coulomb_start;
            //printf("diff, %8.3f\n", diff);
            if (diff < env.relax_e_thr) {
            	//printf("env.relax_e_thr: %.2f, smaller\n", env.relax_e_thr);
                for (i_atom=0; i_atom< n_relax; i_atom++) {
                    float angle;
                    float theta,phi,psi;
                    hyd_pos = prot.res[relax_atoms[i_atom]->i_res_prot].n_hyd_pos;
                    //printf("name: %s, n_hyd_pos: %d, move: %d\n", relax_atoms[i_atom]->atom_p->name, hyd_pos, relax_atoms[i_atom]->move);
                    if (!hyd_pos) continue;
                    if (relax_atoms[i_atom]->move == 1) {
                        angle = torsion_angle(relax_atoms[i_atom]->atom_p->xyz,relax_atoms[i_atom]->torsion_atom[0]->xyz,
                        relax_atoms[i_atom]->torsion_atom[1]->xyz,relax_atoms[i_atom]->torsion_atom[2]->xyz);
                        i_bin = (int) (angle/(2.*env.PI)*(float) hyd_pos + 0.5);
                        if (i_bin == hyd_pos) i_bin = 0;
                        
                        i_subres = prot.res[relax_atoms[i_atom]->i_res_prot].conf[relax_atoms[i_atom]->i_conf_res].i_subres_res;
                        ka = prot.res[relax_atoms[i_atom]->i_res_prot].subres[i_subres].conf[0]->atom[relax_atoms[i_atom]->i_atom_conf].i_atom_prot;
                        all_atoms[ka].counter[i_bin]++;
                        //printf("MOVE1\n");
                        //printf("%s, %8.3f, %d, %d, %d, counter=%d\n",relax_atoms[i_atom]->atom_p->name,angle/env.d2r,i_bin,i_subres,ka,all_atoms[ka].counter[i_bin]);
                    }
                    else if (relax_atoms[i_atom]->move == 2) {
                    	if (strcmp("HOH01", prot.res[relax_atoms[i_atom]->i_res_prot].conf[relax_atoms[i_atom]->i_conf_res].confName)) continue;
                    	//printf("MOVE2\n");
                        //printf("name: %s, n_hyd_pos: %d, move: %d\n", relax_atoms[i_atom]->atom_p->name, hyd_pos, relax_atoms[i_atom]->move);
                        //printf("count: %d, resname: %s\n", relax_atoms[i_atom]->counter[i_bin + j_bin*hyd_pos + k_bin*hyd_pos*hyd_pos], prot.res[relax_atoms[i_atom]->i_res_prot].resName);

                        i_subres = prot.res[relax_atoms[i_atom]->i_res_prot].conf[relax_atoms[i_atom]->i_conf_res].i_subres_res;
                        ka = prot.res[relax_atoms[i_atom]->i_res_prot].subres[i_subres].conf[0]->atom[relax_atoms[i_atom]->i_atom_conf].i_atom_prot;
                        //printf("ka %d, n_slaved: %d, confName: %s\n", ka, all_atoms[ka].n_slaved, prot.res[relax_atoms[i_atom]->i_res_prot].conf[relax_atoms[i_atom]->i_conf_res].confName);
                        //printf("i_atom, xyz: %.2f, %.2f, %.2f\n", relax_atoms[i_atom]->atom_p->xyz.x,
                        //		relax_atoms[i_atom]->atom_p->xyz.y,
                        //		relax_atoms[i_atom]->atom_p->xyz.z);



                        //printf("i_atom, slave0, xyz: %.2f, %.2f, %.2f\n", relax_atoms[i_atom]->slaved[0]->atom_p->xyz.x,
                        //		relax_atoms[i_atom]->slaved[0]->atom_p->xyz.y,
                        //		relax_atoms[i_atom]->slaved[0]->atom_p->xyz.z);

                        //printf("i_atom, slave1, xyz: %.2f, %.2f, %.2f\n", relax_atoms[i_atom]->slaved[1]->atom_p->xyz.x,
                        //                        		relax_atoms[i_atom]->slaved[1]->atom_p->xyz.y,
                        //                        		relax_atoms[i_atom]->slaved[1]->atom_p->xyz.z);

                        water_orient(relax_atoms[i_atom]->atom_p->xyz,
                        relax_atoms[i_atom]->slaved[0]->atom_p->xyz,
                        relax_atoms[i_atom]->slaved[1]->atom_p->xyz,
                        &theta, &phi, &psi);
                        //printf("MIDDLE, ka: %d\n", ka);
                        i_bin = (int) (theta/(2.*env.PI)*(float) hyd_pos + 0.5);
                        if (i_bin == hyd_pos) i_bin = 0;
                        j_bin = (int) (  phi/(2.*env.PI)*(float) hyd_pos + 0.5);
                        if (j_bin == hyd_pos) j_bin = 0;
                        k_bin = (int) (  psi/(2.*env.PI)*(float) hyd_pos + 0.5);
                        if (k_bin == hyd_pos) k_bin = 0;
                        
                        i_subres = prot.res[relax_atoms[i_atom]->i_res_prot].conf[relax_atoms[i_atom]->i_conf_res].i_subres_res;

                        //printf("i_subres %d bin %d %d %d\n",i_subres, i_bin, j_bin, k_bin);
                        ka = prot.res[relax_atoms[i_atom]->i_res_prot].subres[i_subres].conf[0]->atom[relax_atoms[i_atom]->i_atom_conf].i_atom_prot;
                        //printf("Near end, ka %d\n", ka);
                        all_atoms[ka].counter[i_bin + j_bin*hyd_pos + k_bin*hyd_pos*hyd_pos]++;
                        //(relax_atoms[i_atom]->counter[i_bin + j_bin*hyd_pos + k_bin*hyd_pos*hyd_pos])++;
                    }
                }
            }
            /* Recover x,y,z */
            for (i_atom=0; i_atom< n_relax; i_atom++) {
                relax_atoms[i_atom]->atom_p->xyz = relax_atoms[i_atom]->xyz_back;
            }
            free(relax_atoms);
        }
        for (i_ngh=0; i_ngh<prot.res[i_res].n_ngh; i_ngh++) {
            for (i_conf=0; i_conf<prot.res[i_res].ngh[i_ngh]->n_conf; i_conf++) {
                for (i_atom = 0; i_atom < prot.res[i_res].ngh[i_ngh]->conf[i_conf].n_atom; i_atom++) {
                    ATOM *iatom_p;
                    iatom_p = &prot.res[i_res].ngh[i_ngh]->conf[i_conf].atom[i_atom];
                    if (!iatom_p->on) continue;
                    if (!all_atoms[iatom_p->i_atom_prot].move) continue;
                    free(factor_matrix[all_atoms[iatom_p->i_atom_prot].k_relax]);
                }
            }
        }
        free(factor_matrix);
    }

    /* Add conformers */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        int hyd_pos;
        int nconf_old;
        if(!prot.res[i_res].opt_hyd) continue;
        if(!prot.res[i_res].n_hyd_pos) continue;
        hyd_pos = prot.res[i_res].n_hyd_pos;
        nconf_old = prot.res[i_res].n_conf;
        
        for (i_subres=0; i_subres<prot.res[i_res].n_subres; i_subres++) {
            RES added_conf;
            int i_conf_rec,i_subres_rec,i_conf_subres_rec;
            memset(&added_conf,0,sizeof(RES));
            ins_conf(&added_conf,0,prot.res[i_res].subres[i_subres].conf[0]->n_atom);
            cpy_conf(&added_conf.conf[0],prot.res[i_res].subres[i_subres].conf[0]);
            added_conf.conf[0].counter = 0;
            
            for (i_atom=0; i_atom<prot.res[i_res].subres[i_subres].conf[0]->n_atom; i_atom++) {
                ATOM *atom_p;
                if (!prot.res[i_res].subres[i_subres].conf[0]->atom[i_atom].on) continue;
                atom_p = &prot.res[i_res].subres[i_subres].conf[0]->atom[i_atom];
                ka = atom_p->i_atom_prot;
                if (all_atoms[ka].move != 1 && all_atoms[ka].move != 2) continue;
                
                n_conf = added_conf.n_conf;
                //printf("%s, hyd_pos=%d,ka=%d,move=%d\n",added_conf.conf[0].confName,hyd_pos,ka,all_atoms[ka].move);
                for (i_conf=0;i_conf<n_conf;i_conf++) {
                    
                    if (all_atoms[ka].move==1) {
                        float a;
                        VECTOR i,j,k,r21,r23,r10,r10_k,r10_ij;
                        int i_bin;
                        
                        r21 = vector_vminusv(all_atoms[ka].torsion_atom[0]->xyz,all_atoms[ka].torsion_atom[1]->xyz);
                        r23 = vector_vminusv(all_atoms[ka].torsion_atom[2]->xyz,all_atoms[ka].torsion_atom[1]->xyz);
                        k   = vector_normalize(r21);
                        i   = vector_normalize(vector_vminusv(r23, vector_rescale(k, vdotv(k,r23))));
                        j   = vector_vxv(k,i);
                        //printf("ii=%f,jj=%f,kk=%f,ij=%f,jk=%f,ki=%f\n",vdotv(i,i),vdotv(j,j),vdotv(k,k),vdotv(i,j),vdotv(j,k),vdotv(k,i));
                        r10 = vector_vminusv(all_atoms[ka].atom_p->xyz,  all_atoms[ka].torsion_atom[0]->xyz);
                        r10_k  = vector_rescale(k, vdotv(r10,k));
                        r10_ij = vector_vminusv(r10, r10_k);
                        a = sqrt(vdotv(r10_ij,r10_ij));
                        
                        for (i_bin=0;i_bin<hyd_pos;i_bin++) {
                            float phi;
                            CONF *new_conf;
                            //printf("ka = %d, i_bin = %d, counter=%d\n",ka, i_bin, all_atoms[ka].counter[i_bin]);
                            if (!all_atoms[ka].counter[i_bin]) continue;
                            phi = (2.*env.PI/(float) hyd_pos) * (float) i_bin;
                            
                            ins_conf(&added_conf, added_conf.n_conf, added_conf.conf[i_conf].n_atom);
                            new_conf = &added_conf.conf[added_conf.n_conf-1];
                            cpy_conf(new_conf, &added_conf.conf[i_conf]);
                            new_conf->counter += all_atoms[ka].counter[i_bin];
                            
                            new_conf->atom[all_atoms[ka].i_atom_conf].xyz = 
                            vector_vplusv(all_atoms[ka].torsion_atom[0]->xyz,
                            vector_sum3v(vector_rescale(i,a*cos(phi)), vector_rescale(j,a*sin(phi)), r10_k));
                            //printf("O:%8.3f%8.3f%8.3f\n",all_atoms[ka].torsion_atom[0]->xyz.x,all_atoms[ka].torsion_atom[0]->xyz.y,all_atoms[ka].torsion_atom[0]->xyz.z);
                            //printf("H:%8.3f%8.3f%8.3f\n",all_atoms[ka].atom_p->xyz.x,all_atoms[ka].atom_p->xyz.y,all_atoms[ka].atom_p->xyz.z);
                            //printf("D:%8.3f%8.3f%8.3f\n",new_conf->atom[all_atoms[ka].i_atom_conf].xyz.x,new_conf->atom[all_atoms[ka].i_atom_conf].xyz.y,new_conf->atom[all_atoms[ka].i_atom_conf].xyz.z);
                        }
                    }
                    if (all_atoms[ka].move == 2) {
                        float a;
                        int i_bin,j_bin,k_bin;
                        float theta,phi,psi;
                        a = 0.5*dvv(all_atoms[ka].atom_p->xyz,all_atoms[ka].slaved[0]->atom_p->xyz);
                        for (i_bin=0;i_bin<hyd_pos;i_bin++) {
                            for (j_bin=0;j_bin<hyd_pos;j_bin++) {
                                for (k_bin=0;k_bin<hyd_pos;k_bin++) {
                                    CONF *new_conf;
                                    VECTOR i,k;
                                    int j_atom;
                                    if (!all_atoms[ka].counter[i_bin + j_bin*hyd_pos + k_bin*hyd_pos*hyd_pos]) continue;
                                    theta = (2.*env.PI/(float) hyd_pos) * (float) i_bin;
                                    phi   = (2.*env.PI/(float) hyd_pos) * (float) j_bin;
                                    psi   = (2.*env.PI/(float) hyd_pos) * (float) k_bin;
                                    
                                    ins_conf(&added_conf, added_conf.n_conf, added_conf.conf[i_conf].n_atom);
                                    new_conf = &added_conf.conf[added_conf.n_conf-1];
                                    cpy_conf(new_conf, &added_conf.conf[i_conf]);
                                    new_conf->counter += all_atoms[ka].counter[i_bin + j_bin*hyd_pos + k_bin*hyd_pos*hyd_pos];
                                    
                                    i.x = cos(phi)*cos(theta)*cos(psi)-sin(phi)*sin(psi);
                                    i.y = sin(phi)*cos(theta)*cos(psi)+cos(phi)*sin(psi);
                                    i.z = -sin(theta)*cos(psi);
                                    
                                    k.x = cos(phi)*sin(theta);
                                    k.y = sin(phi)*sin(theta);
                                    k.z = cos(theta);
                                    
                                    j_atom = all_atoms[ka].slaved[0]->i_atom_conf;
                                    new_conf->atom[j_atom].xyz = vector_vplusv(all_atoms[ka].atom_p->xyz,vector_vminusv(vector_rescale(k,a), vector_rescale(i,sqrt(3.)*a)));
                                    j_atom = all_atoms[ka].slaved[1]->i_atom_conf;
                                    new_conf->atom[j_atom].xyz = vector_vplusv(all_atoms[ka].atom_p->xyz,vector_vplusv( vector_rescale(k,a), vector_rescale(i,sqrt(3.)*a)));
                                }
                            }
                        }
                    }
                }
                
                for (i_conf=0;i_conf<n_conf;i_conf++) del_conf(&added_conf,0);
            }
            
            for (i_conf=0;i_conf<added_conf.n_conf;i_conf++) {
                added_conf.conf[i_conf].history[6] = 'H';
                sprintf(sbuff,"%03d",i_conf);
                strncpy(added_conf.conf[i_conf].history+7,sbuff,3);
                
                for (j_conf=1; j_conf<prot.res[i_res].n_conf; j_conf++) {
                    if (!cmp_conf(added_conf.conf[i_conf], prot.res[i_res].conf[j_conf], 0.1)) break;
                }
                
                if (j_conf>=prot.res[i_res].n_conf) {
                    ins_conf(&prot.res[i_res], prot.res[i_res].n_conf, added_conf.conf[i_conf].n_atom);
                    for (i_conf_rec = 1; i_conf_rec < nconf_old; i_conf_rec++) {
                        i_subres_rec = prot.res[i_res].conf[i_conf_rec].i_subres_res;
                        i_conf_subres_rec = prot.res[i_res].conf[i_conf_rec].i_conf_subres;
                        prot.res[i_res].subres[i_subres_rec].conf[i_conf_subres_rec] = &prot.res[i_res].conf[i_conf_rec];
                    }
                    j_conf = prot.res[i_res].n_conf - 1;
                    cpy_conf(&prot.res[i_res].conf[j_conf],&added_conf.conf[i_conf]);
                }
                else {
                    prot.res[i_res].conf[j_conf].counter += added_conf.conf[i_conf].counter;
                }
            }
            for (i_conf=added_conf.n_conf-1;i_conf>=0;i_conf--) {
                del_conf(&added_conf,i_conf);
            }
        }
        
        if (prot.res[i_res].nconf_limit) {
            int i_conf_rec,i_subres_rec,i_conf_subres_rec;
            if (nconf_old-1 > prot.res[i_res].nconf_limit) {
                printf(" Input conformer number over the limit, residue \"%3s %c%04d%c\", n_conf = %d (>%3d)\n",
                prot.res[i_res].resName, prot.res[i_res].chainID, prot.res[i_res].resSeq, prot.res[i_res].iCode, nconf_old, prot.res[i_res].nconf_limit);
                while (prot.res[i_res].n_conf > nconf_old) del_conf(&prot.res[i_res],prot.res[i_res].n_conf-1);
            }
            else {
                while (prot.res[i_res].n_conf-1 > prot.res[i_res].nconf_limit) {
                    int counter_min, iconf_min;
                    counter_min = prot.res[i_res].conf[nconf_old].counter;
                    iconf_min = nconf_old;
                    for (i_conf=nconf_old+1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                        if (prot.res[i_res].conf[i_conf].counter < counter_min) {
                            counter_min = prot.res[i_res].conf[i_conf].counter;
                            iconf_min = i_conf;
                        }
                    }
                    del_conf(&prot.res[i_res], iconf_min);
                    for (i_conf_rec = 1; i_conf_rec < nconf_old; i_conf_rec++) {
                        i_subres_rec = prot.res[i_res].conf[i_conf_rec].i_subres_res;
                        i_conf_subres_rec = prot.res[i_res].conf[i_conf_rec].i_conf_subres;
                        prot.res[i_res].subres[i_subres_rec].conf[i_conf_subres_rec] = &prot.res[i_res].conf[i_conf_rec];
                    }
                }
            }
        }
    }
    
    for (i_res=0;i_res<prot.n_res;i_res++) {
      for (i_subres=0;i_subres<prot.res[i_res].n_subres;i_subres++) {
         free(prot.res[i_res].subres[i_subres].conf);
      }
      free(prot.res[i_res].subres);
      if (prot.res[i_res].n_ngh) {
          prot.res[i_res].n_ngh = 0;
          free(prot.res[i_res].ngh);
      }
    }

    for (ka=0;ka<na;ka++) {
        if (all_atoms[ka].counter)
            free(all_atoms[ka].counter);
    }
    for (ic=0;ic<prot.nc;ic++) {
        if (pair_vdw[ic])
            free(pair_vdw[ic]);
    }
    free(pair_vdw);
    /*
    for (i_elem=0;i_elem<n_elem;i_elem++) {
        free(C6_matrix[i_elem]);
        free(C12_matrix[i_elem]);
    }
    free(C6_matrix);
    free(C12_matrix);
    */
    free(all_atoms);

    return 0;
}

int relax(int n_relax, RELAX_ATOM **relax_atoms, PROT prot)
{
    float       phi_step = env.relax_phi;
    int         i_iter;
    float       TOQ_THR2 = env.relax_torq_thr*env.relax_torq_thr;
    float       C6,C12,factor;
    int         i_slaved;
    
    for (i_iter=0;i_iter<env.relax_niter; i_iter++) {
        int moving = 0, i_atom;
        
        /* force and torque on each atom */
        for (i_atom=0; i_atom< n_relax; i_atom++) {
            int j_atom;
            if (!relax_atoms[i_atom]->move) continue;
            /* Get force */
            relax_atoms[i_atom]->frc.x = 0.;
            relax_atoms[i_atom]->frc.y = 0.;
            relax_atoms[i_atom]->frc.z = 0.;
            for (j_atom=0; j_atom < n_relax; j_atom++) {
                VECTOR frc1, frc2, frc;
                if (j_atom == i_atom) continue;
                float sig_min = relax_atoms[i_atom]->atom_p->vdw_rad + relax_atoms[j_atom]->atom_p->vdw_rad;
                float eps = sqrt(relax_atoms[i_atom]->atom_p->vdw_eps * relax_atoms[j_atom]->atom_p->vdw_eps);
                C12 = eps*pow(sig_min,12);
                C6 = 2.*eps*pow(sig_min,6);
            
                //C6 = C6_matrix[relax_atoms[i_atom]->atom_p->i_elem][relax_atoms[j_atom]->atom_p->i_elem];
                //C12 = C12_matrix[relax_atoms[i_atom]->atom_p->i_elem][relax_atoms[j_atom]->atom_p->i_elem];
                factor = factor_matrix[relax_atoms[i_atom]->k_relax][relax_atoms[j_atom]->k_relax];
                
                frc1 = vdw_frc(relax_atoms[i_atom]->atom_p->xyz, relax_atoms[j_atom]->atom_p->xyz, C6, C12);
                frc2 = coulomb_frc(relax_atoms[i_atom]->atom_p->xyz, relax_atoms[j_atom]->atom_p->xyz, relax_atoms[i_atom]->atom_p->crg, relax_atoms[j_atom]->atom_p->crg);
                frc  = vector_vplusv(frc1,frc2);
                frc  = vector_rescale(frc, factor);
                relax_atoms[i_atom]->frc = vector_vplusv(relax_atoms[i_atom]->frc, frc);
                //printf("%s %s %10.1f,%10.1f,%8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f, %8.3f\n",
                //relax_atoms[i_atom]->atom_p->name,relax_atoms[j_atom]->atom_p->name,C6,C12,factor,frc1.x,frc1.y,frc1.z,frc2.x,frc2.y,frc2.z);
            }
            //printf("%s %8.3f,%8.3f,%8.3f\n",relax_atoms[i_atom]->atom_p->name,relax_atoms[i_atom]->frc.x,relax_atoms[i_atom]->frc.y,relax_atoms[i_atom]->frc.z);
            
            /* Calculate torque */
            if (relax_atoms[i_atom]->move==1) {   //hydroxyl hydrogen
                VECTOR r,norm_frc;
                float phi;
                int i_term;
                norm_frc = vector_vminusv(relax_atoms[i_atom]->frc, vector_rescale(relax_atoms[i_atom]->axis.t, vdotv(relax_atoms[i_atom]->frc, relax_atoms[i_atom]->axis.t)));
                r = vector_vminusv(relax_atoms[i_atom]->atom_p->xyz, relax_atoms[i_atom]->axis.p0);
                r = vector_normalize(vector_vminusv(r, vector_rescale(relax_atoms[i_atom]->axis.t, vdotv(r, relax_atoms[i_atom]->axis.t))));
                relax_atoms[i_atom]->torque = vector_vxv(r, norm_frc);
                
                /* add torque from the bond */
                phi = torsion_angle(relax_atoms[i_atom]->atom_p->xyz,relax_atoms[i_atom]->torsion_atom[0]->xyz,
                    relax_atoms[i_atom]->torsion_atom[1]->xyz,relax_atoms[i_atom]->torsion_atom[2]->xyz);
                for (i_term=0;i_term<relax_atoms[i_atom]->tors.n_term; i_term++) {
                    relax_atoms[i_atom]->torque = vector_vplusv(relax_atoms[i_atom]->torque, torsion_torq(phi, relax_atoms[i_atom]->tors.V2[i_term],
                        relax_atoms[i_atom]->tors.n_fold[i_term], relax_atoms[i_atom]->tors.gamma[i_term], relax_atoms[i_atom]->axis.t));
                }
            }
            else if (relax_atoms[i_atom]->move==3) { // water hydrogen
                VECTOR r;
                r = vector_vminusv(relax_atoms[i_atom]->atom_p->xyz, relax_atoms[i_atom]->origin);
                relax_atoms[i_atom]->torque = vector_vxv(r, relax_atoms[i_atom]->frc);
            }
        }
        
        /* get starting position */
        for (i_atom=0; i_atom< n_relax; i_atom++) {
            relax_atoms[i_atom]->xyz_new = relax_atoms[i_atom]->atom_p->xyz;
        }
        
        for (i_atom=0; i_atom< n_relax; i_atom++) {
            if (relax_atoms[i_atom]->move==1) {          //hydroxyl hydrogen
                float torque_sc2 = vdotv(relax_atoms[i_atom]->torque, relax_atoms[i_atom]->torque);
                if ( torque_sc2 > TOQ_THR2) {
                    GEOM  op;
                    LINE  axis = relax_atoms[i_atom]->axis;
                    //float torque_sc = sqrt(torque_sc2);
                    //torque_sc = (int)torque_sc + 1.;
                    //if (torque_sc > MAX_RELAX_SC) torque_sc = MAX_RELAX_SC;
                    
                    axis.t = vector_normalize(relax_atoms[i_atom]->torque);
                    
                    geom_reset(&op);
                    geom_roll(&op, phi_step, axis);
                    geom_apply(op, &relax_atoms[i_atom]->xyz_new);
                    //printf("%f\n",vdotv(torque,torque));
                    moving = 1;
                }
            }
            else if (relax_atoms[i_atom]->move==2) {     //water oxygen
                VECTOR tot_toq;
                float torque_sc2;
                tot_toq.x=0.;tot_toq.y=0.;tot_toq.z=0.;
                for (i_slaved=0; i_slaved<relax_atoms[i_atom]->n_slaved; i_slaved++) {
                    tot_toq = vector_vplusv(tot_toq, relax_atoms[i_atom]->slaved[i_slaved]->torque);
                    //printf("%f,%f,%f\n",tot_toq.x,tot_toq.y,tot_toq.z);
                }
                torque_sc2 = vdotv(tot_toq,tot_toq);
                if ( torque_sc2 > TOQ_THR2) {
                    LINE axis;
                    GEOM  op;
                    //float torque_sc = sqrt(torque_sc2);
                    //torque_sc = (int)torque_sc + 1.;
                    //if (torque_sc > MAX_RELAX_SC) torque_sc = MAX_RELAX_SC;
                    
                    axis.p0 = relax_atoms[i_atom]->atom_p->xyz;
                    axis.t = vector_normalize(tot_toq);
                    
                    geom_reset(&op);
                    geom_roll(&op, phi_step, axis);
                    
                    for (i_slaved=0; i_slaved<relax_atoms[i_atom]->n_slaved; i_slaved++) {
                        geom_apply(op, &relax_atoms[i_atom]->slaved[i_slaved]->xyz_new);
                        //printf("hoh%f\n",vdotv(tot_toq,tot_toq));
                    }
                    moving = 1;
                }
            }
        }
        
        //printf("moving=%d\n",moving);
        if (!moving) break;

        for (i_atom=0; i_atom< n_relax; i_atom++) {
            relax_atoms[i_atom]->atom_p->xyz = relax_atoms[i_atom]->xyz_new;
        }
        
        /* Print out interactions 
        check_vdw = 0; check_coulomb=0;
        printf("\n");
        for (j_ngh=0; j_ngh<res_extra[i_res].n_nghs; j_ngh++) {
            j_res = res_extra[i_res].ngh_res[j_ngh].j_res;
            j_conf = res_extra[j_res].j_conf_on;
            int k_res, k_conf, k_ngh;
            for (k_ngh=0; k_ngh<res_extra[i_res].n_nghs; k_ngh++) {
                k_res = res_extra[i_res].ngh_res[k_ngh].j_res;
                k_conf = res_extra[k_res].j_conf_on;
                printf("%8.3f,%8.3f",vdw_conf(j_res,j_conf,k_res,k_conf,prot),coulomb_conf(j_res,j_conf,k_res,k_conf,prot));
                check_vdw += vdw_conf(j_res,j_conf,k_res,k_conf,prot);
                check_coulomb += coulomb_conf(j_res,j_conf,k_res,k_conf,prot);
            }
            printf("\n");
        }
        */
        /* Print out each state to a PDB 
        //if (strcmp(relax_atoms[0].atom_p->resName,"HOH")) continue;
        if (counter<1000) {
            FILE *temp_fp;
            char  sbuffer[MAXCHAR_LINE];
            sprintf(sbuffer,"%06d",counter);counter++;
            strcat(sbuffer,".pdb");
            temp_fp = fopen(sbuffer,"w");
            for (i_atom=0; i_atom< n_relax; i_atom++) {
                if (!relax_atoms[i_atom]->atom_p->on) continue;
                fprintf(temp_fp, "ATOM   %d    %4s %3s %c%4d%c   %8.3f%8.3f%8.3f %7.3f  %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f \n",
                relax_atoms[i_atom]->move,
                relax_atoms[i_atom]->atom_p->name,
                prot.res[relax_atoms[i_atom]->i_res_prot].resName,
                prot.res[relax_atoms[i_atom]->i_res_prot].chainID,
                prot.res[relax_atoms[i_atom]->i_res_prot].resSeq,
                prot.res[relax_atoms[i_atom]->i_res_prot].iCode,
                relax_atoms[i_atom]->atom_p->xyz.x,
                relax_atoms[i_atom]->atom_p->xyz.y,
                relax_atoms[i_atom]->atom_p->xyz.z,
                relax_atoms[i_atom]->atom_p->crg,
                relax_atoms[i_atom]->torque.x,
                relax_atoms[i_atom]->torque.y,
                relax_atoms[i_atom]->torque.z,
                relax_atoms[i_atom]->origin.x,
                relax_atoms[i_atom]->origin.y,
                relax_atoms[i_atom]->origin.z,
                relax_atoms[i_atom]->frc.x,
                relax_atoms[i_atom]->frc.y,
                relax_atoms[i_atom]->frc.z);
            }
            fclose(temp_fp);
        }
        */
    }
    
    return i_iter;
    //if (i_iter == env.relax_niter) return -1;
    //else return i_iter;
}

int water_orient(VECTOR v0, VECTOR v1, VECTOR v2, float *theta, float *phi, float *psi)
{
    VECTOR i,j,k,i_p,j_p;
    VECTOR v01,v02,v12,bisect102,bisect102_p;
    /*
    printf("v0: %.2f, %.2f, %.2f\n", v0.x, v0.y, v0.z);
    printf("v1: %.2f, %.2f, %.2f\n", v1.x, v1.y, v1.z);
    printf("v2: %.2f, %.2f, %.2f\n", v2.x, v2.y, v2.z);
    */
    i.x = 1.;i.y = 0.;i.z = 0.;
    j.x = 0.;j.y = 1.;j.z = 0.;
    k.x = 0.;k.y = 0.;k.z = 1.;
    
    v01 = vector_vminusv(v1,v0);
    v02 = vector_vminusv(v2,v0);
    v12 = vector_normalize(vector_vminusv(v2,v1));
    
    bisect102 = vector_normalize(vector_vplusv(v01,v02));
    *theta = acos(vdotv(k,bisect102));
    
    bisect102_p = vector_normalize(vector_vminusv(bisect102, vector_rescale(k,vdotv(k,bisect102))));
    *phi = acos(vdotv(i,bisect102_p));
    if (vdotv(j,bisect102_p)<0) *phi = 2.*env.PI - *phi;
    
    i_p.x = cos(*theta)*cos(*phi);
    i_p.y = cos(*theta)*sin(*phi);
    i_p.z = sin(*theta);
    
    j_p.x = -sin(*phi);
    j_p.y =  cos(*phi);
    j_p.z =  0.;
    
    *psi=acos(vdotv(i_p,v12));
    if (vdotv(j_p,v12)<0) *psi = 2.*env.PI - *psi;
    
    return 0;
}

void load_headlst(PROT prot) {
    int  i_res;
    FILE *fp;
    char sbuff[MAXCHAR_LINE], sbuff2[MAXCHAR_LINE];
    char do_rot,do_sw,opt_hyd, resName[4],chainID,iCode;
    int  resSeq;
    int  rotations,n_hyd_pos,nconf_limit;
    float phi_swing;
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        prot.res[i_res].do_rot      = env.pack;
        prot.res[i_res].do_sw       = env.swing;
        prot.res[i_res].opt_hyd     = env.relax_h;
        prot.res[i_res].rotations   = env.rotations;
        prot.res[i_res].phi_swing   = env.phi_swing;
        prot.res[i_res].n_hyd_pos   = env.relax_n_hyd;
        prot.res[i_res].nconf_limit = env.nconf_limit;
        
        if (!prot.res[i_res].rotations)         prot.res[i_res].do_rot = 0;
        if ( prot.res[i_res].phi_swing < 1e-4)  prot.res[i_res].do_sw = 0;
        if (!prot.res[i_res].opt_hyd)           prot.res[i_res].n_hyd_pos = 0;
    }
    

    /* always read the head1.lst */
//    if (env.rot_specif) {
    if ( (fp = fopen(FN_CONFLIST1,"r")) ) {
        while(fgets(sbuff, sizeof(sbuff), fp)) {
            rm_comment(sbuff2, sbuff);
            if (strlen(sbuff2)<10) continue;
            
            sscanf(sbuff2, "%s %c%4d%c R %c %2d S %c %4f H %c %2d M %3d",
            resName,  &chainID, &resSeq, &iCode,
            &do_rot,  &rotations,
            &do_sw,   &phi_swing,
            &opt_hyd, &n_hyd_pos,
            &nconf_limit );
            
            /*
            printf("%s %c%04d%c R %c %02d S %c %4.1f H %c %02d M %03d\n",
            resName,chainID,resSeq,iCode,
            do_rot,  rotations,
            do_sw,   phi_swing,
            opt_hyd, n_hyd_pos,
            nconf_limit);
            */
            
            for (i_res = 0; i_res < prot.n_res; i_res++) {
                if (!strcmp(resName, prot.res[i_res].resName) &&
                    chainID == prot.res[i_res].chainID &&
                    resSeq  == prot.res[i_res].resSeq  &&
                    iCode   == prot.res[i_res].iCode  )
                {
                    break;
                }
            }
            if (i_res < prot.n_res) {
                if (do_rot == 't')  prot.res[i_res].do_rot   = 1; else prot.res[i_res].do_rot   = 0;
                if (do_sw  == 't')   prot.res[i_res].do_sw   = 1; else prot.res[i_res].do_sw    = 0;
                if (opt_hyd == 't') prot.res[i_res].opt_hyd  = 1; else prot.res[i_res].opt_hyd  = 0;
                prot.res[i_res].rotations  = rotations;
                prot.res[i_res].phi_swing  = phi_swing;
                prot.res[i_res].n_hyd_pos  = n_hyd_pos;
                prot.res[i_res].nconf_limit = nconf_limit;
                
                if (!prot.res[i_res].rotations)         prot.res[i_res].do_rot = 0;
                if ( prot.res[i_res].phi_swing < 1e-4)  prot.res[i_res].do_sw = 0;
                if (!prot.res[i_res].opt_hyd)           prot.res[i_res].n_hyd_pos = 0;
            }
        }
        fclose(fp);
    }
}

void write_headlst(PROT prot) {
    FILE *fp;
    int i_res;
    char do_rot,do_sw,opt_hyd;
    
    fp = fopen(FN_CONFLIST1,"w");
    for (i_res=0; i_res<prot.n_res; i_res++) {
        if (prot.res[i_res].do_rot)  do_rot  = 't'; else do_rot  = 'f';
        if (prot.res[i_res].do_sw)   do_sw   = 't'; else do_sw   = 'f';
        if (prot.res[i_res].opt_hyd) opt_hyd = 't'; else opt_hyd = 'f';
        
        fprintf(fp,"%s %c%04d%c R %c %02d S %c %4.1f H %c %02d M %03d\n",
        prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].iCode,
        do_rot,  prot.res[i_res].rotations,
        do_sw,   prot.res[i_res].phi_swing,
        opt_hyd, prot.res[i_res].n_hyd_pos,
        prot.res[i_res].nconf_limit);
    }
    fclose(fp);
}


/*
{
    Delete conformer if not converged 
    if (convg == -1) {
        for (j_ngh=0; j_ngh<res_extra[i_res].n_nghs; j_ngh++) {
            RES *jres_p;
            j_res = res_extra[i_res].ngh_res[j_ngh].j_res;
            if(!res_extra[j_res].hyd_opt && !res_extra[j_res].hoh_opt) continue;
            j_res = res_extra[i_res].ngh_res[j_ngh].j_res;
            jres_p = &prot.res[j_res];
            del_conf(jres_p, jres_p->n_conf-1);
        }
        //printf("Fail trial No. %4d\n",n_trial);
        continue;
    }
    
    
    Print out interactions 
    check_vdw = 0, check_coulomb=0;
    //printf("\n");
    for (j_ngh=0; j_ngh<res_extra[i_res].n_nghs; j_ngh++) {
        j_res = res_extra[i_res].ngh_res[j_ngh].j_res;
        j_conf = res_extra[j_res].j_conf_on;
        int k_res, k_conf, k_ngh;
        for (k_ngh=0; k_ngh<res_extra[i_res].n_nghs; k_ngh++) {
            k_res = res_extra[i_res].ngh_res[k_ngh].j_res;
            k_conf = res_extra[k_res].j_conf_on;
            //printf("%8.3f,%8.3f",vdw_conf(j_res,j_conf,k_res,k_conf,prot),coulomb_conf(j_res,j_conf,k_res,k_conf,prot));
            check_vdw += vdw_conf(j_res,j_conf,k_res,k_conf,prot);
            check_coulomb += coulomb_conf(j_res,j_conf,k_res,k_conf,prot);
        }
        //printf("\n");
    }
    printf("after %3d steps of relaxation vdw %8.3f, coulomb %8.3f\n", convg, check_vdw/2.,check_coulomb/2.);
    
    
    Check energy 
    check_vdw = 0, check_coulomb=0;
    for (i_atom=0; i_atom< n_relax; i_atom++) {
        if (!relax_atoms[i_atom].atom_p->on) continue;
        if (!relax_atoms[i_atom].move) continue;
        for (j_atom=0; j_atom< n_relax; j_atom++) {
            if (i_atom == j_atom) continue;
            if (!relax_atoms[j_atom].atom_p->on) continue;
            if (!relax_atoms[j_atom].move) {
                check_vdw += vdw(*relax_atoms[i_atom].atom_p,*relax_atoms[j_atom].atom_p)*relax_atoms[i_atom].factor[j_atom];
                check_coulomb += coulomb(*relax_atoms[i_atom].atom_p,*relax_atoms[j_atom].atom_p)*relax_atoms[i_atom].factor[j_atom];
            }
            else {
                check_vdw += vdw(*relax_atoms[i_atom].atom_p,*relax_atoms[j_atom].atom_p)*relax_atoms[i_atom].factor[j_atom]/2.;
                check_coulomb += coulomb(*relax_atoms[i_atom].atom_p,*relax_atoms[j_atom].atom_p)*relax_atoms[i_atom].factor[j_atom]/2.;
            }
            
            if (vdw(*relax_atoms[i_atom].atom_p,*relax_atoms[j_atom].atom_p)*relax_atoms[i_atom].factor[j_atom] > 10e2 ||
                coulomb(*relax_atoms[i_atom].atom_p,*relax_atoms[j_atom].atom_p)*relax_atoms[i_atom].factor[j_atom] > 10e2)
            printf("iatom%d,%s,jatom%d,%s,vdw%8.3f,coulomb%8.3f,factor%f\n",i_atom,relax_atoms[i_atom].atom_p->name,j_atom,relax_atoms[j_atom].atom_p->name,check_vdw,check_coulomb,relax_atoms[i_atom].factor[j_atom]);
            
        }
    }
    tot_crg = 0;
    for (i_atom=0; i_atom< n_relax; i_atom++) {
        if (!relax_atoms[i_atom].atom_p->on) continue;
        tot_crg += relax_atoms[i_atom].atom_p->crg;
    }
    //printf("Trial %d, after %3d steps of relaxation vdw %8.3f, coulomb %8.3f,tot_crg %8.3f\n", n_trial,convg, check_vdw,check_coulomb,tot_crg);
    if (env.relax_e_thr>0) relax_e_thr = tot_crg * env.relax_e_thr;
    else relax_e_thr = -tot_crg * env.relax_e_thr;
    if ( check_vdw+check_coulomb > relax_e_thr) {
        for (j_ngh=0; j_ngh<res_extra[i_res].n_nghs; j_ngh++) {
            RES *jres_p;
            j_res = res_extra[i_res].ngh_res[j_ngh].j_res;
            jres_p = &prot.res[j_res];
            if(!res_extra[j_res].hyd_opt && !res_extra[j_res].hoh_opt) continue;
            del_conf(jres_p, jres_p->n_conf-1);
        }
        
        for (i_atom=0; i_atom< n_relax; i_atom++) {
            if (relax_atoms[i_atom].move) {
                free(relax_atoms[i_atom].C6);
                free(relax_atoms[i_atom].C12);
                free(relax_atoms[i_atom].factor);
            }
        }
        continue;
    }
    
    for (j_ngh=0; j_ngh<res_extra[i_res].n_nghs; j_ngh++) {
        RES *jres_p;
        int j_conf;
        j_res = res_extra[i_res].ngh_res[j_ngh].j_res;
        if(!res_extra[j_res].hyd_opt && !res_extra[j_res].hoh_opt) continue;
        jres_p = &prot.res[j_res];
        
        for (j_conf=1; j_conf<jres_p->n_conf-1; j_conf++) {
            if (!cmp_conf(jres_p->conf[j_conf], jres_p->conf[jres_p->n_conf-1],0.5)) break;
        }
        if ( j_conf < jres_p->n_conf-1 ) {
            //for (i_atom=0; i_atom< jres_p->conf[j_conf].n_atom; i_atom++) {if (!jres_p->conf[j_conf].atom[i_atom].on) continue;printf("%s %8.3f,%8.3f,%8.3f,%s %8.3f,%8.3f,%8.3f\n",jres_p->conf[j_conf].atom[i_atom].name,jres_p->conf[j_conf].atom[i_atom].xyz.x,jres_p->conf[j_conf].atom[i_atom].xyz.y,jres_p->conf[j_conf].atom[i_atom].xyz.z,jres_p->conf[jres_p->n_conf-1].atom[i_atom].name,jres_p->conf[jres_p->n_conf-1].atom[i_atom].xyz.x,jres_p->conf[jres_p->n_conf-1].atom[i_atom].xyz.y,jres_p->conf[jres_p->n_conf-1].atom[i_atom].xyz.z);}
            del_conf(jres_p, jres_p->n_conf-1);
            continue;
        }
    }
    for (i_atom=0; i_atom< n_relax; i_atom++) {
        if (relax_atoms[i_atom].move) {
            free(relax_atoms[i_atom].C6);
            free(relax_atoms[i_atom].C12);
            free(relax_atoms[i_atom].factor);
        }
    }
}
*/



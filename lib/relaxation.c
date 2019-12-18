#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mcce.h"

#define  PRINT_THR        9
#define  contact2        16
#define  CONVERGED_R     1e-5

typedef struct RELAX_STRUCT{
    int    i_relax;
    ATOM   *atom_p;
    
    VECTOR r;
    VECTOR r_p;       /* r' */
    //VECTOR r_1step_back;  /* r of 1 step  back */
    VECTOR r_orig;
    
    int tor_atom1;
    int tor_atom2;
    int tor_atom3;
    TORS tors;
    int n_rotate1, n_rotate2;
    int *rotate1_lst, *rotate2_lst;
    
    int    n_ngh;
    struct RELAX_STRUCT **ngh;
    int    n_ngh14;
    struct RELAX_STRUCT **ngh14;
    int    n_constr;
    struct RELAX_STRUCT **constr_list;
    double  *constr_dsq;
    
    VECTOR lj_frc;
    VECTOR elec_frc;
    VECTOR torsion_frc;
    VECTOR constr_frc;
    VECTOR frc;
    int    movable;
    int    moving;
    int    moved;
} RELAX;

void relaxation_setup(PROT prot);
void complete_constr(int i_res, int i_conf, int j_res, int j_conf);
//int  in_relax_list(int ia);
void setup_nghlst(PROT prot);
void get_frc(float tors_scale, PROT prot);
void get_rp();
int  shake();
void add_conf2relax(CONF *conf_p, int fix);
void add_atom2relax(int ia, int fix);
int res_in_relax(int i_res);
int pick_sidechain(int i_res, PROT prot);
void collect_ngh(int i_res, PROT prot);
int closer_than(int i_res, int i_conf, int j_res, int j_conf, PROT prot, float crg_thr, float dist_thr);
int write_relax_pdb(char *relax_pdb_filename, PROT prot);
extern int ionization(PROT prot);
extern int rm_dupconf_res(PROT prot, int i_res, float prune_thr);
extern int rm_dupconf(PROT prot, float prune_thr);
extern int rm_dupconf_hv(PROT prot);
extern int delete_h(PROT prot);

int     na;
RELAX   *all_atoms;
int     n_relax;
RELAX   *relax_atoms;
int     n_relax_res;
int     relax_res[10000];
int     relax_conf[10000];
VECTOR  zero;
int     n_pair_trial=0, n_pair_converged=0;

//int     n_elem;
//float   C6_matrix[N_ELEM_MAX][N_ELEM_MAX];
//float   C12_matrix[N_ELEM_MAX][N_ELEM_MAX];
float   DTSQ, CONSTRAINT2, CONSTRAINT_FRC;
extern  long    idum;

//float   stepwise2_max;

int relaxation(PROT prot)
{   
    int n_pair_relaxed, i_cycle, done;
    int i_res, i_conf, j_res, j_conf;
    int i_relax, ia, i_relax_step;
    float pair_vdw, pair_vdw_hv;
    RES  *ires_p, *jres_p;
    CONF *iconf_p, *jconf_p;
    FILE *progress_fp;
    float dsq,dsq_max,sum_dsq;
    int   n_movable;
    time_t   timer_start, timer_end;
    
    zero.x = 0.; zero.y = 0.; zero.z = 0.;
    //stepwise2_max = 0.;
    //DTSQ = env.hv_relax_dt*env.hv_relax_dt*0.4187e-6;
    CONSTRAINT2 = env.hv_relax_constraint * env.hv_relax_constraint;
    CONSTRAINT_FRC = env.hv_relax_constraint_frc;
    
    /* setup relaxation */
    printf("   Start setting up for relaxation.\n"); fflush(stdout);
    relaxation_setup(prot);
    id_conf(prot);
    printf("   Setup for relaxation done.\n"); fflush(stdout);
    
    n_pair_relaxed = 1;
    i_cycle = 0;
    //int outSeq = 0;
    printf("   Start relaxation.\n"); fflush(stdout);
    while (n_pair_relaxed) {
        i_cycle++;
        if (i_cycle > env.hv_relax_ncycle) break;
        
        timer_start = time(NULL);
        setup_vdw_fast(prot);
        
        n_pair_relaxed = 0;
        n_pair_converged = 0;
        n_pair_trial = 0;
        
        int n_pair_failed = 0;
        float sum_pair_vdw_failed = 0.;
        float sum_pair_vdw_hv_failed = 0.;
        float max_pair_vdw_failed = -9999.;
        float max_pair_vdw_hv_failed = -9999.;
        float sumsq_pair_vdw_failed = 0.;
        float sumsq_pair_vdw_hv_failed = 0.;

        int n_pair_success = 0;
        float sum_pair_vdw_success = 0.;
        float sum_pair_vdw_hv_success = 0.;
        float max_pair_vdw_success = -9999.;
        float max_pair_vdw_hv_success = -9999.;
        float sumsq_pair_vdw_success = 0.;
        float sumsq_pair_vdw_hv_success = 0.;
        
        for (i_res=0; i_res<prot.n_res; i_res++) {
            ires_p = &prot.res[i_res];
            for (j_res=0; j_res<prot.n_res; j_res++) {
                jres_p = &prot.res[j_res];
                if (out_of_range(ires_p->r_min, ires_p->r_max, jres_p->r_min, jres_p->r_max, contact2)) continue;
                
                for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
                    if (!prot.res[i_res].conf[i_conf].n_atom) continue;
                    iconf_p = &prot.res[i_res].conf[i_conf];
                    for (j_conf=0; j_conf<prot.res[j_res].n_conf; j_conf++) {
                        char relax_flag = 't';
                        if (!prot.res[j_res].conf[j_conf].n_atom) continue;
                        jconf_p = &prot.res[j_res].conf[j_conf];
                        if (i_res == j_res && i_conf == j_conf) continue; /* same conformer */
                        if (out_of_range(iconf_p->r_min, iconf_p->r_max, jconf_p->r_min, jconf_p->r_max, contact2)) continue;
                        
                        if ((!prot.res[i_res].relax||!i_conf) && (!prot.res[j_res].relax||!j_conf)) continue; /* Skip if both conformers are not relaxable */
                        if (i_res == j_res && (i_conf*j_conf != 0 && i_conf != j_conf)) continue; /* Skip if same residue, two different sidechain */
                        
                        if (!strncmp(prot.res[i_res].conf[i_conf].history+2, "O000", 4))
                            if (!strncmp(prot.res[j_res].conf[j_conf].history+2, "O000", 4))
                                continue;
                        
                        pair_vdw    = vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 0)
                        +vdw_conf_fast(i_res, i_conf, i_res, i_conf, prot, 0)
                        +torsion_conf(&prot.res[i_res].conf[i_conf])
                        +vdw_conf_fast(j_res, j_conf, j_res, j_conf, prot, 0)
                        +torsion_conf(&prot.res[j_res].conf[j_conf]);
                        
                        pair_vdw_hv = vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 1)
                        +vdw_conf_fast(i_res, i_conf, i_res, i_conf, prot, 1)
                        +torsion_conf(&prot.res[i_res].conf[i_conf])
                        +vdw_conf_fast(j_res, j_conf, j_res, j_conf, prot, 1)
                        +torsion_conf(&prot.res[j_res].conf[j_conf]);
                        
                        if (i_conf*j_conf != 0) {/* both are sidechain (if one of them is backbone then always relax */
                            if (vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 1) > env.hv_relax_hv_vdw_thr) {
                                /* for sidechains, do not relax if heavy atom clashes */
                                relax_flag = 'f';
                            }
                        }
                        
                        if (coulomb_conf(i_res, i_conf, j_res, j_conf, prot) < env.hv_relax_elec_thr) {
                            /* If there is a favorable electrostatic interaction,
                            check the distance between two conf, but do not check vdw crash */
                            
                            /* if the charge region of two conf are too far away, no relaxation*/
                            if ( !closer_than(i_res, i_conf, j_res, j_conf, prot, env.hv_relax_elec_crg_thr, env.hv_relax_elec_dist_thr) ) {
                                relax_flag = 'f';
                            }
                        }
                        else {
                            if (vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 0) < env.hv_relax_vdw_thr) {
                                relax_flag = 'f';
                            }
                        }
                        
                        if (i_res == j_res) {
                            if (i_conf != 0 && prot.res[i_res].conf[i_conf].history[2] == 'E') relax_flag = 't';
                            if (j_conf != 0 && prot.res[j_res].conf[j_conf].history[2] == 'E') relax_flag = 't';
                        }
                        
                        //printf("%d, %d, %d, %d\n",i_res,i_conf,j_res,j_conf);
                        
                        /*
                        if (pair_vdw >= 100.) {
                            progress_fp = fopen(env.progress_log,"a");
                            if (progress_fp) {
                                fprintf(progress_fp, "   Relaxation: this pair %s-%s has a large vdw pairwise, pair_vdw = %6.3f, pair_vdw_hv = %6.3f\n", iconf_p->uniqID, jconf_p->uniqID, pair_vdw, pair_vdw_hv);
                                fclose(progress_fp);
                            }
                        }
                        int pr_toggle = 0;
                        if (vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 0) > 5.) {
                            pr_toggle = 1;
                        }
                        */
                        
                        if (relax_flag == 'f') continue;
                        
                        /* setup relax_atoms */
                        n_relax=0;
                        relax_atoms = NULL;
                        n_relax_res = 0;
                        
                        add_conf2relax(&prot.res[i_res].conf[0], 0);
                        if (i_conf) add_conf2relax(&prot.res[i_res].conf[i_conf], 0);
                        
                        if (i_res != j_res) {
                            add_conf2relax(&prot.res[j_res].conf[0], 0);
                            if (j_conf) add_conf2relax(&prot.res[j_res].conf[j_conf], 0);
                        }
                        else if (j_conf && i_conf != j_conf) {
                            add_conf2relax(&prot.res[j_res].conf[j_conf], 0);
                        }
                        if (env.hv_relax_include_ngh) {
                            /* complete i_res and j_res if their sidechain is not in relax */
                            if (!res_in_relax(i_res)) {
                                if (prot.res[i_res].n_conf > 1) {
                                    i_conf = pick_sidechain(i_res, prot);
                                    add_conf2relax(&prot.res[i_res].conf[i_conf], 1);
                                }
                                else {
                                    n_relax_res++;
                                    relax_res[n_relax_res-1] = i_res;
                                    relax_conf[n_relax_res-1] = 0;
                                }
                            }
                            if (!res_in_relax(j_res)) {
                                if (prot.res[j_res].n_conf > 1) {
                                    j_conf = pick_sidechain(j_res, prot);
                                    add_conf2relax(&prot.res[j_res].conf[j_conf], 1);
                                }
                                else {
                                    n_relax_res++;
                                    relax_res[n_relax_res-1] = j_res;
                                    relax_conf[n_relax_res-1] = 0;
                                }
                            }
                            /* collect ngh */
                            collect_ngh(i_res, prot);
                            if (i_res!=j_res) collect_ngh(j_res, prot);
                        }
                        complete_constr(i_res, i_conf, j_res, j_conf);
                        setup_nghlst(prot);
                        
                        for (i_relax =0; i_relax<n_relax; i_relax++) {
                            int j_relax;
                            for (j_relax =0; j_relax<n_relax; j_relax++) {
                                if (i_relax == j_relax) continue;
                                if (relax_atoms[i_relax].atom_p == relax_atoms[j_relax].atom_p) {
                                    printf("\nError! same atoms in relaxation.\n");
                                    printf("%d,%d\n",relax_atoms[i_relax].atom_p->i_res_prot,relax_atoms[i_relax].atom_p->i_conf_res);
                                    printf("%d,%d,%d,%d\n",i_res,i_conf,j_res,j_conf);
                                }
                                if (relax_atoms[i_relax].atom_p->i_res_prot == relax_atoms[j_relax].atom_p->i_res_prot) {
                                    if (relax_atoms[i_relax].atom_p->i_conf_res != relax_atoms[j_relax].atom_p->i_conf_res) {
                                        if (relax_atoms[i_relax].atom_p->i_conf_res * relax_atoms[j_relax].atom_p->i_conf_res) {
                                            printf("\nError! different side chain.\n");
                                        }
                                    }
                                }
                            }
                        }
                        
                        /* start relaxation */
                        DTSQ = 10.*env.hv_relax_dt*env.hv_relax_dt*0.4187e-6;
                        done = 0;
                        while (!done) {
                            n_pair_trial++;
                            i_relax_step = env.hv_relax_niter;
                            while (i_relax_step) {
                                float tors_scale = 1.+(env.hv_tors_scale - 1.)*(float)(i_relax_step-1)/(float)env.hv_relax_niter;
                                
                                get_frc(tors_scale, prot);
                                get_rp();
                                if (shake()) break;
                                
                                for (i_relax =0; i_relax<n_relax; i_relax++) {
                                    /*
                                    if (ddvv(relax_atoms[i_relax].r,relax_atoms[i_relax].r_p) > stepwise2_max) {
                                        stepwise2_max = ddvv(relax_atoms[i_relax].r,relax_atoms[i_relax].r_p);
                                        printf("stepwise %6.3f\n",sqrt(stepwise2_max));
                                    }
                                    */
                                    //printf("%d, %s, %13.10f\n",i_relax_step, relax_atoms[i_relax].atom_p->name, dvv(relax_atoms[i_relax].r, relax_atoms[i_relax].r_p));
                                    relax_atoms[i_relax].r = relax_atoms[i_relax].r_p;
                                }
                                
                                int converged = 1;
                                for (i_relax =0; i_relax<n_relax; i_relax++) {
                                    if (converged) {
                                        if (dvv(relax_atoms[i_relax].r, relax_atoms[i_relax].atom_p->xyz) > CONVERGED_R*env.hv_relax_dt*(env.hv_relax_niter-i_relax_step)) converged = 0;
                                    }
                                }                                
                                if (converged) {
                                    //printf("%d\n",i_relax_step);
                                    i_relax_step = 0;
                                }
                                else {
                                    i_relax_step--;
                                }
                            }
                            if (i_relax_step) {
                                /*
                                progress_fp = fopen(env.progress_log,"a");
                                if (progress_fp) {
                                    fprintf(progress_fp, "   Relaxation: this pair %s-%s failed in shake convergence, pair_vdw = %6.3f, pair_vdw_hv = %6.3f\n", iconf_p->uniqID, jconf_p->uniqID, pair_vdw, pair_vdw_hv);
                                    fclose(progress_fp);
                                }

                                char pair_pr[160];
                                sprintf(pair_pr,"%s-%s.pdb",iconf_p->uniqID, jconf_p->uniqID);
                                FILE *pair_fp = fopen(pair_pr,"w");
                                if (pair_fp) {
                                    for (i_relax =0; i_relax<n_relax; i_relax++) {
                                        int i_res_pr, i_conf_pr, i_atom_pr;
                                        i_res_pr = relax_atoms[i_relax].atom_p->i_res_prot;
                                        i_conf_pr = relax_atoms[i_relax].atom_p->i_conf_res;
                                        i_atom_pr = relax_atoms[i_relax].atom_p->i_atom_conf;
                                        fprintf(pair_fp,"ATOM        %4s%c%3s %c%04d%c%03d%8.3f%8.3f%8.3f %8.3f%8.3f%8.3f\n",
                                        prot.res[i_res_pr].conf[i_conf_pr].atom[i_atom_pr].name,
                                        prot.res[i_res_pr].conf[i_conf_pr].altLoc,
                                        prot.res[i_res_pr].resName,
                                        prot.res[i_res_pr].chainID,
                                        prot.res[i_res_pr].resSeq,
                                        prot.res[i_res_pr].iCode,
                                        i_conf_pr,
                                        prot.res[i_res_pr].conf[i_conf_pr].atom[i_atom_pr].xyz.x,
                                        prot.res[i_res_pr].conf[i_conf_pr].atom[i_atom_pr].xyz.y,
                                        prot.res[i_res_pr].conf[i_conf_pr].atom[i_atom_pr].xyz.z,
                                        relax_atoms[i_relax].r.x,
                                        relax_atoms[i_relax].r.y,
                                        relax_atoms[i_relax].r.z);
                                    }
                                    fclose(pair_fp);
                                }
                                */

                                if (DTSQ < 1e-10) {
                                    /* relaxation failed */
                                    n_pair_failed++;
                                    
                                    sum_pair_vdw_failed += pair_vdw;
                                    sum_pair_vdw_hv_failed += pair_vdw_hv;

                                    sumsq_pair_vdw_failed += pair_vdw*pair_vdw;
                                    sumsq_pair_vdw_hv_failed += pair_vdw_hv*pair_vdw_hv;
                                    
                                    if (pair_vdw > max_pair_vdw_failed) max_pair_vdw_failed = pair_vdw;
                                    if (pair_vdw_hv > max_pair_vdw_hv_failed) max_pair_vdw_hv_failed = pair_vdw_hv;

                                    for (i_relax =0; i_relax<n_relax; i_relax++) {
                                        relax_atoms[i_relax].r = relax_atoms[i_relax].atom_p->xyz;
                                    }

                                    break;
                                }
                                DTSQ = DTSQ /2.;
                                for (i_relax =0; i_relax<n_relax; i_relax++) {
                                    relax_atoms[i_relax].r = relax_atoms[i_relax].atom_p->xyz;
                                }
                            }
                            else
                                done = 1;
                        }
                        
                        /* success */
                        n_pair_success++;
                        
                        sum_pair_vdw_success += pair_vdw;
                        sum_pair_vdw_hv_success += pair_vdw_hv;
                        
                        sumsq_pair_vdw_success += pair_vdw*pair_vdw;
                        sumsq_pair_vdw_hv_success += pair_vdw_hv*pair_vdw_hv;
                        
                        if (pair_vdw > max_pair_vdw_success) max_pair_vdw_success = pair_vdw;
                        if (pair_vdw_hv > max_pair_vdw_hv_success) max_pair_vdw_hv_success = pair_vdw_hv;
                        
                        for (i_relax =0; i_relax<n_relax; i_relax++) {
                            relax_atoms[i_relax].atom_p->xyz = relax_atoms[i_relax].r;
                        }
                        for (i_relax =0; i_relax<n_relax; i_relax++) {
                            ia = relax_atoms[i_relax].atom_p->i_atom_prot;
                            all_atoms[ia].i_relax = -1;
                        }

                        float relaxed_pair_vdw    = vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 0)
                        +vdw_conf_fast(i_res, i_conf, i_res, i_conf, prot, 0)
                        +torsion_conf(&prot.res[i_res].conf[i_conf])
                        +vdw_conf_fast(j_res, j_conf, j_res, j_conf, prot, 0)
                        +torsion_conf(&prot.res[j_res].conf[j_conf]);
                        
                        float relaxed_pair_vdw_hv = vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 1)
                        +vdw_conf_fast(i_res, i_conf, i_res, i_conf, prot, 1)
                        +torsion_conf(&prot.res[i_res].conf[i_conf])
                        +vdw_conf_fast(j_res, j_conf, j_res, j_conf, prot, 1)
                        +torsion_conf(&prot.res[j_res].conf[j_conf]);

                        if (pair_vdw > 100.) {
                            progress_fp = fopen(env.progress_log,"a");
                            fprintf(progress_fp,"   Relaxation: pair %s-%s before relaxation pair_vdw = %8.1f, pair_vdw_hv = %8.1f\n",iconf_p->uniqID, jconf_p->uniqID, pair_vdw, pair_vdw_hv);
                            fprintf(progress_fp,"   Relaxation: pair %s-%s after             pair_vdw = %8.1f, pair_vdw_hv = %8.1f\n",iconf_p->uniqID, jconf_p->uniqID, relaxed_pair_vdw, relaxed_pair_vdw_hv);
                            fclose(progress_fp);
                        }
                        
                        
                        /*
                        if (pr_toggle) {
                            progress_fp = fopen(env.progress_log,"a");
                            if (progress_fp) {
                                fprintf(progress_fp, "   Relaxation: this pair %s-%s before relaxation pair_vdw = %6.3f, pair_vdw_hv = %6.3f, after pair_vdw = %6.3f, pair_vdw_hv = %6.3f\n", iconf_p->uniqID, jconf_p->uniqID, pair_vdw, pair_vdw_hv, relaxed_pair_vdw, relaxed_pair_vdw_hv);
                                fclose(progress_fp);
                            }
                        }
                        */
                        for (i_relax =0; i_relax<n_relax; i_relax++) {
                            free(relax_atoms[i_relax].ngh);
                            free(relax_atoms[i_relax].ngh14);
                            free(relax_atoms[i_relax].constr_list);
                            free(relax_atoms[i_relax].constr_dsq);
                            free(relax_atoms[i_relax].rotate1_lst);
                            free(relax_atoms[i_relax].rotate2_lst);
                        }
                        free(relax_atoms);
                        n_pair_relaxed++;
                        n_pair_converged++;
                        
                        /* done with this pair of conformers */
                    }
                }
            }
        }
        
        dsq_max = 0;
        sum_dsq = 0.;
        n_movable = 0;
        for (ia=0;ia<na;ia++) {
            if (!all_atoms[ia].movable) continue;
            dsq = ddvv(all_atoms[ia].atom_p->xyz, all_atoms[ia].r_orig);
            sum_dsq += dsq;
            n_movable++;
            if (dsq > dsq_max) dsq_max = dsq;
        }
        float avg_pair_vdw_failed = sum_pair_vdw_failed/(float)n_pair_failed;
        float dev_pair_vdw_failed = sqrt(sumsq_pair_vdw_failed/(float)n_pair_failed - avg_pair_vdw_failed*avg_pair_vdw_failed);
        float avg_pair_vdw_success = sum_pair_vdw_success/(float)n_pair_success;
        float dev_pair_vdw_success = sqrt(sumsq_pair_vdw_success/(float)n_pair_success - avg_pair_vdw_success*avg_pair_vdw_success);
        
        timer_end = time(NULL);
        progress_fp = fopen(env.progress_log, "a");
        if (progress_fp) {
            fprintf(progress_fp, "   Relaxation Cycle = %3d, n_relaxed = %5d, shake convergence rate = %6.2f %%, success rate = %6.2f %%\n",i_cycle,n_pair_relaxed, 100.*(float)n_pair_converged/(float)n_pair_trial, 100.*(float)n_pair_converged/(float)n_pair_relaxed); fflush(stdout);
            fprintf(progress_fp, "   Relaxation Cycle = %3d, rmsd = %6.3f, max displacement = %6.3f, time used %ld seconds\n",i_cycle,sqrt(sum_dsq/(float)n_movable),sqrt(dsq_max),timer_end-timer_start); fflush(stdout);
            fprintf(progress_fp, "   Relaxation Cycle = %3d, # of failed  = %5d, avg_failed_vdw  = %6.3f, deviation = %6.3f, maximum = %6.3f\n",i_cycle,n_pair_failed,avg_pair_vdw_failed,dev_pair_vdw_failed,max_pair_vdw_failed); fflush(stdout);
            fprintf(progress_fp, "   Relaxation Cycle = %3d, # of success = %5d, avg_success_vdw = %6.3f, deviation = %6.3f, maximum = %6.3f\n",i_cycle,n_pair_success,avg_pair_vdw_success,dev_pair_vdw_success,max_pair_vdw_success); fflush(stdout);
            fclose(progress_fp);
        }
    }

    for (i_res=0; i_res<prot.n_res; i_res++) {
        free_connect_res(prot, i_res);
    }

    return 0;
}

int initial_relaxation(PROT prot)
{
    FILE *fp;
    int i_cycle, done;
    int i_res, i_conf;
    int i_relax, ia, i_relax_step;
    RES  *ires_p;
    time_t   timer_start, timer_end;

    timer_start = time(NULL);
    
    zero.x = 0.; zero.y = 0.; zero.z = 0.;
    
    //stepwise2_max = 0.;
    //DTSQ = env.hv_relax_dt*env.hv_relax_dt*0.4187e-6;
    
    CONSTRAINT2 = env.hv_relax_constraint * env.hv_relax_constraint;
    CONSTRAINT_FRC = env.hv_relax_constraint_frc;
    
    /* setup relaxation */
    //printf("   Start setting up for relaxation.\n"); fflush(stdout);

    /* rebuild sidechain */
    for (i_res=0; i_res<prot.n_res; ++i_res) {
        char copy_atoms[MAXCHAR_LINE];
        //int n_conf = prot.res[i_res].n_conf;
        
        if ( param_get("REBUILD", prot.res[i_res].resName, "", copy_atoms) ) {
            /* skip if no REBUILD parameter found for this residue */
            continue;
        }
        /* pad copy_atoms string with a few spaces */
        strcpy(copy_atoms+strlen(copy_atoms), "    ");
        
        if (prot.res[i_res].n_conf > 1) {
            int i_atom;
            /* make a copy of the first sidechain */
            int k_conf = prot.res[i_res].n_conf;
            ins_conf(&prot.res[i_res], k_conf, prot.res[i_res].conf[1].n_atom);
            cpy_conf(&prot.res[i_res].conf[k_conf], &prot.res[i_res].conf[1]);
            
            prot.res[i_res].conf[k_conf].history[2] = 'I'; /* Label as initial */

            for (i_atom=0; i_atom<prot.res[i_res].conf[k_conf].n_atom; i_atom++) {
                if (!strstr(copy_atoms, prot.res[i_res].conf[k_conf].atom[i_atom].name)) {
                    /* if atom name does not match those defined in parameter, then do not keep */
                    prot.res[i_res].conf[k_conf].atom[i_atom].on = 0;
                }
            }
            
            while(place_missing(prot, 1) > 0); /* handle_addconf=2 adds more conformers than suggested by torsion parameter */
            //printf("%s%4d n_conf=%6d\n",prot.res[i_res].resName, prot.res[i_res].resSeq, prot.res[i_res].n_conf);
            //rm_dupconf_res(prot, i_res, 0.001);
            //printf("%s%4d n_conf=%6d\n",prot.res[i_res].resName, prot.res[i_res].resSeq, prot.res[i_res].n_conf);
            
            /* find which one is closest to the original conf */
            int   closest_kconf = 0;
            float closest_dist = 999.;
            for (i_conf = 1; i_conf < prot.res[i_res].n_conf; i_conf++) {
                if (prot.res[i_res].conf[i_conf].history[2] != 'I') continue;
                if (closest_dist > rmsd_conf_hv(prot.res[i_res].conf[i_conf], prot.res[i_res].conf[1])) {
                    closest_kconf = i_conf;
                    closest_dist = rmsd_conf_hv(prot.res[i_res].conf[i_conf],prot.res[i_res].conf[1]);
                }
            }
            //printf("%s%4d n_conf=%6d %4d %8.3f\n",prot.res[i_res].resName, prot.res[i_res].resSeq, prot.res[i_res].n_conf, closest_kconf, closest_dist);
            
            /* keep the closest one and delete the rest */
            for (i_conf = prot.res[i_res].n_conf-1; i_conf >= 1; i_conf--) {
                if (prot.res[i_res].conf[i_conf].history[2] != 'I') continue;
                if (i_conf == closest_kconf) continue;
                del_conf(&prot.res[i_res], i_conf);
            }
        }
    }
    //printf("   01-Current time is %ds\n",(int)(time(NULL)-timer_start));

    for (i_res=0; i_res<prot.n_res; ++i_res) {
        char copy_atoms[MAXCHAR_LINE];
        if ( param_get("REBUILD", prot.res[i_res].resName, "", copy_atoms) ) {
            /* skip if no REBUILD parameter found for this residue */
            continue;
        }
        if (prot.res[i_res].n_conf <= 2) continue;

    }
    delete_h(prot);
    //printf("   02-Current time is %ds\n",(int)(time(NULL)-timer_start));

    /* Make ionization states */
    if (ionization(prot)) {
        printf("   FATAL: Fatal error reported by ionization()\n"); fflush(stdout);
        return USERERR;
    }
    while(place_missing(prot,1) > 0); rm_dupconf(prot, 0.0001);
    //printf("   Current time is %ds\n",(int)(time(NULL)-timer_start));

    /* which residues are allowed to relax: */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        if (param_get("RELAX", prot.res[i_res].resName, "", &prot.res[i_res].relax)) {
            prot.res[i_res].relax = 0;
        }
    }
    //printf("   Current time is %ds\n",(int)(time(NULL)-timer_start));

    /* add conformers for relaxtion */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        int n_conf;
        char def_conf_list[MAXCHAR_LINE];
        
        ires_p = &prot.res[i_res];
        n_conf = prot.res[i_res].n_conf;
        
        int with_rebuild = 0;
        for (i_conf=1; i_conf<n_conf; i_conf++) {
            if (prot.res[i_res].conf[i_conf].history[2] == 'I') with_rebuild = 1;
        }
        
        for (i_conf=1; i_conf<n_conf; i_conf++) {
            int k_conf;
            if (!param_get("DEF_CONF", prot.res[i_res].resName, "", def_conf_list)) {
                if (!strstr(def_conf_list, prot.res[i_res].conf[i_conf].confName)) {
                    continue;
                }
            }
            if (with_rebuild) {
                if (prot.res[i_res].conf[i_conf].history[2] != 'I') continue;
            }
            
            k_conf = prot.res[i_res].n_conf;
            ins_conf(&prot.res[i_res], k_conf, prot.res[i_res].conf[i_conf].n_atom);
            cpy_conf(&prot.res[i_res].conf[k_conf], &prot.res[i_res].conf[i_conf]);
            prot.res[i_res].conf[k_conf].history[2] = 'X';
        }
        
        /* If no 'X' conformer made, make one from 'O' conformer */
        /* This only happens occasionally when 'I' is very close to 'O' conf and being removed by rm_dupconf() func. */
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            if (prot.res[i_res].conf[i_conf].history[2] == 'X') break;
        }
        if (i_conf>=prot.res[i_res].n_conf && prot.res[i_res].n_conf!=1) {
            int j_conf;
            for (j_conf=1; j_conf<n_conf; j_conf++) {
                int k_conf;
                if (!param_get("DEF_CONF", prot.res[i_res].resName, "", def_conf_list)) {
                    if (!strstr(def_conf_list, prot.res[i_res].conf[j_conf].confName)) {
                        continue;
                    }
                }
                if (prot.res[i_res].conf[j_conf].history[2] != 'O') continue;
                
                k_conf = prot.res[i_res].n_conf;
                ins_conf(&prot.res[i_res], k_conf, prot.res[i_res].conf[j_conf].n_atom);
                cpy_conf(&prot.res[i_res].conf[k_conf], &prot.res[i_res].conf[j_conf]);
                prot.res[i_res].conf[k_conf].history[2] = 'X';
            }
        }
        
        /* check once more */
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            if (prot.res[i_res].conf[i_conf].history[2] == 'X') break;
        }
        if (i_conf>=prot.res[i_res].n_conf && prot.res[i_res].n_conf!=1) {
            printf("   No conformer made for relaxation for residue %s. It is likely a bug in the program or pdb file.\n",prot.res[i_res].resName);
            return -1;
        }
    }
    //printf("   03-Current time is %ds\n",(int)(time(NULL)-timer_start));

    /* setup commonly used parameter (radii, charge connectivity etc.) */
    id_conf(prot);
    assign_rad(prot);
    assign_crg(prot);
    assign_vdw_param(prot);
    get_connect12(prot);
    setup_vdw_fast(prot);
    for (i_res = 0; i_res< prot.n_res; i_res++) {
        setup_connect_res(prot, i_res);
    }
    //printf("   04-Current time is %ds\n",(int)(time(NULL)-timer_start));

    /*
    FILE *checking;
    checking = fopen("checking.pdb","w");
    write_pdb(checking,prot);
    fclose(checking);
    */
    
    /* make all_atom array */
    na = 0;
    all_atoms = NULL;
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            int i_atom;
            prot.res[i_res].conf[i_conf].i_conf_res  = i_conf;
            prot.res[i_res].conf[i_conf].i_res_prot  = i_res;
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                ATOM *atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                na++;
                atom_p->i_atom_prot = na-1;
                atom_p->i_res_prot  = i_res;
                atom_p->i_conf_res  = i_conf;
                atom_p->i_atom_conf = i_atom;
                
                all_atoms = realloc(all_atoms, na*sizeof(RELAX));
                memset(&all_atoms[na-1],0,sizeof(RELAX));
                all_atoms[na-1].r_orig = atom_p->xyz;
                all_atoms[na-1].atom_p = atom_p;
                all_atoms[na-1].i_relax = -1;
                
                /* which atoms are allowed to move */
                if (i_conf == 0) { /* not moving backbone */
                    all_atoms[na-1].movable = 0;
                }
                else if (!strncmp(prot.res[i_res].conf[i_conf].history+2, "O000", 4)) {
                    if (atom_p->name[1] != 'H') /* not moving heavy atom if it's crystal rotamer */
                        all_atoms[na-1].movable = 0;
                    else
                        all_atoms[na-1].movable = prot.res[i_res].relax;
                }
                else {
                    all_atoms[na-1].movable = prot.res[i_res].relax;
                }
            }
        }
    }
    //printf("   05-Current time is %ds\n",(int)(time(NULL)-timer_start));

    if (env.test_seed < 0) idum = time(NULL); //allows random numbers to be fixed for testing
    else idum = env.test_seed;
    
    for (i_conf=0;i_conf<2000;i_conf++) ran2(&idum);
    //printf("   06-Current time is %ds\n",(int)(time(NULL)-timer_start));

    for (i_cycle = 0; i_cycle < env.n_initial_relax; i_cycle++) {
        FILE *progress_fp;
        float Total_e, Total_e_back;
        int j_res, j_conf;
        time_t timer_cycle_start, timer_cycle_end;
        timer_cycle_start = time(NULL);

        /* choose a conformer for each residue to be relaxed */
        n_relax = 0;
        relax_atoms = NULL;
        for (i_res=0; i_res<prot.n_res; i_res++) {
            add_conf2relax(&prot.res[i_res].conf[0], 0);
            if (prot.res[i_res].n_conf == 1) continue;
            while(1) {
                for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                    if (prot.res[i_res].conf[i_conf].history[2] != 'X') continue;
                    //printf("%d\n",i_conf);
                    
                    if (ran2(&idum) * (float) (prot.res[i_res].n_conf-1) < 0.1) {
                        add_conf2relax(&prot.res[i_res].conf[i_conf], 0);
                        prot.res[i_res].k_conf_on = i_conf;
                        break;
                    }
                }
                if (i_conf<prot.res[i_res].n_conf) break;
            }
        }
        //printf("   07-Current time is %ds\n",(int)(time(NULL)-timer_start));

        /* Setup constraint, 14, ngh list */
        for (i_relax=0; i_relax<n_relax; i_relax++) {
            int n_connect12, n_connect13, n_connect14;
            ATOM **connect12, **connect13, **connect14, *atom_p;
            int i_atom,i_connect, i_const, j_relax;
            
            atom_p = relax_atoms[i_relax].atom_p;
            i_res = atom_p->i_res_prot;
            i_conf = atom_p->i_conf_res;
            i_atom = atom_p->i_atom_conf;
            
            n_connect12 = prot.res[i_res].n_connect12[i_conf][i_atom];
            connect12 = prot.res[i_res].connect12[i_conf][i_atom];
            n_connect13 = prot.res[i_res].n_connect13[i_conf][i_atom];
            connect13 = prot.res[i_res].connect13[i_conf][i_atom];
            
            relax_atoms[i_relax].n_constr = n_connect12 + n_connect13;
            relax_atoms[i_relax].constr_list = malloc(relax_atoms[i_relax].n_constr * sizeof(RELAX *));
            for (i_connect=0; i_connect<n_connect12; i_connect++) {
                int j_relax = all_atoms[connect12[i_connect]->i_atom_prot].i_relax;
                if (j_relax == -1) {
                    relax_atoms[i_relax].constr_list[i_connect] = NULL;
                }
                else {
                    relax_atoms[i_relax].constr_list[i_connect] = &relax_atoms[j_relax];
                }
            }
            for (i_connect=0; i_connect<n_connect13; i_connect++) {
                int j_relax = all_atoms[connect13[i_connect]->i_atom_prot].i_relax;
                if (j_relax == -1) {
                    relax_atoms[i_relax].constr_list[i_connect+n_connect12] = NULL;
                }
                else {
                    relax_atoms[i_relax].constr_list[i_connect+n_connect12] = &relax_atoms[j_relax];
                }
            }
            for (i_const=relax_atoms[i_relax].n_constr-1; i_const>=0; i_const--) {
                if (relax_atoms[i_relax].constr_list[i_const] == NULL) {
                    memmove(relax_atoms[i_relax].constr_list+i_const,
                        relax_atoms[i_relax].constr_list+i_const+1,
                        (relax_atoms[i_relax].n_constr-i_const-1)*sizeof(RELAX *));
                    relax_atoms[i_relax].n_constr--;
                }
            }
            
            relax_atoms[i_relax].constr_dsq  = malloc(relax_atoms[i_relax].n_constr * sizeof(double));
            for (i_const=0; i_const<relax_atoms[i_relax].n_constr; i_const++) {
                relax_atoms[i_relax].constr_dsq[i_const] = ddvv(atom_p->xyz, relax_atoms[i_relax].constr_list[i_const]->atom_p->xyz);
            }
            
            /* 14 list */
            n_connect14 = prot.res[i_res].n_connect14[i_conf][i_atom];
            connect14 = prot.res[i_res].connect14[i_conf][i_atom];
            for (i_connect=0; i_connect < n_connect14; i_connect++) {
                j_relax = all_atoms[connect14[i_connect]->i_atom_prot].i_relax;
                if (j_relax != -1) {
                    int ja;
                    relax_atoms[i_relax].n_ngh14++;
                    relax_atoms[i_relax].ngh14 = realloc(relax_atoms[i_relax].ngh14, relax_atoms[i_relax].n_ngh14*sizeof(RELAX *));
                    ja = prot.res[i_res].connect14[i_conf][i_atom][i_connect]->i_atom_prot;
                    relax_atoms[i_relax].ngh14[relax_atoms[i_relax].n_ngh14-1] = &relax_atoms[all_atoms[ja].i_relax];
                }
            }
            
            /* ngh list */
            for (j_relax =0; j_relax<n_relax; j_relax++) {
                int i_constr, i_ngh;
                float cutoff_far2  = 36.;
                if (i_relax == j_relax) continue;
                if (ddvv(relax_atoms[i_relax].r, relax_atoms[j_relax].r) > cutoff_far2) continue;
                
                /* exclude members in constraint list */
                for (i_constr=0; i_constr<relax_atoms[i_relax].n_constr; i_constr++) {
                    if (relax_atoms[i_relax].constr_list[i_constr] == &relax_atoms[j_relax]) break;
                }
                if (i_constr<relax_atoms[i_relax].n_constr) continue;
                
                /* exclude members in ngh14 list */
                for (i_ngh=0; i_ngh<relax_atoms[i_relax].n_ngh14; i_ngh++) {
                    if (relax_atoms[i_relax].ngh14[i_ngh] == &relax_atoms[j_relax]) break;
                }
                if (i_ngh<relax_atoms[i_relax].n_ngh14) continue;
                
                relax_atoms[i_relax].n_ngh++;
                relax_atoms[i_relax].ngh = realloc(relax_atoms[i_relax].ngh, relax_atoms[i_relax].n_ngh*sizeof(RELAX *));
                relax_atoms[i_relax].ngh[relax_atoms[i_relax].n_ngh-1] = &relax_atoms[j_relax];
            }
        }
        //printf("   08-Current time is %ds\n",(int)(time(NULL)-timer_start));

        /* Setup torsion */
        for (i_relax=0; i_relax<n_relax; i_relax++) {
            ATOM *atom_p, *atom0_p, *atom1_p, *atom2_p, *atom3_p, **connect12;
            int i_atom, n_connect12,i_connect;
            TORS tors;
            
            atom_p = relax_atoms[i_relax].atom_p;
            i_res = atom_p->i_res_prot;
            i_conf = atom_p->i_conf_res;
            i_atom = atom_p->i_atom_conf;
            
            /* get torsion atoms */
            if (torsion_atoms(&prot.res[i_res].conf[i_conf], i_atom, &atom0_p, &atom1_p, &atom2_p, &atom3_p, &tors, 1) == -1) {
                relax_atoms[i_relax].tor_atom1 = -1;
                continue;
            }
            relax_atoms[i_relax].tor_atom1 = all_atoms[atom1_p->i_atom_prot].i_relax;
            relax_atoms[i_relax].tor_atom2 = all_atoms[atom2_p->i_atom_prot].i_relax;
            relax_atoms[i_relax].tor_atom3 = all_atoms[atom3_p->i_atom_prot].i_relax;
            relax_atoms[i_relax].tors      = tors;
            
            /* collect rotatable atoms on the left side of the bond (atom1-atom2) */
            relax_atoms[i_relax].n_rotate1 = 0;
            relax_atoms[i_relax].rotate1_lst = NULL;
            
            /* use connectivity of atom1 */
            n_connect12 = prot.res[atom1_p->i_res_prot].n_connect12[atom1_p->i_conf_res][atom1_p->i_atom_conf];
            connect12 = prot.res[atom1_p->i_res_prot].connect12[atom1_p->i_conf_res][atom1_p->i_atom_conf];
            for (i_connect=0; i_connect<n_connect12; i_connect++) {
                /* exclude atom2 */
                if (connect12[i_connect] == atom2_p) continue;
                if (connect12[i_connect]->i_res_prot == atom2_p->i_res_prot) {
                    if (!strcmp(connect12[i_connect]->name, atom2_p->name)) {
                        continue;
                    }
                }
                
                /* exclude atom not in relax list */
                if (all_atoms[connect12[i_connect]->i_atom_prot].i_relax == -1) continue;
                
                relax_atoms[i_relax].n_rotate1++;
                relax_atoms[i_relax].rotate1_lst = realloc(relax_atoms[i_relax].rotate1_lst, relax_atoms[i_relax].n_rotate1*sizeof(int));
                relax_atoms[i_relax].rotate1_lst[relax_atoms[i_relax].n_rotate1-1] = all_atoms[connect12[i_connect]->i_atom_prot].i_relax;
            }
            
            /* collect rotated atoms on the right side of the bond (atom1-atom2) */
            relax_atoms[i_relax].n_rotate2 = 0;
            relax_atoms[i_relax].rotate2_lst = NULL;
            
            /* use connectivity of atom2 */
            n_connect12 = prot.res[atom2_p->i_res_prot].n_connect12[atom2_p->i_conf_res][atom2_p->i_atom_conf];
            connect12 = prot.res[atom2_p->i_res_prot].connect12[atom2_p->i_conf_res][atom2_p->i_atom_conf];
            for (i_connect=0; i_connect<n_connect12; i_connect++) {
                /* exclude atom1 */
                if (connect12[i_connect] == atom1_p) continue;
                if (connect12[i_connect]->i_res_prot == atom1_p->i_res_prot) {
                    if (!strcmp(connect12[i_connect]->name, atom1_p->name)) {
                        continue;
                    }
                }
                
                /* exclude atom not in relax list */
                if (all_atoms[connect12[i_connect]->i_atom_prot].i_relax == -1) continue;
                
                relax_atoms[i_relax].n_rotate2++;
                relax_atoms[i_relax].rotate2_lst = realloc(relax_atoms[i_relax].rotate2_lst, relax_atoms[i_relax].n_rotate2*sizeof(int));
                relax_atoms[i_relax].rotate2_lst[relax_atoms[i_relax].n_rotate2-1] = all_atoms[connect12[i_connect]->i_atom_prot].i_relax;
            }
        }
        //printf("   09-Current time is %ds\n",(int)(time(NULL)-timer_start));
        //printf("   Setup for relaxation done.\n"); fflush(stdout);
        //printf("   Start relaxation.\n"); fflush(stdout);
        
        //write_relax_pdb("initial.pdb", prot);
        
        /* start relaxation */
        DTSQ = env.hv_relax_dt*env.hv_relax_dt*0.4187e-6;

        Total_e = 0.0;
        for (i_res=0; i_res<prot.n_res; i_res++) {
            if (prot.res[i_res].n_conf == 1) continue;
            i_conf = prot.res[i_res].k_conf_on;
            
            Total_e += torsion_conf(&prot.res[i_res].conf[i_conf]);
            Total_e += vdw_conf_fast(i_res, i_conf, i_res, i_conf, prot, 0);

            for (j_res=0; j_res<prot.n_res; j_res++) {
                float vdw_temp;
                vdw_temp = vdw_conf_fast(i_res, i_conf, j_res, 0, prot, 0);
                Total_e += vdw_temp;
                if (vdw_temp >=100.) {
                    progress_fp = fopen(env.progress_log, "a");
                    if (progress_fp) {
                        fprintf(progress_fp, "   Initial Relaxation Cycle = %3d, vdw between %s and %s = %11.3f\n",
                            i_cycle+1,prot.res[i_res].conf[i_conf].uniqID,prot.res[j_res].conf[0].uniqID,vdw_temp);
                        fclose(progress_fp);
                    }
                }

                if (i_res >= j_res) continue;
                if (prot.res[j_res].n_conf == 1) continue;
                j_conf = prot.res[j_res].k_conf_on;

                vdw_temp = vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 0);
                Total_e += vdw_temp;
                if (vdw_temp >=100.) {
                    progress_fp = fopen(env.progress_log, "a");
                    if (progress_fp) {
                        fprintf(progress_fp, "   Initial Relaxation Cycle = %3d, vdw between %s and %s = %11.3f\n",
                            i_cycle+1,prot.res[i_res].conf[i_conf].uniqID,prot.res[j_res].conf[j_conf].uniqID,vdw_temp);
                        fclose(progress_fp);
                    }
                }
            }
        }
        
        for (i_res=0; i_res<prot.n_res; i_res++) {
            if (prot.res[i_res].n_conf == 1) continue;
            i_conf = prot.res[i_res].k_conf_on;
            
            for (j_res=0; j_res<prot.n_res; j_res++) {
                
                Total_e += coulomb_conf_fast(i_res, i_conf, j_res, 0, prot);

                if (i_res >= j_res) continue;
                if (prot.res[j_res].n_conf == 1) continue;
                j_conf = prot.res[j_res].k_conf_on;

                Total_e += coulomb_conf_fast(i_res, i_conf, j_res, j_conf, prot);
            }
        }
        //printf("   Initial relaxation: Total_e = %11.3f, DTSQ = %.10f\n", Total_e, DTSQ);
        
        progress_fp = fopen(env.progress_log, "a");
        if (progress_fp) {
            fprintf(progress_fp, "   Initial Relaxation Cycle = %3d, initial energy = %11.3f\n",
                i_cycle+1,Total_e);
            fclose(progress_fp);
        }
        
        done = 0;
        i_relax_step = -1;
        while (!done) {
            n_pair_trial++;
            while (i_relax_step) {
                float tors_scale = 1.;
                
                get_frc(tors_scale, prot);
                get_rp();
                if (shake()) break;
                
                for (i_relax =0; i_relax<n_relax; i_relax++) {
                    relax_atoms[i_relax].r = relax_atoms[i_relax].r_p;
                }

                /*
                if (fmod(-i_relax_step, 500) == 0) {
                    char relax_out[160];
                    sprintf(relax_out,"relax%05d_%05d.pdb",i_cycle,(-i_relax_step)/500);
                    write_relax_pdb(relax_out, prot);
                }
                */
                
                /* check convergence */
                if (fmod(-i_relax_step, 500) == 0) {
                    //float hv_relax_dt = sqrt(DTSQ);
                    int converged;
                    //printf("   10a-Current time is %ds\n",(int)(time(NULL)-timer_start));
                    
                    Total_e_back = Total_e;
                    Total_e = 0.0;
                    /* temporarily return xyz to atom */
                    for (i_relax =0; i_relax<n_relax; i_relax++) {
                        relax_atoms[i_relax].atom_p->r_orig = relax_atoms[i_relax].atom_p->xyz;
                        relax_atoms[i_relax].atom_p->xyz = relax_atoms[i_relax].r;
                    }
                    
                    for (i_res=0; i_res<prot.n_res; i_res++) {
                        if (prot.res[i_res].n_conf == 1) continue;
                        i_conf = prot.res[i_res].k_conf_on;
                        
                        Total_e += torsion_conf(&prot.res[i_res].conf[i_conf]);
                        Total_e += vdw_conf_fast(i_res, i_conf, i_res, i_conf, prot, 0);

                        for (j_res=0; j_res<prot.n_res; j_res++) {
                            
                            Total_e += vdw_conf_fast(i_res, i_conf, j_res, 0, prot, 0);
                            
                            if (i_res >= j_res) continue;
                            if (prot.res[j_res].n_conf == 1) continue;
                            j_conf = prot.res[j_res].k_conf_on;
                            
                            Total_e += vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 0);
                        }
                    }
                    //printf("   Initial relaxation: Total_e = %11.3f, DTSQ = %.10f, i_step = %d\n", Total_e, DTSQ, i_relax_step);
                    for (i_res=0; i_res<prot.n_res; i_res++) {
                        if (prot.res[i_res].n_conf == 1) continue;
                        i_conf = prot.res[i_res].k_conf_on;
                        
                        for (j_res=0; j_res<prot.n_res; j_res++) {
                            
                            Total_e += coulomb_conf_fast(i_res, i_conf, j_res, 0, prot);
                            
                            if (i_res >= j_res) continue;
                            if (prot.res[j_res].n_conf == 1) continue;
                            j_conf = prot.res[j_res].k_conf_on;
                            
                            Total_e += coulomb_conf_fast(i_res, i_conf, j_res, j_conf, prot);
                        }
                    }
                    
                    for (i_relax =0; i_relax<n_relax; i_relax++) {
                        relax_atoms[i_relax].atom_p->xyz = relax_atoms[i_relax].atom_p->r_orig;
                    }
                    //printf("   Initial relaxation: Total_e = %11.3f, DTSQ = %.10f, i_step = %d\n", Total_e, DTSQ, i_relax_step);
                    //printf("   10b-Current time is %ds\n",(int)(time(NULL)-timer_start));

                    /*
                    for (i_relax =0; i_relax<n_relax; i_relax++) {
                        if (converged) {
                            if (dvv(relax_atoms[i_relax].r, relax_atoms[i_relax].atom_p->xyz)
                                > CONVERGED_R*hv_relax_dt*(-i_relax_step)) {
                            converged = 0;
                                }
                        }
                    }
                    */
                    
                    /* energy drop by more than 10 Kcal/mol, then not converged */
                    if (Total_e - Total_e_back < -10.) converged = 0;
                    else converged = 1;
                    
                    if (converged) {
                        //printf("%d\n",-i_relax_step);
                        i_relax_step = 0;
                    }
                    else {
                        DTSQ = DTSQ * 2.;
                        i_relax_step = i_relax_step / 2.;
                        i_relax_step--;
                    }
                    
                    //printf("   10c-Current time is %ds\n",(int)(time(NULL)-timer_start));
                    
                }
                else {
                    i_relax_step--;
                }
            }
            if (i_relax_step) {
                if (DTSQ < 1e-10) {
                    /* relaxation failed */
                    for (i_relax =0; i_relax<n_relax; i_relax++) {
                        relax_atoms[i_relax].r = relax_atoms[i_relax].atom_p->xyz;
                    }
                    
                    break;
                }
                
                /* reduce step size and continue */
                DTSQ = DTSQ /2.;
                i_relax_step = i_relax_step * 2. + 1;
                //printf("   Initial relaxation: DTSQ = %.10f, i_step = %d\n", DTSQ, i_relax_step);
            }
            else {
                done = 1;
            }
        }
        //printf("   10-Current time is %ds\n",(int)(time(NULL)-timer_start));
    
        /* return final xyz to atom */
        for (i_relax =0; i_relax<n_relax; i_relax++) {
            relax_atoms[i_relax].atom_p->xyz = relax_atoms[i_relax].r;
        }
        //printf("   11-Current time is %ds\n",(int)(time(NULL)-timer_start));

        /* free memory */
        for (i_relax =0; i_relax<n_relax; i_relax++) {
            ia = relax_atoms[i_relax].atom_p->i_atom_prot;
            all_atoms[ia].i_relax = -1;
        }
        
        for (i_relax =0; i_relax<n_relax; i_relax++) {
            free(relax_atoms[i_relax].ngh);
            free(relax_atoms[i_relax].ngh14);
            free(relax_atoms[i_relax].constr_list);
            free(relax_atoms[i_relax].constr_dsq);
        }
        free(relax_atoms);
        
        timer_cycle_end = time(NULL);
        progress_fp = fopen(env.progress_log, "a");
        if (progress_fp) {
            fprintf(progress_fp, "   Initial Relaxation Cycle = %3d, final energy = %11.3f, time used %ld seconds\n",
                i_cycle+1,Total_e, timer_cycle_end-timer_cycle_start);
            fclose(progress_fp);
        }
        //printf("   12-Current time is %ds\n",(int)(time(NULL)-timer_start));
    }
    //printf("   13-Current time is %ds\n",(int)(time(NULL)-timer_start));

    fp = fopen("initial_relax.pdb", "w");
    write_pdb(fp, prot);
    fclose(fp);
    
    delete_h(prot);
    float orig_prune_thr = env.prune_thr;
    env.prune_thr = 0.005;
    rm_dupconf_hv(prot);
    env.prune_thr = orig_prune_thr;
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            if (prot.res[i_res].conf[i_conf].history[2] != 'X') continue;
            prot.res[i_res].conf[i_conf].history[2] = 'Y'; /* relabel 'X' to 'Y' and avoid relaxation later */
        }
    }
    
    timer_end = time(NULL);
    printf("   Total time used for initial relaxation is %ds\n",(int)(timer_end-timer_start));
    return 0;
}

void relaxation_setup(PROT prot)
{
    int i_res, i_conf, i_atom;
    int ia, ja, i_connect, i_const, n_connect12, n_connect13;
    ATOM *atom_p, *atom0_p, *atom1_p, *atom2_p, *atom3_p, **connect12, **connect13;
    TORS tors;
    
    /* define relax parameter for each res */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        if (param_get("RELAX", prot.res[i_res].resName, "", &prot.res[i_res].relax)) {
            prot.res[i_res].relax = 0;
        }
        /* if it's original rotamer, make a copy 'X' and do relaxation on 'X' */
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            if (!strncmp(prot.res[i_res].conf[i_conf].history+2, "O000", 4)) {
                if (prot.res[i_res].relax) {
                    ins_conf(&prot.res[i_res], prot.res[i_res].n_conf, prot.res[i_res].conf[i_conf].n_atom);
                    cpy_conf(&prot.res[i_res].conf[prot.res[i_res].n_conf-1], &prot.res[i_res].conf[i_conf]);
                    prot.res[i_res].conf[prot.res[i_res].n_conf-1].history[2] = 'X';
                }
            }
        }
    }
    
    /* setup commonly used parameter (radii, charge connectivity etc.) */
    assign_rad(prot);
    assign_crg(prot);
    get_connect12(prot);
    //setup_C6_C12(prot);
    setup_vdw_fast(prot);
    for (i_res = 0; i_res< prot.n_res; i_res++) {
        setup_connect_res(prot, i_res);
    }
    
    /* make all_atoms array */
    na = 0;
    all_atoms = NULL;
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            prot.res[i_res].conf[i_conf].i_conf_res  = i_conf;
            prot.res[i_res].conf[i_conf].i_res_prot  = i_res;
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                na++;
                atom_p->i_atom_prot = na-1;
                atom_p->i_res_prot  = i_res;
                atom_p->i_conf_res  = i_conf;
                atom_p->i_atom_conf = i_atom;
                
                all_atoms = realloc(all_atoms, na*sizeof(RELAX));
                memset(&all_atoms[na-1],0,sizeof(RELAX));
                all_atoms[na-1].r_orig = atom_p->xyz;
                all_atoms[na-1].atom_p = atom_p;
                all_atoms[na-1].i_relax = -1;
                
                if (i_conf == 0) { /* not moving backbone */
                    all_atoms[na-1].movable = 0;
                }
                else if (!strncmp(prot.res[i_res].conf[i_conf].history+2, "O000", 4)
                    || prot.res[i_res].conf[i_conf].history[2] == 'Y') {
                    if (atom_p->name[1] != 'H') /* not moving heavy atom if it's crystal rotamer */
                        all_atoms[na-1].movable = 0;
                    else
                        all_atoms[na-1].movable = prot.res[i_res].relax;
                }
                else {
                    all_atoms[na-1].movable = prot.res[i_res].relax;
                }
            }
        }
    }
    
    /* Constraint list */
    for (ia=0;ia<na;ia++) {
        atom_p = all_atoms[ia].atom_p;
        i_res = atom_p->i_res_prot;
        i_conf = atom_p->i_conf_res;
        i_atom = atom_p->i_atom_conf;
        n_connect12 = prot.res[i_res].n_connect12[i_conf][i_atom];
        connect12 = prot.res[i_res].connect12[i_conf][i_atom];
        n_connect13 = prot.res[i_res].n_connect13[i_conf][i_atom];
        connect13 = prot.res[i_res].connect13[i_conf][i_atom];
        
        all_atoms[ia].n_constr = n_connect12 + n_connect13;
        all_atoms[ia].constr_list = malloc(all_atoms[ia].n_constr * sizeof(RELAX *));
        for (i_connect=0; i_connect<n_connect12; i_connect++) {
            ja = connect12[i_connect]->i_atom_prot;
            all_atoms[ia].constr_list[i_connect] = &all_atoms[ja];
        }
        for (i_connect=0; i_connect<n_connect13; i_connect++) {
            ja = connect13[i_connect]->i_atom_prot;
            all_atoms[ia].constr_list[i_connect+n_connect12] = &all_atoms[ja];
        }
        
        all_atoms[ia].constr_dsq  = malloc(all_atoms[ia].n_constr * sizeof(double));
        for (i_const=0; i_const<all_atoms[ia].n_constr; i_const++) {
            all_atoms[ia].constr_dsq[i_const] = ddvv(atom_p->xyz, all_atoms[ia].constr_list[i_const]->atom_p->xyz);
        }
    }
    
    /* torsion parameter */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                RELAX *relax_atom_p;
                if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                ia = prot.res[i_res].conf[i_conf].atom[i_atom].i_atom_prot;
                all_atoms[ia].tor_atom1 = -1;
                
                /* get torsion atoms */
                if (torsion_atoms(&prot.res[i_res].conf[i_conf], i_atom, &atom0_p, &atom1_p, &atom2_p, &atom3_p, &tors, 1) == -1) continue;
                all_atoms[ia].tor_atom1 = atom1_p->i_atom_prot;
                all_atoms[ia].tor_atom2 = atom2_p->i_atom_prot;
                all_atoms[ia].tor_atom3 = atom3_p->i_atom_prot;
                all_atoms[ia].tors      = tors;
                
                /* collect rotated atoms on the left side of the bond (atom1-atom2) */
                all_atoms[ia].n_rotate1 = 0;
                all_atoms[ia].rotate1_lst = NULL;
                relax_atom_p = &all_atoms[atom1_p->i_atom_prot];
                /* use connectivity of atom1 */
                n_connect12 = prot.res[atom1_p->i_res_prot].n_connect12[atom1_p->i_conf_res][atom1_p->i_atom_conf];
                connect12 = prot.res[atom1_p->i_res_prot].connect12[atom1_p->i_conf_res][atom1_p->i_atom_conf];
                for (i_connect=0; i_connect<n_connect12; i_connect++) {
                    /* exclude atom2 */
                    if (connect12[i_connect] == atom2_p) continue;
                    if (connect12[i_connect]->i_res_prot == atom2_p->i_res_prot) {
                        if (!strcmp(connect12[i_connect]->name, atom2_p->name)) {
                            continue;
                        }
                    }
                    /* exclude atom from other sidechain */
                    if (connect12[i_connect]->i_res_prot == i_res) /* same residue */
                        if (connect12[i_connect]->i_conf_res != i_conf) /* different conformer */
                            if (connect12[i_connect]->i_conf_res) /* it's another side chain */
                                continue;
                    
                    //printf("%s-%s-%s\n",atom2_p->name,atom1_p->name,relax_atom_p->connect12[i_connect]->name);
                    all_atoms[ia].n_rotate1++;
                    all_atoms[ia].rotate1_lst = realloc(all_atoms[ia].rotate1_lst, all_atoms[ia].n_rotate1*sizeof(int));
                    all_atoms[ia].rotate1_lst[all_atoms[ia].n_rotate1-1] = connect12[i_connect]->i_atom_prot;
                }
                
                /* collect rotated atoms on the right side of the bond (atom1-atom2) */
                all_atoms[ia].n_rotate2 = 0;
                all_atoms[ia].rotate2_lst = NULL;
                relax_atom_p = &all_atoms[atom2_p->i_atom_prot];
                /* use connectivity of atom2 */
                n_connect12 = prot.res[atom2_p->i_res_prot].n_connect12[atom2_p->i_conf_res][atom2_p->i_atom_conf];
                connect12 = prot.res[atom2_p->i_res_prot].connect12[atom2_p->i_conf_res][atom2_p->i_atom_conf];
                for (i_connect=0; i_connect<n_connect12; i_connect++) {
                    /* exclude atom1 */
                    if (connect12[i_connect] == atom1_p) continue;
                    if (connect12[i_connect]->i_res_prot == atom1_p->i_res_prot) {
                        if (!strcmp(connect12[i_connect]->name, atom1_p->name)) {
                            continue;
                        }
                    }
                    /* exclude atom from other sidechain */
                    if (connect12[i_connect]->i_res_prot == i_res) /* same residue */
                        if (connect12[i_connect]->i_conf_res != i_conf) /* different conformer */
                            if (connect12[i_connect]->i_conf_res) /* it's another side chain */
                                continue;
                    
                    //printf("%s-%s-%s\n",atom1_p->name,atom2_p->name,relax_atom_p->connect12[i_connect]->name);
                    all_atoms[ia].n_rotate2++;
                    all_atoms[ia].rotate2_lst = realloc(all_atoms[ia].rotate2_lst, all_atoms[ia].n_rotate2*sizeof(int));
                    all_atoms[ia].rotate2_lst[all_atoms[ia].n_rotate2-1] = connect12[i_connect]->i_atom_prot;
                }
                /*
                printf("%d,%s,%s,%s\n",ia,all_atoms[ia].atom_p->name,atom1_p->name,atom2_p->name);
                for (i_rotate = 0; i_rotate< all_atoms[ia].n_rotate1; i_rotate++)
                    printf(" %s ",all_atoms[all_atoms[ia].rotate1_lst[i_rotate]].atom_p->name);
                printf("\n");
                for (i_rotate = 0; i_rotate< all_atoms[ia].n_rotate2; i_rotate++)
                    printf(" %s ",all_atoms[all_atoms[ia].rotate2_lst[i_rotate]].atom_p->name);
                printf("\n");
                */
            }
        }
    }
    
    /* setup residue ngh list */
    if (env.hv_relax_include_ngh) {
        float ngh_thr2 = env.hv_relax_ngh_thr * env.hv_relax_ngh_thr;
        float pair_vdw_hv;
        int j_res, j_conf, j_atom, found;
        for (i_res = 0; i_res< prot.n_res; i_res++) {
            prot.res[i_res].n_ngh = 0;
            prot.res[i_res].ngh = NULL;
            if (!prot.res[i_res].cal_vdw) continue;
            
            for (j_res = 0; j_res< prot.n_res; j_res++) {
                if (i_res == j_res) continue;
                if (!prot.res[j_res].cal_vdw) continue;
                
                /* check if distance within the threshold */
                if (out_of_range(prot.res[i_res].r_min,prot.res[i_res].r_max,prot.res[j_res].r_min,prot.res[j_res].r_max,ngh_thr2)) continue;
                found = 0;
                for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
                    for (j_conf=0; j_conf<prot.res[j_res].n_conf; j_conf++) {
                        
                        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                            if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                            for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
                                if (!prot.res[j_res].conf[j_conf].atom[j_atom].on) continue;
                                
                                if (ddvv(prot.res[i_res].conf[i_conf].atom[i_atom].xyz,prot.res[j_res].conf[j_conf].atom[j_atom].xyz) < ngh_thr2) {
                                    found = 1;
                                    break;
                                }
                            }
                            if (found) break;
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (!found) continue;
                
                /* exclude jres if two residues always clash with each other (for example two water close to each other */
                if (!prot.res[j_res].conf[0].n_atom) {  /* only check res with no backbone */
                    found = 0;
                    for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                        if (!prot.res[i_res].conf[i_conf].n_atom) continue;
                        for (j_conf=1; j_conf<prot.res[j_res].n_conf; j_conf++) {
                            if (!prot.res[j_res].conf[j_conf].n_atom) continue;
                            
                            pair_vdw_hv = vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 1);
                            if (pair_vdw_hv > env.hv_relax_hv_vdw_thr) continue;
                            
                            found = 1;
                            break;
                        }
                        if (found) break;
                    }
                    if (!found) continue;
                }
                
                prot.res[i_res].n_ngh++;
                prot.res[i_res].ngh = realloc(prot.res[i_res].ngh, prot.res[i_res].n_ngh * sizeof(RES *));
                prot.res[i_res].ngh[prot.res[i_res].n_ngh-1] = &prot.res[j_res];
            }
        }
    }
}

void get_frc(float tors_scale, PROT prot) {
    int i_relax, i_ngh;
    int i_rotate;
    float C6,C12, dsq, phi;
    VECTOR atom0_r, atom1_r, atom2_r, atom3_r, k, torq, r, r_p;
    
    for (i_relax=0;i_relax<n_relax;i_relax++) {
        relax_atoms[i_relax].lj_frc = zero;
        relax_atoms[i_relax].elec_frc = zero;
        relax_atoms[i_relax].torsion_frc = zero;
    }
    
    for (i_relax=0;i_relax<n_relax;i_relax++) {
        if (relax_atoms[i_relax].movable) {
            
            /* LJ interactions */
            //printf("%s\n",relax_atoms[i_relax].atom_p->name);
            for (i_ngh=0; i_ngh<relax_atoms[i_relax].n_ngh; i_ngh++) {
                float sig_min = relax_atoms[i_relax].atom_p->vdw_rad + relax_atoms[i_relax].ngh[i_ngh]->atom_p->vdw_rad;
                float eps = sqrt(relax_atoms[i_relax].atom_p->vdw_eps*relax_atoms[i_relax].ngh[i_ngh]->atom_p->vdw_eps);
                C12 = eps*pow(sig_min,12);
                C6 = 2.*eps*pow(sig_min,6);
                //printf("%11.6f, %11.6f, %11.6f, %11.6f\n",eps,sig_min,C12,C6);
                
                relax_atoms[i_relax].lj_frc = vector_vplusv(relax_atoms[i_relax].lj_frc,
                    vdw_frc(relax_atoms[i_relax].r, relax_atoms[i_relax].ngh[i_ngh]->r, C6, C12));
                //VECTOR frc = vdw_frc(relax_atoms[i_relax].r, relax_atoms[i_relax].ngh[i_ngh]->r, C6, C12);
                //printf("   %s %11.6f%11.6f%11.6f\n",relax_atoms[i_relax].ngh[i_ngh]->atom_p->name, frc.x, frc.y, frc.z);
                /* scaling using the parameter for later monte carlo */
                //relax_atoms[i_relax].lj_frc = vector_rescale(relax_atoms[i_relax].lj_frc, env.scale_vdw);
                /*
                if (!strcmp(relax_atoms[i_relax].atom_p->name, "2HG ")) {
                    printf("%2d, %s, %11.6f, %11.6f, %11.6f\n",i_ngh,relax_atoms[i_relax].ngh[i_ngh]->atom_p->name,relax_atoms[i_relax].lj_frc.x,relax_atoms[i_relax].lj_frc.y,relax_atoms[i_relax].lj_frc.z);
                }
                */
            }
            for (i_ngh=0; i_ngh<relax_atoms[i_relax].n_ngh14; i_ngh++) {
                float sig_min = relax_atoms[i_relax].atom_p->vdw_rad + relax_atoms[i_relax].ngh14[i_ngh]->atom_p->vdw_rad;
                float eps = sqrt(relax_atoms[i_relax].atom_p->vdw_eps*relax_atoms[i_relax].ngh14[i_ngh]->atom_p->vdw_eps);
                C12 = eps*pow(sig_min,12);
                C6 = 2.*eps*pow(sig_min,6);
                
                //C6 = C6_matrix[relax_atoms[i_relax].atom_p->i_elem][relax_atoms[i_relax].ngh14[i_ngh]->atom_p->i_elem];
                //C12 = C12_matrix[relax_atoms[i_relax].atom_p->i_elem][relax_atoms[i_relax].ngh14[i_ngh]->atom_p->i_elem];
                relax_atoms[i_relax].lj_frc = vector_vplusv(relax_atoms[i_relax].lj_frc,
                    vector_rescale(vdw_frc(relax_atoms[i_relax].r, relax_atoms[i_relax].ngh14[i_ngh]->r, C6, C12), env.factor_14lj));
                /* scaling using the parameter for later monte carlo */
                //relax_atoms[i_relax].lj_frc = vector_rescale(relax_atoms[i_relax].lj_frc, env.scale_vdw);
                /*
                if (!strcmp(relax_atoms[i_relax].atom_p->name, "2HG ")) {
                    printf("%2d, %s, %11.6f, %11.6f, %11.6f\n",i_ngh,relax_atoms[i_relax].ngh14[i_ngh]->atom_p->name,relax_atoms[i_relax].lj_frc.x,relax_atoms[i_relax].lj_frc.y,relax_atoms[i_relax].lj_frc.z);
                }
                */
            }
            /* scaling using the parameter for later monte carlo */
            relax_atoms[i_relax].lj_frc = vector_rescale(relax_atoms[i_relax].lj_frc, env.scale_vdw);
            
            /* Electrostatic interactions */
            if (relax_atoms[i_relax].atom_p->crg > 1e-6 || relax_atoms[i_relax].atom_p->crg < -1e-6) {
                for (i_ngh=0; i_ngh<relax_atoms[i_relax].n_ngh; i_ngh++) {
                    if (fabs(relax_atoms[i_relax].ngh[i_ngh]->atom_p->crg) < 1e-6) continue;
                    relax_atoms[i_relax].elec_frc = vector_vplusv(relax_atoms[i_relax].elec_frc,
                        coulomb_frc(relax_atoms[i_relax].r,
                            relax_atoms[i_relax].ngh[i_ngh]->r,
                            relax_atoms[i_relax].atom_p->crg,
                            relax_atoms[i_relax].ngh[i_ngh]->atom_p->crg));
                }
                for (i_ngh=0; i_ngh<relax_atoms[i_relax].n_ngh14; i_ngh++) {
                    if (fabs(relax_atoms[i_relax].ngh14[i_ngh]->atom_p->crg) < 1e-6) continue;
                    relax_atoms[i_relax].elec_frc = vector_vplusv(relax_atoms[i_relax].elec_frc,
                        coulomb_frc(relax_atoms[i_relax].r,
                            relax_atoms[i_relax].ngh14[i_ngh]->r,
                            relax_atoms[i_relax].atom_p->crg,
                            relax_atoms[i_relax].ngh14[i_ngh]->atom_p->crg));
                }
            }
            
            /* Extra constraint from its original position */
            dsq = ddvv(relax_atoms[i_relax].r,relax_atoms[i_relax].r_orig);
            if (dsq < CONSTRAINT2 || relax_atoms[i_relax].atom_p->name[1] == 'H') relax_atoms[i_relax].constr_frc = zero;
            else {
                relax_atoms[i_relax].constr_frc = 
                vector_rescale(vector_vminusv(relax_atoms[i_relax].r_orig, relax_atoms[i_relax].r), CONSTRAINT_FRC*(dsq-CONSTRAINT2));
            }
        }
        
        /* Torsion */
        if (relax_atoms[i_relax].tor_atom1 != -1) {
            int i_term;
            atom0_r = relax_atoms[i_relax].r;
            atom1_r = relax_atoms[relax_atoms[i_relax].tor_atom1].r;
            atom2_r = relax_atoms[relax_atoms[i_relax].tor_atom2].r;
            atom3_r = relax_atoms[relax_atoms[i_relax].tor_atom3].r;
            k = vector_normalize(vector_vminusv(atom1_r, atom2_r));
            phi = torsion_angle(atom0_r,atom1_r,atom2_r,atom3_r);
            torq = zero;
            for (i_term=0;i_term<relax_atoms[i_relax].tors.n_term; i_term++) {
                torq = vector_vplusv(torq,
                torsion_torq(phi, relax_atoms[i_relax].tors.V2[i_term],
                relax_atoms[i_relax].tors.n_fold[i_term],
                relax_atoms[i_relax].tors.gamma[i_term], k));
            }
            
            for (i_rotate = 0; i_rotate < relax_atoms[i_relax].n_rotate1; i_rotate++) {
                int k_relax = relax_atoms[i_relax].rotate1_lst[i_rotate];
                if (k_relax == -1) continue;
                r = vector_vminusv(relax_atoms[k_relax].r, atom1_r);
                r_p = vector_vminusv(r, vector_rescale(k,vdotv(r,k)));
                
                relax_atoms[k_relax].torsion_frc = 
                vector_vplusv(relax_atoms[k_relax].torsion_frc,
                vector_rescale(vector_vxv(torq, r_p), tors_scale/vdotv(r_p,r_p)/(float)(relax_atoms[i_relax].n_rotate1) ));
            }
            for (i_rotate = 0; i_rotate < relax_atoms[i_relax].n_rotate2; i_rotate++) {
                int k_relax = relax_atoms[i_relax].rotate2_lst[i_rotate];
                if (k_relax == -1) continue;
                r = vector_vminusv(relax_atoms[k_relax].r, atom1_r);
                r_p = vector_vminusv(r, vector_rescale(k,vdotv(r,k)));
                
                relax_atoms[k_relax].torsion_frc = 
                vector_vplusv(relax_atoms[k_relax].torsion_frc,
                vector_rescale(vector_vxv(torq, r_p), -tors_scale/vdotv(r_p,r_p)/(float)(relax_atoms[i_relax].n_rotate2) ));
            }
        }
        
    }
    
    /* Total force */
    for (i_relax=0;i_relax<n_relax;i_relax++) {
        if (!relax_atoms[i_relax].movable) continue;
        
        relax_atoms[i_relax].frc = vector_vplusv(
        vector_vplusv(
        relax_atoms[i_relax].lj_frc,
        relax_atoms[i_relax].elec_frc),
        vector_vplusv(
        relax_atoms[i_relax].torsion_frc,
        relax_atoms[i_relax].constr_frc));
        /*
        if (!strcmp(relax_atoms[i_relax].atom_p->name, " OG1")) {
            printf("lj %11.6f, %11.6f, %11.6f\n",relax_atoms[i_relax].lj_frc.x,relax_atoms[i_relax].lj_frc.y,relax_atoms[i_relax].lj_frc.z);
            printf("to %11.6f, %11.6f, %11.6f\n",relax_atoms[i_relax].frc.x,relax_atoms[i_relax].frc.y,relax_atoms[i_relax].frc.z);
        }
        */
    }
}

/*
DT is in pico sec;
mass in a.u; (all atoms are asumed to be 1 a.u. in this relaxation function)
force in kcal/(mol*anstrom);
0.4187 makes the conversion and the displacement would be in anstrom
*/
void get_rp() {
    int i_relax;
    
    for (i_relax=0;i_relax<n_relax;i_relax++) {
        if (!relax_atoms[i_relax].movable) continue;
        relax_atoms[i_relax].r_p = vector_vplusv(relax_atoms[i_relax].r, vector_rescale(relax_atoms[i_relax].frc,DTSQ));
        //printf("%s %11.6f%11.6f%11.6f%11.6f%11.6f%11.6f\n",relax_atoms[i_relax].atom_p->name,relax_atoms[i_relax].r_p.x,relax_atoms[i_relax].r_p.y,relax_atoms[i_relax].r_p.z,relax_atoms[i_relax].r.x,relax_atoms[i_relax].r.y,relax_atoms[i_relax].r.z);
    }
}

int shake() {
    int i_relax, i_constr, iter;
    VECTOR  r12, rp12, d;
    double  r12sq;
    double rp12sq, r_rp, diffsq, g;
    int done;
    double TOL2 = env.hv_relax_shake_tol * 2.;
    
    for (i_relax =0; i_relax<n_relax; i_relax++) {
        relax_atoms[i_relax].moving = 0;
        if (!relax_atoms[i_relax].movable) relax_atoms[i_relax].moved = 0;
        else relax_atoms[i_relax].moved = 1;
    }
    
    done = 0;
    iter = env.hv_relax_n_shake;
    while (iter) {
        if (done) break;
        done = 1;
        
        for (i_relax =0; i_relax<n_relax; i_relax++) {
            if (!relax_atoms[i_relax].moved) continue;
            for (i_constr=0; i_constr<relax_atoms[i_relax].n_constr; i_constr++) {
                //if (!relax_atoms[i_relax].constr_list[i_constr]->moved) continue;
                
                /* calculate the distance between 1 and 2 before the force is applied */
                r12 = vector_vminusv(relax_atoms[i_relax].r, relax_atoms[i_relax].constr_list[i_constr]->r);
                r12sq = relax_atoms[i_relax].constr_dsq[i_constr];
                
                /* calculate the distance after */
                rp12 = vector_vminusv(relax_atoms[i_relax].r_p, relax_atoms[i_relax].constr_list[i_constr]->r_p);
                rp12sq = vdotv(rp12,rp12);
                
                /* dot product of two vectors */
                r_rp = vdotv(r12, rp12);
                diffsq = (r12sq-rp12sq);
                
                if (fabs(diffsq) < (r12sq * TOL2)) continue;
                
                //if (relax_atoms[i_relax].movable && relax_atoms[i_relax].constr_list[i_constr]->movable) {
                if (relax_atoms[i_relax].constr_list[i_constr]->movable) {
                    /* move both */
                    g = diffsq / (2. * r_rp);
                    d = vector_rescale(r12, 0.5 * g);
                    
                    relax_atoms[i_relax].r_p = vector_vplusv(relax_atoms[i_relax].r_p, d);
                    relax_atoms[i_relax].moving = 1;
                    relax_atoms[i_relax].constr_list[i_constr]->r_p = vector_vminusv(relax_atoms[i_relax].constr_list[i_constr]->r_p, d);
                    relax_atoms[i_relax].constr_list[i_constr]->moving = 1;
                }
                else {
                    /* move one */
                    g = diffsq / (2. * r_rp);
                    d = vector_rescale(r12, g);
                    
                    //if (relax_atoms[i_relax].movable) {
                        relax_atoms[i_relax].r_p = vector_vplusv(relax_atoms[i_relax].r_p, d);
                        relax_atoms[i_relax].moving = 1;
                    //}
                    /*
                    if (relax_atoms[i_relax].constr_list[i_constr]->movable) {
                        relax_atoms[i_relax].constr_list[i_constr]->r_p = vector_vminusv(relax_atoms[i_relax].constr_list[i_constr]->r_p, d);
                        relax_atoms[i_relax].constr_list[i_constr]->moving = 1;
                    }
                    */
                }
                done = 0;
            }
        }
        
        for (i_relax =0; i_relax<n_relax; i_relax++) {
            relax_atoms[i_relax].moved = relax_atoms[i_relax].moving;
            relax_atoms[i_relax].moving = 0;
        }
        iter--;
    }
    
    if (!done) {
        return -1;
    }
    else return 0;
}

void add_conf2relax(CONF *conf_p, int fix) {
    int i_atom, ia;
    
    if (conf_p->i_conf_res) {
        n_relax_res++;
        relax_res[n_relax_res-1] = conf_p->i_res_prot;
        relax_conf[n_relax_res-1] = conf_p->i_conf_res;
    }
    for (i_atom=0; i_atom < conf_p->n_atom; i_atom++) {
        if (!conf_p->atom[i_atom].on) continue;
        ia = conf_p->atom[i_atom].i_atom_prot;
        add_atom2relax(ia, fix);
    }
}

void add_atom2relax(int ia, int fix) {
    n_relax++;
    all_atoms[ia].i_relax = n_relax-1;
    relax_atoms = realloc(relax_atoms, n_relax * sizeof(RELAX));
    relax_atoms[n_relax-1] = all_atoms[ia];
    relax_atoms[n_relax-1].r = relax_atoms[n_relax-1].atom_p->xyz;
    relax_atoms[n_relax-1].r_p = relax_atoms[n_relax-1].atom_p->xyz;

    if (fix)
        relax_atoms[n_relax-1].movable = 0;
}

void complete_constr(int i_res, int i_conf, int j_res, int j_conf) {
    int i_relax, i_constr;
    
    /* check if all atoms in constraint list are in the relaxation list */
    for (i_relax =0; i_relax<n_relax; i_relax++) {
        int ia;
        if (!relax_atoms[i_relax].movable) continue;
        ia = relax_atoms[i_relax].atom_p->i_atom_prot;
        
        for (i_constr=0; i_constr<all_atoms[ia].n_constr; i_constr++) {
            int ja = all_atoms[ia].constr_list[i_constr]->atom_p->i_atom_prot;
            int j_relax = all_atoms[ia].constr_list[i_constr]->i_relax;
            
            if (j_relax == -1) {
                /* if atom j in the constraint list is not being relaxed */
                
                if (!res_in_relax(all_atoms[ja].atom_p->i_res_prot)) {
                    /* if atom j is from a residues being relaxed (different sidechain),
                    do nothing, contraint doesn't apply.
                    otherwise, if atom j is from another residue:*/
                    
                    if (all_atoms[ja].atom_p->i_conf_res == 0) {
                        /* if it's a backbone atom, add atom j and force it to be fixed */
                        add_atom2relax(ja, 1);
                    }
                    else {
                        /* if it's a sidechain atom, then do not relax atom i */
                        relax_atoms[i_relax].movable = 0;
                        break;
                    }
                }
            }
        }
    }
    
    /* update index numbers from all_atom list to be relax_atom list 
    for (i_relax =0; i_relax<n_relax; i_relax++) {
        int i_rotate;
        if (relax_atoms[i_relax].tor_atom1 == -1) continue;
        
        if (all_atoms[relax_atoms[i_relax].tor_atom1].i_relax == -1) {
            relax_atoms[i_relax].tor_atom1 = -1; continue;
        }
        else {
            relax_atoms[i_relax].tor_atom1 = all_atoms[relax_atoms[i_relax].tor_atom1].i_relax;
        }
        
        if (all_atoms[relax_atoms[i_relax].tor_atom2].i_relax == -1) {
            relax_atoms[i_relax].tor_atom1 = -1; continue;
        }
        else {
            relax_atoms[i_relax].tor_atom2 = all_atoms[relax_atoms[i_relax].tor_atom2].i_relax;
        }
        
        if (all_atoms[relax_atoms[i_relax].tor_atom3].i_relax == -1) {
            relax_atoms[i_relax].tor_atom1 = -1; continue;
        }
        else {
            relax_atoms[i_relax].tor_atom3 = all_atoms[relax_atoms[i_relax].tor_atom3].i_relax;
        }
        
        for (i_rotate = 0; i_rotate < relax_atoms[i_relax].n_rotate1; i_rotate++) {
            relax_atoms[i_relax].rotate1_lst[i_rotate]
            = all_atoms[relax_atoms[i_relax].rotate1_lst[i_rotate]].i_relax;
        }
        for (i_rotate = 0; i_rotate < relax_atoms[i_relax].n_rotate2; i_rotate++) {
            relax_atoms[i_relax].rotate2_lst[i_rotate]
            = all_atoms[relax_atoms[i_relax].rotate2_lst[i_rotate]].i_relax;
        }
    }
    */
}

int in_relax_list(int ia) {
     /* check if ia is in relax_atoms list */
    int i_relax;
    
    if (all_atoms[ia].i_relax != -1) return 1;
    
     /* If it's a sidechain atom, make sure there's no another sidechain in the relaxed atoms */
     
    if (all_atoms[ia].atom_p->i_conf_res) {
        for (i_relax = 0; i_relax < n_relax; i_relax++) {
             /* If it's in the same residue with an already relaxed atom */ 
            if (relax_atoms[i_relax].atom_p->i_res_prot != all_atoms[ia].atom_p->i_res_prot) continue;
             /* If that already relaxed atom is in a sidechain */ 
            if (!relax_atoms[i_relax].atom_p->i_conf_res) continue;
             /* If they belong to different sidechain */ 
            if (relax_atoms[i_relax].atom_p->i_conf_res != all_atoms[ia].atom_p->i_conf_res) {
                return 1;
            }
        }
    }
    return 0;
}

void setup_nghlst(PROT prot) {
    int i_relax, j_relax, ia, ja, i_constr, i_connect, i_ngh;
    int i_res, i_conf, i_atom, n_connect14;
    //float cutoff_far2  = VDW_CUTOFF_FAR  * VDW_CUTOFF_FAR;
    
    for (i_relax =0; i_relax<n_relax; i_relax++) {
        ia = relax_atoms[i_relax].atom_p->i_atom_prot;
        relax_atoms[i_relax].n_ngh = 0;
        relax_atoms[i_relax].ngh = NULL;
        relax_atoms[i_relax].n_ngh14 = 0;
        relax_atoms[i_relax].ngh14 = NULL;
        relax_atoms[i_relax].n_constr = 0;
        relax_atoms[i_relax].constr_list = NULL;
        relax_atoms[i_relax].constr_dsq = NULL;
        relax_atoms[i_relax].rotate1_lst =  NULL;
        relax_atoms[i_relax].rotate2_lst =  NULL;
        
        i_res = relax_atoms[i_relax].atom_p->i_res_prot;
        i_conf = relax_atoms[i_relax].atom_p->i_conf_res;
        i_atom = relax_atoms[i_relax].atom_p->i_atom_conf;
        
        /* no need to set up if this atom can't move */
        if (relax_atoms[i_relax].movable) {
            
            /* constr list */
            for (i_constr=0; i_constr<all_atoms[ia].n_constr; i_constr++) {
                j_relax = all_atoms[ia].constr_list[i_constr]->i_relax;
                
                if (j_relax != -1) {
                    relax_atoms[i_relax].n_constr++;
                    relax_atoms[i_relax].constr_list = realloc(relax_atoms[i_relax].constr_list, relax_atoms[i_relax].n_constr*sizeof(RELAX *));
                    relax_atoms[i_relax].constr_list[relax_atoms[i_relax].n_constr-1] = &relax_atoms[all_atoms[ia].constr_list[i_constr]->i_relax];
                    relax_atoms[i_relax].constr_dsq = realloc(relax_atoms[i_relax].constr_dsq, relax_atoms[i_relax].n_constr*sizeof(double));
                    relax_atoms[i_relax].constr_dsq[relax_atoms[i_relax].n_constr-1] = all_atoms[ia].constr_dsq[i_constr];
                }
            }
            
            /* 14 list */
            n_connect14 = prot.res[i_res].n_connect14[i_conf][i_atom];
            for (i_connect=0; i_connect < n_connect14; i_connect++) {
                ja = prot.res[i_res].connect14[i_conf][i_atom][i_connect]->i_atom_prot;
                j_relax = all_atoms[ja].i_relax;
                if (j_relax == -1) continue;
                relax_atoms[i_relax].n_ngh14++;
                relax_atoms[i_relax].ngh14 = realloc(relax_atoms[i_relax].ngh14, relax_atoms[i_relax].n_ngh14*sizeof(RELAX *));
                relax_atoms[i_relax].ngh14[relax_atoms[i_relax].n_ngh14-1] = &relax_atoms[all_atoms[ja].i_relax];
            }
            
            /* ngh list */
            for (j_relax =0; j_relax<n_relax; j_relax++) {
                if (i_relax == j_relax) continue;
                //if (ddvv(relax_atoms[i_relax].r, relax_atoms[j_relax].r) > cutoff_far2) continue;
                
                /* exclude members in constraint list */
                for (i_constr=0; i_constr<relax_atoms[i_relax].n_constr; i_constr++) {
                    if (relax_atoms[i_relax].constr_list[i_constr] == &relax_atoms[j_relax]) break;
                }
                if (i_constr<relax_atoms[i_relax].n_constr) continue;
                
                /* exclude members in ngh14 list */
                for (i_ngh=0; i_ngh<relax_atoms[i_relax].n_ngh14; i_ngh++) {
                    if (relax_atoms[i_relax].ngh14[i_ngh] == &relax_atoms[j_relax]) break;
                }
                if (i_ngh<relax_atoms[i_relax].n_ngh14) continue;
                
                relax_atoms[i_relax].n_ngh++;
                relax_atoms[i_relax].ngh = realloc(relax_atoms[i_relax].ngh, relax_atoms[i_relax].n_ngh*sizeof(RELAX *));
                relax_atoms[i_relax].ngh[relax_atoms[i_relax].n_ngh-1] = &relax_atoms[j_relax];
            }
        }
        
        /* torsion list */
        ja = all_atoms[ia].tor_atom1;
        if (ja != -1) {
            relax_atoms[i_relax].tor_atom1 = all_atoms[ja].i_relax;
        }
        
        ja = relax_atoms[i_relax].tor_atom2;
        if (ja != -1) {
            relax_atoms[i_relax].tor_atom2 = all_atoms[ja].i_relax;
            /* disable torsion calculation */
            if (relax_atoms[i_relax].tor_atom2 == -1) {
                relax_atoms[i_relax].tor_atom1 = -1;
            }
        }

        ja = relax_atoms[i_relax].tor_atom3;
        if (ja != -1) {
            relax_atoms[i_relax].tor_atom3 = all_atoms[ja].i_relax;
            /* disable torsion calculation */
            if (relax_atoms[i_relax].tor_atom3 == -1) {
                relax_atoms[i_relax].tor_atom1 = -1;
            }
        }
        
        if (relax_atoms[i_relax].tor_atom1 != -1) {
            int i_rotate;
            relax_atoms[i_relax].rotate1_lst = malloc(relax_atoms[i_relax].n_rotate1*sizeof(int));
            relax_atoms[i_relax].rotate2_lst = malloc(relax_atoms[i_relax].n_rotate2*sizeof(int));
            
            for (i_rotate = 0; i_rotate < relax_atoms[i_relax].n_rotate1; i_rotate++) {
                relax_atoms[i_relax].rotate1_lst[i_rotate]
                = all_atoms[all_atoms[ia].rotate1_lst[i_rotate]].i_relax;
            }
            for (i_rotate = 0; i_rotate < relax_atoms[i_relax].n_rotate2; i_rotate++) {
                relax_atoms[i_relax].rotate2_lst[i_rotate]
                = all_atoms[all_atoms[ia].rotate2_lst[i_rotate]].i_relax;
            }
        }
    }
}

#define N_TRIAL_MAX     10
int pick_sidechain(int i_res, PROT prot) {
    int i_conf, j_relax_res, j_res, j_conf, i_conf_min;
    int clash = 1, n_trial = 0;
    float pair_vdw_hv, e, e_min;
    
    while (clash) {
        if (n_trial > N_TRIAL_MAX) break; /* after trying n times, stop */
        
        i_conf = (int)(ran2(&idum) * (float)(prot.res[i_res].n_conf-1)) + 1;
        
        /* see if i_conf clash to any others */
        clash = 0;
        for (j_relax_res=0; j_relax_res<n_relax_res; j_relax_res++) {
            j_res = relax_res[j_relax_res];
            j_conf = relax_conf[j_relax_res];
            if (!j_conf) continue;
            
            pair_vdw_hv = vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 1);
            if (pair_vdw_hv > env.hv_relax_hv_vdw_thr) {
                break;
            }
        }
        if (j_relax_res<n_relax_res) {
            clash = 1;
            n_trial++;
        }
    }
    
    /* pick the lowest energy one */
    if (n_trial > N_TRIAL_MAX) {
        for (i_conf=1;i_conf<prot.res[i_res].n_conf;i_conf++) {
            /* calc. total interaction */
            e=0;
            for (j_relax_res=0; j_relax_res<n_relax_res; j_relax_res++) {
                j_res = relax_res[j_relax_res];
                j_conf = relax_conf[j_relax_res];
                if (!j_conf) continue;
                e += vdw_conf_fast(i_res, i_conf, j_res, j_conf, prot, 1);
            }
            
            if (e < e_min || i_conf == 1) {
                e_min = e;
                i_conf_min = i_conf;
            }
        }
        i_conf = i_conf_min;
    }
    return i_conf;
}

void collect_ngh(int i_res, PROT prot) {
    int i_ngh, k_res, k_conf;
    /* collect ngh list of i_res */
    for (i_ngh=0;i_ngh<prot.res[i_res].n_ngh;i_ngh++) {
        k_res = prot.res[i_res].ngh[i_ngh]->i_res_prot;
        if (!res_in_relax(k_res)) {
            add_conf2relax(&prot.res[k_res].conf[0], 1);
            if (prot.res[k_res].n_conf > 1) {
                k_conf = pick_sidechain(k_res, prot);
                add_conf2relax(&prot.res[k_res].conf[k_conf], 1);
            }
            else {
                n_relax_res++;
                relax_res[n_relax_res-1] = k_res;
                relax_conf[n_relax_res-1] = 0;
            }
        }
    }
}

int res_in_relax(int i_res) {
    int j_relax_res, j_res;
    for (j_relax_res=0; j_relax_res<n_relax_res; j_relax_res++) {
        j_res = relax_res[j_relax_res];
        if (i_res == j_res) break;
    }
    if (j_relax_res<n_relax_res) return 1;
    else return 0;
}

int closer_than(int i_res, int i_conf, int j_res, int j_conf, PROT prot, float crg_thr, float dist_thr)
{
    float d2 = dist_thr*dist_thr;
    int i_atom, j_atom;
    for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
        ATOM *iatom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
        if (!iatom_p->on) continue;
        if (fabs(iatom_p->crg) < crg_thr) continue;
                
        for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
            ATOM *jatom_p = &prot.res[j_res].conf[j_conf].atom[j_atom];
            if (!jatom_p->on) continue;
            if (fabs(jatom_p->crg) < crg_thr) continue;
                
                if (ddvv(iatom_p->xyz, jatom_p->xyz) < d2) return 1;
        }
    }
    
    return 0; /* can't find two charged atoms within dist_thr */
}

int write_relax_pdb(char *relax_pdb_filename, PROT prot) {
    int i_relax;
    if (relax_pdb_filename) {
        FILE *relax_fp = fopen(relax_pdb_filename,"w");
        if (relax_fp) {
            for (i_relax =0; i_relax<n_relax; i_relax++) {
                int i_res_pr, i_conf_pr, i_atom_pr;
                i_res_pr = relax_atoms[i_relax].atom_p->i_res_prot;
                i_conf_pr = relax_atoms[i_relax].atom_p->i_conf_res;
                i_atom_pr = relax_atoms[i_relax].atom_p->i_atom_conf;
                fprintf(relax_fp, "ATOM        %4s%c%3s %c%04d%c%03d%8.3f%8.3f%8.3f %7.3f      %6.3f      %-11s\n",
                    prot.res[i_res_pr].conf[i_conf_pr].atom[i_atom_pr].name,
                    prot.res[i_res_pr].conf[i_conf_pr].altLoc,
                    prot.res[i_res_pr].resName,
                    prot.res[i_res_pr].chainID,
                    prot.res[i_res_pr].resSeq,
                    prot.res[i_res_pr].iCode,
                    i_conf_pr,
                    relax_atoms[i_relax].r.x,
                    relax_atoms[i_relax].r.y,
                    relax_atoms[i_relax].r.z,
                    prot.res[i_res_pr].conf[i_conf_pr].atom[i_atom_pr].rad,
                    prot.res[i_res_pr].conf[i_conf_pr].atom[i_atom_pr].crg,
                    prot.res[i_res_pr].conf[i_conf_pr].history);
            }
        }
        fclose(relax_fp);
    }
    return 0;
}


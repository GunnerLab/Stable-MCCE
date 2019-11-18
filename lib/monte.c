#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include "mcce.h"

/* normal pw is garanteed to be smaller than 2000. When it is is bigger than 5000, it is
 * represents a value refenced to a clashed pair of conformers.
 */

/* self energy in prot, prot_free, prot_fixed and copied to conflist
 * before sampling
 * pairwise in pairwise
 * counter in conflist
 */
/* This Monte Carlo Program adds entropy correction to each conformer type different by
 * proton and electron, so that the ionized and neutral will get equal chance regardless of
 * the population of conformers.
 * It has been found that the drift of calculated pKa with rotamer number is caused by this
 * entropy effect
 */
typedef struct {
    float vdw0;
    float vdw1;
    float tors;
    float ebkb;
    float dsol;
    float offset;
    float pHpK0;
    float EhEm0;
    float TS;
    float residues;
    float total;
    float *mfePair;
    float *crg;
    float *vdw;
    float *ele;
} MFE;  // all the mfe energy terms

typedef struct  {
    int n;      /* Number of free confs on this residue */
    int on;
    int *conf;  /* idices of free confs */
    MFE mfe;   // let residue include mfe
} RESIDUE;

typedef struct  {
    int n;
    int *res;
} BIGLIST;

typedef struct {
    char  resName[4];
    char  chainID;
    char  iCode;
    int   resSeq;
    int   e;
    int   H;
    char  DM;
    float E_TS;
} CONFTYPE;

typedef struct {
   int n;
   CONFTYPE *conftype;
}  TYPES;

/* for the microstate ---By Cai*/
typedef struct {
    unsigned short *conf_id;
    double H;  //to store acumulate energy value
//    double Hsq;
    int counter; //to store counter for MC sampling
    double occ;  // to store occ for enumeration
    double Hav; // to store average state energy
} MSRECORD;

int enum_flag = 0; //flag to judge if run enumerate(): analytical solution;

int   load_ms_gold(STRINGS *str);
int   update_conf_id(unsigned short *conf_id, int *state);
int   write_ms(MSRECORD *ms_state);
void MC_smp(int n);

STRINGS ms_spe_lst;
FILE   *ms_fp;
FILE   *re_ms_fp;  //readable ms.dat: re_ms.dat
/* for the microstate */


/* public variables */
PROT     prot;
RES      conflist, conflist_bak;
float    **pairwise, **pairwise_vdw, **pairwise_ele;
float    E_base, E_state;
float    ph, eh;
int      *state, *state_bak;

RESIDUE  *free_res, *fixed_res, *all_res, *mfe_res; // mfe_res is all the mfe residue
int      n_free, n_fixed, n_all, n_mfe;
BIGLIST  *biglist;   /* same length as free residues */
TYPES    Sconverge, SconvergeBak;
float    E_minimum;
FILE     *fp;
time_t   timerA, timerB, timerC;

float    **MC_occ;    /* occ of 3 parallel MC */
float    **occ_table;   /* occ of conformers at various pH/Eh */
float    **entropy_table;

float fround3(float x);
int   print_mfe(int i_res, float mfeP, FILE *pK_fp,  FILE *res_fp);
int   get_mfe(int i_res, int t_point);
int   verify_flag();
int   load_pairwise();
int   load_pairwise_fround3();
int   load_pairwise_vdw();
void  group_confs();
float get_base();
float get_E();
void  mk_neighbors();
void MC(int n);
int reduce_conflist();
int fitit();
int enumerate(int i_ph_eh);
int enumerate_new(int i_ph_eh); //Cai: new enumerate subroutine to ouput microstate 
int load_conflist();
int group_conftype();
int cmp_conftype(CONFTYPE t1, CONFTYPE t2);
CONFTYPE get_conftype(CONF conf);
void update_Sconvergence();
float s_stat();
float get_entropy(int iconf, int start, int end);
char **shead;
float get_totalTS();

#define mev2Kcal 0.0235  // to keep consistent with mfe.py, use this constant instead of PH2KCAL/58, just for mfe part
int monte()
{  int i, j, k, counter, k_run;
    float *sigma, sigma_max, t;
    int N_smp;
    float S_max;

    timerA = time(NULL);
    strcpy(env.entropy_converge_error, "");

    /* Load conformer list from FN_CONFLIST3 */
    printf("   Load conformer list from file \"%s\" ...\n", FN_CONFLIST3);
    fflush(stdout);
    prot = new_prot();
    load_conflist();
    if (conflist.n_conf == 0) {
        printf("   FATAL: error in reading conformer list %s", FN_CONFLIST3);
        return USERERR;
    }
    printf("   Scaling factors:\n");
    printf("   VDW0  = %.3f:\n", env.scale_vdw0);
    printf("   VDW1  = %.3f:\n", env.scale_vdw1);
    printf("   VDW   = %.3f:\n", env.scale_vdw);
    printf("   TORS  = %.3f:\n", env.scale_tor);
    printf("   ELE   = %.3f:\n", env.scale_ele);
    printf("   DSOLV = %.3f:\n", env.scale_dsolv);
    printf("   Done\n\n"); fflush(stdout);


    /* load pairwise */
    printf("   Load pairwise interactions ...\n"); fflush(stdout);
    if (load_pairwise()) {
        printf("   FATAL: pairwise interaction not loaded\n");
        return USERERR;
    }
    printf("   Done\n\n");
    fflush(stdout);

    /* Verify self-consistancy between occupancy and flag */
    printf("   Verify self consistancy of conformer flag ...\n");
    if (verify_flag()) {
        printf("   Conformer flags updated due to self consistancy.\n\n");
    }
    else {
        printf("   Done\n\n");
    }
    fflush(stdout);

    /* backup */
    conflist_bak.n_conf = 0; conflist_bak.conf = NULL;
    cpy_res(&conflist_bak, &conflist);

    all_res = NULL;
    free_res = NULL;
    fixed_res = NULL;
    biglist  = NULL;
    state = NULL;


    /* titration */
    printf("   Do titration at %d points...\n", env.titr_steps);
    printf("   Detailed progress is in file \"%s\"\n", MC_OUT);
    fflush(stdout);
    if (env.monte_seed < 0) srand(time(NULL));
    else srand(env.monte_seed);

    MC_occ = (float **) malloc(env.monte_runs * sizeof(float *));
    for (i=0; i<env.monte_runs; i++) {
        if (!(MC_occ[i] = (float *) malloc(conflist.n_conf * sizeof(float)))) {
            printf("   FATAL: memory error.\n");
            return USERERR;
        }

    }

    occ_table = (float **) malloc(conflist.n_conf*sizeof(float *));
    for (i=0; i<conflist.n_conf; i++) occ_table[i] = (float *) malloc(env.titr_steps * sizeof(float));

    entropy_table = (float **) malloc(conflist.n_conf*sizeof(float *));
    for (i=0; i<conflist.n_conf; i++) entropy_table[i] = (float *) malloc(env.titr_steps * sizeof(float));

    sigma = (float *) malloc(conflist.n_conf*sizeof(float *));

    if (!env.minimize_size) {
       if (!(fp = fopen(MC_OUT, "w"))) {
          printf("   FATAL: Can not write to file \"%s\"\n", MC_OUT);
          return USERERR;
       }
    }
    else {
       fp = tmpfile();
    }

    /* Microstate ---By Cai: inititalize ms.dat writing */
    if (env.ms_out) {
        memset(&ms_spe_lst, 0, sizeof(STRINGS));
        ms_spe_lst.n = 0;
        if (load_ms_gold(&ms_spe_lst)) {
            printf("   Load file ms_gold error: ms_gold file for Output Microstate (MONTE_MS) is missing. If you do not want microstates reported, turn (MS_OUT) to f. If you do want microstate reported, include ms_gold file.\n");
            printf("   Formate of ms_gold:\n   GLUA0286\n   TYRA0288\n   THRA0288\n   Or formate of pdb file.\n");
            return USERERR;
        }
        if (ms_spe_lst.n == 0) {
            printf("Load 0 ms gold residue, don't output ms_out\n");
            env.ms_out = 0;
        }
        else {
            if(env.re_ms_out) { //output microstate in binary file and readable file
                ms_fp = fopen("ms.dat", "wb");
                re_ms_fp=fopen("re_ms.dat", "w");
                int i_spe;
                fwrite(&ms_spe_lst.n, 1, sizeof(int), ms_fp);
//checkpoint
                fprintf(re_ms_fp, "residue number: %d\n", ms_spe_lst.n);
                for (i_spe=0; i_spe<ms_spe_lst.n; i_spe++) {
//checkpoint
                    fprintf(re_ms_fp,"%s\t", ms_spe_lst.strings[i_spe]);
                    fwrite(ms_spe_lst.strings[i_spe], 8, sizeof(char), ms_fp);
                }
                fprintf(re_ms_fp,"\n");
            }
            else { //only output microstate in binary formate
                ms_fp = fopen("ms.dat", "wb");
                int i_spe;
                fwrite(&ms_spe_lst.n, 1, sizeof(int), ms_fp);
//checkpoint
                //printf("       writing ms.dat...\nresidue number: %d\n", ms_spe_lst.n);
                for (i_spe=0; i_spe<ms_spe_lst.n; i_spe++) {
//checkpoint
                    //printf("       %s\n", ms_spe_lst.strings[i_spe]);
                    fwrite(ms_spe_lst.strings[i_spe], 8, sizeof(char), ms_fp);
                }
             }
        }    
    }   


    timerB = time(NULL);
    printf("   Monte Carlo set up time: %ld seconds.\n", timerB-timerA); fflush(stdout);


    /* entropy sampling and correction */
    /* group conformers */
    group_conftype();  /* This ininitializes Sconverge and SconvergeBak */

    for (i=0; i<env.titr_steps; i++) {
        timerB = time(NULL);

        ph = env.titr_ph0;
        eh = env.titr_eh0;
        if (env.titr_type == 'p') ph = env.titr_ph0 + i*env.titr_phd;
        else eh = env.titr_eh0 + i*env.titr_ehd;
        fprintf(fp, "Titration (%dth) at pH = %6.2f Eh = %6.2f:\n", i+1, ph, eh); fflush(fp);


        /* reduce prot to free residues, conformers fixed as 0 occ are excluded from this free residues */
        free(conflist.conf);
        conflist.n_conf = 0; conflist.conf = NULL;
        cpy_res(&conflist, &conflist_bak);

        group_confs();
        fprintf(fp, "Identified %d free residues from %d residues.\n", n_free, n_all);
        fflush(fp);

        /* get big list */
        mk_neighbors();

        /* Compute base energy, self and mfe */
        E_base = get_base();

        /* get a microstate */
        state = realloc(state, n_free*sizeof(int));
        for (j=0; j<n_free; j++)
            state[j] = free_res[j].conf[rand() / (RAND_MAX/free_res[j].n + 1)];

        /* DEBUG
        for (j=0; j<n_free; j++) printf("%03d ", state[j]);
        printf("\n");
        */

        /* Do annealing */
        /*
        fprintf(fp, "Conformer list before reduction: E_base = %10.3f\n", E_base);
        for (j=0; j<conflist.n_conf; j++) {
            fprintf(fp, "%s %c %4.2f self = %8.3f mfe = %8.3f\n", conflist.conf[j].uniqID,
            conflist.conf[j].on,
            conflist.conf[j].occ,
            conflist.conf[j].E_self0,
            conflist.conf[j].E_mfe);
        }
        fprintf(fp, "\n"); fflush(fp);
        */

        counter = 0;
        for (j=0; j<n_free; j++) counter+=free_res[j].n;
        N_smp = env.monte_nstart * counter;
        fprintf(fp, "Doing annealing... \n"); fflush(fp);
        if (N_smp) MC(N_smp);
        fprintf(fp, "Done\n\n");

        /* memcpy(state_check, state, sizeof(int)*n_free);  DEBUG */

        /* Do equalibriation */
        counter = 0;
        for (j=0; j<n_free; j++) counter+=free_res[j].n;
        N_smp = env.monte_neq * counter;
        fprintf(fp, "Doing equalibration ... \n"); fflush(fp);
        if (N_smp) MC(N_smp);
        fprintf(fp, "Done\n\n");

        /* do reduction */
        fprintf(fp, "Reducing conflist ... %d conformers are marked as fixed.\n\n", reduce_conflist());
        group_confs();
        fprintf(fp, "%d free residues from %d residues, after reduction.\n\n", n_free, n_all);
        fflush(fp);

        mk_neighbors();


        j = 0; S_max = 999.9;
        for (k=0; k<Sconverge.n; k++) SconvergeBak.conftype[k].E_TS = 0.0;
   
        if (env.monte_tsx) {
           while (S_max > 0.5 &&  j < 10) {
             /* a new state */
             for (k=0; k<n_free; k++)
                 state[k] = free_res[k].conf[rand() / (RAND_MAX/free_res[k].n + 1)];

             E_base = get_base();
             
             if (enumerate(i) == -1) { // Non-ANALYTICAL SOLOTION 
                counter = 0;
                for (k=0; k<n_free; k++) counter+=free_res[k].n;
                N_smp = env.monte_niter * counter;

                fprintf(fp, "Doing Entropy sampling cycle %d...\n", j+1); fflush(fp);
                if (env.monte_nstart * counter) MC(env.monte_nstart * counter);
                if (N_smp) MC(N_smp);
             }

             update_Sconvergence(); // calculate entropy from occupancy 
             S_max = s_stat();
             fprintf(fp, "Max delta = %8.3f\n", S_max);

             // backup 
             for (k=0; k<Sconverge.n; k++) SconvergeBak.conftype[k] = Sconverge.conftype[k];

             j++;
           }
            if (j > 10) {
               fprintf(fp, "WARNING: entropy correction didn't converge.\n");
               strcpy(env.entropy_converge_error, "WARNING: entropy correction didn't converge at one point, see mc_out for details.");
            }
           fprintf(fp, "Done, exit at max entropy convergence %.3f \n\n", S_max);
        }



        E_base = get_base();
        /* ANALYTICAL SOLOTION */
        //if (enumerate(i) == -1) {
        if (enumerate_new(i) == -1) { // use new enumerate subroutine to output microstate
            if (env.ms_out){
                fwrite("MONTERUNS", 9, sizeof(char), ms_fp); //The third line of ms.dat: method
                if(env.re_ms_out) { //write out the method
                    fprintf(re_ms_fp, "METHOD: %s\n", "MONTERUNS");
                }
            }

            /*
            fprintf(fp, "Conformer list after reduction: E_base = %10.3f\n", E_base);
            for (j=0; j<conflist.n_conf; j++) {
                fprintf(fp, "%s %c %4.2f self = %8.3f mfe = %8.3f\n", conflist.conf[j].uniqID,
                conflist.conf[j].on,
                conflist.conf[j].occ,
                conflist.conf[j].E_self0,
                conflist.conf[j].E_mfe);
            }
            fprintf(fp, "\n"); fflush(fp);
            */
            counter = 0;
            for (j=0; j<n_free; j++) counter+=free_res[j].n;
            N_smp = env.monte_niter * counter;

            for (j=0; j<env.monte_runs; j++) {
                /* a new state */
                for (k=0; k<n_free; k++)
                    state[k] = free_res[k].conf[rand() / (RAND_MAX/free_res[k].n + 1)];

                fprintf(fp, "Doing annealing of MC %2d ...\n", j+1); fflush(fp);
                if (env.monte_nstart * counter) MC(env.monte_nstart * counter);

                fprintf(fp, "Doing MC %2d ... \n", j+1); fflush(fp);
                //if (N_smp) MC(N_smp);
                if (N_smp) {            //Cai: microstate output or not
                    if (env.ms_out) MC_smp(N_smp);   // Using MC_smp to write out microstate
                    else MC(N_smp);  //initial MC without writing out microstate  
                }

                for (k=0; k<conflist.n_conf; k++) {
                    MC_occ[j][k] = conflist.conf[k].occ;
                }
                fprintf(fp, "Done\n\n");
                 /* average -Yifan */
                for (k=0; k<conflist.n_conf; k++) conflist.conf[k].occ = 0.0;
                for (k_run=0; k_run<j+1; k_run++) {
                    for (k=0; k<conflist.n_conf; k++) {
                        conflist.conf[k].occ +=  MC_occ[k_run][k]/(j+1);
                    }
                }
                if (env.monte_tsx) {
                    update_Sconvergence(); /* calculate entropy from occupancy */
                    E_base = get_base();
                }
            }

            /* average */
            for (k=0; k<conflist.n_conf; k++) conflist.conf[k].occ = 0.0;
            for (j=0; j<env.monte_runs; j++) {
                for (k=0; k<conflist.n_conf; k++) {
                    conflist.conf[k].occ +=  MC_occ[j][k]/env.monte_runs;
                }
            }

            /* standard deviation */
            sigma_max = 0.0;
            for (j=0; j<conflist.n_conf; j++) {
                t = 0.0;
                for (k=0; k<env.monte_runs; k++) {
                    t += (MC_occ[k][j] - conflist.conf[j].occ) * (MC_occ[k][j] - conflist.conf[j].occ) ;
                }
                if (env.monte_runs > 1) sigma[j] = sqrt(t/(env.monte_runs-1));
                else sigma[j] = 999.00;

                if (sigma_max < sigma[j]) sigma_max = sigma[j];
            }
        }
        /* write statstics */
        fprintf(fp, "Conformer     flag   E_self");
        for (j=0; j<env.monte_runs; j++) fprintf(fp, "   mc%02d", j+1);
        fprintf(fp, "      occ Sgm(n-1)\n");
        for (j=0; j<conflist.n_conf; j++) {
            occ_table[j][i] = conflist.conf[j].occ;
            entropy_table[j][i] = conflist.conf[j].E_TS;
            fprintf(fp, "%s   %c %8.3f", conflist.conf[j].uniqID,
            conflist.conf[j].on,
            conflist.conf[j].E_self0);
            for (k=0; k<env.monte_runs; k++) fprintf(fp, " %6.3f", MC_occ[k][j]);
            fprintf(fp, " Av=%5.3f Sg=%5.3f\n", conflist.conf[j].occ, sigma[j]);
        }
        fprintf(fp, "\n"); fflush(fp);

        timerC = time(NULL);
        printf("   Titration %2d: %5ld seconds, biggest stdev of conformer occ = %5.3f\n",
        i+1, timerC-timerB, sigma_max);
        fflush(stdout);

    }

    //fclose(fp);
    if (env.ms_out && env.re_ms_out){ // close the microstate output file if true
    fclose(fp);
    fclose(ms_fp);
    fclose(re_ms_fp);}
    else if (env.ms_out && !env.re_ms_out){
    fclose(fp);
    fclose(ms_fp);}      // output microstate -- By Cai
    else if (!env.ms_out && env.re_ms_out)
    printf("   Warning: need to turn on (MS_OUT) in run.prm to get readable microstate file.\n");
    else
    fclose(fp);

    printf("   Done\n\n"); fflush(stdout);


    /* writing occ table */
    if (!(fp = fopen(OCC_TABLE, "w"))) {
        printf("   FATAL: Can not write occupancy to file \"%s\"\n", OCC_TABLE);
        return USERERR;
    }

   //if (env.entropy_converge_error) {
  //    fprintf(fp, "%s\n", env.entropy_converge_error);
  // }

    if (env.titr_type == 'p') {
        fprintf(fp, " ph           ");
        for (i=0; i<env.titr_steps; i++) fprintf(fp, " %5.1f", env.titr_ph0+i*env.titr_phd);
    }
    else {
        fprintf(fp, " eh           ");
        for (i=0; i<env.titr_steps; i++) fprintf(fp, " %5.0f", env.titr_eh0+i*env.titr_ehd);
    }
    fprintf(fp, "\n");
    for (i=0; i<conflist.n_conf; i++) {
        fprintf(fp, "%s", conflist.conf[i].uniqID);
        for (j=0; j<env.titr_steps; j++) {
            fprintf(fp, " %5.3f", occ_table[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    /* writing entropy table */
    if (!(fp = fopen("entropy.out", "w"))) {
        printf("   FATAL: Can not write occupancy to file \"%s\"\n", "entropy.out");
        return USERERR;
    }
    if (env.titr_type == 'p') {
        fprintf(fp, " ph           ");
        for (i=0; i<env.titr_steps; i++) fprintf(fp, " %5.1f", env.titr_ph0+i*env.titr_phd);
    }
    else {
        fprintf(fp, " eh           ");
        for (i=0; i<env.titr_steps; i++) fprintf(fp, " %5.0f", env.titr_eh0+i*env.titr_ehd);
    }
    fprintf(fp, "\n");
    for (i=0; i<conflist.n_conf; i++) {
        fprintf(fp, "%s", conflist.conf[i].uniqID);
        for (j=0; j<env.titr_steps; j++) {
            fprintf(fp, " %5.3f", entropy_table[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);




    /* write self energy terms in right unit */


    /* curve fitting */
    printf("   Fit titration curves to get pKa/Em ...\n");
    fflush(stdout);
    if (fitit()) {
        printf("   Fatal error detected in fitting program.\n");
        return USERERR;
    }
    printf("   Done\n\n"); fflush(stdout);

    timerC=time(NULL);
    printf("   Total time on MC: %ld seconds\n\n", timerC-timerA);

    printf("   Output files:\n");
    printf("      %-16s: MC progress and convergence.\n", MC_OUT);
    printf("      %-16s: Occupancy table.\n", OCC_TABLE);
    printf("      %-16s: Entropy correction table\n", "entropy.out");
    printf("      %-16s: pKa or Em from titration curve fitting.\n", CURVE_FITTING);
    printf("      %-16s: Summary of residue charges.\n", "sum_crg.out");


    /* free MC stat */
    for (i=0; i<env.monte_runs; i++) {
        free(MC_occ[i]);
    }
    free(MC_occ);
    free(sigma);

    /* free conflist */
    free(conflist.conf);
    free(conflist_bak.conf);

    /* free pairwise */
    for (i=0; i<conflist.n_conf; i++)
        free(pairwise[i]);
    free(pairwise);

    /* free residues */
    for (i=0; i<n_free; i++)
        free(free_res[i].conf);
    free(free_res);

    for (i=0; i<n_fixed; i++)
        free(fixed_res[i].conf);
    free(fixed_res);

    for (i=0; i<n_all; i++)
        free(all_res[i].conf);
    free(all_res);

    free(Sconverge.conftype);
    //free(SconvergeBak.conftype);


    return 0;
}

int load_conflist()
{  FILE *fp;
    char sbuff[MAXCHAR_LINE];
    char stemp[MAXCHAR_LINE];
    CONF conf_temp;
    char notfound;
    int iconf, ires;
    int kr;
    int counter;

    conflist.n_conf = 0;
    conflist.conf   = NULL;

    if (!(fp=fopen(FN_CONFLIST3, "r"))) {
        printf("   FATAL: Can't open file %s\n", FN_CONFLIST3);
        return USERERR;
    }
    fgets(sbuff, sizeof(sbuff), fp); /* skip the first line */
    counter = 0;
    while(fgets(sbuff, sizeof(sbuff), fp)) {
        /* load this line to a conf template */
        if (strlen(sbuff) < 20) continue;
        sscanf(sbuff, "%d %s %c %f %f %f %f %d %d %f %f %f %f %f %f %s", &conf_temp.iConf,
        conf_temp.uniqID,
        &conf_temp.on,
        &conf_temp.occ,
        &conf_temp.netcrg,
        &conf_temp.Em,
        &conf_temp.pKa,
        &conf_temp.e,
        &conf_temp.H,
        &conf_temp.E_vdw0,
        &conf_temp.E_vdw1,
        &conf_temp.E_tors,
        &conf_temp.E_epol,
        &conf_temp.E_dsolv,
        &conf_temp.E_extra,
        conf_temp.history);

        conf_temp.E_TS = 0.0; /* initialize entropy effect at the time of loading conflist */

        /* rescale */
        conf_temp.E_vdw0 *= env.scale_vdw0;
        conf_temp.E_vdw1 *= env.scale_vdw1;
        conf_temp.E_epol *= env.scale_ele;
        conf_temp.E_tors *= env.scale_tor;
        conf_temp.E_dsolv*= env.scale_dsolv;

        strncpy(conf_temp.resName, conf_temp.uniqID, 3); conf_temp.resName[3] = '\0';
        strncpy(conf_temp.confName, conf_temp.uniqID, 5); conf_temp.confName[5] = '\0';
        conf_temp.chainID = conf_temp.uniqID[5];
        strncpy(stemp, conf_temp.uniqID+6, 4); stemp[4] = '\0';
        conf_temp.resSeq = atoi(stemp);
        conf_temp.iCode = conf_temp.uniqID[10];
        conf_temp.n_atom = 0;
        if (conf_temp.on == 't' || conf_temp.on == 'T') conf_temp.on = 't';
        else conf_temp.on = 'f';
        conf_temp.iConf = counter;
        /* creating conflist */
        iconf = ins_conf(&conflist, conflist.n_conf, 0);
        cpy_conf(&conflist.conf[iconf], &conf_temp);
        counter++;

        notfound = 1;
        for (kr=0; kr<prot.n_res; kr++) {
            if (conf_temp.chainID != prot.res[kr].chainID ||
                strncmp(conf_temp.resName, prot.res[kr].resName, 3) ||
            conf_temp.resSeq != prot.res[kr].resSeq ||
            conf_temp.iCode != prot.res[kr].iCode)
            continue;
            /* residue found */
            iconf = ins_conf(&prot.res[kr], prot.res[kr].n_conf, 0);
            cpy_conf(&prot.res[kr].conf[iconf], &conf_temp);

            notfound = 0;
            break;
        }
        if (notfound) { /* belongs to new residue */
            ires  = ins_res(&prot, prot.n_res);
            iconf = ins_conf(&prot.res[kr], prot.res[kr].n_conf, 0);

            cpy_conf(&prot.res[kr].conf[iconf], &conf_temp);

            strncpy(prot.res[ires].resName, conf_temp.resName,3);
            prot.res[ires].chainID = conf_temp.chainID;
            prot.res[ires].resSeq = conf_temp.resSeq;
            prot.res[ires].iCode = conf_temp.iCode;
        }
    }
    return 0;
}

int verify_flag()
/* this program will update conflist, the backup is in conflist_bak */
{  PROT old_prot;
    int kr, kc;
    float socc;
    int changed = 0;
    int n_free;


    memset(&old_prot, 0, sizeof(PROT));

    cpy_prot(&old_prot, &prot);

    for (kr=0; kr<prot.n_res; kr++) {
        /* set occ of free conformer to be 0 */
        for (kc=0; kc<prot.res[kr].n_conf; kc++) {
            if (prot.res[kr].conf[kc].on == 'f') {
                prot.res[kr].conf[kc].occ = 0.0;
                if (fabs(old_prot.res[kr].conf[kc].occ) > 0.001) {
                    printf("   %s %c %4.2f -> %c %4.2f (free conformer has 0.0 occ)\n",
                    old_prot.res[kr].conf[kc].uniqID,
                    old_prot.res[kr].conf[kc].on,
                    old_prot.res[kr].conf[kc].occ,
                    prot.res[kr].conf[kc].on,
                    prot.res[kr].conf[kc].occ);
                    changed++;
                }
            }
        }

        /* set t for single free conformer */
        socc = 0.0;
        n_free = prot.res[kr].n_conf;
        for (kc=0; kc<prot.res[kr].n_conf; kc++) {
            if (prot.res[kr].conf[kc].on == 't') {
                n_free --;
                socc += prot.res[kr].conf[kc].occ;
            }
        }
        if (n_free == 1) {         /* single free conf, convert it to t */
            for (kc=0; kc<prot.res[kr].n_conf; kc++) {
                if (prot.res[kr].conf[kc].on == 'f') {
                    prot.res[kr].conf[kc].on = 't';
                    prot.res[kr].conf[kc].occ = 1.0-socc;
                    if (prot.res[kr].conf[kc].occ < 0.0) prot.res[kr].conf[kc].occ = 0.0;
                    if (prot.res[kr].conf[kc].occ > 1.0) prot.res[kr].conf[kc].occ = 1.0;
                    printf("   %s %c %4.2f -> %c %4.2f (single conformer treated as fixed)\n",
                    old_prot.res[kr].conf[kc].uniqID,
                    old_prot.res[kr].conf[kc].on,
                    old_prot.res[kr].conf[kc].occ,
                    prot.res[kr].conf[kc].on,
                    prot.res[kr].conf[kc].occ);
                    changed++;
                }
            }
        }
        else if (n_free>1) {      /* multiple free conf, free confs take 1.00 occ */
            for (kc=0; kc<prot.res[kr].n_conf; kc++) {
                if (prot.res[kr].conf[kc].on == 't' && prot.res[kr].conf[kc].occ > 0.001) {
                    prot.res[kr].conf[kc].occ = 0.0;
                    printf("   %s %c %4.2f -> %c %4.2f (with multiple free conformers)\n",
                    old_prot.res[kr].conf[kc].uniqID,
                    old_prot.res[kr].conf[kc].on,
                    old_prot.res[kr].conf[kc].occ,
                    prot.res[kr].conf[kc].on,
                    prot.res[kr].conf[kc].occ);
                    changed++;
                }
            }
        }
    }

    del_prot(&old_prot);

    /* update conflist */
    for (kr=0; kr<prot.n_res; kr++) {
        for (kc=0; kc<prot.res[kr].n_conf; kc++) {
            conflist.conf[prot.res[kr].conf[kc].iConf].occ = prot.res[kr].conf[kc].occ;
            conflist.conf[prot.res[kr].conf[kc].iConf].on  = prot.res[kr].conf[kc].on;
        }
    }

    return changed;
}

int load_pairwise()
{   int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;
    if (load_energies(&ematrix, ".")<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    /* scan conformers to see if all pw were calculated, DISABLED TEMPERORY
    for (i=0; i<ematrix.n; i++) {
        if (!ematrix.conf[i].on && (ematrix.conf[i].uniqID[3] != 'D' || ematrix.conf[i].uniqID[4] != 'M')) {
           printf("      Incompleted delphi run, the first place detected at %s\n ", ematrix.conf[i].uniqID);
           return USERERR;
        }
    }
    */


    if (!(pairwise = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }


    for (kc=0; kc<ematrix.n; kc++) {
        if (!(pairwise[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }

    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
           /*
           if (ematrix.pw[i][j].vdw > 990. && ematrix.pw[j][i].vdw > 990.) {
               pairwise[i][j] = pairwise[j][i] = 999.;
           }
           else {
               pairwise[i][j] = pairwise[j][i] = ((ematrix.pw[i][j].ori + ematrix.pw[j][i].ori)*env.scale_ele \
                   +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           }
           */
           pairwise[i][j] = pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele
                                            +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           //pairwise[i][j] = fround3(ematrix.pw[i][j].ele) *env.scale_ele + fround3(ematrix.pw[i][j].vdw)*env.scale_vdw;
           /* proprocessing
           if ((ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw) > 999.0) {
              pairwise[i][j] = pairwise[j][i] = 999.0;
           }
           else {
              pairwise[i][j] = pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele \
                                                +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
           }
           */
        }
    }

    /* DEBUG print pairwise table
    int kr;
    printf("                ");
    for (kc=0; kc<conflist.n_conf; kc++) {
        printf("%16s", conflist.conf[kc].uniqID);
    }
    printf("\n");
    for (kr=0; kr<conflist.n_conf; kr++) {
        printf("%s ", conflist.conf[kr].uniqID);
        for (kc=0; kc<conflist.n_conf; kc++) {
            printf("%16.3f", pairwise[kr][kc]);
        }
        printf("\n");
    }
    */


    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

int load_pairwise_vdw()
{
    int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;
    if (load_energies(&ematrix, ".")<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    if (!(pairwise_vdw = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }
    for (kc=0; kc<ematrix.n; kc++) {
        if (!(pairwise_vdw[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }

    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
           pairwise_vdw[i][j] = pairwise_vdw[j][i] = (ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw/2.0 ;
        }
    }

    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

void group_confs()
{  int  kr, kc;
    int free_conf;
    int ir, ic;

    /* delete existing big list when regrouping */
    if (biglist != NULL) {
        for (ir=0; ir<n_free; ir++) {
            if (biglist[ir].res != NULL)
                free(biglist[ir].res);
        }
        free(biglist);
    }


    /* create all_res */
    if (all_res != NULL) {
        for (ir=0; ir<n_all; ir++) {
            if (all_res[ir].conf != NULL)
                free(all_res[ir].conf);
        }
        free(all_res);
    }
    all_res = (RESIDUE *) malloc(prot.n_res * sizeof(RESIDUE));
    n_all = prot.n_res;
    for (ir=0; ir<n_all; ir++) {
        all_res[ir].n    = prot.res[ir].n_conf;
        all_res[ir].conf = (int *) malloc(prot.res[ir].n_conf * sizeof(int));
        all_res[ir].mfe.crg = (float *) calloc(env.titr_steps, sizeof(float));
    }
    for (ir=0; ir<n_all; ir++) {
        for (ic=0; ic<all_res[ir].n; ic++) {
            all_res[ir].conf[ic] = prot.res[ir].conf[ic].iConf;
        }
    }


    /* create free residues */
    if (free_res != NULL) {
        for (ir=0; ir<n_free; ir++) {
            free(free_res[ir].conf);
        }
        free(free_res);
        free_res = NULL;
    }

    if (fixed_res != NULL) {
        for (ir=0; ir<n_fixed; ir++) {
            free(fixed_res[ir].conf);
        }
        free(fixed_res);
        fixed_res = NULL;
    }

    n_free = 0; n_fixed = 0;
    for (kr=0; kr<n_all; kr++) {
        free_conf = 0;
        for (kc=0; kc<all_res[kr].n; kc++) {
            if (conflist.conf[all_res[kr].conf[kc]].on == 'f') free_conf ++;
            /*
            printf("%s, %c, %f\n", conflist.conf[all_res[kr].conf[kc]].uniqID,
            conflist.conf[all_res[kr].conf[kc]].on,
            conflist.conf[all_res[kr].conf[kc]].occ);
            */
        }
        if (free_conf) {       /* free_conf stores number of free conformers in this residue */
            n_free++;
            free_res = realloc(free_res, n_free * sizeof(RESIDUE));
            free_res[n_free-1].n = free_conf;
            free_res[n_free-1].conf = (int *) malloc(free_conf*sizeof(int));

            ic = 0;
            for (kc=0; kc<all_res[kr].n; kc++) {
                if (conflist.conf[all_res[kr].conf[kc]].on == 'f') {
                    free_res[n_free-1].conf[ic] = conflist.conf[all_res[kr].conf[kc]].iConf;
                    ic++;
                }
                else { /* 't' must be with occ = 0.0, ignored */
                    if (fabs(conflist.conf[all_res[kr].conf[kc]].occ) > 0.000000001)  {
                        printf("   WARNING: %s has a fixed occ %.3f in a free residue\n",
                        conflist.conf[all_res[kr].conf[kc]].uniqID, conflist.conf[all_res[kr].conf[kc]].occ);
                        fflush(stdout);
                    }
                }
            }
        }
        else {
            n_fixed++;
            fixed_res = realloc(fixed_res, n_fixed * sizeof(RESIDUE));
            fixed_res[n_fixed-1].n = all_res[kr].n;
            fixed_res[n_fixed-1].conf = (int *) malloc(all_res[kr].n*sizeof(int));

            ic = 0;
            for (kc=0; kc<all_res[kr].n; kc++) {
                fixed_res[n_fixed-1].conf[ic] = conflist.conf[all_res[kr].conf[kc]].iConf;
                ic++;
            }
        }
    }

    return;
}

float get_base()
{  float E;
    int kr, kc, ir, ic;
    float mfe;

    E = 0.0;

    /* compute pH and Eh effect for each conformer */
    for (kc=0; kc<conflist.n_conf; kc++) {
        conflist.conf[kc].E_ph =
        env.monte_temp/ROOMT * conflist.conf[kc].H * (ph-conflist.conf[kc].pKa) * PH2KCAL;
        conflist.conf[kc].E_eh =
        env.monte_temp/ROOMT * conflist.conf[kc].e * (eh-conflist.conf[kc].Em) * PH2KCAL/58.0;
    }

    /* self without mfe */
    for (kc=0; kc<conflist.n_conf; kc++) {
        conflist.conf[kc].E_self0 = conflist.conf[kc].E_vdw0
        + conflist.conf[kc].E_vdw1
        + conflist.conf[kc].E_epol
        + conflist.conf[kc].E_tors
        + conflist.conf[kc].E_dsolv
        + conflist.conf[kc].E_extra
        + conflist.conf[kc].E_ph
        + conflist.conf[kc].E_eh
        + conflist.conf[kc].E_TS; /* Monte Carlo "knows" entropy, we want to "correct" it, so the signe is "+" */
/*        if (conflist.conf[kc].E_self0 > 900.0)
           conflist.conf[kc].E_self0 = 999.0;
*/
    }

    /* mfe on all conformers from fixed conformers */
    for (kr=0; kr<n_all; kr++) {
        for (kc=0; kc<all_res[kr].n; kc++) {
            mfe = 0.0;
            for (ir=0; ir<n_fixed; ir++) {
                for (ic=0; ic<fixed_res[ir].n; ic++) {
                    mfe += pairwise[all_res[kr].conf[kc]][fixed_res[ir].conf[ic]]
                    *conflist.conf[fixed_res[ir].conf[ic]].occ; /* self - self was already set to 0 */
                }
            }
            conflist.conf[all_res[kr].conf[kc]].E_mfe = mfe;
            conflist.conf[all_res[kr].conf[kc]].E_self = mfe + conflist.conf[all_res[kr].conf[kc]].E_self0;
        }
    }

    /* base energy of fixed conformers */
    for (ir=0; ir<n_fixed; ir++) {
        for (ic=0; ic<fixed_res[ir].n; ic++) {
            E+=conflist.conf[fixed_res[ir].conf[ic]].E_self * conflist.conf[fixed_res[ir].conf[ic]].occ;
        }
    }
    /* substract double counted mfe on fixed conformers */
    for (ir=0; ir<n_fixed; ir++) {
        for (kr=ir+1; kr<n_fixed; kr++) {
            for (ic=0; ic<fixed_res[ir].n; ic++) {
                for (kc=0; kc<fixed_res[kr].n; kc++) {
                    E-= pairwise[fixed_res[ir].conf[ic]][fixed_res[kr].conf[kc]]
                    *conflist.conf[fixed_res[ir].conf[ic]].occ
                    *conflist.conf[fixed_res[kr].conf[kc]].occ;
                }
            }
        }
    }

    return E;
}


float get_E()
{  float E = 0.0;
    int kr, ir;

    /* Triangl calculation, no bias but expensive */
    for (kr=0; kr<n_free; kr++) E += conflist.conf[state[kr]].E_self;
    for (kr=0; kr<n_free; kr++) {
        for (ir=0; ir<kr; ir++)
            E += pairwise[state[kr]][state[ir]];
    }

    return E;
}

void mk_neighbors()
{  int kr, ir, kc, ic;
    char big;

    /* here */
    biglist = (BIGLIST *) malloc(n_free * sizeof(BIGLIST));

    for (kr=0; kr<n_free; kr++) {
        biglist[kr].n = 0;
        biglist[kr].res = NULL;
        for (ir=0; ir<n_free; ir++) {
            if (kr==ir) continue;
            big = 0;
            for (kc=0; kc<free_res[kr].n; kc++) {
                if (big) break;
                for (ic=0; ic<free_res[ir].n; ic++) {
                    if (big) break;
                    if (fabs(pairwise[free_res[kr].conf[kc]][free_res[ir].conf[ic]])>env.big_pairwise)
                        big = 1;
                }
            }
            if (big) {
                biglist[kr].n++;
                biglist[kr].res = realloc(biglist[kr].res, biglist[kr].n*sizeof(int));
                biglist[kr].res[biglist[kr].n-1] = ir;
            }
        }
    }

    return;
}


void MC(int n)
{  int cycles, n_total, n_cycle;
    int i, j, k;
    register int iters;
    int mem;
    int *old_state;
    float  old_E;
    float dE;
    float b;
    int nflips;
    double H_average;
    double H_noTS;
    float E_entropy;

    int iflip, ires, iconf; /* iconf is 0 to n of the conf in a res */
    int old_conf, new_conf; /* old_conf and new_conf are from 0 to n_conf in conflist */


    b = -KCAL2KT/(env.monte_temp/ROOMT);

    mem =n_free * sizeof(int);
    old_state = (int *) malloc(mem);
    E_minimum = E_state = get_E();

    /* number of cycles and iters in each cycle */
    if (env.monte_trace > 0) {
        cycles = (n-1)/env.monte_trace + 1; /* round up */
        n_total = cycles*env.monte_trace;
        n_cycle = env.monte_trace;
    }
    else {
        cycles = 1;
        n_total = n_cycle = n;
    }

    /* clear counters */
    for (i=0; i<conflist.n_conf; i++) conflist.conf[i].counter = 0;
    H_average = 0.0;
    H_noTS = 0.0;

    for (i=0; i<cycles; i++) {
        /*
        fprintf(fp, "Step %10d, E_minimum = %10.2f, E_running = %10.2f, E_reset = %10.2f\n",
        i*n_cycle, E_minimum+E_base, E_state+E_base, get_E()+E_base);
        */
        E_entropy = get_totalTS();
        fprintf(fp, "Step %10d, E_minimum= %10.2f, E_running = %10.2f ,  E_running without entropy = %10.2f\n",
        i*n_cycle, E_minimum+E_base, E_state+E_base, E_state+E_base - E_entropy);
        fflush(fp);
        iters = n_cycle;
        while (iters) {
            /*  save state */
            old_E = E_state;
            memcpy(old_state, state, mem);

            /* 1st flip */
            ires  = rand()/(RAND_MAX/n_free + 1);
            while (1) {
                iconf = rand()/(RAND_MAX/free_res[ires].n + 1);
                old_conf = state[ires];
                new_conf = free_res[ires].conf[iconf];
                if (old_conf != new_conf) break;
            }
            state[ires] = new_conf;
            E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
            for (j=0; j<n_free; j++) {
                E_state += pairwise[new_conf][state[j]] - pairwise[old_conf][state[j]];
            }

            /* now multiple flip */
            /*   1st flip -> No (50% probablity)
             *     | Yes (50% probablity)
             *   2nd flip (any res in big list)
             *     |
             *   3rd flip (any res in big list)
             *     |
             *   4th flip (any res in big list)
             */
            if (rand() & 1) {   /* do multiple flip if odd number */
                if (biglist[ires].n) {
                    nflips = env.monte_flips > (biglist[ires].n+1) ? biglist[ires].n+1: env.monte_flips;
                    for (k=1; k<nflips; k++) {

                        iflip = biglist[ires].res[rand()/(RAND_MAX/biglist[ires].n + 1)];
                        iconf = rand()/(RAND_MAX/free_res[iflip].n + 1);
                        old_conf = state[iflip];
                        new_conf = free_res[iflip].conf[iconf];

                        state[iflip] = new_conf;
                        E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
                        for (j=0; j<n_free; j++) {
                            E_state += pairwise[new_conf][state[j]] - pairwise[old_conf][state[j]];
                        }
                    }
                }
            }


            /* DEBUG
            for (j=0; j<n_free; j++) printf("%03d ", state[j]);
            printf("\n");
            */

            if (E_minimum > E_state) E_minimum = E_state;

            dE = E_state - old_E;
            if (dE < 0.0) {                                 /* go to new low */
            }
            /*<<< Boltman distribution >>>*/
            else if ((float) rand()/RAND_MAX < exp(b*dE)) { /* Go */
            }
            else {                                                    /* stay, restore the state */
                memcpy(state, old_state, mem);
                E_state = old_E;
            }

            /*
            if (!memcmp(state, state_check, sizeof(int)*n)) {
                printf("E=%12.5f\n", E_state);
            }  DEBUG */

            /* count this state energy */
            H_average += E_state/n_total;
            H_noTS += (E_state-get_totalTS())/n_total;
            /*<<< Do statistics >>>*/
            for (j=0; j<n_free; j++) {
                conflist.conf[state[j]].counter++;
                //+= pairwise_vdw[j]
            }

            iters --;
        }
    }

    E_entropy = get_totalTS();
    fprintf(fp, "Exit %10d, E_minimum = %10.2f, E_running = %10.2f E_running without TS = %10.2f\n", n_total, E_minimum+E_base, E_state+E_base, E_state+E_base-E_entropy);
    fprintf(fp, "The average running energy, corresponding to H = %8.3f kCal/mol, H without TS = %8.3f\n", H_average+E_base, H_noTS + E_base);
    fflush(fp);

    /* compute occ */
    for (i=0; i<conflist.n_conf; i++) {
        if (conflist.conf[i].on == 't') continue;
        conflist.conf[i].occ = (float) conflist.conf[i].counter / n_total;
    }

    free(old_state);

    return;
}

int reduce_conflist()
{  int i, j, t;
    int counter = 0;

    for (i=0; i<n_free; i++) {
        t = 0;
        for (j=0; j<free_res[i].n; j++) {
            if (conflist.conf[free_res[i].conf[j]].on == 'f'
            && conflist.conf[free_res[i].conf[j]].occ < env.monte_reduce) {
                conflist.conf[free_res[i].conf[j]].on = 't';
                conflist.conf[free_res[i].conf[j]].occ = 0.0;
                counter++;
                t++;
            }
        }
        if ((free_res[i].n-t) == 1) {  /* leave only one free conf, got to be 't' 1.00 */
            for (j=0; j<free_res[i].n; j++) {
                if (conflist.conf[free_res[i].conf[j]].on == 'f') {
                    conflist.conf[free_res[i].conf[j]].on = 't';
                    conflist.conf[free_res[i].conf[j]].occ = 1.0;
                    counter++;
                }
            }
        }
    }

    return counter;
}

float get_totalTS() {
    int kr, kc, ir, ic;
    float E = 0.0;

    /* TS from fixed conformers */
    for (ir=0; ir<n_fixed; ir++) {
        for (ic=0; ic<fixed_res[ir].n; ic++) {
            E+=conflist.conf[fixed_res[ir].conf[ic]].E_TS * conflist.conf[fixed_res[ir].conf[ic]].occ;
        }
    }

    /* TS from free conformers */
    for (kr=0; kr<n_free; kr++) E += conflist.conf[state[kr]].E_TS;

   return E;
}

int Nx;       /* number of titration points */
float *xp, *yp;   /* titration points */

struct STAT {
    float a;
    float b;
    float chi2;
};

struct STAT fit(float a, float b);
float score(float v[]);
void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk);

int fitit()
{  int i, j;
    int N_crg;
    int N_res;
    int Counter;
    char line[20];
    struct STAT stat;	/* a structure of statistcs */
    char **head, **mhead;  // add mhead for mfe header
    float *netcrg;
    float **ypp, **ysp;
    float a, b, mid;
    float *crg, *protons, *electrons;    /* total net charge at pHs */
    float n;
    char sbuff[21];
    int n_protons, n_electrons, n_crg;
    int n_protons_grnd, n_electrons_grnd, n_crg_grnd;
    int k; // just for cycle
    int tifound; // whether do fitit
    FILE *blist_fp;
    float mfePoint;

    /*<<< Get titration points >>>*/
    /*--- Determine the titration points, x values ---*/
    Nx = env.titr_steps;

    /*--- Assign titration points to an array ---*/
    xp = (float *) malloc(Nx * sizeof(float));		/* store x points */
    yp = (float *) malloc(Nx * sizeof(float));		/* store y points */
    ypp = (float **) malloc(conflist.n_conf * sizeof(float *));
    for (i=0; i<conflist.n_conf; i++) ypp[i] = (float *) malloc(Nx * sizeof(float));
    ysp = (float **) malloc(conflist.n_conf * sizeof(float *));
    for (i=0; i<conflist.n_conf; i++) ysp[i] = (float *) malloc(Nx * sizeof(float));

    /*--- Convert x points to float ---*/
    for (i=0; i<Nx; i++) {
        if (env.titr_type == 'p') xp[i] = env.titr_ph0 + i*env.titr_phd;
        else  xp[i] = env.titr_eh0 + i*env.titr_ehd;
    }

    head = (char **) malloc(conflist.n_conf * sizeof(char *));
    for (i=0; i<conflist.n_conf; i++) head[i] = (char *) malloc(20 * sizeof(char));
    shead = (char **) malloc(conflist.n_conf * sizeof(char *));
    for (i=0; i<conflist.n_conf; i++) shead[i] = (char *) malloc(20 * sizeof(char));
    mhead = (char **) malloc(conflist.n_conf * sizeof(char *));
    for (i=0; i<conflist.n_conf; i++) mhead[i] = (char *) malloc(20 * sizeof(char));
    netcrg = (float *) malloc(conflist.n_conf*sizeof(float));

    /*<<< Group into residues >>>*/
    /* keep only charged conformers */
    Counter = 0;
    for (i=0; i <conflist.n_conf; i++) {
        if (strchr(conflist.conf[i].uniqID, '+') || strchr(conflist.conf[i].uniqID, '-')) {
            strncpy(head[Counter], conflist.conf[i].uniqID, 4); head[Counter][4] = '\0';
            strncat(head[Counter], conflist.conf[i].uniqID+5, 6);
            head[Counter][10] = '\0';
            for (j=0; j<Nx; j++) {
                ypp[Counter][j] = occ_table[i][j];
            }
            Counter++;
        }
    }
    N_crg = Counter;


    /* group residues */
    Counter = 0;
    strncpy(line, head[0], 10);
    for (j=0; j<Nx; j++) ysp[0][j] = ypp[0][j];
    for (i=1; i<N_crg; i++) {
        if (strncmp(head[i], line, 10)) {      /* not equal, a new residue */
            strncpy(shead[Counter], line, 10); shead[Counter][10] = '\0';
            Counter++;
            strncpy(line, head[i], 10);
            for (j=0; j<Nx; j++) ysp[Counter][j] = ypp[i][j];
        }
        else {                                 /* same residue */
            for (j=0; j<Nx; j++) ysp[Counter][j] += ypp[i][j];
        }
    }
    strncpy(shead[Counter], line, 10); shead[Counter][10] = '\0';
    N_res = Counter + 1;

    n_mfe = N_res;
    for (i=0; i<n_all; i++) {
        for (j=0; j<Nx; j++) all_res[i].mfe.crg[j] = 0.0;
    }
    // N_res are the number we do mfe
    for (i=0; i<n_mfe; i++) {
        strncpy(mhead[i], shead[i], 3); mhead[i][3] = '\0';
        strncat(mhead[i], shead[i]+4, 6); mhead[i][9] = '\0';
    }
    if (env.mfe_flag) {
        if (env.titr_type == 'p') printf("       Doing mfe at pH %.3f for all the residues\n", env.mfe_point);
        else printf("       Doing mfe at Eh %.3f for all the residues\n", env.mfe_point);
        if ((env.mfe_point>xp[0] && env.mfe_point>xp[Nx-1]) || (env.mfe_point<xp[0] && env.mfe_point<xp[Nx-1]))
            printf("         MFE: mfe point not in the titration range, do mfe at pKa or Em\n");
    }
    else printf("   MFE: didn't specify mfe point, do mfe at pKa or Em \n");

    /* free pairwise */
    for (i=0; i<conflist.n_conf; i++)
        free(pairwise[i]);
    free(pairwise);
   // load the pairwise interaction again, round the ele and vdw to keep consistent with mfe.py
    if (load_pairwise_fround3()) {
        printf("   FATAL: mfe pairwise interaction not loaded\n");
        return USERERR;
    }

    // load all the conformers of mfe residues
    mfe_res = (RESIDUE *) calloc(n_mfe, sizeof(RESIDUE));
    for (k=0; k<n_mfe; k++) {
        mfe_res[k].n = 0;
        for (i=0; i<conflist.n_conf; i++) {
            strncpy(sbuff, conflist.conf[i].uniqID, 3); sbuff[3] = '\0';
            strncat(sbuff, conflist.conf[i].uniqID+5, 6); sbuff[9] = '\0';
            if (!strcmp(mhead[k], sbuff)) mfe_res[k].n++;
        }
    }
    for (i=0; i<n_mfe; i++) {
        mfe_res[i].conf = (int *) calloc(mfe_res[i].n, sizeof(int));
        mfe_res[i].mfe.mfePair = (float *) calloc(n_all, sizeof(float));
        mfe_res[i].mfe.vdw = (float *) calloc(n_all, sizeof(float));
        mfe_res[i].mfe.ele = (float *) calloc(n_all, sizeof(float));
        mfe_res[i].mfe.crg = (float *) calloc(env.titr_steps, sizeof(float));
    }

    for (i=0; i<n_mfe; i++) {
        mfe_res[i].n = 0;
        for (j=0; j<conflist.n_conf; j++) {
            strncpy(sbuff, conflist.conf[j].uniqID, 3); sbuff[3] = '\0';
            strncat(sbuff, conflist.conf[j].uniqID+5, 6); sbuff[9] = '\0';
            if (!strcmp(mhead[i], sbuff)) {
                mfe_res[i].n++;
                mfe_res[i].conf[mfe_res[i].n-1] = conflist.conf[j].iConf;
            }
        }
    }
    for (i=0; i<n_mfe; i++) {
        for (j=0; j<mfe_res[i].n; j++) {
            conflist.conf[mfe_res[i].conf[j]].E_self0 = conflist.conf[mfe_res[i].conf[j]].E_vdw0
                                                      + conflist.conf[mfe_res[i].conf[j]].E_vdw1
                                                      + conflist.conf[mfe_res[i].conf[j]].E_epol
                                                      + conflist.conf[mfe_res[i].conf[j]].E_tors
                                                      + conflist.conf[mfe_res[i].conf[j]].E_dsolv
                                                      + conflist.conf[mfe_res[i].conf[j]].E_extra;
        }
    }

    /* Write net charge ,protons, electrons to sum_crg.out, to be implemented */
    crg       = (float *) malloc(Nx * sizeof(float));
    protons   = (float *) malloc(Nx * sizeof(float));
    electrons = (float *) malloc(Nx * sizeof(float));

    memset(crg, 0, Nx * sizeof(float));
    memset(protons, 0, Nx * sizeof(float));
    memset(electrons, 0, Nx * sizeof(float));

    if ((fp = fopen("sum_crg.out", "w")) == NULL) {
        printf("   Can not open \"sum_crg.out\" to write, abort ...\n");
        return USERERR;
    }
    if (env.titr_type == 'p') {   /* pH titration */
        fprintf(fp, "  pH      ");
    }
    else {      /* Eh titration assumed */
        fprintf(fp, "  Eh      ");
    }
    for (i=0; i<Nx; i++) fprintf(fp, " %5d", (int) xp[i]);
    fprintf(fp, "\n");

    for(i=0; i<N_res; i++) {
        fprintf(fp, "%s", shead[i]);
        strncpy(sbuff, shead[i], 4); sbuff[4] = '1'; sbuff[5] = '\0';
        if (param_get( "PROTON", sbuff, "", &n_protons)) n_protons = 0;
        if (param_get( "ELECTRON", sbuff, "", &n_electrons)) n_electrons = 0;
        /* n_crg = n_protons-n_electrons; */
        n_crg = n_protons-n_electrons;

        /* number of protons on ground conformer type */
        strncpy(sbuff, shead[i], 3); sbuff[3] = '0'; sbuff[4] = '1'; sbuff[5] = '\0';
        if (param_get( "PROTON", sbuff, "", &n_protons_grnd)) n_protons_grnd = 0;
        if (param_get( "ELECTRON", sbuff, "", &n_electrons_grnd)) n_electrons_grnd = 0;
        n_crg_grnd = n_protons_grnd - n_electrons_grnd;

        for(j=0; j<Nx; j++) {
            fprintf(fp, " %5.2f", n_crg*ysp[i][j]);
            mfe_res[i].mfe.crg[j] = n_crg*ysp[i][j] + n_crg_grnd*(1.0-ysp[i][j]);
            for (k=0; k<n_all; k++) {
                strncpy(sbuff, conflist.conf[all_res[k].conf[0]].uniqID, 3); sbuff[3] = '\0';
                strncat(sbuff, conflist.conf[all_res[k].conf[0]].uniqID+5, 6); sbuff[9] = '\0';
                // get all the charges of all the residues
                if (!strcmp(sbuff, mhead[i])) all_res[k].mfe.crg[j] = n_crg*ysp[i][j] + n_crg_grnd*(1.0-ysp[i][j]);
            }
            crg[j] += n_crg*ysp[i][j] + n_crg_grnd*(1.0-ysp[i][j]);
            protons[j] += n_protons*ysp[i][j] + n_protons_grnd*(1.0-ysp[i][j]);
            electrons[j] += n_electrons*ysp[i][j] + n_electrons_grnd*(1.0-ysp[i][j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "----------\n");

    fprintf(fp, "Net_Charge");
    for(j=0; j<Nx; j++) {
        fprintf(fp, " %5.2f", crg[j]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "Protons   ");
    for(j=0; j<Nx; j++) {
        fprintf(fp, " %5.2f", protons[j]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "Electrons ");
    for(j=0; j<Nx; j++) {
        fprintf(fp, " %5.2f", electrons[j]);
    }
    fprintf(fp, "\n");

    fclose(fp);


    /*<<< Loop over y values >>>*/
    if (!(fp = fopen(CURVE_FITTING, "w"))) {
        printf("   FATAL: can not write to file \"%s\".", CURVE_FITTING);
        return USERERR;
    }
    if (env.titr_type == 'p') {   /* pH titration */
        fprintf(fp, "  pH      ");
    }
    else {      /* Eh titration assumed */
        fprintf(fp, "  Eh      ");
    }
    fprintf(fp, "       pKa/Em  n(slope) 1000*chi2      vdw0    vdw1    tors    ebkb    dsol   offset  pHpK0   EhEm0    -TS   residues   total\n");

    if (!(blist_fp = fopen("respair.lst", "w"))) {
        printf("can't write file respair.lst\n");
    }
    fprintf(blist_fp, " residue    partner         vdw     ele  pairwise  charge\n");



    for (i=0; i <N_res; i++) {
        /*--- Convert y points to float ---*/

        /* a reasonable guess makes optimization easier */
        a = 0.0;
        mid = 0.6;
        for (j=0; j<Nx; j++) {
            yp[j] = ysp[i][j];
            if (fabs(yp[j] - 0.5) < mid) {
                mid = fabs(yp[j] - 0.5);
                b = xp[j];
            }
        }
        tifound=1;
        if (mid >= 0.485) {
            tifound=0;
            if (fabs(yp[0]-yp[Nx-1])>0.5) {  /* jumps form <0.015 to >0.985 */
               fprintf(fp, "%s        %-25s", shead[i],"titration curve too sharp");
               mfePoint=xp[Nx/2];
            }
            /*else {                           /* all < 0.015 or all >0.985 */
               //fprintf(fp, "%s          pKa or Em out of range   \n", shead[i]);
	   /* 	if (yp[0] > yp[Nx-1]) {
                    fprintf(fp, "%s        pKa/Em more than %-8.1f", shead[i], xp[Nx-1]);
                    i_low=i_high=Nx-1;
                }
                else if (yp[0] < yp[Nx-1]) {
                    fprintf(fp, "%s        pKa/Em less than %-8.1f", shead[i], xp[0]);
                    i_low=i_high=0;
                }
                else {
                    strncpy(sbuff, shead[i], 3); sbuff[3] = '\0';
                    if (!strcmp(sbuff, "TYR") || !strcmp(sbuff, "ASP") || !strcmp(sbuff, "GLU") || !strcmp(sbuff, "CTR")) {
                        if (yp[0] == 1.0) {
                            fprintf(fp, "%s        pKa/Em less than %-8.1f", shead[i], xp[0]);
                            i_low=i_high=0;
                        }
                        else {
                            fprintf(fp, "%s        pKa/Em more than %-8.1f", shead[i], xp[Nx-1]);
                            i_low=i_high=Nx-1;
                        }
                    }
                    else {
                        if (yp[0] == 1.0) {
                            fprintf(fp, "%s        pKa/Em more than %-8.1f", shead[i], xp[Nx-1]);
                            i_low=i_high=Nx-1;
                        }
                        else {
                             fprintf(fp, "%s        pKa/Em less than %-8.1f", shead[i], xp[0]);
                             i_low=i_high=0;
                        }
                    }
	        }
	    }*/

            else {    // assume either all the yp >0.985 or all of them <0.015
                if (yp[0] > 0.985) {
                    if (strchr(shead[i], '+')) {
                        fprintf(fp, "%s        >%-24.1f", shead[i], xp[Nx-1]);
                        mfePoint=xp[Nx-1];
                    }
                    else {
                        fprintf(fp, "%s        <%-24.1f", shead[i], xp[0]);
                        mfePoint=xp[0];
                    }
                }
                else {
                    if (strchr(shead[i], '+')) {
                        fprintf(fp, "%s        <%-24.1f", shead[i], xp[0]);
                        mfePoint=xp[0];
                    }
                    else {
                        fprintf(fp, "%s        >%-24.1f", shead[i], xp[Nx-1]);
                        mfePoint=xp[Nx-1];
                    }
                }
            }
        }
        if (tifound == 0) {
            print_mfe(i, mfePoint, fp, blist_fp);
            continue;
        }
         //   fprintf(fp, "
        /* printf("a=%.3f; b=%.3f\n",a,b); */
        stat = fit(a, b);

        if (env.titr_type == 'p') n = fabs(stat.a * 8.617342E-2 * env.monte_temp / 58.0);
        else if (env.titr_type == 'e') n = fabs(stat.a * 8.617342E-2 * env.monte_temp);
        else n = stat.a;


        //if (stat.b < xp[0] || stat.b > xp[Nx-1]) fprintf(fp, "%s        pKa or Em out of range   \n", shead[i]);
        if (stat.b < xp[0]) {
            fprintf(fp, "%s        <%-24.1f", shead[i], xp[0]);
            mfePoint = xp[0];
        }
 	else if (stat.b > xp[Nx-1]) {
            fprintf(fp, "%s        >%-24.1f", shead[i], xp[Nx-1]);
            mfePoint = xp[Nx-1];
        }
        else {
            fprintf(fp, "%s    %9.3f %9.3f %9.3f", shead[i], stat.b, n, 1000*stat.chi2);
            mfePoint = stat.b;
        }
        print_mfe(i, mfePoint, fp, blist_fp);
    }

    free(crg);
    free(protons);
    free(electrons);
    free(xp);
    free(yp);
    for (i=0; i<conflist.n_conf; i++) free(ypp[i]);
    free(ypp);
    for (i=0; i<conflist.n_conf; i++) free(ysp[i]);
    free(ysp);
    for (i=0; i<conflist.n_conf; i++) free(head[i]);
    free(head);
    for (i=0; i<conflist.n_conf; i++) free(shead[i]);
    free(shead);

    return 0;
}

float fround3(float x)
{
    // round off the float point number
    char sbuff[128];
    float y;

    sprintf(sbuff, "%.3f", x);
    y=atof(sbuff);
//    if (x>0.0000001) return (int)(x*1000+0.5)/1000.0;
//    if (x<-0.0000001) return -(int)(-x*1000+0.5)/1000.0;
//    else return 0.000;
    return y;

}

int load_pairwise_fround3()
{   int i, j, kc;
    EMATRIX ematrix;

    ematrix.n = 0;
    if (load_energies(&ematrix, ".")<0) {
        printf("   File %s not found\n", ENERGY_TABLE);
        return USERERR;
    }

    if (!(pairwise = (float **) malloc(ematrix.n * sizeof(float *)))) {
        printf("   FATAL: memory error in make_matrices()\n");
        return USERERR;
    }
    pairwise_vdw = (float **) malloc(ematrix.n * sizeof(float *));
    pairwise_ele = (float **) malloc(ematrix.n * sizeof(float *));

    for (kc=0; kc<ematrix.n; kc++) {
        if (!(pairwise[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
        if (!(pairwise_vdw[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
        if (!(pairwise_ele[kc] = (float *) malloc(ematrix.n * sizeof(float)))) {
            printf("   FATAL: memory error in make_matrices()\n");
            return USERERR;
        }
    }

    for (i=0; i<ematrix.n; i++) {
       for (j=0; j<ematrix.n; j++) {
          // pairwise[i][j] = pairwise[j][i] = ((ematrix.pw[i][j].ele + ematrix.pw[j][i].ele)*env.scale_ele
          //                                 +(ematrix.pw[i][j].vdw + ematrix.pw[j][i].vdw)*env.scale_vdw)/2.0 ;
          // if (ematrix.pw[i][j].vdw > 998.0) pairwise[i][j] = 999.0;
            pairwise[i][j] = fround3(ematrix.pw[i][j].ele) * env.scale_ele + fround3(ematrix.pw[i][j].vdw) * env.scale_vdw;
            pairwise_vdw[i][j] = fround3(ematrix.pw[i][j].vdw) * env.scale_vdw;
            pairwise_ele[i][j] = fround3(ematrix.pw[i][j].ele) * env.scale_ele;
        }
    }

    /* free memory */
    free_ematrix(&ematrix);

    return 0;
}

int print_mfe(int i_res, float mfeP, FILE *pK_fp, FILE *res_fp)
{
    int i_low, i_high;
    int j, k;
    float cut_off, ratio;
    char head1[20], head2[20];
    RESIDUE old_mfe;

    if (env.mfe_flag) {
        if (! ((env.mfe_point>xp[0] && env.mfe_point>xp[Nx-1]) || (env.mfe_point<xp[0] && env.mfe_point<xp[Nx-1]))) mfeP = env.mfe_point;
    }

    cut_off = env.mfe_cutoff;
    // if (cut_off < 0.0000) cut_off = 0.5;
    // xuyu commented it, 8/11

    for (k=0; k<=Nx-1; k++) {
        if (xp[k] >= mfeP) {
            i_high=k; break;
        }
    }
    for (k=Nx-1; k>=0; k--) {
        if (xp[k] <= mfeP) {
            i_low=k; break;
        }
    }

    if (abs(i_high-i_low) == 0) { // do mfe at only one point
        get_mfe(i_res, i_low);
        if (env.titr_type == 'p')
            fprintf(pK_fp, "  %8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f\n",
                   mfe_res[i_res].mfe.vdw0/PH2KCAL,
                   mfe_res[i_res].mfe.vdw1/PH2KCAL,
                   mfe_res[i_res].mfe.tors/PH2KCAL,
                   mfe_res[i_res].mfe.ebkb/PH2KCAL,
                   mfe_res[i_res].mfe.dsol/PH2KCAL,
                   mfe_res[i_res].mfe.offset/PH2KCAL,
                   mfe_res[i_res].mfe.pHpK0/PH2KCAL,
                   mfe_res[i_res].mfe.EhEm0/PH2KCAL,
                   mfe_res[i_res].mfe.TS/PH2KCAL,
                   mfe_res[i_res].mfe.residues/PH2KCAL,
                   mfe_res[i_res].mfe.total/PH2KCAL);
        else
            fprintf(pK_fp, "  %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%10.1f\n",
                   mfe_res[i_res].mfe.vdw0/mev2Kcal,
                   mfe_res[i_res].mfe.vdw1/mev2Kcal,
                   mfe_res[i_res].mfe.tors/mev2Kcal,
                   mfe_res[i_res].mfe.ebkb/mev2Kcal,
                   mfe_res[i_res].mfe.dsol/mev2Kcal,
                   mfe_res[i_res].mfe.offset/mev2Kcal,
                   mfe_res[i_res].mfe.pHpK0/mev2Kcal,
                   mfe_res[i_res].mfe.EhEm0/mev2Kcal,
                   mfe_res[i_res].mfe.TS/mev2Kcal,
                   mfe_res[i_res].mfe.residues/mev2Kcal,
                   mfe_res[i_res].mfe.total/mev2Kcal);

        for (j=0; j<n_all; j++) {
            head1[0]='\0'; head2[0]='\0';
            if (fabs(mfe_res[i_res].mfe.mfePair[j]) < cut_off) continue;
            strcpy(head1, shead[i_res]);
            strncpy(head2, conflist.conf[all_res[j].conf[0]].uniqID, 3); head2[3] = '\0';
            strncat(head2, conflist.conf[all_res[j].conf[0]].uniqID+5, 6); head2[9] = '\0';
            if (env.titr_type == 'p')
                fprintf(res_fp, "%s  %s   %8.2f%8.2f%8.2f%8.2f\n", head1, head2,
                mfe_res[i_res].mfe.vdw[j]/PH2KCAL, mfe_res[i_res].mfe.ele[j]/PH2KCAL, mfe_res[i_res].mfe.mfePair[j]/PH2KCAL, all_res[j].mfe.crg[i_low]);
            else
                fprintf(res_fp, "%s  %s   %8.2f%8.2f%8.2f%8.2f\n", head1, head2,
                mfe_res[i_res].mfe.vdw[j]/mev2Kcal, mfe_res[i_res].mfe.ele[j]/mev2Kcal, mfe_res[i_res].mfe.mfePair[j]/mev2Kcal, all_res[j].mfe.crg[i_low]);
        }
    }
    else {    /* get mfe from two titration points */
        ratio = (mfeP - xp[i_low])/(xp[i_high]-xp[i_low]);
        get_mfe(i_res, i_low);
        old_mfe.n = mfe_res[i_res].n;
        old_mfe.conf = (int *) calloc(old_mfe.n, sizeof(int));
        old_mfe.mfe.mfePair = (float *) calloc(n_all, sizeof(float));
        old_mfe.mfe.vdw = (float *) calloc(n_all, sizeof(float));
        old_mfe.mfe.ele = (float *) calloc(n_all, sizeof(float));
        old_mfe.mfe.vdw0 = mfe_res[i_res].mfe.vdw0;
        old_mfe.mfe.vdw1 = mfe_res[i_res].mfe.vdw1;
        old_mfe.mfe.ebkb = mfe_res[i_res].mfe.ebkb;
        old_mfe.mfe.tors = mfe_res[i_res].mfe.tors;
        old_mfe.mfe.dsol = mfe_res[i_res].mfe.dsol;
        old_mfe.mfe.offset = mfe_res[i_res].mfe.offset;
        old_mfe.mfe.pHpK0 = mfe_res[i_res].mfe.pHpK0;
        old_mfe.mfe.EhEm0 = mfe_res[i_res].mfe.EhEm0;
        old_mfe.mfe.TS = mfe_res[i_res].mfe.TS;
        old_mfe.mfe.residues = mfe_res[i_res].mfe.residues;
        old_mfe.mfe.total = mfe_res[i_res].mfe.total;

        for (k=0; k<n_all; k++) {
            old_mfe.mfe.mfePair[k] = mfe_res[i_res].mfe.mfePair[k];
            old_mfe.mfe.vdw[k] = mfe_res[i_res].mfe.vdw[k];
            old_mfe.mfe.ele[k] = mfe_res[i_res].mfe.ele[k];
        }

        get_mfe(i_res, i_high);
        mfe_res[i_res].mfe.vdw0 = old_mfe.mfe.vdw0*(1-ratio) + mfe_res[i_res].mfe.vdw0*ratio;
        mfe_res[i_res].mfe.vdw1 = old_mfe.mfe.vdw1*(1-ratio) + mfe_res[i_res].mfe.vdw1*ratio;
        mfe_res[i_res].mfe.ebkb = old_mfe.mfe.ebkb*(1-ratio) + mfe_res[i_res].mfe.ebkb*ratio;
        mfe_res[i_res].mfe.tors = old_mfe.mfe.tors*(1-ratio) + mfe_res[i_res].mfe.tors*ratio;
        mfe_res[i_res].mfe.dsol = old_mfe.mfe.dsol*(1-ratio) + mfe_res[i_res].mfe.dsol*ratio;
        mfe_res[i_res].mfe.offset = old_mfe.mfe.offset*(1-ratio) + mfe_res[i_res].mfe.offset*ratio;
        mfe_res[i_res].mfe.pHpK0 = old_mfe.mfe.pHpK0*(1-ratio) + mfe_res[i_res].mfe.pHpK0*ratio;
        mfe_res[i_res].mfe.EhEm0 = old_mfe.mfe.EhEm0*(1-ratio) + mfe_res[i_res].mfe.EhEm0*ratio;
        mfe_res[i_res].mfe.TS = old_mfe.mfe.TS*(1-ratio) + mfe_res[i_res].mfe.TS*ratio;
        mfe_res[i_res].mfe.residues = old_mfe.mfe.residues*(1-ratio) + mfe_res[i_res].mfe.residues*ratio;
        mfe_res[i_res].mfe.total = old_mfe.mfe.total*(1-ratio) + mfe_res[i_res].mfe.total*ratio;

        for (k=0; k<n_all; k++) {
            mfe_res[i_res].mfe.mfePair[k] = old_mfe.mfe.mfePair[k]*(1-ratio) + mfe_res[i_res].mfe.mfePair[k]*ratio;
            mfe_res[i_res].mfe.vdw[k] = old_mfe.mfe.vdw[k]*(1-ratio) + mfe_res[i_res].mfe.vdw[k]*ratio;
            mfe_res[i_res].mfe.ele[k] = old_mfe.mfe.ele[k]*(1-ratio) + mfe_res[i_res].mfe.ele[k]*ratio;
        }

        if (env.titr_type == 'p')
            fprintf(pK_fp, "  %8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f\n",
                   mfe_res[i_res].mfe.vdw0/PH2KCAL,
                   mfe_res[i_res].mfe.vdw1/PH2KCAL,
                   mfe_res[i_res].mfe.tors/PH2KCAL,
                   mfe_res[i_res].mfe.ebkb/PH2KCAL,
                   mfe_res[i_res].mfe.dsol/PH2KCAL,
                   mfe_res[i_res].mfe.offset/PH2KCAL,
                   mfe_res[i_res].mfe.pHpK0/PH2KCAL,
                   mfe_res[i_res].mfe.EhEm0/PH2KCAL,
                   mfe_res[i_res].mfe.TS/PH2KCAL,
                   mfe_res[i_res].mfe.residues/PH2KCAL,
                   mfe_res[i_res].mfe.total/PH2KCAL);
        else
            fprintf(pK_fp, "  %8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%8.1f%10.1f\n",
                   mfe_res[i_res].mfe.vdw0/mev2Kcal,
                   mfe_res[i_res].mfe.vdw1/mev2Kcal,
                   mfe_res[i_res].mfe.tors/mev2Kcal,
                   mfe_res[i_res].mfe.ebkb/mev2Kcal,
                   mfe_res[i_res].mfe.dsol/mev2Kcal,
                   mfe_res[i_res].mfe.offset/mev2Kcal,
                   mfe_res[i_res].mfe.pHpK0/mev2Kcal,
                   mfe_res[i_res].mfe.EhEm0/mev2Kcal,
                   mfe_res[i_res].mfe.TS/mev2Kcal,
                   mfe_res[i_res].mfe.residues/mev2Kcal,
                   mfe_res[i_res].mfe.total/mev2Kcal);

        for (j=0; j<n_all; j++) {
            head1[0]='\0'; head2[0]='\0';
            if (fabs(mfe_res[i_res].mfe.mfePair[j]) < cut_off) continue;
            strcpy(head1, shead[i_res]);
            strncpy(head2, conflist.conf[all_res[j].conf[0]].uniqID, 3); head2[3] = '\0';
            strncat(head2, conflist.conf[all_res[j].conf[0]].uniqID+5, 6); head2[9] = '\0';
            if (env.titr_type == 'p')
                fprintf(res_fp, "%s  %s   %8.2f%8.2f%8.2f%8.2f\n", head1, head2,
                       mfe_res[i_res].mfe.vdw[j]/PH2KCAL, mfe_res[i_res].mfe.ele[j]/PH2KCAL, mfe_res[i_res].mfe.mfePair[j]/PH2KCAL,
                       (all_res[j].mfe.crg[i_low]*(1-ratio) + all_res[j].mfe.crg[i_high]*ratio));
            else
                fprintf(res_fp, "%s  %s   %8.1f%8.1f%8.1f%8.2f\n", head1, head2,
                       mfe_res[i_res].mfe.vdw[j]/mev2Kcal, mfe_res[i_res].mfe.ele[j]/mev2Kcal, mfe_res[i_res].mfe.mfePair[j]/mev2Kcal,
                       (all_res[j].mfe.crg[i_low]*(1-ratio) + all_res[j].mfe.crg[i_high]*ratio));
        }
     }

     return 0;
}

int get_mfe(int i_res, int t_point)
{
    // get the mfe of the ith mfe_res at titration point t_point
    int j, k, q;
    float ph, eh;
    double *nocc, *rocc, socc[2]={0.0, 0.0};
    double Eref[2];
    MFE E_ground, E_ionize;
    float *mfe;

    memset(&E_ground, 0, sizeof(MFE));
    memset(&E_ionize, 0, sizeof(MFE));
    memset(&socc, 0, 2*sizeof(float));
    E_ground.mfePair = (float *) calloc(n_all, sizeof(float));
    E_ground.vdw = (float *) calloc(n_all, sizeof(float));
    E_ground.ele = (float *) calloc(n_all, sizeof(float));
    E_ionize.mfePair = (float *) calloc(n_all, sizeof(float));
    E_ionize.vdw = (float *) calloc(n_all, sizeof(float));
    E_ionize.ele = (float *) calloc(n_all, sizeof(float));
    mfe = (float *) calloc(mfe_res[i_res].n, sizeof(float));
    nocc = (double *) calloc(mfe_res[i_res].n, sizeof(double));
    rocc = (double *) calloc(mfe_res[i_res].n, sizeof(double));

    ph = env.titr_ph0;
    eh = env.titr_eh0;
    if (env.titr_type == 'p') ph = env.titr_ph0 + t_point*env.titr_phd;
    else eh = env.titr_eh0 + t_point*env.titr_ehd;

    for (j=0; j<mfe_res[i_res].n; j++) {
        if (!(strncmp(conflist.conf[mfe_res[i_res].conf[j]].uniqID+3, "DM", 2))) {  // dummy conformer
            conflist.conf[mfe_res[i_res].conf[j]].E_self = conflist.conf[mfe_res[i_res].conf[j]].E_self0;
            continue;
        }

        conflist.conf[mfe_res[i_res].conf[j]].E_ph =  conflist.conf[mfe_res[i_res].conf[j]].H * (ph-conflist.conf[mfe_res[i_res].conf[j]].pKa) * PH2KCAL;
        conflist.conf[mfe_res[i_res].conf[j]].E_eh =  conflist.conf[mfe_res[i_res].conf[j]].e * (eh-conflist.conf[mfe_res[i_res].conf[j]].Em) * mev2Kcal;

        for (k=0; k<conflist.n_conf; k++)
                mfe[j] += pairwise[mfe_res[i_res].conf[j]][conflist.conf[k].iConf] * fround3(occ_table[k][t_point]);

        conflist.conf[mfe_res[i_res].conf[j]].E_self = conflist.conf[mfe_res[i_res].conf[j]].E_self0
                                                     + conflist.conf[mfe_res[i_res].conf[j]].E_ph
                                                     + conflist.conf[mfe_res[i_res].conf[j]].E_eh
                                                     + mfe[j];
     }

    /* initial Eref */
    // '0' and 'DM' go to neutral form
    for (j=0; j<mfe_res[i_res].n; j++) {
        if (conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == '0' || conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == 'D') {
            Eref[0] = conflist.conf[mfe_res[i_res].conf[j]].E_self;
            break;
        }
    }
    for (j=0; j<mfe_res[i_res].n; j++) {
        if (conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] != '0' && conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] != 'D') {
            Eref[1] = conflist.conf[mfe_res[i_res].conf[j]].E_self;
            break;
        }
    }

      /* get Eref for ground and ionized state */
      /* if there is only one kind of conformer, this may be a problem */
    for (j=0; j<mfe_res[i_res].n; j++) {
        if (conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == '0' || conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == 'D') {
            if (Eref[0] > conflist.conf[mfe_res[i_res].conf[j]].E_self) Eref[0] = conflist.conf[mfe_res[i_res].conf[j]].E_self;
        }
        else {
            if (Eref[1] > conflist.conf[mfe_res[i_res].conf[j]].E_self) Eref[1] = conflist.conf[mfe_res[i_res].conf[j]].E_self;
        }
    }

    for (j=0; j<mfe_res[i_res].n; j++) {
        if (conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == '0' || conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == 'D') {
            rocc[j] = exp(-(conflist.conf[mfe_res[i_res].conf[j]].E_self-Eref[0]) * KCAL2KT);
        }
        else {
            rocc[j] = exp(-(conflist.conf[mfe_res[i_res].conf[j]].E_self-Eref[1]) * KCAL2KT);
        }
    }

    for (j=0; j<mfe_res[i_res].n; j++) {
        if (conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == '0' || conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == 'D') socc[0] += rocc[j];
        else socc[1] += rocc[j];
    }

    for (j=0; j<mfe_res[i_res].n; j++) {
        if (conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == '0' || conflist.conf[mfe_res[i_res].conf[j]].uniqID[3] == 'D') {
            nocc[j] = rocc[j]/socc[0];
        }
        else {
            nocc[j] = rocc[j]/socc[1];
        }
    }


    for (j=0; j<mfe_res[i_res].n; j++) {
        if (strchr(conflist.conf[mfe_res[i_res].conf[j]].uniqID, '+') || strchr(conflist.conf[mfe_res[i_res].conf[j]].uniqID, '-')) {
            E_ionize.vdw0 += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_vdw0;
            E_ionize.vdw1 += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_vdw1;
            E_ionize.ebkb += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_epol;
            E_ionize.tors += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_tors;
            E_ionize.dsol += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_dsolv;
            E_ionize.offset += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_extra;
            E_ionize.pHpK0 += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_ph;
            E_ionize.EhEm0 += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_eh;
            if (nocc[j] > 0.000001) E_ionize.TS += nocc[j] * log(nocc[j])/KCAL2KT;
            E_ionize.residues += nocc[j] * mfe[j];

            E_ionize.total = E_ionize.vdw0 + E_ionize.vdw1 + E_ionize.ebkb + E_ionize.tors + E_ionize.dsol + E_ionize.offset + E_ionize.pHpK0 + E_ionize.EhEm0 + E_ionize.TS + E_ionize.residues;
            for (k=0; k<n_all; k++) {
                for (q=0; q<all_res[k].n; q++) {
                    E_ionize.mfePair[k] += nocc[j] * (pairwise[mfe_res[i_res].conf[j]][all_res[k].conf[q]] * fround3(occ_table[all_res[k].conf[q]][t_point])) ;
                    E_ionize.vdw[k] += nocc[j] * (pairwise_vdw[mfe_res[i_res].conf[j]][all_res[k].conf[q]] * fround3(occ_table[all_res[k].conf[q]][t_point])) ;
                    E_ionize.ele[k] += nocc[j] * (pairwise_ele[mfe_res[i_res].conf[j]][all_res[k].conf[q]] * fround3(occ_table[all_res[k].conf[q]][t_point])) ;
                }
            }
        }
        else {
            E_ground.vdw0 += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_vdw0;
            E_ground.vdw1 += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_vdw1;
            E_ground.ebkb += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_epol;
            E_ground.tors += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_tors;
            E_ground.dsol += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_dsolv;
            E_ground.offset += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_extra;
            E_ground.pHpK0 += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_ph;
            E_ground.EhEm0 += nocc[j] * conflist.conf[mfe_res[i_res].conf[j]].E_eh;
            if (nocc[j] > 0.000001) E_ground.TS += nocc[j] * log(nocc[j])/KCAL2KT;
            E_ground.residues += nocc[j] * mfe[j];
            E_ground.total = E_ground.vdw0 + E_ground.vdw1 + E_ground.ebkb + E_ground.tors + E_ground.dsol + E_ground.offset + E_ground.pHpK0 + E_ground.EhEm0 + E_ground.TS + E_ground.residues;
            for (k=0; k<n_all; k++) {
                for (q=0; q<all_res[k].n; q++) {
                   E_ground.mfePair[k] += nocc[j] * pairwise[mfe_res[i_res].conf[j]][all_res[k].conf[q]] * fround3(occ_table[all_res[k].conf[q]][t_point]);
                   E_ground.vdw[k] += nocc[j] * (pairwise_vdw[mfe_res[i_res].conf[j]][all_res[k].conf[q]] * fround3(occ_table[all_res[k].conf[q]][t_point])) ;
                   E_ground.ele[k] += nocc[j] * (pairwise_ele[mfe_res[i_res].conf[j]][all_res[k].conf[q]] * fround3(occ_table[all_res[k].conf[q]][t_point])) ;
                }
            }
        }
    }

    mfe_res[i_res].mfe.vdw0 = E_ionize.vdw0 - E_ground.vdw0;
    mfe_res[i_res].mfe.vdw1 = E_ionize.vdw1 - E_ground.vdw1;
    mfe_res[i_res].mfe.ebkb = E_ionize.ebkb - E_ground.ebkb;
    mfe_res[i_res].mfe.tors = E_ionize.tors - E_ground.tors;
    mfe_res[i_res].mfe.dsol = E_ionize.dsol - E_ground.dsol;
    mfe_res[i_res].mfe.offset = E_ionize.offset - E_ground.offset;
    mfe_res[i_res].mfe.pHpK0 = E_ionize.pHpK0 - E_ground.pHpK0;
    mfe_res[i_res].mfe.EhEm0 = E_ionize.EhEm0 - E_ground.EhEm0;
    mfe_res[i_res].mfe.TS = E_ionize.TS - E_ground.TS;
    mfe_res[i_res].mfe.residues = E_ionize.residues - E_ground.residues;
    mfe_res[i_res].mfe.total = E_ionize.total - E_ground.total;
    for (j=0; j<n_all; j++) {
        mfe_res[i_res].mfe.mfePair[j] = E_ionize.mfePair[j] - E_ground.mfePair[j];
        mfe_res[i_res].mfe.vdw[j] = E_ionize.vdw[j] - E_ground.vdw[j];
        mfe_res[i_res].mfe.ele[j] = E_ionize.ele[j] - E_ground.ele[j];
    }

    free(E_ground.mfePair); free(E_ground.vdw); free(E_ground.ele);
    free(E_ionize.mfePair); free(E_ionize.vdw); free(E_ionize.ele);
    free(mfe); free(nocc); free(rocc);

    return 0;
}

struct STAT fit(float a, float b)
/* initialize the simplex and set up optimization */
{  float **p;
    float y[3];
    int neval;
    struct STAT result;

    p = (float **) malloc(3 * sizeof(float *));
    p[0] = (float *) malloc(2 * sizeof(float));
    p[1] = (float *) malloc(2 * sizeof(float));
    p[2] = (float *) malloc(2 * sizeof(float));

    p[0][0] = a;      p[0][1] = b;      y[0] = score(p[0]);
    p[1][0] = a+0.01; p[1][1] = b;      y[1] = score(p[1]);
    p[2][0] = a;      p[2][1] = b+1.0;  y[2] = score(p[2]);

    dhill(p, y, 2, 0.0001, &score, &neval);

    /* DEBUG
    printf("(%8.3f%8.3f)=%8.3f exit at %d\n", p[0][0], p[0][1], y[0], neval);
    */

    /* fit this eq. y = exp(a(x-b))/(1+exp(a(x-b))) */

    result.a = p[0][0];
    result.b = p[0][1];
    result.chi2 = y[0];

    free(p[0]);
    free(p[1]);
    free(p[2]);
    free(p);

    return result;
}


float score(float v[])
{  float S = 0.0;
    float yt;
    float T;
    int i;

    for (i=0; i<Nx; i++) {
        T = exp(v[0]*(xp[i]-v[1]));
        yt = T/(1.0+T);
        S += (yt-yp[i])*(yt-yp[i]);
    }

    return S;
}


#define TINY 1.0E-10
#define NMAX 5000
#define GET_PSUM for (j=0; j<ndim; j++) {\
                        for (sum=0.0, i=0; i<mpts; i++) sum += p[i][j];\
		                 psum[j] = sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk)
{
    float dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac);
    int i, ihi, ilo, inhi,j, mpts = ndim+1;
    float rtol,sum,swap,ysave,ytry,*psum;


    psum = (float *)malloc(ndim * sizeof(float));
    *nfunk = 0;
    GET_PSUM

    for (;;) {
        ilo = 0;
        ihi = y[0] > y[1] ? (inhi = 1,0):(inhi = 0,1);
        for(i=0; i<mpts;i++) {
            if(y[i] <= y[ilo]) ilo=i;
            if(y[i] > y[ihi]) {
                inhi = ihi;
                ihi  = i;
            }
            else if (y[i] > y[inhi] && i != ihi) inhi = i;
        }

        /* DEBUG
        for (i=0; i<mpts; i++) {
            printf("%8.3f at (%8.3f, %8.3f)\n", y[i], p[i][0], p[i][1]);
        }
        printf("\n");
        */

        rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);

        if(rtol < ftol || *nfunk >= NMAX) {
            SWAP(y[0],y[ilo])
            for (i=0; i<ndim; i++) SWAP(p[0][i],p[ilo][i])
                break;
        }

        *nfunk += 2;

        ytry = dhtry(p,y,psum,ndim,funk,ihi,-1.0);

        if (ytry <= y[ilo]) ytry=dhtry(p,y,psum,ndim,funk,ihi,2.0);
        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry=dhtry(p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave) {
                for (i=0; i<mpts; i++) {
                    if (i !=ilo) {
                        for (j=0; j<ndim; j++)
                            p[i][j] = psum[j] = 0.5*(p[i][j]+p[ilo][j]);
                        y[i] = (*funk)(psum);
                    }
                }
                *nfunk += ndim;
                GET_PSUM
            }
        }
        else --(*nfunk);
    }
    free(psum);
}

float dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac)
{
    int j;
    float fac1, fac2, ytry, *ptry;

    ptry = (float *)malloc(ndim*sizeof(float));
    fac1 = (1.0-fac)/ndim;
    fac2 = fac1 - fac;
    for (j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry = (*funk)(ptry);   /* evaluate the function at the trial point */
    if (ytry < y[ihi]) {    /* if it's better than the highest, then replace the highest */
        y[ihi] = ytry;
        for (j=0; j<ndim; j++) {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
    free(ptry);
    return ytry;
}

int enumerate(int i_ph_eh)
{
    long    nstate,istate;
    float   b;
    int     ires,jres,iconf,old_conf,new_conf;
    float   *E_states,*occ_states;
    float   E_min, tot_occ;
    b = -KCAL2KT/(env.monte_temp/ROOMT);

    if (!n_free) {
        for (iconf=0; iconf<conflist.n_conf; iconf++) {
            occ_table[iconf][i_ph_eh] = conflist.conf[iconf].occ;
        }
        return 0;
    }

    /* each residue state with first conformer */
    nstate = 1;
    for (ires=0;ires<n_free;ires++) {
        free_res[ires].on=0;
        state[ires]=free_res[ires].conf[0];
        nstate=nstate*free_res[ires].n;
        if (nstate>env.nstate_max) return -1;
    }
    if (nstate>env.nstate_max) return -1;


    fprintf(fp, "%ld <= %d microstates, will use analytical solution\n\n", nstate, env.nstate_max);
    fflush(fp);

    E_states = malloc(nstate*sizeof(float));
    occ_states = malloc(nstate*sizeof(float));
    enum_flag = 1; //Cai


    istate = 0;
    E_min = E_state = get_E();  /* get microstate energy */
    E_states[istate]=E_state;

    while (1) {
        istate++;
        ires = 0;
        old_conf = state[ires];
        free_res[ires].on++;
        while (free_res[ires].on>=free_res[ires].n) {
            free_res[ires].on=0;
            new_conf = state[ires] = free_res[ires].conf[0];

            E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
            for (jres=0; jres<n_free; jres++) {
                E_state += pairwise[new_conf][state[jres]] - pairwise[old_conf][state[jres]];
            }

            ires++;
            if (ires==n_free) break;
            old_conf = state[ires];
            free_res[ires].on++;
        }
        if (ires==n_free) break;

        new_conf = state[ires] = free_res[ires].conf[free_res[ires].on];
        E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
        for (jres=0; jres<n_free; jres++) {
            E_state += pairwise[new_conf][state[jres]] - pairwise[old_conf][state[jres]];
        }

        E_states[istate]=E_state;
        if (E_state<E_min) E_min = E_state;
    }

    tot_occ=0.;
    for (istate=0;istate<nstate;istate++) {
        occ_states[istate] = exp(b*(E_states[istate]-E_min));
        tot_occ += occ_states[istate];
    }
    for (istate=0;istate<nstate;istate++) {
        occ_states[istate] = occ_states[istate]/tot_occ;
    }

    /* get conformer occ */
    /* reset occ to 0 */
    for (iconf=0;iconf<conflist.n_conf;iconf++) {
        if (conflist.conf[iconf].on == 't') continue;
        conflist.conf[iconf].occ = 0.;
    }

    /* first microstate */
    istate = 0;
    for (ires=0;ires<n_free;ires++) {
        free_res[ires].on=0;
        state[ires]=free_res[ires].conf[0];

        iconf = free_res[ires].conf[free_res[ires].on];
        conflist.conf[iconf].occ += occ_states[istate];
    }

    while (1) {
        istate++;
        /* flipping (looping over) microstate */
        ires = 0;
        free_res[ires].on++;
        while (free_res[ires].on>=free_res[ires].n) {
            free_res[ires].on=0;
            new_conf = free_res[ires].conf[0];

            ires++;
            if (ires==n_free) break;
            free_res[ires].on++;
        }
        if (ires==n_free) break;

        /* collect occ for each working conformer */
        for (ires=0;ires<n_free;ires++) {
            iconf = free_res[ires].conf[free_res[ires].on];
            conflist.conf[iconf].occ += occ_states[istate];
        }
    }

    for (iconf=0; iconf<conflist.n_conf; iconf++) {
        occ_table[iconf][i_ph_eh] = conflist.conf[iconf].occ;
    }

    free(E_states);
    free(occ_states);
    return 0;
}

int enumerate_new(int i_ph_eh)  // new eneumerate subroutine to output microstate if true
{
    long    nstate,istate;
    float   b;
    int     ires,jres,iconf,old_conf,new_conf;
    float   *E_states,*occ_states;
    float   E_min, tot_occ;
    b = -KCAL2KT/(env.monte_temp/ROOMT);

    if (!n_free) {
        for (iconf=0; iconf<conflist.n_conf; iconf++) {
            occ_table[iconf][i_ph_eh] = conflist.conf[iconf].occ;
        }
        return 0;
    }

    /* each residue state with first conformer */
    nstate = 1;
    for (ires=0;ires<n_free;ires++) {
        free_res[ires].on=0;
        state[ires]=free_res[ires].conf[0];
        nstate=nstate*free_res[ires].n;
        if (nstate>env.nstate_max) return -1;
    }
    if (nstate>env.nstate_max) return -1;


    fprintf(fp, "%ld <= %d microstates, will use analytical solution\n\n", nstate, env.nstate_max);
    fflush(fp);

    enum_flag = 1; //turn on the enumeration flag

    E_states = malloc(nstate*sizeof(float));
    occ_states = malloc(nstate*sizeof(float));

    istate = 0;
    E_min = E_state = get_E();  /* get microstate energy */
    E_states[istate]=E_state;


    //initialize microstate
    MSRECORD ms_state;
    ms_state.conf_id = (unsigned short *) calloc(ms_spe_lst.n,  sizeof(unsigned short));
    ms_state.H       = 0.0;
    ms_state.counter = 0;
    ms_state.occ = 0.0;


    while (1) {
        istate++;
        ires = 0;
        old_conf = state[ires];
        free_res[ires].on++;
        while (free_res[ires].on>=free_res[ires].n) {
            free_res[ires].on=0;
            new_conf = state[ires] = free_res[ires].conf[0];

            E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
            for (jres=0; jres<n_free; jres++) {
                E_state += pairwise[new_conf][state[jres]] - pairwise[old_conf][state[jres]];
            }

            ires++;
            if (ires==n_free) break;
            old_conf = state[ires];
            free_res[ires].on++;
        }
        if (ires==n_free) break;

        new_conf = state[ires] = free_res[ires].conf[free_res[ires].on];
        E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
        for (jres=0; jres<n_free; jres++) {
            E_state += pairwise[new_conf][state[jres]] - pairwise[old_conf][state[jres]];
        }

        E_states[istate]=E_state;
        if (E_state<E_min) E_min = E_state;
    } //calculate energy for each microstates and find out the minimum energy state

    tot_occ=0.;
    for (istate=0;istate<nstate;istate++) {
        occ_states[istate] = exp(b*(E_states[istate]-E_min));
        tot_occ += occ_states[istate];
    }
    for (istate=0;istate<nstate;istate++) {
        occ_states[istate] = occ_states[istate]/tot_occ;
    }
/*   
    //write out occupancy for each microstate at mc_out
    if (env.ms_out){
        fprintf(fp, "microstate_id       OCC\n");
        for (istate=0; istate<nstate; istate++){
        fprintf(fp, "%ld        %5.3f\n", istate, occ_states[istate]);
        }
        fprintf(fp, "\n");
        fflush(fp);
    }
*/
    /* get conformer occ */
    /* reset occ to 0 */
    for (iconf=0;iconf<conflist.n_conf;iconf++) {
        if (conflist.conf[iconf].on == 't') continue;
        conflist.conf[iconf].occ = 0.;
    }

    /* first microstate */
    istate = 0;
    for (ires=0;ires<n_free;ires++) {
        free_res[ires].on=0;
        state[ires]=free_res[ires].conf[0];

        iconf = free_res[ires].conf[free_res[ires].on];
        conflist.conf[iconf].occ += occ_states[istate];
    }

    if (env.ms_out) {
        fwrite("ENUMERATE", 9, sizeof(char), ms_fp); //The third line of ms.dat: method
        if(env.re_ms_out) { //write out the method
            fprintf(re_ms_fp, "METHOD: %s\n", "ENUMERATE");
        }

        //write each microstate: write first microstate
        update_conf_id(ms_state.conf_id, state);
        ms_state.occ = occ_states[istate];
        ms_state.H = E_states[istate] + E_base;
        // ms_state.Hsq = (E_states[istate] + E_base) * (E_states[istate] + E_base);
        write_ms(&ms_state);
    }


    while (1) {
        istate++;
        /* flipping (looping over) microstate */
        ires = 0;
        free_res[ires].on++;
        while (free_res[ires].on>=free_res[ires].n) {
            free_res[ires].on=0;
            new_conf = state[ires] = free_res[ires].conf[0];

            ires++;
            if (ires==n_free) break;
            free_res[ires].on++;
        }
        if (ires==n_free) break;

        new_conf = state[ires] = free_res[ires].conf[free_res[ires].on];

        if (env.ms_out) {
        //write each microstate
        update_conf_id(ms_state.conf_id, state);
        ms_state.occ = occ_states[istate];
        ms_state.H = E_states[istate] + E_base;
        //ms_state.Hsq = (E_states[istate] + E_base) * (E_states[istate] + E_base);
        write_ms(&ms_state);
        }


        /* collect occ for each working conformer */
        for (ires=0;ires<n_free;ires++) {
            iconf = free_res[ires].conf[free_res[ires].on];
            conflist.conf[iconf].occ += occ_states[istate];
        }
    }

    for (iconf=0; iconf<conflist.n_conf; iconf++) {
        occ_table[iconf][i_ph_eh] = conflist.conf[iconf].occ;
    }

    free(E_states);
    free(occ_states);
    return 0;
}



int group_conftype()
{  int i;
   int n = 0;
   CONFTYPE oldID, newID;

   oldID = get_conftype(conflist.conf[0]);
   for (i=1; i<conflist.n_conf; i++) {
      newID = get_conftype(conflist.conf[i]);
      if (cmp_conftype(oldID, newID) == 0); /* same group of conformers */
      else {
         n++;
         oldID = newID;
      }
   }
   n++;

   /* now we have number of conformer groups */
   Sconverge.n = SconvergeBak.n = n;
   if (!(Sconverge.conftype = (CONFTYPE *) malloc(n*sizeof(CONFTYPE)))) {
      printf("   Memory allocation error in \"group_conftype\"\n");
      return USERERR;
   }
   if (!(SconvergeBak.conftype = (CONFTYPE *) malloc(n*sizeof(CONFTYPE)))) {
      printf("   Memory allocation error in \"group_conftype\"\n");
      return USERERR;
   }

   /* real grouping */
   oldID = get_conftype(conflist.conf[0]);
   n = 0;
   Sconverge.conftype[0] = oldID;
   for (i=1; i<conflist.n_conf; i++) {
      newID = get_conftype(conflist.conf[i]);
      if (cmp_conftype(oldID, newID) != 0) {
         n++;
         oldID = newID;
         Sconverge.conftype[n] = oldID;
      }
   }

   /* backup */
   for (i=0; i<Sconverge.n; i++) SconvergeBak.conftype[i] = Sconverge.conftype[i];

   return 0;
}

int cmp_conftype(CONFTYPE t1, CONFTYPE t2)
{  if (!strcmp(t1.resName, t2.resName) && \
        t1.chainID==t2.chainID && \
        t1.iCode==t2.iCode && \
        t1.resSeq==t2.resSeq && \
        t1.e == t2.e && \
        t1.H == t2.H &&\
        t1.DM == t2.DM)
      return 0;
   return -1;
}

CONFTYPE get_conftype(CONF conf)
{  CONFTYPE typeid;
   strcpy(typeid.resName, conf.resName);
   typeid.chainID = conf.chainID;
   typeid.iCode = conf.iCode;
   typeid.resSeq = conf.resSeq;
   typeid.e = conf.e;
   typeid.H = conf.H;
   if (conf.uniqID[3] == 'D' && conf.uniqID[4] == 'M') typeid.DM = 't';
   else typeid.DM = '\0';
   return typeid;
}

void update_Sconvergence()
/* This function actually computes the entropy of each comformer type */
{  int i, pstart, pend;
   CONFTYPE newID;
   pstart = 0; /* pointer to the line of conflist of the first conformer */
   pend = 1;   /* pointer to the line of the next new conf type */

   for (i=0; i<Sconverge.n; i++) {
      while (pend < conflist.n_conf) {
         newID = get_conftype(conflist.conf[pend]);
         if (cmp_conftype(Sconverge.conftype[i], newID) == 0) {
            pend ++;
            continue; /* same group of conformers */
         }
         else break;
      }
      Sconverge.conftype[i].E_TS = get_entropy(i, pstart, pend);
      /* Have done that in get_entropy */
      /* for (j=pstart; j<pend; j++) conflist.conf[j].E_TS = Sconverge.conftype[i].E_TS; */

      /* done with is type, go to the next */
      pstart = pend;
      pend++;
   }

   return;
}

float s_stat()
/* This function returns the maximum difference of Sconvergence and the backup */
{  int i;
   float max = 0.0;
   float diff;

   for (i=0; i<Sconverge.n; i++) {
      diff = fabs(Sconverge.conftype[i].E_TS - SconvergeBak.conftype[i].E_TS);
      if (diff > max) max = diff;
   }
   return max;
}

float get_entropy(int iconf, int start, int end)
{  int i;
   float TS = 0.0;
   float sum = 0.0;

   for (i=start; i<end; i++) {
      sum += conflist.conf[i].occ;
   }

   if (sum<0.0001) {
      for (i=start; i<end; i++) {
         conflist.conf[i].E_TS = 0.0;
      }
      return 0.0;
   }

   for (i=start; i<end; i++) {
      if (conflist.conf[i].occ/sum > 0.00001)
         TS -= conflist.conf[i].occ/sum*log(conflist.conf[i].occ/sum)/1.688;
   }

   /* assign the calculated entropy to conformers */
   for (i=start; i<end; i++) {
      conflist.conf[i].E_TS = TS;
   }


   return TS;
}


int load_ms_gold(STRINGS *str)  //Cai
{
    FILE *ms_gold_fp;
    char line[MAXCHAR_LINE], sbuff[MAXCHAR_LINE], cbuff[MAXCHAR_LINE];
    ATOM atom;
    int i_spe, j;

    if ((ms_gold_fp=fopen("ms_gold", "r"))) {
        printf("\n   Loading file \"ms_gold\"\n"); fflush(stdout);

        while (fgets(line, sizeof(line), ms_gold_fp)) {
            if (strlen(line) == 9) {
                strncpy(sbuff, line, 8);
                sbuff[8] = '\0';
            }
            else {
                if (strlen(line) < 24) continue;
                atom = pdbline2atom(line);
                if (atom.chainID == ' ') atom.chainID = '_';  // change the chain ID,if it's space         
                sprintf(sbuff, "%s%c%04d", atom.resName, atom.chainID, atom.resSeq);
            }
            for (j=0; j<conflist.n_conf; j++) {
                sprintf(cbuff, "%s%c%04d", conflist.conf[j].resName, conflist.conf[j].chainID, conflist.conf[j].resSeq);
                if (!strcmp(sbuff, cbuff)) break;
            }
            if (j==conflist.n_conf) {
     //           printf("      no conformer for residue %s, ignore it\n", sbuff);
                continue;
            }

            for (i_spe=str->n-1; i_spe>=0; i_spe--) {
                if (!strcmp(str->strings[i_spe], sbuff)) break;
            }
            if (i_spe == -1) {
                str->n++;
                str->strings = (char **) realloc(str->strings, str->n * sizeof(char *));
                str->strings[str->n-1] = (char *) malloc(strlen(sbuff) * sizeof(char));
                strcpy(str->strings[str->n-1], sbuff);
            }
        }
        printf("      Done, %d residue(s) marked.\n\n", str->n); fflush(stdout);
    }
    else {
        return USERERR;
    }

    return 0;
}


void MC_smp(int n)
{  int cycles, n_total, n_cycle;
    int i, j, k;
    register int iters;
    int mem;
    int *old_state;
    float  old_E;
    float dE;
    float b;
    int nflips;
    double H_average;
    double H_noTS;
    float E_entropy;

    int iflip, ires, iconf; /* iconf is 0 to n of the conf in a res */
    int old_conf, new_conf; /* old_conf and new_conf are from 0 to n_conf in conflist */


    b = -KCAL2KT/(env.monte_temp/ROOMT);

    mem =n_free * sizeof(int);
    old_state = (int *) malloc(mem);
    E_minimum = E_state = get_E();

    /* number of cycles and iters in each cycle */
    if (env.monte_trace > 0) {
        cycles = (n-1)/env.monte_trace + 1; /* round up */
        n_total = cycles*env.monte_trace;
        n_cycle = env.monte_trace;
    }
    else {
        cycles = 1;
        n_total = n_cycle = n;
    }

    /* clear counters */
    for (i=0; i<conflist.n_conf; i++) conflist.conf[i].counter = 0;
    H_average = 0.0;
    H_noTS = 0.0;

    /* Cai */
    MSRECORD ms_state;
    ms_state.conf_id = (unsigned short *) calloc(ms_spe_lst.n,  sizeof(unsigned short));
    ms_state.H       = 0.0;
    ms_state.counter = 0;
    ms_state.occ = 0.0;


    for (i=0; i<cycles; i++) {
        /*
        fprintf(fp, "Step %10d, E_minimum = %10.2f, E_running = %10.2f, E_reset = %10.2f\n",
        i*n_cycle, E_minimum+E_base, E_state+E_base, get_E()+E_base);
        */
        E_entropy = get_totalTS();
        fprintf(fp, "Step %10d, E_minimum= %10.2f, E_running = %10.2f ,  E_running without entropy = %10.2f\n",
        i*n_cycle, E_minimum+E_base, E_state+E_base, E_state+E_base - E_entropy);
        fflush(fp);
        iters = n_cycle;
        while (iters) {
            /*  save state */
            old_E = E_state;
            memcpy(old_state, state, mem);

            /* 1st flip */
            ires  = rand()/(RAND_MAX/n_free + 1);
            while (1) {
                iconf = rand()/(RAND_MAX/free_res[ires].n + 1);
                old_conf = state[ires];
                new_conf = free_res[ires].conf[iconf];
                if (old_conf != new_conf) break;
            }
            state[ires] = new_conf;
            E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
            for (j=0; j<n_free; j++) {
                E_state += pairwise[new_conf][state[j]] - pairwise[old_conf][state[j]];
            }

            /* now multiple flip */
            /*   1st flip -> No (50% probablity)
             *     | Yes (50% probablity)
             *   2nd flip (any res in big list)
             *     |
             *   3rd flip (any res in big list)
             *     |
             *   4th flip (any res in big list)
             */
            if (rand() & 1) {   /* do multiple flip if odd number */
                if (biglist[ires].n) {
                    nflips = env.monte_flips > (biglist[ires].n+1) ? biglist[ires].n+1: env.monte_flips;
                    for (k=1; k<nflips; k++) {

                        iflip = biglist[ires].res[rand()/(RAND_MAX/biglist[ires].n + 1)];
                        iconf = rand()/(RAND_MAX/free_res[iflip].n + 1);
                        old_conf = state[iflip];
                        new_conf = free_res[iflip].conf[iconf];

                        state[iflip] = new_conf;
                        E_state += conflist.conf[new_conf].E_self - conflist.conf[old_conf].E_self;
                        for (j=0; j<n_free; j++) {
                            E_state += pairwise[new_conf][state[j]] - pairwise[old_conf][state[j]];
                        }
                    }
                }
            }


            /* DEBUG
            for (j=0; j<n_free; j++) printf("%03d ", state[j]);
            printf("\n");
            */

            if (E_minimum > E_state) E_minimum = E_state;

            dE = E_state - old_E;
            /*if (dE < 0.0) {                                 // go to new low //
            }
            //<<< Boltman distribution >>>//
            else if ((float) rand()/RAND_MAX < exp(b*dE)) { // Go //
            }
            else {                                                    // stay, restore the state //
                memcpy(state, old_state, mem);
                E_state = old_E;
            }  */

            if (dE < 0.0 || (float) rand()/RAND_MAX < exp(b*dE)) {                                 /* go to new low */
                if (ms_state.counter != 0) write_ms(&ms_state);
                update_conf_id(ms_state.conf_id, state);
                ms_state.counter = 1;
                ms_state.H       = E_state + E_base;
                //ms_state.Hsq     = (E_state + E_base) * (E_state+E_base);
            }
            else {                                                    /* stay, restore the state */
                memcpy(state, old_state, mem);
                E_state = old_E;
                if (ms_state.counter != 0) {
                    ms_state.counter++;
                    ms_state.H    += E_state + E_base;
                    //ms_state.Hsq  += (E_state + E_base) * (E_state+E_base);
                }
                else {
                    update_conf_id(ms_state.conf_id, state);
                    ms_state.counter = 1;
                    ms_state.H       = E_state + E_base;
                    //ms_state.Hsq     = (E_state + E_base) * (E_state+E_base);
                }
            }


            /*
            if (!memcmp(state, state_check, sizeof(int)*n)) {
                printf("E=%12.5f\n", E_state);
            }  DEBUG */

            /* count this state energy */
            H_average += E_state/n_total;
            H_noTS += (E_state-get_totalTS())/n_total;
            /*<<< Do statistics >>>*/
            for (j=0; j<n_free; j++) {
                conflist.conf[state[j]].counter++;
                //+= pairwise_vdw[j]
            }

            iters --;
        }
    }

    E_entropy = get_totalTS();
    fprintf(fp, "Exit %10d, E_minimum = %10.2f, E_running = %10.2f E_running without TS = %10.2f\n", n_total, E_minimum+E_base, E_state+E_base, E_state+E_base-E_entropy);
    fprintf(fp, "The average running energy, corresponding to H = %8.3f kCal/mol, H without TS = %8.3f\n", H_average+E_base, H_noTS + E_base);
    fflush(fp);

    /* compute occ */
    for (i=0; i<conflist.n_conf; i++) {
        if (conflist.conf[i].on == 't') continue;
        conflist.conf[i].occ = (float) conflist.conf[i].counter / n_total;
    }

    free(old_state);

    return;
}


int write_ms(MSRECORD *ms_state)
{
    int i_spe;

    //printf("writing ms ...\n");
    if(env.re_ms_out){
    for (i_spe=0; i_spe<ms_spe_lst.n; i_spe++) {
        fwrite(&ms_state->conf_id[i_spe], 1, sizeof(unsigned short), ms_fp);
        fprintf(re_ms_fp,"%d\t", ms_state->conf_id[i_spe]);
    }
    }
    else{
    for (i_spe=0; i_spe<ms_spe_lst.n; i_spe++) {
        fwrite(&ms_state->conf_id[i_spe], 1, sizeof(unsigned short), ms_fp);
        //printf("%d\t", ms_state->conf_id[i_spe]);
    }
    }

    //write microstate at ms.dat, binary file
    fwrite(&ms_state->H, 1, sizeof(double), ms_fp);
    //fwrite(&ms_state->Hsq, 1, sizeof(double), ms_fp);
    if (enum_flag == 0) { //for MC sampling
        ms_state->Hav = ms_state->H/ms_state->counter;
        fwrite(&ms_state->Hav, 1, sizeof(double), ms_fp);
        fwrite(&ms_state->counter, 1, sizeof(int), ms_fp);
    }
    else {  //for enumerate, ms_state->counter is the occ of the microstate
        ms_state->Hav = ms_state->H;
        fwrite(&ms_state->Hav, 1, sizeof(double), ms_fp);
        fwrite(&ms_state->occ, 1, sizeof(double), ms_fp);
    }


    //for readable ms.dat (re_ms.dat)
    if (env.re_ms_out){
        fprintf(re_ms_fp,"\ncumulative energy: %lf\t", ms_state->H);
        fprintf(re_ms_fp,"state energy: %lf\t", ms_state->Hav);

        //fprintf(re_ms_fp,"cumulative energy^2: %lf\t", ms_state->Hsq);
        if (enum_flag == 0 ){   //for MC sampling, ms_state->counter is the times the microstate stays
        fprintf(re_ms_fp,"count: %d\n", ms_state->counter);}
    else{ //for enumerate, ms_state->counter is the occ of the microstate
        fprintf(re_ms_fp,"occ: %5.3f\n", ms_state->occ);}
    }

    return 0;
}

int update_conf_id(unsigned short *conf_id, int *state)
{
    bool found;
    int i_spe, i_free, i_fix, ic;

    for (i_spe=0; i_spe<ms_spe_lst.n; i_spe++) {
        found = false;
        for (i_free=0; i_free<n_free; i_free++) {
            //if (!strncmp(conflist.conf[state[i_free]].uniqID+5, ms_spe_lst.strings[i_spe]+3, 5)) {
            if ((!strncmp(conflist.conf[state[i_free]].uniqID+5, ms_spe_lst.strings[i_spe]+3, 5))
&& (!strncmp(conflist.conf[state[i_free]].uniqID, ms_spe_lst.strings[i_spe], 3))) {
                found = true;
                conf_id[i_spe]  = state[i_free];
                break;
            }
        }
        if (!found) {
            for (i_fix=0; i_fix<n_fixed; i_fix++) {
                //if (!strncmp(conflist.conf[fixed_res[i_fix].conf[0]].uniqID+5, ms_spe_lst.strings[i_spe]+3, 5)) {
                if ((!strncmp(conflist.conf[fixed_res[i_fix].conf[0]].uniqID+5, ms_spe_lst.strings[i_spe]+3, 5))  && (!strncmp(conflist.conf[fixed_res[i_fix].conf[0]].uniqID, ms_spe_lst.strings[i_spe], 3))) {
                    for (ic=0; ic<fixed_res[i_fix].n; ic++) {
                        // this will be a problem, if partial occ is assigned
                        if (conflist.conf[fixed_res[i_fix].conf[ic]].occ > 0.99) {
                            conf_id[i_spe] = fixed_res[i_fix].conf[ic];
                            found = true;
                            break;
                        }
                    }
                }
                if (found) break;
            }
        }
        
        if (!found) {
            //printf("special list can't find micro states\n");
        }   
    }   
    
    return 0;
}   



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <syscall.h>
#include <math.h>
#include "mcce.h"
#include <sys/stat.h>

#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "ga_engine.h"

typedef struct {
   int *nc_start;
   int *nc_rmdup;
   int *nc_swap;
   int *nc_rotate;
   int *nc_vdw;
   int *nc_hdir;
   int *nc_repack;
   int *nc_ioniz;
   int *nc_h;
   int *nc_oh;
   int *nc_elect;
} CONFSTAT;

extern void wait ( int seconds );

time_t nowA, nowB, nowStart, nowEnd;
long    idum;

extern int place_rot(PROT prot);
extern int place_rot_rule(int i_res, ROTAMER rule, int n, PROT prot);
extern int swing_rot(PROT prot);
extern int swing_rot_rule(int i_res, ROTAMER rule, float phi, PROT prot);
extern void write_step2stat(FILE *fp, PROT prot, CONFSTAT stat);
extern int prune_by_vdw(PROT prot, float delta_E);
extern int prune_by_vdw_res(int kr, PROT prot, float delta_E);
extern int rot_pack(PROT prot, int n);
extern int extra_rot(PROT prot);
extern int rot_refine(PROT prot, MICROSTATE state, float ****pariwise);
extern int ionization(PROT prot);
extern int rm_dupconf(PROT prot, float prune_thr);
extern int rm_dupconf_hv(PROT prot);
extern int rm_dupconf_res(PROT prot, int i_res, float prune_thr);
extern int write_conflist(FILE *fp, PROT prot);
extern int load_listrot(char *fname, PROT prot);
extern int get_demetri_out(void);
extern int prune_hdirected(PROT prot);
extern float relax_this_pair(PROT prot, RES *res_a, int conf_a, RES *res_b, int conf_b);
extern int swing_conf_rot_rule(RES *res, int j, ROTAMER rule, float phi);
extern int rebuild_sc(PROT prot);
extern void del_non_common_h(PROT prot);
extern int rot_swap(PROT prot);
extern int get_resnbrs(PROT prot, float dlimit);
extern int relaxation(PROT prot);
/*extern int self_relaxation(PROT prot);*/
extern int relax_h(PROT prot);
extern int relax_water(PROT prot);
extern void collect_all_connect(int i_res, int i_conf, int i_atom, PROT prot, int *n_connect12, ATOM ***connect12, int *n_connect13, ATOM ***connect13, int *n_connect14, ATOM ***connect14);
extern int print_vdw_map(PROT prot);
extern int initial_relaxation(PROT prot);
extern int prune_pv(PROT prot, float c1, float c2, float c3);
extern int over_geo(CONF conf1, CONF conf2, float cutoff);
extern int over_ele(PROT prot, int ir, int ic, int jc, float cutoff);
extern int over_vdw(PROT prot, int ir, int ic, int jc, float cutoff);
extern float max_confSAS(PROT prot, int ir, int ic);
extern int label_exposed(PROT prot);
extern int rand_conf_prune(PROT prot);
extern int place_missing_res(PROT prot, int i_res, int handle_addconf);

/* Pascal MSC*/
PROT load_pdb_binary_no_param(FILE *fp);
int atom2pdbline(char *line, ATOM atom);
int write_pdb_binary(FILE *stream, PROT prot);
extern PATCH compute_patches(PROT prot, int mode, int N, float eps, float b1_param, float b2_param, int , int);
/* Pascal MSC */

void prep_for_step3_and_step4(PROT prot) {
   int c, i, j, k;
   FILE *fp;
   char sbuff[MAXCHAR_LINE];

   /* Relabel rotamer history */
   for (i=0; i<prot.n_res; i++) {
     for (j=0; j<prot.res[i].n_conf; j++) { 
     //for (j=prot.res[i].n_conf_ori; j<prot.res[i].n_conf; j++) {
         if (prot.res[i].conf[j].history[2] == 'H') {
            sprintf(sbuff, "H%03X", j-1);
            strncpy(prot.res[i].conf[j].history+2, sbuff, 4);
         }
         else if (prot.res[i].conf[j].history[2] == 'E') {
            sprintf(sbuff, "E%03X", j-1);
            strncpy(prot.res[i].conf[j].history+2, sbuff, 4);
         }
         else if (prot.res[i].conf[j].history[2] == 'I') {
            sprintf(sbuff, "I%03X", j-1);
            strncpy(prot.res[i].conf[j].history+2, sbuff, 4);
         }
         else if (prot.res[i].conf[j].history[2] == 'B') {
            sprintf(sbuff, "B%03X", j-1);
            strncpy(prot.res[i].conf[j].history+2, sbuff, 4);
         }
         else if (prot.res[i].conf[j].history[2] == 'Y') {
            sprintf(sbuff, "Y%03X", j-1);
            strncpy(prot.res[i].conf[j].history+2, sbuff, 4);
         }
         else {
            sprintf(sbuff, "R%03X", j-1);
            strncpy(prot.res[i].conf[j].history+2, sbuff, 4);
         }
      }
   }

   printf("   Making ionization conformers...\n"); fflush(stdout);
   if (ionization(prot)) {
      printf("   FATAL: Fatal error reported by ionization()\n");
      return exit(0);
   }
   printf("   Done\n\n"); fflush(stdout);

   /*need to fix PDB read/write for more than 1000 conformers or remove those above count of 1000?*/ 
   for(i=0; i<prot.n_res; i++) {
   	if(prot.res[i].n_conf>1000) {
		printf("Res %s(%i) has more than 1000 conformers. Those above count of 1000 are removed.\n",prot.res[i].resName,prot.res[i].resSeq);
		check:
		for(j=1000; j<prot.res[i].n_conf; j++) {
			del_conf(&prot.res[i],j);
			goto check;
		}
	}
   }
   
   // add h 
   printf("   Add H atoms...\n"); fflush(stdout);
   while(place_missing(prot,1) > 0);
   rm_dupconf(prot, 0.001);

   assign_rad(prot);
   assign_vdw_param(prot);
   assign_crg(prot); // reassign parameters before rm_dupconf(), Yifan -02/07/07
   printf("   Done\n\n"); fflush(stdout);

   /* IPECE add membrane */
   if (env.ipece.add_mem) {
        printf("   Add membrane.\n");   fflush(stdout);
        add_membrane(&prot, &env.ipece);
        printf("   Done\n\n");   fflush(stdout);
   }

   /* write out */
   printf("   Write output...\n"); fflush(stdout);
   /* mcce pdb */
   printf("Writing Step2_out.pdb...\n");
   fp = fopen(STEP2_OUT, "w");
   write_pdb(fp, prot);
   fclose(fp);
  
   printf("Writing headlist2...\n");
   /* conformer summary */
   fp = fopen(FN_CONFLIST2, "w");
   write_conflist(fp, prot);
   fclose(fp);
   
   /* Do step 3. energies */
   if (env.do_energies) {
      printf("Step 3. Compute energy lookup table\n"); fflush(stdout);
      if (energies()) {db_close(); exit(0);}
      else printf("Step 3 Done.\n\n");
   }
   else printf("Not doing \"Step 3. Compute energy lookup table\"\n\n");

   /* Do step 4. Monte Carlo */
   if (env.do_monte) {
      printf("Step 4. Monte Carlo Sampling\n"); fflush(stdout);
      if (!env.monte_adv_opt) {
      if (monte()) {db_close();} 
           else printf("Step 4 Done.\n\n");
       }
       else {
           if (monte2()) {db_close();}
           else printf("Step 4 Done.\n\n");
       }
   }
   else printf("Not doing \"Step 4. Monte Carlo Sampling\"\n\n");
};

int rotamers_GA(int argc, char *argv[]) {
   FILE *fp;
   PROT prot;
   CONFSTAT confstat;
   int c, i, j, kr, kc,ic,k,m;
   char sbuff[MAXCHAR_LINE];

   nowStart = time(NULL);
   idum = time(NULL);

   FILE *pdb_fp;
   int i_res, i_conf, i_atom;
   char line[MAXCHAR_LINE];
   
   int one_patch=1;
   int **result;
   GA_STRUCTURE ga;

   ga.dE = env.ga_deltaE;
   ga.pop_size = env.pop_size;
   ga.mutation_rate = env.mutation_rate;
   ga.migration_rate = env.migration_rate;
   ga.xover_rate = env.xover_rate;
   ga.elitism = 0.01;/*not currently used, only THE best is kept at each generation*/
   ga.seed = env.ga_seed;/* -1 for random one*/
   ga.k_select = 3;//env.k_select;/*not currently used; using Boltzmann selection pressure*/
   ga.generations = env.generations;

   ga.phase = env.ga_phase;
   ga.shift = env.ga_shift;
   ga.sphere_focus_resid = env.ga_focus_resid;
   ga.sphere_focus_probe_radius = env.ga_focus_probe_radius;

   ga.dist_center = env.ga_dist_center;
   ga.pop_bucket = env.pop_bucket;/*final unique rotamers in final population bucket*/
   ga.dist_center_eps = env.ga_dist_center_eps;/*how far within KT*KT_EPS from the min do you want your solutions to be*/
   ga.residue_check = env.residue_check;/*up to how much from each residue lowest occupied minimal energy do you want to keep rotamers for*/
   ga.rand_cut_points = env.rand_cut_points;
   ga.occupancy = env.occupancy;

   /*check for resume_patches.txt file which tells you which .pdb files to open, 1 at a time, open them all.*/
   GA:
   if ((pdb_fp=fopen("./pdb_patches/patches_resume.txt", "rt"))) {
        int nb_p=0;
	int tot_ga_conf;
	char fname[32];
	char dir[64];
        fscanf(pdb_fp,"%i\n",&nb_p);
        if(nb_p>0) {
                printf("Binary PDB file found...\n");
                for(i=0; i<nb_p; i++) {fscanf(pdb_fp,"%s\n",fname);}
                fclose(pdb_fp);

                pdb_fp = 0;
                sprintf(dir,"./pdb_patches/"); strcat(dir,fname);
                pdb_fp = fopen(dir,"rb");
                prot = load_pdb_binary_no_param(pdb_fp);//called _no_param but it does load params.
                fclose(pdb_fp);
		/*check if file for this ga_patch already exists before computing GA over again*/
		/*if it doesn't exist, then we need to compute GA*/
                sprintf(dir,"./pdb_patches/ga_output_patch%04i", i); fp = fopen(dir,"rt");
		if(fp==NULL) {
			printf("No GA Resume output file found. Starting new GA...\n");
			GA_SETUP(&prot, ga, nb_p);
		}
		else {
			printf("Loading results from GA...\n");
			/*load results saved from GA*/
			result = calloc(prot.n_res, sizeof(int*));
			for(j=0; j<prot.n_res; j++) {
				result[j] = calloc(prot.res[j].n_conf,sizeof(int));
				for(k=0; k<prot.res[j].n_conf; k++) {fscanf(fp,"%i\n",&result[j][k]);}
			}
			fclose(fp);
		
			printf("Flagging rotamers to keep...\n");
			/*flag which rotamers to keep based on ga output*/
		        for(j=0; j<prot.n_res; j++) {
		                prot.res[j].conf[0].tmp_flag = 1;/*always keep the backbone*/
                		if(prot.res[j].n_conf>1) {prot.res[j].conf[1].tmp_flag =1;}
                		for(k=2; k<prot.res[j].n_conf; k++) {
                       			if(result[j][k] > 0) {prot.res[j].conf[k].tmp_flag = 1;}
                       			else {prot.res[j].conf[k].tmp_flag = 0;}
                		}
        		}

			printf("Removing unused rotamers...\n");
			/*delete conformers that are flagged 0*/
		        for(j=0; j<prot.n_res; j++) {
                		rem_conf:
                		for(k=2; k<prot.res[j].n_conf; k++) {/*don't delete backbone or initial conformation*/
                       			if(prot.res[j].conf[k].tmp_flag==0) {del_conf(&prot.res[j], k); goto rem_conf;}
                		}
        		}
			/*reset the flag of all conformers to 0*/
                        for(j=0; j<prot.n_res; j++) {
                        	for(k=0; k<prot.res[j].n_conf; k++) {prot.res[j].conf[k].tmp_flag =0;}
                        }
		}
		
		tot_ga_conf=0;
                for(j=0; j<prot.n_res; j++) {
                	tot_ga_conf += (prot.res[j].n_conf-1);
                }
                printf("%i conformers from Genetic Algorithm\nResuming MCCE Step3 and Step4...\n", tot_ga_conf);
                prep_for_step3_and_step4(prot);
                exit(0);
        }
	else {
		goto step1_out;
	}
        if (prot.n_res == 0) {
            printf("   Fail to load binary resume pdb file for GA Run. Starting regular MCCE Step2...\n");
        }
    }
    else {
        printf("   Fail to open binary resume pdb file for GA Run. Starting regular MCCE Step2...\n");
    }

  step1_out:
  if (!(fp=fopen(STEP1_OUT, "r"))) {
      printf("   FATAL: rotamers(): \"No step 1 output \"%s\".\n", STEP1_OUT);
      return USERERR;
   }
  prot = load_pdb(fp);

   if (prot.n_res == 0) {
      printf("   There are errors in pdb file, quiting ...\n");
      return USERERR;
   }
   printf("   Done loading %i residues\n\n", prot.n_res); fflush(stdout);
   fclose(fp);

   /* load CONFLIST1 */
   load_headlst(prot);

   /* prepare arrys to store number of conformers after r tatyy swing, prune, repack,
    * place_missing, and optimize OH, last slot for total */
   confstat.nc_start  = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_rmdup  = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_swap   = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_rotate = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_vdw    = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_hdir   = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_repack = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_ioniz  = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_h      = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_oh     = (int *) calloc((prot.n_res+1) , sizeof(int));
   confstat.nc_elect  = (int *) calloc((prot.n_res+1) , sizeof(int));
   get_resnbrs(prot, 5.0);

   /* count conformers */
   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_start[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_start[prot.n_res] = c;

   printf("   Rotamer statistics is dynamically updated in file \"%s\"\n\n", ROTSTAT);

   /* Delete extra confs */
   printf("   Remove redundant heavy atom conformers in %s\n", STEP1_OUT);
   printf("   %d conformers were deleted.\n   Done\n\n", rm_dupconf_hv(prot));
   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_rmdup[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_rmdup[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   if (fp) {
       fprintf(fp,"   Rotamer making statitics:\n");
       write_step2stat(fp, prot, confstat);
       printf("\n");
       fclose(fp);
   }

   /* Relabel rotamer history 
   - disable in mcce2.4 version to avoid a conformer explosion in the relaxation
   step when there are multiple conformers in step1_out - Yifan 10/18/2007 */

   /*BIG BUG FIXED by PASCAL 2009:
 * this loop needs to count the number of original conformers based on the history
 * DO NOT assume the number of original conformers to come from 'n_conf'
 * 'n_conf_ori' is used in REPACK phase and is critical to remove low occupancy conformers (rotamers)
 * as it will NOT remove the original conformers */   
   for(i=0; i<prot.n_res; i++) {
      prot.res[i].n_conf_ori=0;
      prot.res[i].original_index = i;
      for(j=0; j<prot.res[i].n_conf; j++) {
        if ( strncmp(&prot.res[i].conf[j].history[2],"R",1) )  prot.res[i].n_conf_ori++;
      }
   }


/* Spherical Focus Code */
   int ga_focus_resid;
   float ga_focus_probe_radius;

   PATCH pat;
   int *list;
   if(env.ga_focus_resid!=-1) {
        int iso_res=-1;
        for(i=0; i<prot.n_res; i++) {
                if(prot.res[i].resSeq==env.ga_focus_resid) {
                	printf("Creating spherical focus GA centered at RES[%s](Seq:%i)\n",prot.res[i].resName,env.ga_focus_resid);
                        iso_res = i;
                        break;
                }
        }
	if(iso_res!=-1) {
		pat = compute_patches(prot,3,0,0,0,0,env.ga_focus_resid,iso_res);
		goto process;
	}
	printf("Residue sequence %i not found. Exiting...\n",env.ga_focus_resid);
	exit(1);
   }
   else {
        printf("Regular sidechain packing\n");
        one_patch = 1;
        goto patch_resume;
   }
   process:
   one_patch=1;
   list = calloc(prot.n_res,sizeof(int));
   for(i=0; i<prot.n_res; i++) {list[i] = 0;}
   for(i=0; i<pat.pat[0].res_count; i++) {list[pat.pat[0].res_ids[i]] = 1;}
   for(i=0; i<pat.pat[0].b1_count; i++) {list[pat.pat[0].boundary1[i]] = 1;}

   if (env.n_initial_relax) {
       printf("   Initial relaxation...\n"); fflush(stdout);
       initial_relaxation(prot);
       printf("   Done\n\n"); fflush(stdout);
   }

   /* Relax water */
   if (env.relax_wat) {
       printf("   Relax crystal water..\n"); fflush(stdout);
       printf("   Done\n\n"); fflush(stdout);
   }

   if (env.rebuild_sc) {
       nowA = time(NULL);
       printf("   Rebuild sidechains...\n"); fflush(stdout);
       rebuild_sc(prot);
       nowB= time(NULL);
       printf("   Done. Time spent = %d\n\n", (int)(nowB - nowA)); fflush(stdout);
   }

   get_connect12(prot); delete_h(prot); assign_vdw_param(prot); assign_rad(prot); surfw(prot, 1.4);

   /*if res is not flagged as in the list, then remove its conformers from protein structure but leave the backbone there*/
   for(i=0; i<prot.n_res; i++) {
        if(list[i]==0 && (strncmp(prot.res[i].resName,"CTR",3)!=0) && (strncmp(prot.res[i].resName,"NTR",3)!=0) && (strncmp(prot.res[i].resName,"NTG",3)!=0)) {
                del:
                for(j=1; j<prot.res[i].n_conf; j++) {
                        del_conf(&prot.res[i],j);
                        goto del;
                }
        }
   }
   goto patch_continue;
/* Spherical Focus Code */


   patch_resume:
   /* Relax original conformer */
   if (env.n_initial_relax) {
       printf("   Initial relaxation...\n"); fflush(stdout);
       initial_relaxation(prot);
       printf("   Done\n\n"); fflush(stdout);
   }
   
   /* Relax water */
   if (env.relax_wat) {
       printf("   Relax crystal water..\n"); fflush(stdout);
       printf("   Done\n\n"); fflush(stdout);
   }
   
   if (env.rebuild_sc) {
       nowA = time(NULL);
       printf("   Rebuild sidechains...\n"); fflush(stdout);
       rebuild_sc(prot);
       nowB= time(NULL);
       printf("   Done. Time spent = %d\n\n", (int)(nowB - nowA)); fflush(stdout);
   }

   printf("   Prepare for rotamer making ...\n");
   get_connect12(prot);
   printf("   Deleting H atoms...%d H atoms were deleted.\n", delete_h(prot)); fflush(stdout);
   assign_vdw_param(prot);
   printf("   Assigning radii.\n"); fflush(stdout); assign_rad(prot);
   printf("   Estimating Solvent Accessible Surface (SAS).\n"); fflush(stdout); surfw(prot, 1.4);
   
   patch_continue:
   write_headlst(prot);
   printf("   Done.\n\n");

   /* Swap */
   if (env.rot_swap) {
       printf("   Swap atoms...\n"); fflush(stdout);
       rot_swap(prot);
       printf("   Done\n\n"); fflush(stdout);
   }
   /* count conformers */
   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_swap[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_swap[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   if (fp) {
       fprintf(fp,"   Rotamer making statitics:\n");
       write_step2stat(fp, prot, confstat);
       printf("\n");
       fclose(fp);
   }

   /* Place conformers */
   if (env.pack) {
      printf("   Place rotamers...\n"); fflush(stdout);
      place_rot(prot);
      printf("   Done\n\n"); fflush(stdout);
   }
   else {
      printf("NOT MAKING ROTAMERS!!!\n");
   }
   /* Do swing */
   if (env.swing) {
      printf("   Swing rotamers...\n"); fflush(stdout);
      swing_rot(prot);
      printf("   Done\n\n"); fflush(stdout);
   }
   /* Extra rotamers (translation etc.) */
   printf("   Extra rotamers...\n"); fflush(stdout);
   extra_rot(prot);
   printf("   Done\n\n"); fflush(stdout);

   /* count conformers */
   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_rotate[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_rotate[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   if (fp) {
       fprintf(fp,"   Rotamer making statitics:\n");
       write_step2stat(fp, prot, confstat);
       printf("\n");
       fclose(fp);
   }

   /* delete duplicate conformers */
   nowA = time(NULL);
   printf("   Delete duplicate conformers ..."); fflush(stdout);
   int counter_total_deleted = 0;
   for (kr=0; kr<prot.n_res; kr++) {
       /* not using rm_dupconf_hv() function because it saves time to skip residues if it starts with only 1 conformer */
       if (prot.res[kr].n_conf_ori <=2) continue;
       for (kc=1; kc<prot.res[kr].n_conf; kc++) {
           for (ic=kc+1; ic<prot.res[kr].n_conf; ic++) {
               if(!cmp_conf_hv(prot.res[kr].conf[kc], prot.res[kr].conf[ic], env.prune_thr)) {
                   del_conf(&prot.res[kr], ic);
                   ic--;
                   counter_total_deleted++;
               }
           }
       }
   }
   printf(" %d conformers deleted.\n", counter_total_deleted);
   nowB= time(NULL);
   i=nowB - nowA;
   printf("   Done. Time spent = %d\n\n", i); fflush(stdout);

   /* prune by self vdw */
   nowA = time(NULL);
   printf("   Prune rotamers by self VDW potential...\n");
   
   /* add protons for self energy pruning */
   for (kr=0; kr<prot.n_res; kr++) {
       while(place_missing_res(prot,kr,1) > 0); rm_dupconf_res(prot, kr, 0.005);
   }
   del_non_common_h(prot);
   for (kr=0; kr<prot.n_res; kr++) {
       rm_dupconf_res(prot, kr, 0.005);
   }
   
   printf("   Creating connectivity table...\n"); fflush(stdout);
   get_connect12(prot);
   
   printf("   Computing self VDW potential. It may take a while...\n"); fflush(stdout);
   assign_vdw_param(prot);
   get_vdw0_no_sas(prot);
   get_vdw1(prot);
   for (kr=0; kr<prot.n_res; kr++) {
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
          prot.res[kr].conf[kc].E_torsion = torsion_conf(&prot.res[kr].conf[kc]);
          prot.res[kr].conf[kc].E_self = prot.res[kr].conf[kc].E_vdw0 + prot.res[kr].conf[kc].E_vdw1 + prot.res[kr].conf[kc].E_torsion;
      }
   }
   printf("   Pruning rotamers...");fflush(stdout);
   printf("%d rotamers deleted.\n", prune_by_vdw(prot, env.vdw_cutoff));
   /* remove protons again */
   delete_h(prot);
   rm_dupconf_hv(prot);
   
   nowB= time(NULL); i = nowB - nowA;
   printf("   Done. Time spent = %d\n\n", i); fflush(stdout);

   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_vdw[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_vdw[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   if (fp) {
       fprintf(fp,"   Rotamer making statitics:\n");
       write_step2stat(fp, prot, confstat);
       printf("\n");
       fclose(fp);
   }

   /* find the most exposed conformer */
   /* using swing subroutine to optimize surface sidechains */
   nowA = time(NULL);
   printf("   Tuning exposed rotamers...\n");fflush(stdout);
   label_exposed(prot); /* history[2] = 'E' for the most exposed */

   //get_connect12(prot);
   //get_vdw0(prot);
   //get_vdw1(prot);
   nowB= time(NULL); i = nowB - nowA;
   printf("   Done. Time spent = %d\n\n", i); fflush(stdout);

   /* Hydrogen bond directed rotamer making */
   //get_resnbrs(prot, 5.0);
   if (env.hdirected) {
      nowA = time(NULL);
      printf("   Hydrogen bond directed rotamer making...\n"); fflush(stdout);

      if (prune_hdirected(prot)) {
         printf("   Fatal error.\n");
         return USERERR;
      }

      /* refine and update */
      printf("   Updating self vdw\n"); fflush(stdout);
      get_connect12(prot);
      get_vdw0_no_sas(prot);
      get_vdw1(prot);
      for (kr=0; kr<prot.n_res; kr++) {
         for (kc=1; kc<prot.res[kr].n_conf; kc++) {
             prot.res[kr].conf[kc].E_torsion = torsion_conf(&prot.res[kr].conf[kc]);
             prot.res[kr].conf[kc].E_self = prot.res[kr].conf[kc].E_vdw0 + prot.res[kr].conf[kc].E_vdw1 + prot.res[kr].conf[kc].E_torsion;
         }
      }
      printf("   Pruning rotamers by new self energy...");fflush(stdout);
      printf("%d rotamers deleted.\n", prune_by_vdw(prot, env.vdw_cutoff));
      rm_dupconf_hv(prot);

      nowB= time(NULL);
      i=nowB - nowA;
      printf("   Done. Time spent = %d\n\n", i); fflush(stdout);
   }
   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_hdir[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_hdir[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   if (fp) {
       fprintf(fp,"   Rotamer making statitics:\n");
       write_step2stat(fp, prot, confstat);
       printf("\n");
       fclose(fp);
   }

   FILE* output_resume;
   char filepatch[32];

   mkdir("./pdb_patches",S_IRWXU);
   sprintf(filepatch,"./pdb_patches/output_resume_patch%04i.pdb",one_patch);
   FILE* pf = fopen("./pdb_patches/patches_resume.txt","wt");
   fprintf(pf,"%i\n",one_patch);
   fprintf(pf,"output_resume_patch%04i.pdb\n",one_patch);
   fclose(pf);

   output_resume = fopen(filepatch,"wb");
   write_pdb_binary(output_resume,prot);
   fclose(output_resume);

   del_prot(&prot);
   goto GA;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mcce.h"

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

time_t nowA, nowB, nowStart, nowEnd;
long    idum;

int place_rot(PROT prot);
int place_rot_rule(int i_res, ROTAMER rule, int n, PROT prot);
int swing_rot(PROT prot);
int swing_rot_rule(int i_res, ROTAMER rule, float phi, PROT prot);
void write_step2stat(FILE *fp, PROT prot, CONFSTAT stat);
int prune_by_vdw(PROT prot, float delta_E);
int prune_by_vdw_res(int kr, PROT prot, float delta_E);
int rot_pack(PROT prot, int n);
int extra_rot(PROT prot);
int rot_refine(PROT prot, MICROSTATE state, float ****pariwise);
int ionization(PROT prot);
int rm_dupconf(PROT prot, float prune_thr);
int rm_dupconf_hv(PROT prot);
int rm_dupconf_res(PROT prot, int i_res, float prune_thr);
int write_conflist(FILE *fp, PROT prot);
int load_listrot(char *fname, PROT prot);
int get_demetri_out(void);
int prune_hdirected(PROT prot);
float relax_this_pair(PROT prot, RES *res_a, int conf_a, RES *res_b, int conf_b);
int swing_conf_rot_rule(RES *res, int j, ROTAMER rule, float phi);
int rebuild_sc(PROT prot);
void del_non_common_h(PROT prot);
extern int rot_swap(PROT prot);
extern int get_resnbrs(PROT prot, float dlimit);
extern int relaxation(PROT prot);
//extern int self_relaxation(PROT prot);
extern int relax_h(PROT prot);
extern int relax_water(PROT prot);
extern void collect_all_connect(int i_res, int i_conf, int i_atom, PROT prot, int *n_connect12, ATOM ***connect12, int *n_connect13, ATOM ***connect13, int *n_connect14, ATOM ***connect14);
extern int print_vdw_map(PROT prot);
extern int initial_relaxation(PROT prot);
int prune_pv(PROT prot, float c1, float c2, float c3);
int over_geo(CONF conf1, CONF conf2, float cutoff);
int over_ele(PROT prot, int ir, int ic, int jc, float cutoff);
int over_vdw(PROT prot, int ir, int ic, int jc, float cutoff);
int over_vdw_mhd(PROT prot, int ir, int ic, int jc, float cutoff);
float max_confSAS(PROT prot, int ir, int ic);
int label_exposed(PROT prot);
int rand_conf_prune(PROT prot);
int rand_conf_prune_mhd(PROT prot);
int place_missing_res(PROT prot, int i_res, int handle_addconf);

int rotamers()
{
    FILE *fp;
    PROT prot;
    CONFSTAT confstat;
    int c, i, j, kr, kc,ic;
    char sbuff[MAXCHAR_LINE];

   nowStart = time(NULL);
   if (env.test_seed < 0) idum = time(NULL); //allows random numbers to be fixed for testing
   else idum = env.test_seed;

   /* Load step 1 output pdb file */
   printf("   Load step 1 output file %s...\n", STEP1_OUT);
   if (!(fp=fopen(STEP1_OUT, "r"))) {
      printf("   FATAL: rotamers(): \"No step 1 output \"%s\".\n", STEP1_OUT);
      return USERERR;
   }
   prot = load_pdb(fp);
   if (prot.n_res == 0) {
      printf("   There are errors in pdb file, quiting ...\n");
      return USERERR;
   }
   printf("   Done\n\n"); fflush(stdout);
   fclose(fp);

   /* write out ms_gold --by Cai*/
   if (env.ms_gold_out) {
       printf("   Writing ms_gold ...\n");
       fflush(stdout);
       fp = fopen("ms_gold", "w");
       write_ms_gold(fp, prot);
       fclose(fp);
       printf("   Done. ms_gold is created.\n\n", i); fflush(stdout);
   }
   else {
      printf("   NOT output ms_gold file.\n");
   }


   /* load CONFLIST1 */
   load_headlst(prot);

   /* prepare arrys to store number of conformers after rotate, swing, prune, repack,
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

   /* DEBUG
   for (i=0; i<prot.n_res; i++) {
      for (j=0; j<prot.res[i].n_ngh; j++) {
         printf("RES %s %04d: %s %04d\n", prot.res[i].resName, prot.res[i].resSeq, prot.res[i].ngh[j]->resName, prot.res[i].ngh[j]->resSeq);
      }
      printf("\n");
   }
   return 0;
   */

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
   
   for (i=0; i<prot.n_res; i++) {
      prot.res[i].n_conf_ori = prot.res[i].n_conf;
      /*
      for (j=1; j<prot.res[i].n_conf; j++) {
         sprintf(sbuff, "O%03X", j-1);
         strncpy(prot.res[i].conf[j].history+2, sbuff, 4);
      }
      */
   }
   
   /* Relax original conformer */
   if (env.n_initial_relax) {
       printf("   Initial relaxation...\n"); fflush(stdout);
       initial_relaxation(prot);
       printf("   Done\n\n"); fflush(stdout);
   }
   
   /* Relax water */
   if (env.relax_wat) {
       printf("   Relax crystal water..\n"); fflush(stdout);
       relax_water(prot);
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
   
   /* repack */
   nowA = time(NULL);
   printf("   Repack side chains %d times, remove inaccessible conformers...\n", env.repacks);
   fflush(stdout);
   rot_pack(prot, env.repacks);
   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_repack[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_repack[prot.n_res] = c;
   nowB= time(NULL); i = nowB - nowA;
   printf("   Done. Time spent = %d\n\n", i); fflush(stdout);
   fflush(stdout);
   fp = fopen(ROTSTAT,"w");
   if (fp) {
       fprintf(fp,"   Rotamer making statitics:\n");
       write_step2stat(fp, prot, confstat);
       printf("\n");
       fclose(fp);
   }

   /* Relabel rotamer history */
   for (i=0; i<prot.n_res; i++) {
      for (j=prot.res[i].n_conf_ori; j<prot.res[i].n_conf; j++) {
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

   /* Write heavy atom rotamer pdb file, for recursive use of step 2 */
   if (!env.minimize_size) {
      if (!(fp=fopen(FN_HVROT, "w"))) {
         printf("   WARNING: Error in writing file %s\n\n", FN_HVROT);
      }
      else {
         write_pdb(fp, prot);
         fclose(fp);
      }
   }
   else {
      printf("   Skip writing heavy atome rotamer pdb file %s\n", FN_HVROT);
   }
   
   if (env.n_hv_conf_limit > 0) {
       printf("   Randomly prune conformers...\n"); fflush(stdout);
       if (env.rot_mhd_prune == 1) { //execution path that allows for MHD's pruning function
       		rand_conf_prune_mhd(prot);
       }
       else {
       		rand_conf_prune(prot);
       }
       printf("   Done\n\n"); fflush(stdout);
   }
   
   /* Make ionization states */
   printf("   Making ionization conformers...\n"); fflush(stdout);
   if (ionization(prot)) {
      printf("   FATAL: Fatal error reported by ionization()\n");
      return USERERR;
   }
   printf("   Done\n\n"); fflush(stdout);

   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_ioniz[i] = prot.res[i].n_conf-1;
      c += (prot.res[i].n_conf-1);
   }
   confstat.nc_ioniz[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   fprintf(fp,"   Rotamer making statitics:\n");
   write_step2stat(fp, prot, confstat);
   printf("\n");
   fclose(fp);
   
   /* add h */
   printf("   Add H atoms...\n"); fflush(stdout);
   while(place_missing(prot,1) > 0); rm_dupconf(prot, 0.001);
   assign_rad(prot);
   assign_vdw_param(prot);
   assign_crg(prot); /* reassign parameters before rm_dupconf(), Yifan -02/07/07 */
   printf("   Done\n\n"); fflush(stdout);

   c = prune_pv(prot, env.prune_rmsd, env.prune_ele, env.prune_vdw);
   if (c) printf("    %5d conformers deleted in this cycle at %8.3f %8.3f %8.3f\n", c, env.prune_rmsd, env.prune_ele, env.prune_vdw);
   
   //assign_vdw_param(prot);

   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_h[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_h[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   fprintf(fp,"   Rotamer making statitics:\n");
   write_step2stat(fp, prot, confstat);
   printf("\n");
   fclose(fp);

   /* heavy atom relaxation */
   if (env.hv_relax_ncycle > 0) {
       int i_res;
       nowA = time(NULL);

       printf("   Relaxation...\n"); fflush(stdout);
       relaxation(prot);
       /* add h to torsion minima again, this is due to protons are easily moved in relaxation program */
       for (i_res=0; i_res<prot.n_res; ++i_res) {
           int n_conf = prot.res[i_res].n_conf;
           int i_conf;
           for (i_conf=1; i_conf<n_conf; ++i_conf) {
               int i_atom;
               int k_conf = prot.res[i_res].n_conf;
               ins_conf(&prot.res[i_res], k_conf, prot.res[i_res].conf[i_conf].n_atom);
               cpy_conf(&prot.res[i_res].conf[k_conf], &prot.res[i_res].conf[i_conf]);
               prot.res[i_res].conf[k_conf].history[2] = 'T';
               
               for (i_atom=0; i_atom<prot.res[i_res].conf[k_conf].n_atom; i_atom++) {
                   if (prot.res[i_res].conf[k_conf].atom[i_atom].on) {
                       if (prot.res[i_res].conf[k_conf].atom[i_atom].name[1] == 'H') {
                           prot.res[i_res].conf[k_conf].atom[i_atom].on = 0;
                       }
                   }
               }
           }
       }
       
       /* remaking ionization due to different relaxation results on neutral and ionized forms, Yifan -01/25/07 */
       if (env.rot_mhd_prune == 1) { //execution path that allows for MHD's pruning function
       		if (ionization(prot)) {
       			printf("   FATAL: Fatal error reported by ionization()\n");
       			return USERERR;
       		}
       }
       
       assign_crg(prot); /* reassign parameters before rm_dupconf(), Yifan -02/07/07 */
       assign_rad(prot);
       rm_dupconf(prot, 0.001); /* identical coordinates */

       while(place_missing(prot,1) > 0); rm_dupconf(prot, 0.001);

       assign_crg(prot); /* reassign parameters before rm_dupconf(), Yifan -02/07/07 */
       assign_rad(prot);
       rm_dupconf(prot, 0.001); /* identical coordinates */
       
       nowB= time(NULL);
       i=nowB - nowA;
       printf("   Done. Time spent = %d\n\n", i); fflush(stdout);
   }

   /* optimize hydroxyl and water */
   if (env.relax_h) {
       nowA = time(NULL);
       printf("   Optimizing hydroxyl and water...\n");
       fflush(stdout);
       relax_h(prot);
       nowB= time(NULL);
       i=nowB - nowA;
       printf("   Done. Time spent = %d\n\n", i); fflush(stdout);
   }
   else {
      printf("   NOT optimizing hydroxyl and water.\n");
   }
   /*
   else {
       opt_h(prot);
   }
   */

   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_oh[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_oh[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   fprintf(fp,"   Rotamer making statitics:\n");
   write_step2stat(fp, prot, confstat);
   printf("\n");
   fclose(fp);

   printf("   Sorting conformers...\n"); fflush(stdout);
   sort_conf(prot);
   printf("   Done\n\n"); fflush(stdout);
   
   /* write full step2_out */
   assign_crg(prot);
   assign_rad(prot);
   if (!env.minimize_size) {
      fp = fopen("step2_out.full", "w");
      write_pdb(fp, prot);
      fclose(fp);
   }
   else {
      printf("   Skip writing full rotamer file step2_out.full.\n");
   }
   
   nowA = time(NULL);
   printf("   Delete duplicate conformers ...\n");
   fflush(stdout);

   c = prune_pv(prot, env.prune_rmsd, env.prune_ele, env.prune_vdw);
   if (c) printf("      %d conformers deleted in this cycle at %8.3f %8.3f %8.3f\n", c, env.prune_rmsd, env.prune_ele, env.prune_vdw);
   
   nowB= time(NULL);
   i=nowB - nowA;
   printf("   Done. Time spent = %d\n\n", i); fflush(stdout);
   c = 0;
   for (i=0; i<prot.n_res; i++) {
      confstat.nc_elect[i] = prot.res[i].n_conf-1;
      c+=(prot.res[i].n_conf-1);
   }
   confstat.nc_elect[prot.n_res] = c;
   fp = fopen(ROTSTAT,"w");
   fprintf(fp,"   Rotamer making statitics:\n");
   write_step2stat(fp, prot, confstat);
   printf("\n");
   fclose(fp);
   
   /* IPECE add membrane */
   if (env.ipece.add_mem) {
        printf("   Add membrane.\n");   fflush(stdout);
        add_membrane(&prot, &env.ipece);
        printf("   Done\n\n");   fflush(stdout);
   }
   
   /* write out */
   printf("   Write output...\n"); fflush(stdout);
   /* mcce pdb */
   fp = fopen(STEP2_OUT, "w");
   write_pdb(fp, prot);
   fclose(fp);
   /* conformer summary */
   fp = fopen(FN_CONFLIST2, "w");
   write_conflist(fp, prot);
   fclose(fp);

   /*
   get_demetri_out();
   */
   
   printf("   Done\n\n"); fflush(stdout);

   nowEnd = time(NULL);
   printf("   Total time for step2 (rotamer making) is %ld seconds.\n\n", nowEnd - nowStart);


   printf("   Output files:\n");
   printf("      %-16s: mcce pdb file, the input of step 3\n", STEP2_OUT);
   printf("      %-16s: conformer summary\n", FN_CONFLIST2);
   printf("      %-16s: rotamers without pairwise pruning.\n", "step2_out.full");
   printf("      %-16s: heavy atom rotamers, can be used recursively by step 2\n", FN_HVROT);
   printf("      %-16s: residue name list, can be used for step 4 to output microstate\n", "ms_gold");        // for ms_gold file, Cai
   printf("\n"); fflush(stdout);

   /* clean up memory */
   free(confstat.nc_start);
   free(confstat.nc_rmdup);
   free(confstat.nc_rotate);
   free(confstat.nc_vdw);
   free(confstat.nc_repack);
   free(confstat.nc_ioniz);
   free(confstat.nc_h);
   free(confstat.nc_oh);
   free(confstat.nc_elect);
   del_prot(&prot);

   return 0;
}

int delete_h(PROT prot)
{  int n = 0;
   int i, j, k;

   for (i=0; i<prot.n_res; i++) {
       char toggle[MAXCHAR_LINE];
       /* check DEL_HYDR parameter */
       if ( !param_get("DEL_HYDR", prot.res[i].resName, "", toggle) ) {
           if (strchr(toggle,'f')) continue;
           if (strchr(toggle,'F')) continue;
       }
       
       for (j=0; j<prot.res[i].n_conf; j++)
           for (k=0; k<prot.res[i].conf[j].n_atom; k++)
               if (prot.res[i].conf[j].atom[k].name[1] == 'H' && prot.res[i].conf[j].atom[k].on) {
                   prot.res[i].conf[j].atom[k].on = 0;
                   n++;
               }
   }
   return n;
}

int place_rot(PROT prot)
{  int i, i_conf;
   int C;
   char C_str[5];
   ROTAMER rule;
   char sbuff[MAXCHAR_LINE];

   /* construct rotamers */
   for (i=0; i<prot.n_res; i++) {
      if (!prot.res[i].do_rot) {
         printf("   Skip rotamer making for residue \"%s%04d%c\"\n", prot.res[i].resName,
                                                                                 prot.res[i].resSeq,
                                                                                 prot.res[i].chainID);
         continue;
      }
      if (prot.res[i].sas > env.sas_cutoff) {
         C = (int) (prot.res[i].rotations/2.0+0.5);
         if (C <= 0) C = 1;
         printf("   Reduce rotation steps from %2d to %2d for residue \"%s %c%04d\" as SAS=%3.f%%\n",  prot.res[i].rotations,
                                                                                 C,                                 
                                                                                 prot.res[i].resName,
                                                                                 prot.res[i].chainID,
                                                                                 prot.res[i].resSeq,
                                                                                 prot.res[i].sas*100.0);
         prot.res[i].rotations = C;
      }
      
      for (i_conf = 1; i_conf<prot.res[i].n_conf; i_conf++) {
          if (prot.res[i].conf[i_conf].history[2] == 'I' || prot.res[i].conf[i_conf].history[2] == 'B') {
              break;
          }
      }
      if (i_conf >= prot.res[i].n_conf) {
          /* if no rebuilt conformer is found */
          for (i_conf=1; i_conf<prot.res[i].n_conf; i_conf++) {
              if (prot.res[i].conf[i_conf].history[2] == 'O') {
                  break;
              }
          }
      }
      if (i_conf < prot.res[i].n_conf) {
          int ins;
          /* add a copy of the current conformer */
          ins = ins_conf(&prot.res[i], prot.res[i].n_conf, prot.res[i].conf[i_conf].n_atom);
          if (ins == USERERR) return USERERR;
          cpy_conf(&prot.res[i].conf[ins], &prot.res[i].conf[i_conf]);
          prot.res[i].conf[ins].history[2] = 'R';
          get_connect12_conf(i,ins,prot);
      }
      
      C = 0;
      while (1) {
         sprintf(C_str, "%d", C);
         if (param_get("ROTAMER", prot.res[i].resName, C_str, &rule)) break;
         sprintf(sbuff, " %s  ", rule.affected);
         sprintf(rule.affected, sbuff);
         if (place_rot_rule(i, rule, prot.res[i].rotations,prot)) {
            printf("   WARNING: place_rot(): \"failed placing rotamers for residue \"%s %d %c\"\"\n",
                   prot.res[i].resName, prot.res[i].resSeq, prot.res[i].chainID);
            printf("            No rotamers were made for this residue.\n");
         }
         C++;
      }
   }

   return 0;
}

int place_rot_rule(int i_res, ROTAMER rule, int n, PROT prot)
{  VECTOR v1, v2, v3;
   GEOM op;
   LINE axis;
   int n_conf;
   int ins;
   int i, j, k, l, k_conf;
   float phi;
   ATOM * atom1_p;
   ATOM * atom2_p;
   char found;
   int n_connected;
   ATOM **connected;
   RES *res = &prot.res[i_res];

   /* rotate */
   phi = 3.1415926*2.0/(float) n;
   n_conf = res->n_conf;

   for (j=1; j<n_conf; j++) { /* apply the rotation on passed in confs */
       /* search for the bond: atom2 is the second atom of the rotatable. This atom must
       * be in this residue. atom2->v2
       */
       
       if (prot.res[i_res].conf[j].history[2] != 'R') continue; /* only rotate 'R' conformers */
       
       if ((k=iatom(res->conf[j].confName, rule.atom2)) == -1) {
           
           /* search in the backbone */
           if ((k=iatom(res->conf[0].confName, rule.atom2)) == -1) {
               printf("   Error: place_rot_rule(): can't find atom \"%s\" in residue \"%s\"\n",
               rule.atom2, res->resName);
               //return USERERR;
               continue; /* jump to the next conformer, instead of exiting the subroutine - Yifan */
           }
           else {
               /* found in the backbone */
               atom2_p = &res->conf[0].atom[k];
           }
       }
       else {
           /* found in the sidechain */
           atom2_p = &res->conf[j].atom[k];
       }
       v2 = atom2_p->xyz;
       
       /* search for the bond: atom1 is the first atom of the bond, it is one of the connected
       * atoms of atom2, but not necessary in this conformer or this residue */
       found = 0;
       if ((k=iatom(res->conf[j].confName, rule.atom1)) == -1) {
           l=0;
           while (atom2_p->connect12[l]!=NULL) {
               if (!strcmp(atom2_p->connect12[l]->name, rule.atom1)) {
                   atom1_p = atom2_p->connect12[l];
                   found = 1;
                   break;  /* atom1 found */
               }
               l++;
           }
       }
       else {
           atom1_p = &res->conf[j].atom[k];
           found = 1;
       }
       
       if (found)  v1 = atom1_p->xyz;
       else {
           printf("   FATAL: place_rot_rule(): can't find atom \"%s\" connected to \"%s\" in residue \"%s\"\n",
           rule.atom2, rule.atom1, res->resName);
           //return USERERR;
           continue; /* jump to the next conformer, instead of exiting the subroutine - Yifan */
       }
       
       /* If "ALL_CONNECTED" is in the parameter, then search for the connected atoms */
       if (strstr(rule.affected, "ALL_CONNECTED")) {
           int i_connect;
           n_connected = 1;
           connected = malloc(sizeof(ATOM*));
           connected[0] = atom2_p;
           
           for (i_connect=0; i_connect<n_connected; i_connect++) {
               int j_connect;
               for (j_connect=0; j_connect<MAX_CONNECTED; j_connect++) {
                   ATOM *atom_in_network = connected[i_connect]->connect12[j_connect];
                   int k_connect;
                   if (!atom_in_network) break; /* end of connect12 array */
                   if (!atom_in_network->on) continue; /* empty atom slot skipped */
                   if (atom1_p == atom_in_network) continue; /* back connect to atom1 is disallowed */
                   
                   /* search the array to see if this atom is already in the array */
                   for (k_connect=0; k_connect<n_connected; k_connect++) {
                       if (connected[k_connect] == atom_in_network) {
                           break;
                       }
                   }
                   if (k_connect<n_connected) continue;
                   
                   /* add this connected atom into the array */
                   n_connected++;  /* n_connect changes in the loop!! The search will continue  */
                   connected = realloc(connected, n_connected*sizeof(ATOM *));
                   connected[n_connected-1] = atom_in_network;
               }
           }
       }
       
       /* rotatable bond is found, now use the coordinates to get the rotation axis */
       axis = line_2v(v1, v2);
       
       geom_reset(&op);
       for (i=0; i<n-1; i++) { /* rotate n-1 times */
           /* get rotation operator */
           geom_roll(&op, phi, axis);
           
           /* add a copy of the current conformer */
           ins = ins_conf(res, res->n_conf, res->conf[j].n_atom);
           if (ins == USERERR) return USERERR;
           cpy_conf(&res->conf[ins], &res->conf[j]);
           strncpy(res->conf[ins].history+2, "Ro", 2);
           get_connect12_conf(i_res,ins,prot);

           /* apply rotation operator to the atoms that match the parameter */
           if (strstr(rule.affected, "ALL_CONNECTED")) {
               int i_connect;
               for (i_connect=1; i_connect<n_connected; i_connect++) {
                   /* loop over all connected atoms, slot 0 is atom2 so no need to rotate */
                   v3 = connected[i_connect]->xyz;
                   geom_apply(op, &v3);
                   connected[i_connect]->xyz = v3;
               }
               /* a problem here: connected array points to the original copy,
               so the rotation is applied to the origianal conf */
               for (k=0; k<res->conf[j].n_atom; k++) {
                   /* swap the coordinates of conf j and the new one, conf ins */
                   v3 = res->conf[ins].atom[k].xyz;
                   res->conf[ins].atom[k].xyz = res->conf[j].atom[k].xyz;
                   res->conf[j].atom[k].xyz = v3;
               }
           }
           else if (strstr(rule.affected, "WHOLE_CONF")) {
               for (k=0; k<res->conf[j].n_atom; k++) {
                   v3 = res->conf[j].atom[k].xyz;
                   geom_apply(op, &v3);
                   res->conf[ins].atom[k].xyz = v3;
               }
           }
           else {
               for (k=0; k<res->conf[j].n_atom; k++) {
                   //printf("%s %s\n", rule.affected, res->conf[j].atom[k].name);
                   if (strstr(rule.affected, res->conf[j].atom[k].name)) {
                       //printf("%s\n", res->conf[j].atom[k].name);
                       v3 = res->conf[j].atom[k].xyz;
                       geom_apply(op, &v3);
                       res->conf[ins].atom[k].xyz = v3;
                   }
               }
           }
           
           /* check for duplication */
           for (k_conf=1; k_conf<res->n_conf-1;k_conf++) {
               if(!cmp_conf_hv(res->conf[res->n_conf-1], res->conf[k_conf], 0.2)) {
                   del_conf(res, res->n_conf-1);
                   break;
               }
           }
       }
       
       if (strstr(rule.affected, "ALL_CONNECTED")) {
           n_connected=0;
           free(connected);
       }
   }

   return 0;
}

int swing_rot(PROT prot)
{  int i, i_conf;
   int C;
   char C_str[5];
   ROTAMER rule;
   char sbuff[MAXCHAR_LINE];

   assign_vdw_param(prot);
   
   /* construct rotamers */
   for (i=0; i<prot.n_res; i++) {
      C = 0;
      if (!prot.res[i].do_sw) {
         printf("   Skip rotamer making for residue \"%s%04d%c\"\n", prot.res[i].resName,
                                                                                 prot.res[i].resSeq,
                                                                                 prot.res[i].chainID);
         continue;
      }
      if (prot.res[i].sas > env.sas_cutoff) {
         printf("   Skip rotamer making for residue \"%s%03d%c\" as SAS=%3.f%%\n", prot.res[i].resName,
                                                                                 prot.res[i].resSeq,
                                                                                 prot.res[i].chainID,
                                                                                 prot.res[i].sas*100.0);
         continue;
      }
      
      for (i_conf = 1; i_conf<prot.res[i].n_conf; i_conf++) {
          if (prot.res[i].conf[i_conf].history[2] == 'I' ||
              prot.res[i].conf[i_conf].history[2] == 'B' ||
              prot.res[i].conf[i_conf].history[2] == 'O') {
          
          int ins;
          /* add a copy of the current conformer */
          ins = ins_conf(&prot.res[i], prot.res[i].n_conf, prot.res[i].conf[i_conf].n_atom);
          if (ins == USERERR) return USERERR;
          cpy_conf(&prot.res[i].conf[ins], &prot.res[i].conf[i_conf]);
          prot.res[i].conf[ins].history[2] = 'S';
          get_connect12_conf(i,ins,prot);
              }
      }

      while (1) {
         sprintf(C_str, "%d", C);
         if (param_get("ROTAMER", prot.res[i].resName, C_str, &rule)) break;
         sprintf(sbuff, " %s  ", rule.affected);
         sprintf(rule.affected, sbuff);
         if (swing_rot_rule(i, rule, prot.res[i].phi_swing, prot)) {
            printf("   WARNING: swing_rot(): \"failed swinging rotamers for residue \"%s %d %c\"\"\n",
                   prot.res[i].resName, prot.res[i].resSeq, prot.res[i].chainID);
            printf("            No rotamers were made for this residue.\n");
         }
         C++;
      }
      
      /* prune by self energy here to save memory */
      if (prot.res[i].n_conf > 1000) {
          int kc;
          rm_dupconf_res(prot, i, 0.005);
          for (kc=0;kc<prot.res[i].n_conf;kc++) {
              get_connect12_conf(i,kc,prot);
          }
          get_vdw0_res_no_sas(i, prot);
          get_vdw1_res(i, prot);
          
          prune_by_vdw_res(i,prot, env.vdw_cutoff);
      }
      
   }
   return 0;
}

int swing_rot_rule(int i_res, ROTAMER rule, float phi, PROT prot)
{  VECTOR v1, v2, v3;
   GEOM op;
   LINE axis;
   int n_conf;
   int ins;
   int j, k, l;
   ATOM *atom1_p, *atom2_p;
   char found;
   int n_connected;
   ATOM **connected;
   RES *res = &prot.res[i_res];

   /* rotate */
   n_conf = res->n_conf;

   for (j=1; j<n_conf; j++) { /* apply the rotation on passed in confs */
       /* search for the bond: atom2 is the second atom of the rotatable. This atom must
       * be in this residue. atom2->v2
       */
       if (prot.res[i_res].conf[j].history[2] != 'S') continue; /* only rotate 'S' conformers */
       
       if ((k=iatom(res->conf[j].confName, rule.atom2)) == -1) {
           
           /* search in the backbone */
           if ((k=iatom(res->conf[0].confName, rule.atom2)) == -1) {
               printf("   Error: place_rot_rule(): can't find atom \"%s\" in residue \"%s\"\n",
               rule.atom2, res->resName);
               //return USERERR;
               continue; /* jump to the next conformer, instead of exiting the subroutine - Yifan */
           }
           else {
               /* found in the backbone */
               atom2_p = &res->conf[0].atom[k];
           }
       }
       else {
           /* found in the sidechain */
           atom2_p = &res->conf[j].atom[k];
       }
       v2 = atom2_p->xyz;

       /* search for the bond: atom1 is the first atom of the bond, it is one of the connected
       * atoms of atom2, but not necessary in this conformer or this residue */
       found = 0;
       if ((k=iatom(res->conf[j].confName, rule.atom1)) == -1) {
           l=0;
           while (atom2_p->connect12[l]!=NULL) {
               if (!strcmp(atom2_p->connect12[l]->name, rule.atom1)) {
                   atom1_p = atom2_p->connect12[l];
                   found = 1;
                   break;  /* atom1 found */
               }
               l++;
           }
       }
       else {
           atom1_p = &res->conf[j].atom[k];
           found = 1;
       }
       
       if (found)  v1 = atom1_p->xyz;
       else {
           printf("   FATAL: place_rot_rule(): can't find atom \"%s\" connected to \"%s\" in residue \"%s\"\n",
           rule.atom2, rule.atom1, res->resName);
           //return USERERR;
           continue; /* jump to the next conformer, instead of exiting the subroutine - Yifan */
       }

       /* If "ALL_CONNECTED" is in the parameter, then search for the connected atoms */
       if (strstr(rule.affected, "ALL_CONNECTED")) {
           int i_connect;
           n_connected = 1;
           connected = malloc(sizeof(ATOM*));
           connected[0] = atom2_p;
           
           for (i_connect=0; i_connect<n_connected; i_connect++) {
               int j_connect;
               for (j_connect=0; j_connect<MAX_CONNECTED; j_connect++) {
                   ATOM *atom_in_network = connected[i_connect]->connect12[j_connect];
                   int k_connect;
                   if (!atom_in_network) break; /* end of connect12 array */
                   if (!atom_in_network->on) continue; /* empty atom slot skipped */
                   if (atom1_p == atom_in_network) continue; /* back connect to atom1 is disallowed */
                   
                   /* search the array to see if this atom is already in the array */
                   for (k_connect=0; k_connect<n_connected; k_connect++) {
                       if (connected[k_connect] == atom_in_network) {
                           break;
                       }
                   }
                   if (k_connect<n_connected) continue;
                   
                   /* add this connected atom into the array */
                   n_connected++;  /* n_connect changes in the loop!! The search will continue  */
                   connected = realloc(connected, n_connected*sizeof(ATOM *));
                   connected[n_connected-1] = atom_in_network;
               }
           }
       }
       
       /* rotatable bond is found, now use the coordinates to get the rotation axis */
       axis = line_2v(v1, v2);

      /* swing left */
      geom_reset(&op);
      geom_roll(&op, -phi, axis);
      
      /* add a copy of the current conformer */
      ins = ins_conf(res, res->n_conf, res->conf[j].n_atom);
      if (ins == USERERR) return USERERR;
      cpy_conf(&res->conf[ins], &res->conf[j]);
      strncpy(res->conf[ins].history+4, "Sw", 2);
      get_connect12_conf(i_res,ins,prot);

      /* apply rotation operator to the atoms that match the parameter */
      if (strstr(rule.affected, "ALL_CONNECTED")) {
          int i_connect;
          for (i_connect=1; i_connect<n_connected; i_connect++) {
              /* loop over all connected atoms, slot 0 is atom2 so no need to rotate */
              v3 = connected[i_connect]->xyz;
              geom_apply(op, &v3);
              connected[i_connect]->xyz = v3;
          }
          /* a problem here: connected array points to the original copy,
          so the rotation is applied to the origianal conf */
          for (k=0; k<res->conf[j].n_atom; k++) {
              /* swap conf j and the new one, conf ins */
              v3 = res->conf[ins].atom[k].xyz;
              res->conf[ins].atom[k].xyz = res->conf[j].atom[k].xyz;
              res->conf[j].atom[k].xyz = v3;
          }
      }
      else if (strstr(rule.affected, "WHOLE_CONF")) {
          for (k=0; k<res->conf[j].n_atom; k++) {
              v3 = res->conf[j].atom[k].xyz;
              geom_apply(op, &v3);
              res->conf[ins].atom[k].xyz = v3;
          }
      }
      else {
          for (k=0; k<res->conf[j].n_atom; k++) {
              if (strstr(rule.affected, res->conf[j].atom[k].name)) {
                  v3 = res->conf[j].atom[k].xyz;
                  geom_apply(op, &v3);
                  res->conf[ins].atom[k].xyz = v3;
              }
          }
      }

      /* swing right */
      geom_reset(&op);
      geom_roll(&op, phi, axis);
      
      /* add a copy of the current conformer */
      ins = ins_conf(res, res->n_conf, res->conf[j].n_atom);
      if (ins == USERERR) return USERERR;
      cpy_conf(&res->conf[ins], &res->conf[j]);
      strncpy(res->conf[ins].history+4, "Sw", 2);
      get_connect12_conf(i_res,ins,prot);
      
      /* apply rotation operator to the atoms that match the parameter */
      if (strstr(rule.affected, "ALL_CONNECTED")) {
          int i_connect;
          for (i_connect=1; i_connect<n_connected; i_connect++) {
              /* loop over all connected atoms, slot 0 is atom2 so no need to rotate */
              v3 = connected[i_connect]->xyz;
              geom_apply(op, &v3);
              connected[i_connect]->xyz = v3;
          }
          /* a problem here: connected array points to the original copy,
          so the rotation is applied to the origianal conf */
          for (k=0; k<res->conf[j].n_atom; k++) {
              /* swap conf j and the new one, conf ins */
              v3 = res->conf[ins].atom[k].xyz;
              res->conf[ins].atom[k].xyz = res->conf[j].atom[k].xyz;
              res->conf[j].atom[k].xyz = v3;
          }
      }
      else if (strstr(rule.affected, "WHOLE_CONF")) {
          for (k=0; k<res->conf[j].n_atom; k++) {
              v3 = res->conf[j].atom[k].xyz;
              geom_apply(op, &v3);
              res->conf[ins].atom[k].xyz = v3;
          }
      }
      else {
          for (k=0; k<res->conf[j].n_atom; k++) {
              if (strstr(rule.affected, res->conf[j].atom[k].name)) {
                  v3 = res->conf[j].atom[k].xyz;
                  geom_apply(op, &v3);
                  res->conf[ins].atom[k].xyz = v3;
              }
          }
      }
      if (strstr(rule.affected, "ALL_CONNECTED")) {
          n_connected=0;
          free(connected);
      }

   }

   return 0;
}

int extra_rot(PROT prot)
{
    /* this subroutine creates additional rotamers by translation and rotation,
    designed for ligand binding, such as quinones */
    char toggle[MAXCHAR_LINE];
    int i_res;
    
    /* translations */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        RES *res_p = &prot.res[i_res];
        int i_trans;
        
        /* check if this residue will be translated */
        if ( !param_get("TRANS", prot.res[i_res].resName, "", toggle) ) {
            if (strchr(toggle,'f')) continue;
            if (strchr(toggle,'F')) continue;
        }
        else {
            continue;
        }
        
        /* loop over number of translation steps, defined in run.prm */
        for (i_trans = 0; i_trans < env.n_trans; i_trans++) {
            int n_conf, i_conf;
            
            /* loop over all existing conformers */
            n_conf = prot.res[i_res].n_conf;
            for (i_conf = 1; i_conf < n_conf; i_conf++) {
                int i_direc, j_direc, k_direc;
                
                /* move in three directions */
                for (i_direc = -1; i_direc <= 1; i_direc++) {
                    for (j_direc = -1; j_direc <= 1; j_direc++) {
                        for (k_direc = -1; k_direc <= 1; k_direc++) {
                            int ins, i_atom, k_conf;
                            
                            if (!(i_direc||j_direc||k_direc)) continue; /* i,j,k all 0, no movement */
                            
                            /* add a copy of the current conformer */
                            ins = ins_conf(res_p, res_p->n_conf, res_p->conf[i_conf].n_atom);
                            if (ins == USERERR) {
                                printf("   FATAL! Can't add conformer. Probably out of memory.\n");
                                return USERERR;
                            }
                            cpy_conf(&res_p->conf[ins], &res_p->conf[i_conf]);
                            strncpy(res_p->conf[ins].history+2, "Tr", 2);
                            
                            /* move the new conformer */
                            for (i_atom=0; i_atom<res_p->conf[ins].n_atom; i_atom++) {
                                if (!res_p->conf[ins].atom[i_atom].on) continue;
                                res_p->conf[ins].atom[i_atom].xyz.x += env.trans_dist*(float)i_direc;
                                res_p->conf[ins].atom[i_atom].xyz.y += env.trans_dist*(float)j_direc;
                                res_p->conf[ins].atom[i_atom].xyz.z += env.trans_dist*(float)k_direc;
                            }
                            
                            /* check for duplication */
                            for (k_conf=1; k_conf<ins; k_conf++) {
                                if(!cmp_conf_hv(res_p->conf[ins], res_p->conf[k_conf], 0.05)) {
                                    del_conf(res_p, ins);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* rotation (using different parameter from rotamer making subroutine) */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        int  counter;
        char counter_str[5];
        char rule[MAXCHAR_LINE];
        RES * res_p = &prot.res[i_res];
        
        counter = 0;
        while (1) {
            char atom1[5], atom2[5], atom3[5];
            int i_conf;
            int n_conf = prot.res[i_res].n_conf;
            
            sprintf(counter_str, "%d", counter);
            if (param_get("SLIDE", prot.res[i_res].resName, counter_str, rule)) break;
            counter++;
            
            while (strlen(rule)<15) strcat(rule, " ");
            
            strncpy(atom1, rule,    4); atom1[4] = '\0';
            strncpy(atom2, rule+5,  4); atom2[4] = '\0';
            strncpy(atom3, rule+10, 4); atom3[4] = '\0';
            
            for (i_conf = 1; i_conf<n_conf; i_conf++) {
                int ins, i_atom, k_atom1, k_atom2, k_atom3;
                PLANE plane;
                VECTOR v0, v1;
                LINE axis;
                GEOM op;
                
                if ((k_atom1 = iatom(res_p->conf[i_conf].confName, atom1)) == -1) continue;
                if ((k_atom2 = iatom(res_p->conf[i_conf].confName, atom2)) == -1) continue;
                if ((k_atom3 = iatom(res_p->conf[i_conf].confName, atom3)) == -1) continue;
                
                plane = plane_3v(res_p->conf[i_conf].atom[k_atom1].xyz,res_p->conf[i_conf].atom[k_atom2].xyz,res_p->conf[i_conf].atom[k_atom3].xyz);
                
                /* v0 is the position of atom1 */
                v0 = res_p->conf[i_conf].atom[k_atom1].xyz;
                v1 = vector_vplusv(v0,plane.t);
                
                axis = line_2v(v0, v1);
                
                if (res_p->do_rot) {
                    int i_rot;
                    float phi = env.PI*2.0 /(float) res_p->rotations;
                    geom_reset(&op);
                    
                    for (i_rot=0; i_rot<res_p->rotations-1; i_rot++) { /* rotate n-1 times */
                        /* get rotation operator */
                        /* op is not reset after each rotation step, therefore kept rotating */
                        geom_roll(&op, phi, axis);
                        
                        /* add a copy of the current conformer */
                        ins = ins_conf(res_p, res_p->n_conf, res_p->conf[i_conf].n_atom);
                        if (ins == USERERR) return USERERR;
                        cpy_conf(&res_p->conf[ins], &res_p->conf[i_conf]);
                        strncpy(res_p->conf[ins].history+4, "Ro", 2);
                        
                        /* apply rotation to each atom */
                        for (i_atom = 0; i_atom < res_p->conf[i_conf].n_atom; i_atom++) {
                            if (!res_p->conf[ins].atom[i_atom].on) continue;
                            geom_apply(op, &res_p->conf[ins].atom[i_atom].xyz);
                        }
                    }
                }
                if (res_p->do_sw) {
                    /* rotate in one direction */
                    geom_reset(&op);
                    geom_roll(&op, res_p->phi_swing, axis);
                    
                    /* add a copy of the current conformer */
                    ins = ins_conf(res_p, res_p->n_conf, res_p->conf[i_conf].n_atom);
                    if (ins == USERERR) return USERERR;
                    cpy_conf(&res_p->conf[ins], &res_p->conf[i_conf]);
                    strncpy(res_p->conf[ins].history+4, "Sw", 2);
                    
                    /* apply rotation to each atom */
                    for (i_atom = 0; i_atom < res_p->conf[i_conf].n_atom; i_atom++) {
                        if (!res_p->conf[ins].atom[i_atom].on) continue;
                        geom_apply(op, &res_p->conf[ins].atom[i_atom].xyz);
                    }
                    
                    /* rotate in the other direction */
                    geom_reset(&op);
                    geom_roll(&op, -res_p->phi_swing, axis);
                    
                    /* add a copy of the current conformer */
                    ins = ins_conf(res_p, res_p->n_conf, res_p->conf[i_conf].n_atom);
                    if (ins == USERERR) return USERERR;
                    cpy_conf(&res_p->conf[ins], &res_p->conf[i_conf]);
                    strncpy(res_p->conf[ins].history+4, "Sw", 2);
                    
                    /* apply rotation to each atom */
                    for (i_atom = 0; i_atom < res_p->conf[i_conf].n_atom; i_atom++) {
                        if (!res_p->conf[ins].atom[i_atom].on) continue;
                        geom_apply(op, &res_p->conf[ins].atom[i_atom].xyz);
                    }
                    
                }
            }
            
        }
        /* same as SLIDE, except using the midpoint of first two atoms for axis */
        while (1) {
            char atom1[5], atom2[5], atom3[5];
            int i_conf;
            int n_conf = prot.res[i_res].n_conf;
            
            sprintf(counter_str, "%d", counter);
            if (param_get("SPIN", prot.res[i_res].resName, counter_str, rule)) break;
            counter++;
            
            while (strlen(rule)<15) strcat(rule, " ");
            
            strncpy(atom1, rule,    4); atom1[4] = '\0';
            strncpy(atom2, rule+5,  4); atom2[4] = '\0';
            strncpy(atom3, rule+10, 4); atom3[4] = '\0';
            
            for (i_conf = 1; i_conf<n_conf; i_conf++) {
                int ins, i_atom, k_atom1, k_atom2, k_atom3;
                PLANE plane;
                VECTOR v0, v1;
                LINE axis;
                GEOM op;
                
                if ((k_atom1 = iatom(res_p->conf[i_conf].confName, atom1)) == -1) continue;
                if ((k_atom2 = iatom(res_p->conf[i_conf].confName, atom2)) == -1) continue;
                if ((k_atom3 = iatom(res_p->conf[i_conf].confName, atom3)) == -1) continue;
                
                plane = plane_3v(res_p->conf[i_conf].atom[k_atom1].xyz,res_p->conf[i_conf].atom[k_atom2].xyz,res_p->conf[i_conf].atom[k_atom3].xyz);
                
                /* v0 is the midpoint if atom1 and atom2 */
                v0 = vector_rescale(vector_vplusv(res_p->conf[i_conf].atom[k_atom1].xyz, res_p->conf[i_conf].atom[k_atom2].xyz), 0.5);
                v1 = vector_vplusv(v0,plane.t);
                
                axis = line_2v(v0, v1);
                
                if (res_p->do_rot) {
                    int i_rot;
                    float phi = env.PI*2.0 /(float) res_p->rotations;
                    geom_reset(&op);
                    
                    for (i_rot=0; i_rot<res_p->rotations-1; i_rot++) { /* rotate n-1 times */
                        /* get rotation operator */
                        /* op is not reset after each rotation step, therefore kept rotating */
                        geom_roll(&op, phi, axis);
                        
                        /* add a copy of the current conformer */
                        ins = ins_conf(res_p, res_p->n_conf, res_p->conf[i_conf].n_atom);
                        if (ins == USERERR) return USERERR;
                        cpy_conf(&res_p->conf[ins], &res_p->conf[i_conf]);
                        strncpy(res_p->conf[ins].history+4, "Ro", 2);
                        
                        /* apply rotation to each atom */
                        for (i_atom = 0; i_atom < res_p->conf[i_conf].n_atom; i_atom++) {
                            if (!res_p->conf[ins].atom[i_atom].on) continue;
                            geom_apply(op, &res_p->conf[ins].atom[i_atom].xyz);
                        }
                    }
                }
                if (res_p->do_sw) {
                    /* rotate in one direction */
                    geom_reset(&op);
                    geom_roll(&op, res_p->phi_swing, axis);
                    
                    /* add a copy of the current conformer */
                    ins = ins_conf(res_p, res_p->n_conf, res_p->conf[i_conf].n_atom);
                    if (ins == USERERR) return USERERR;
                    cpy_conf(&res_p->conf[ins], &res_p->conf[i_conf]);
                    strncpy(res_p->conf[ins].history+4, "Sw", 2);
                    
                    /* apply rotation to each atom */
                    for (i_atom = 0; i_atom < res_p->conf[i_conf].n_atom; i_atom++) {
                        if (!res_p->conf[ins].atom[i_atom].on) continue;
                        geom_apply(op, &res_p->conf[ins].atom[i_atom].xyz);
                    }
                    
                    /* rotate in the other direction */
                    geom_reset(&op);
                    geom_roll(&op, -res_p->phi_swing, axis);
                    
                    /* add a copy of the current conformer */
                    ins = ins_conf(res_p, res_p->n_conf, res_p->conf[i_conf].n_atom);
                    if (ins == USERERR) return USERERR;
                    cpy_conf(&res_p->conf[ins], &res_p->conf[i_conf]);
                    strncpy(res_p->conf[ins].history+4, "Sw", 2);
                    
                    /* apply rotation to each atom */
                    for (i_atom = 0; i_atom < res_p->conf[i_conf].n_atom; i_atom++) {
                        if (!res_p->conf[ins].atom[i_atom].on) continue;
                        geom_apply(op, &res_p->conf[ins].atom[i_atom].xyz);
                    }
                    
                }
                
            }
            
        }
    }
    return 0;
}

void write_step2stat(FILE *fp, PROT prot, CONFSTAT stat)
{  int i;

   fprintf(fp, "   Residue   Start  Clean   Swap Rotate   Self  Hbond Repack  Ioni.   TorH     OH  Elect\n");
   for (i=0; i<prot.n_res; i++) {
      fprintf(fp, "   %3s%c%04d%7d%7d%7d%7d%7d%7d%7d%7d%7d%7d%7d\n",
                  prot.res[i].resName,
                  prot.res[i].chainID,
                  prot.res[i].resSeq,
                  stat.nc_start[i],
                  stat.nc_rmdup[i],
                  stat.nc_swap[i],
                  stat.nc_rotate[i],
                  stat.nc_vdw[i],
                  stat.nc_hdir[i],
                  stat.nc_repack[i],
                  stat.nc_ioniz[i],
                  stat.nc_h[i],
                  stat.nc_oh[i],
                  stat.nc_elect[i]);
   }
   fprintf(fp, "   Total   %7d%7d%7d%7d%7d%7d%7d%7d%7d%7d%7d\n",
               stat.nc_start[i],
               stat.nc_rmdup[i],
               stat.nc_swap[i],
               stat.nc_rotate[i],
               stat.nc_vdw[i],
               stat.nc_hdir[i],
               stat.nc_repack[i],
               stat.nc_ioniz[i],
               stat.nc_h[i],
               stat.nc_oh[i],
               stat.nc_elect[i]);

   return;
}

int prune_by_vdw(PROT prot, float delta_E)
{
    int n=0;
    int kr;
    for (kr=0; kr<prot.n_res; kr++) {
        n += prune_by_vdw_res(kr, prot, delta_E);
    }
    return n;
}

int prune_by_vdw_res(int kr, PROT prot, float delta_E)
{
    int n = 0;
    int kc;
    float E_low;
    
    if (prot.res[kr].n_conf > 1) E_low = prot.res[kr].conf[1].E_self;
    
    for (kc=1; kc<prot.res[kr].n_conf; kc++) {
        if (E_low>prot.res[kr].conf[kc].E_self) E_low = prot.res[kr].conf[kc].E_self;
    }
    
    /* keep all conformers that have favorable self energy */
    if (E_low < 0) E_low = 0;
    
    for (kc=2; kc<prot.res[kr].n_conf; kc++) {
        if (prot.res[kr].conf[kc].E_self - E_low > delta_E) {
            del_conf(&prot.res[kr], kc);
            kc--;
            n++;
        }
    }
    
    return n;
}

int rot_pack(PROT prot, int n)
{
    FILE *fp;
    int C, i, j;
    MICROSTATE state;
    int flips;
    float ****pairwise;
    int i_res,i_conf,j_res,j_conf;
    float pair_vdw;
    int counter;
    char pipe;
    float cutoff_far2  = VDW_CUTOFF_FAR  * VDW_CUTOFF_FAR;

    if (n ==0) return 0;
    for (i_res=0; i_res< prot.n_res; i_res++) {
        STRINGS confs;
        int n_atom, ins;
        if (param_get("CONFLIST", prot.res[i_res].resName, "", &confs)) {
            continue;
        }
        for (i_conf=0; i_conf<confs.n; i_conf++) {
            if (!strncmp(confs.strings[i_conf]+3,"BK",2)) continue;
            if (param_get("NATOM", confs.strings[i_conf], "", &n_atom)) {
                continue;
            }
            if (!n_atom) {
               ins = ins_conf(&prot.res[i_res], prot.res[i_res].n_conf, n_atom);
               strcpy(prot.res[i_res].conf[ins].confName, confs.strings[i_conf]);
               strcpy(prot.res[i_res].conf[ins].history,  confs.strings[i_conf]+3);
               strcpy(prot.res[i_res].conf[ins].history+2,  "________");
            }
        }
    }

    /* preparing energy lookup table */
    while(place_missing(prot,1) > 0); rm_dupconf(prot, 0.001);
    del_non_common_h(prot);
    rm_dupconf(prot, 0.001);
    assign_vdw_param(prot);
    assign_rad(prot);
    counter = 0;
    for (i=0; i<prot.n_res; i++)  counter+=(prot.res[i].n_conf-1);
    printf("   Computing pairwise LJ potential. This may take a while.\n");
    i= (int) counter*counter*3.5E-5;
    printf("      Estimated time on AMD 1.6GHz is %d seconds.\n",i);
    fflush(stdout);

    /* timing begins, estimate T from counter */
    nowA = time(NULL);
    
    /* increase vdw radii of carbon atoms 
    for (i_res = 0; i_res < prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                if (atom_p->name[1] != 'C') continue;
                atom_p->vdw_rad += 0.5;
            }
        }
    }
    */
    
    get_vdw0(prot);
    get_vdw1(prot);
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            prot.res[i_res].conf[i_conf].E_torsion = torsion_conf(&prot.res[i_res].conf[i_conf]);
            prot.res[i_res].conf[i_conf].E_self = prot.res[i_res].conf[i_conf].E_vdw0
            + prot.res[i_res].conf[i_conf].E_vdw1 + prot.res[i_res].conf[i_conf].E_torsion;
            
            /* ignore favorable energies */
            if (prot.res[i_res].conf[i_conf].E_self < 0) prot.res[i_res].conf[i_conf].E_self = 0.;
        }
    }
    
    /* Setup for fast vdw calculations */
    setup_vdw_fast(prot);

    /* setup connectivity table */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        setup_connect_res(prot,i_res);
    }

    pairwise = calloc(prot.n_res, sizeof(void *));
    /* initialize ngh list */
    for (i_res = 0; i_res < prot.n_res; i_res++) {
        if (prot.res[i_res].n_ngh) free(prot.res[i_res].ngh);
        prot.res[i_res].n_ngh = 0;
        prot.res[i_res].ngh = NULL;
    }

    for (i_res=0; i_res<prot.n_res; i_res++) {
        /* lower half of the matrix has already been calculated */
        for (j_res=0; j_res<i_res; j_res++) {
            int j_ngh;
            /* search if i_res is already in ngh list of j_res */
            for (j_ngh=0; j_ngh<prot.res[j_res].n_ngh; j_ngh++) {
                if (i_res == prot.res[j_res].ngh[j_ngh]->i_res_prot) {
                    int n_ngh;
                    /* add j_res into ngh list */
                    prot.res[i_res].n_ngh++;
                    n_ngh = prot.res[i_res].n_ngh;
                    prot.res[i_res].ngh = realloc(prot.res[i_res].ngh, n_ngh*sizeof(void *));
                    prot.res[i_res].ngh[n_ngh-1] = &prot.res[j_res];
                    /* pointing pairwise(i,j) to pairwise(j,i), which would be transposed for i_res */
                    pairwise[i_res] = realloc(pairwise[i_res], n_ngh*sizeof(void *));
                    pairwise[i_res][n_ngh-1] = pairwise[j_res][j_ngh];
                    break;
                }
            }
        }

        /* upper half of the matrix */
        for (j_res=i_res+1; j_res<prot.n_res; j_res++) {
            int n_ngh,all_small;
            /* check if i_res and j_res are close enough. */
            if (out_of_range(prot.res[i_res].r_min, prot.res[i_res].r_max,
                prot.res[j_res].r_min, prot.res[j_res].r_max, cutoff_far2)) continue;

            /* add j_res into ngh list */
            prot.res[i_res].n_ngh++;
            n_ngh = prot.res[i_res].n_ngh;
            prot.res[i_res].ngh = realloc(prot.res[i_res].ngh, n_ngh*sizeof(void *));
            prot.res[i_res].ngh[n_ngh-1] = &prot.res[j_res];
            
            /* memory for vdw matrix */
            pairwise[i_res] = realloc(pairwise[i_res], n_ngh*sizeof(void *));
            pairwise[i_res][n_ngh-1] = malloc(prot.res[i_res].n_conf * sizeof(void *));
            if (!pairwise[i_res][n_ngh-1]) {printf("Memory Error\n"); return USERERR;}
            for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                pairwise[i_res][n_ngh-1][i_conf] = malloc(prot.res[j_res].n_conf * sizeof(float));
                if (!pairwise[i_res][n_ngh-1][i_conf]) {printf("Memory Error\n"); return USERERR;}
            }

            /* calculate vdw matrix for (i_res,j_res) */
            all_small = 1;
            for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                if (!prot.res[i_res].conf[i_conf].n_atom) continue;
                for (j_conf=1; j_conf<prot.res[j_res].n_conf; j_conf++) {
                    if (!prot.res[j_res].conf[j_conf].n_atom) continue;
                    
                    pair_vdw = vdw_conf_fast(i_res,i_conf,j_res,j_conf,prot,0);

                    /* ignore favorable energies */
                    if (env.repack_fav_vdw_off == 1) {
                        if (pair_vdw < 0.) pair_vdw = 0.;
                    }
                    
                    /* H bond energy correction */
                    /* turned off, low energy from h bond traps conformer in certain positions, and loses occupancy for other (exposed) conformers -Yifan */
                    //pair_vdw += hbond_extra(prot.res[i_res].conf[i_conf], prot.res[j_res].conf[j_conf]);

                    pairwise[i_res][n_ngh-1][i_conf][j_conf] = pair_vdw;
                    
                    if (fabs(pair_vdw) > env.ngh_vdw_thr) {
                        all_small = 0;
                    }
                }
            }
            /* if all pairwise are smaller than threshold, remove j_res from neighbor list */
            if (all_small) {
                for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) free(pairwise[i_res][n_ngh-1][i_conf]);
                free(pairwise[i_res][n_ngh-1]);
                prot.res[i_res].n_ngh--;
                n_ngh = prot.res[i_res].n_ngh;
                prot.res[i_res].ngh = realloc(prot.res[i_res].ngh, n_ngh*sizeof(void *));
                pairwise[i_res] = realloc(pairwise[i_res], n_ngh*sizeof(void *));
            }
        }

        /*
        printf("%s %c%4d, n_ngh=%3d, sas=%8.2f%%\n",
        prot.res[i_res].resName,prot.res[i_res].chainID,prot.res[i_res].resSeq,prot.res[i_res].n_ngh,100.*prot.res[i_res].sas);
        */
        
    }
    
    /* free the memory used by connectivity table */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        free_connect_res(prot,i_res);
    }

    /*
    memset(vdw_counter,0,21*sizeof(int));
    vdw_stat = fopen("vdw_stat.dat","w");
    for (i_res=0; i_res<prot.n_res; i_res++) {
        if (!(prot.res[i_res].n_conf-1)) continue;
        for (j_res=0; j_res<prot.n_res; j_res++) {
            if (!(prot.res[j_res].n_conf-1)) continue;
            i_counter = 20.*(float)stat[i_res][j_res]/(float)((prot.res[i_res].n_conf-1)*(prot.res[j_res].n_conf-1));
            vdw_counter[i_counter]++;
            fprintf(vdw_stat, "%4d %4d %8d %8d %8.2f\n",
            i_res, j_res,
            stat[i_res][j_res],
            (prot.res[i_res].n_conf-1)*(prot.res[j_res].n_conf-1),
            100.*(float)stat[i_res][j_res]/(float)((prot.res[i_res].n_conf-1)*(prot.res[j_res].n_conf-1)));
        }
    }
    for (i_counter=0;i_counter<21;i_counter++) printf("%d\n", vdw_counter[i_counter]);
    */
    
    nowB= time(NULL);
    i=nowB - nowA;
    printf("      Actual time is %d\n", i);
    /* timing end  */

    /* prepare a radom microstate */
    state.n = prot.n_res;
    state.res = (int *) malloc(prot.n_res * sizeof(int));

    for (i=0; i<prot.n_res; i++) {
        for (j=1; j<prot.res[i].n_conf; j++) {
            prot.res[i].conf[j].counter = 0;
            prot.res[i].conf[j].on      = 0;
        }
    }


    printf("   Repacking in progress, see %s for details...\n", env.progress_log);fflush(stdout);

    pipe = 1;
    if (!(fp=(fopen(env.progress_log, "a")))) {
       printf("   WARNING: can not open file %s to dump messages, use stdout instead\n", env.progress_log);
       fp = stdout;
       pipe = 0;
    }
    fprintf(fp, "   Repacking:\n");fflush(fp);
    for (i=0; i<n; i++) {
        /* initialize the microstate */
        for (j=0; j<prot.n_res; j++) {
            if (prot.res[j].n_conf == 0) {
                printf("   FATAL: rot_pack(): empty residue %s%c%04d.\n", prot.res[j].resName,prot.res[j].chainID,prot.res[j].resSeq);
                return USERERR;
            }
            else if (prot.res[j].n_conf <= 2) state.res[j] = prot.res[j].n_conf - 1;
            else state.res[j] = rand() % (prot.res[j].n_conf-1) + 1;
            if (state.res[j])
                prot.res[j].conf[state.res[j]].on = 1;
        }

        fprintf(fp, "   Repacking %d:", i); fflush(fp);
        for (j=0; j<100; j++) { /* maximum converge trials */
            /* refine the microstate */
            flips = rot_refine(prot, state, pairwise);
            fprintf(fp, " %d", flips); fflush(fp);
            if (flips == 0) break;
        }
        fprintf(fp, "\n"); fflush(fp);

        /* write out the microstate as a pdbfile, for marilyn's presentation
        {  int mi, mk, mc;
            mc = 0;
            FILE *stream;
            char sbuff[256];

            sprintf(sbuff,"state%03d.pdb", i);
            stream = fopen(sbuff, "w");

            for (mi=0; mi<prot.n_res; mi++) {
                for (mk=0; mk<prot.res[mi].conf[0].n_atom; mk++) {
                    if (!prot.res[mi].conf[0].atom[mk].on) continue;
                    mc++;
                    fprintf(stream, "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f            %5s\n",
                    mc, prot.res[mi].conf[0].atom[mk].name,
                    prot.res[mi].conf[0].atom[mk].altLoc,
                    prot.res[mi].conf[0].atom[mk].resName,
                    prot.res[mi].conf[0].atom[mk].chainID,
                    prot.res[mi].conf[0].atom[mk].resSeq,
                    prot.res[mi].conf[0].atom[mk].iCode,
                    prot.res[mi].conf[0].atom[mk].xyz.x,
                    prot.res[mi].conf[0].atom[mk].xyz.y,
                    prot.res[mi].conf[0].atom[mk].xyz.z,
                    prot.res[mi].conf[0].confName);
                }

                if (prot.res[mi].n_conf < 2) continue;
                for (mk=0; mk<prot.res[mi].conf[state.res[mi]].n_atom; mk++) {
                    if (!prot.res[mi].conf[state.res[mi]].atom[mk].on) continue;
                    mc++;
                    fprintf(stream, "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f            %5s\n",
                    mc, prot.res[mi].conf[state.res[mi]].atom[mk].name,
                    prot.res[mi].conf[state.res[mi]].atom[mk].altLoc,
                    prot.res[mi].conf[state.res[mi]].atom[mk].resName,
                    prot.res[mi].conf[state.res[mi]].atom[mk].chainID,
                    prot.res[mi].conf[state.res[mi]].atom[mk].resSeq,
                    prot.res[mi].conf[state.res[mi]].atom[mk].iCode,
                    prot.res[mi].conf[state.res[mi]].atom[mk].xyz.x,
                    prot.res[mi].conf[state.res[mi]].atom[mk].xyz.y,
                    prot.res[mi].conf[state.res[mi]].atom[mk].xyz.z,
                    prot.res[mi].conf[state.res[mi]].confName);
                }
            }
            fclose(stream);
        }
        */

        /* increase counters */
        for (j=0; j<prot.n_res; j++) {
            prot.res[j].conf[state.res[j]].counter ++;
        }
    }
    if (pipe) fclose(fp);

    for (i=0; i<prot.n_res; i++) {
        float total_occ = 0;
        for (j=1; j<prot.res[i].n_conf; j++) {
            total_occ += prot.res[i].conf[j].counter;
        }
        for (j=1; j<prot.res[i].n_conf; j++) {
            prot.res[i].conf[j].occ = (float) prot.res[i].conf[j].counter/ (float)n/2.;
            //printf("res %4d, conf %3d, counter %6d, %5d, total counter %6.0f, occ %8.3f\n", i,j,prot.res[i].conf[j].counter, n, total_occ, prot.res[i].conf[j].occ);
            //prot.res[i].conf[j].occ = (float) prot.res[i].conf[j].counter/ total_occ;
        }
    }


    for (i_res=0; i_res<prot.n_res; i_res++) {
        int i_ngh;
        for (i_ngh = 0; i_ngh < prot.res[i_res].n_ngh; i_ngh++) {
            j_res = prot.res[i_res].ngh[i_ngh]->i_res_prot;
            if (j_res < i_res) continue;
            for (i_conf = 1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                free(pairwise[i_res][i_ngh][i_conf]);
            }
            free(pairwise[i_res][i_ngh]);
        }
        free(pairwise[i_res]);
        free(prot.res[i_res].ngh);
        prot.res[i_res].ngh = NULL;
        prot.res[i_res].n_ngh = 0;
    }
    
    /* delete low occpancy conformers */
    C = 0;
    for (i=0; i<prot.n_res; i++) {
        float highest_occ_exposed = -1.;
        int highest_occ_kconf = -1;
        for (j=prot.res[i].n_conf_ori; j<prot.res[i].n_conf; j++) {
            /* if (prot.res[i].conf[j].history[2] != 'H' && prot.res[i].conf[j].occ <= env.repack_cutoff) { */
            if (prot.res[i].conf[j].occ <= env.repack_cutoff &&
                prot.res[i].conf[j].history[2] != 'E') { /* keep exposed */
                del_conf(&prot.res[i], j);
                j--;
                C++;
            }
        }
        
        /* if there are more than 1 exposed conformer, delete the low occupancy ones */
        for (j=prot.res[i].n_conf_ori; j<prot.res[i].n_conf; j++) {
            if (prot.res[i].conf[j].history[2] != 'E') continue;
            if (prot.res[i].conf[j].occ > highest_occ_exposed) {
                highest_occ_exposed = prot.res[i].conf[j].occ;
                highest_occ_kconf = j;
            }
            //printf("max_occ %s %4d %3d %8.3f %8.3f\n",prot.res[i].resName,prot.res[i].resSeq,j,highest_occ_exposed,prot.res[i].conf[highest_occ_kconf].E_rxn);
        }
        for (j=prot.res[i].n_conf-1; j>=prot.res[i].n_conf_ori; j--) {
            if (prot.res[i].conf[j].history[2] != 'E') continue;
            if (prot.res[i].conf[j].occ >= env.repack_cutoff) continue;
            if (j == highest_occ_kconf) continue;
            del_conf(&prot.res[i], j);
            C++;
        }
    }

    /* remove protons again */
    delete_h(prot);
    rm_dupconf_hv(prot);
    
    free(state.res);
    for (i_res = 0; i_res < prot.n_res; i_res++) {
        if (!prot.res[i_res].n_ngh) continue;
        prot.res[i_res].n_ngh = 0;
        free(prot.res[i_res].ngh);
        prot.res[i_res].ngh = NULL;
    }

    return C;
}

int rot_refine(PROT prot, MICROSTATE state, float ****pairwise) {
    int i, i_ngh;
    float E_min;
    char switching;
    int i_res, j_res, i_conf, j_conf, n_candidate, candidate[9999];
    int n;
    int iRes[9999];  /* random residue and conformer indices */
    float E_state[9999];
    n = 0;

    /* remove Purify UMR error */
    memset(E_state, 0, 9999*sizeof(float));

    shuffle_n(iRes, prot.n_res);
    for (i=0; i<prot.n_res; i++) {                  /* loop over residues */
        float repack_e_thr;
        i_res = iRes[i];
        if (prot.res[i_res].n_conf < 1) continue;
        i_conf = state.res[i_res];

        E_min = 99999.;
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {       /* calculate energy for every conformer */
            E_state[i_conf] = prot.res[i_res].conf[i_conf].E_self;
            for (i_ngh = 0; i_ngh < prot.res[i_res].n_ngh; i_ngh++) {
                j_res = prot.res[i_res].ngh[i_ngh]->i_res_prot;
                j_conf = state.res[j_res];
                if (i_res < j_res)
                    E_state[i_conf] += pairwise[i_res][i_ngh][i_conf][j_conf];
                else
                    E_state[i_conf] += pairwise[i_res][i_ngh][j_conf][i_conf];  /* if i_res > j_res, matrix is transposed, look at the place when pairwise are calculated */
            }
            if (E_state[i_conf] > 99999.) {
                E_state[i_conf] = 99999.;
            }
            if (E_state[i_conf] < E_min) {
                E_min = E_state[i_conf];
            }
        }

        //printf("%d,%8.3f\n",i_res,E_min);
        /* if difference btw current conformer and minimum is bigger than threshold, switch */
        switching = 0;
        i_conf = state.res[i_res];
        if (prot.res[i_res].sas > 0.5) repack_e_thr = env.repack_e_thr_exposed;
        else repack_e_thr = 2.*prot.res[i_res].sas*env.repack_e_thr_exposed + (1.-2.*prot.res[i_res].sas)*env.repack_e_thr_buried;
        if (E_state[i_conf] - E_min > repack_e_thr) {
            switching = 1;
            n_candidate = 0;
        }
        
        if (switching) {
            /* collect candidates */
            for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                //printf("%d,%8.3f,%8.3f,%8.3f,%8.3f\n",i_res,E_min,E_state[i_conf],E_state[i_conf] - E_min,env.repack_e_thr);
                if (E_state[i_conf] - E_min < repack_e_thr) {
                    n_candidate++;
                    candidate[n_candidate-1] = i_conf;
                }
            }
            
            i_conf = state.res[i_res];
            prot.res[i_res].conf[i_conf].on = 0;
            i_conf = candidate[(int)(ran2(&idum) * (float)n_candidate)];
            prot.res[i_res].conf[i_conf].on = 1;
            //printf("switch %3d, from %3d to %3d. n_candi=%3d. n_conf=%3d,E_min=%10.3f\n",i_res,state.res[i_res],i_conf,n_candidate,prot.res[i_res].n_conf,E_min);
            state.res[i_res] = i_conf;
            n++;
        }
        
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            if (E_state[i_conf] - E_min < repack_e_thr) {
                prot.res[i_res].conf[i_conf].counter_trial = 1;
            }
            else {
                prot.res[i_res].conf[i_conf].counter_trial = 0;
            }
        }
    }
    if (n == 0) {
        for (i_res=0; i_res<prot.n_res; i_res++) {
            for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                prot.res[i_res].conf[i_conf].counter += prot.res[i_res].conf[i_conf].counter_trial;
            }
        }
    }
    return n;
}

void sp3_3known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR r3, VECTOR *r4, float bond_len_04);
void sp3_2known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, VECTOR *r4, float bond_len_03, float bond_angle_304);
void sp3_1known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, VECTOR *r4, VECTOR *r5, float bond_len_03, float bond_angle_103, float torsion_angle_3012);
void sp2_2known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, float bond_len_03);
void sp2_2known_with_angle(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, float bond_len_03, float bond_angle_103);
void sp2_1known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, float bond_len_03, float bond_angle_103, float torsion_angle_3012);
void sp2d_3known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR r3, VECTOR *r4, float bond_len_04);
void sp3d2_5known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR r3, VECTOR r4, VECTOR r5, VECTOR *r6, float bond_len_06);
float get_bond_length(CONF *conf_p, ATOM *atom1_p, ATOM *atom2_p);
float get_bond_angle(CONF *conf_p, ATOM *atom0_p, ATOM *atom1_p, ATOM *atom2_p, char *orbital);

int place_missing(PROT prot, int handle_addconf) {
    int         i_res, i_conf, i_atom, ins;
    FILE        *debug_fp;
    CONNECT     connect;
    char        orbital[10],sbuffer[5],name[MAXCHAR_LINE];
    RES         *res_p;
    CONF        *conf_p;
    ATOM        *atom_p, *back_atom_p;
    ATOM        *known_atoms[MAX_CONNECTED], *to_complete_atoms[MAX_CONNECTED], dummy_atom[MAX_CONNECTED];
    int         n_known, n_complete;
    int         i_connect, i_dummy, i_fold, i_complete, t_connect;
    int         i_corner, j_corner, k_corner, l_corner, start;
    VECTOR      corners[4], v;
    float       bond_length, bond_angle, torsion_angle, a;
    TORS        tors;
    int         fatal = 0;
    int error;
    int kr,kc,ka;
    char sbuff[MAXCHAR_LINE],sbuff2[MAXCHAR_LINE], siatom[MAXCHAR_LINE];
    char resName[4];
    int  Missing, n_added=0;
    STRINGS     conflist;
    
    memset(dummy_atom,0,MAX_CONNECTED*sizeof(ATOM));
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        res_p = &prot.res[i_res];
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            conf_p = &prot.res[i_res].conf[i_conf];
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (atom_p->on) continue;

                sprintf(sbuffer,"%d",i_atom);
                if( param_get("ATOMNAME", conf_p->confName, sbuffer, name) ) {
                    strcpy(name, "    ");
                    fatal++;
                }
                while (strlen(name)<4) strcat(name, " ");
                strncpy(atom_p->name, name, 4);
                atom_p->name[4] = '\0';
                if (atom_p->name[1] != 'H') { /* heavy atom */
                    if(!param_get("CONNECT", conf_p->confName, atom_p->name, &connect) ) { /* find connectivity */
                        strip(orbital, connect.orbital);
                        if (strcmp(orbital,"ion")) { /* Not ion */
                            if (!strcmp(atom_p->name, " N  ") || !strcmp(atom_p->name, " CA ")) {
                                if (i_res>0) {
                                    if (!strcmp(prot.res[i_res-1].resName, "NTR") || !strcmp(prot.res[i_res-1].resName, "NTG")) {
                                        continue;
                                    }
                                }
                            }
                            if (!strcmp(atom_p->name, " C  ") || !strcmp(atom_p->name, " O  ")) {
                                if (i_res<prot.n_res-1) {
                                    if (!strcmp(prot.res[i_res+1].resName, "CTR")){
                                        continue;
                                    }
                                }
                            }

                            //printf("    Warning! May add heavy atom %s onto residue \"%s %c%04d\"\n",atom_p->name, res_p->resName, res_p->chainID, res_p->resSeq);
                        }
                    }
                }
            }
        }
    }

    get_connect12(prot);

    error = 0;
    for (i_res=0; i_res<prot.n_res; i_res++) {
        res_p = &prot.res[i_res];
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            //int k_conf; /* used for checking duplicates later in the loop */
            conf_p = &prot.res[i_res].conf[i_conf];
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                if( param_get("CONNECT", conf_p->confName, atom_p->name, &connect) ) {
                    debug_fp = fopen(env.debug_log, "a");
                    fprintf(debug_fp, "   Error! place_missing(): Can't find CONNECT parameter of conformer \"%s\" atom \"%s\"\n", conf_p->confName, atom_p->name);
                    fclose(debug_fp);
                    continue;
                }
                //printf("orbital %s\n",connect.orbital);
                strip(orbital, connect.orbital);
                if (!strcmp(orbital,"ion")) continue;
                
                /* Collect known atoms and unknown atoms */
                memset(known_atoms,0,MAX_CONNECTED*sizeof(void *));
                memset(to_complete_atoms,0,MAX_CONNECTED*sizeof(void *));
                n_known = 0;
                n_complete = 0;
                for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                    if (!atom_p->connect12[i_connect]) {
                        if (!strncmp(atom_p->name, " CA ", 4) && prot.res[i_res].n_conf == 1) {
                            for (t_connect=0; t_connect<connect.n; t_connect++) {
                                if (!strncmp(connect.atom[t_connect].name, " CB ", 4)) break;
                            }
                            if (t_connect == connect.n) break;
                            else { 
                                param_get("CONFLIST", prot.res[i_res].resName, "", &conflist);
                                int tem_n_atom;
                                param_get("NATOM", conflist.strings[1], "", &tem_n_atom); 
                                ins = ins_conf(&prot.res[i_res], prot.res[i_res].n_conf, tem_n_atom);
                                strcpy(prot.res[i_res].conf[ins].confName, conflist.strings[1]);
                                strcpy(prot.res[i_res].conf[ins].history, atom_p->history);
                                strncpy(prot.res[i_res].conf[ins].history, "01", 2);
                                prot.res[i_res].conf[ins].altLoc = ' ';
                                get_connect12_conf(i_res, i_conf, prot);
                            }
                        }
                        else break;
                    }
                    if (atom_p->connect12[i_connect]->on) {
                        n_known++;
                        known_atoms[n_known-1] = atom_p->connect12[i_connect];
                    }
                    else {
                        n_complete++;
                        to_complete_atoms[n_complete-1] = atom_p->connect12[i_connect];
                    }
                }
                
                if (!n_complete) continue;
                
                if (!strcmp(orbital, "sp3")) {
                    /* If total number connected atoms is less than 4 for sp3, use dummy atoms to complete 4 slots */
                    if ( (n_known + n_complete) < 4 ) {
                        for (i_dummy = n_complete; i_dummy < (4-n_known); i_dummy++) {
                            to_complete_atoms[i_dummy] = &dummy_atom[i_dummy];
                        }
                    }
                    
                    if (n_known == 3) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        
                        sp3_3known( atom_p->xyz,
                        known_atoms[0]->xyz,
                        known_atoms[1]->xyz,
                        known_atoms[2]->xyz,
                        &to_complete_atoms[0]->xyz,
                        bond_length );
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res,i_conf,prot);
                    }
                    else if (n_known == 2) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        //bond_angle = 109.;
                        bond_angle = get_bond_angle(conf_p, atom_p, known_atoms[0], to_complete_atoms[0], orbital);

                        //printf("   Debugging! Case sp3, n_known = 2\n");
                        
                        sp3_2known( atom_p->xyz,
                        known_atoms[0]->xyz,
                        known_atoms[1]->xyz,
                        &to_complete_atoms[0]->xyz,
                        &to_complete_atoms[1]->xyz,
                        bond_length,
                        bond_angle);
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res,i_conf,prot);
                        
                        if (!i_conf) continue; /* do not add extra conf for backbone */
                        //printf("   Debugging! Case sp3, n_known = 2, i_conf!=0\n");
                        if (!handle_addconf) continue; /* do not add extra conf if the flag is 0 */
                        if (to_complete_atoms[0]->name[1] == 'H') {
                            if (to_complete_atoms[1]->name[1] == 'H') continue; /* do not add extra conf if added atoms are all protons */
                        }
                        //printf("   Debugging! Case sp3, n_known = 2, handle_addconf !=0\n");
                        
                        ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                        if (ins == USERERR) return USERERR;
                        conf_p = &prot.res[i_res].conf[i_conf]; /* reassign conf_p because ins_conf changes the memory position of conf array */
                        if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                        get_connect12_conf(i_res, res_p->n_conf-1, prot);
                        
                        sp3_2known( atom_p->xyz,
                            known_atoms[0]->xyz,
                            known_atoms[1]->xyz,
                            &to_complete_atoms[1]->xyz,
                            &to_complete_atoms[0]->xyz,
                            bond_length,
                            bond_angle);
                        //printf("debug %s %s\n",prot.res[i_res].resName,to_complete_atoms[0]->name);
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res, i_conf, prot);
                    }
                    else if (n_known == 1) {
                        int n_fold;
                        //printf("   Debugging! residue %s%4d,conformer %s\n", res_p->resName,res_p->resSeq,conf_p->confName);
                        /* Get TORSION parameter, if not exit then make one */
                        if ( param_get("TORSION",conf_p->confName, to_complete_atoms[0]->name, &tors) ) {
                            strncpy(resName, conf_p->confName, 3); resName[3] = '\0';
                            if ( param_get("TORSION",resName, to_complete_atoms[0]->name, &tors) ) {
                                memset(&tors,0,sizeof(TORS));
                                strcpy(tors.atom1, atom_p->name);
                                strcpy(tors.atom2, known_atoms[0]->name);
                                for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                                    if (!known_atoms[0]->connect12[i_connect]) break;
                                    if (!known_atoms[0]->connect12[i_connect]->on) continue;
                                    if (!strcmp(known_atoms[0]->connect12[i_connect]->name, atom_p->name)) continue;
                                    strcpy(tors.atom3, known_atoms[0]->connect12[i_connect]->name);
                                    break;
                                }
                                tors.V2[0]     = 0.;
                                tors.n_fold[0] = 3.;
                                tors.gamma[0]  = 0.;
                                tors.opt_hyd = 0;
                                param_sav("TORSION", conf_p->confName, to_complete_atoms[0]->name, &tors, sizeof(TORS));
                                debug_fp = fopen(env.debug_log, "a");
                                fprintf(debug_fp, "TORSION  %s %s %s %s %s  f  %9.1f %9.0f %9.2f\n",
                                conf_p->confName,
                                to_complete_atoms[0]->name,
                                tors.atom1,
                                tors.atom2,
                                tors.atom3,
                                tors.V2[0],
                                tors.n_fold[0],
                                tors.gamma[0]/env.d2r);
                                fclose(debug_fp);
                                error++;
                            }
                        }
                        
                        if (strcmp(atom_p->name, tors.atom1)) {
                            printf("   Error! Atom %s is in torsion parameter of conformer %s atom %s, but connected atom %s was found\n",
                            tors.atom1,conf_p->confName,to_complete_atoms[0]->name,atom_p->name);
                            continue;
                        }
                        if (strcmp(known_atoms[0]->name, tors.atom2)) {
                            printf("   Error! Atom %s is in torsion parameter of conformer %s atom %s, but connected atom %s was found\n",
                            tors.atom2,conf_p->confName,to_complete_atoms[0]->name,known_atoms[0]->name);
                            continue;
                        }
                        back_atom_p = NULL;
                        for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                            if (!known_atoms[0]->connect12[i_connect]) break;
                            if (!known_atoms[0]->connect12[i_connect]->on) continue;
                            if (!strcmp(known_atoms[0]->connect12[i_connect]->name, tors.atom3)) {back_atom_p = known_atoms[0]->connect12[i_connect]; break;}
                        }
                        if (!back_atom_p) {
                            printf("   Error! place_missing(): cannot put hydrogen on atom %s in residue %s %d, because %s-%s- N/A is unexpected.\n",
                            atom_p->name, res_p->resName, res_p->resSeq, atom_p->name, known_atoms[0]->name);
                            continue;
                        }
                        
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        bond_angle  = get_bond_angle(conf_p, atom_p, known_atoms[0], to_complete_atoms[0], orbital);
                        
                        n_fold = tors.n_fold[0];
                        /*
                        if (tors.n_fold[0] != 1 && handle_addconf ==2) {
                            printf("ires=%d, nconf=%d\n",i_res,prot.res[i_res].n_conf);
                            if (prot.res[i_res].n_conf <500) {
                                n_fold = 72;
                            }
                            else if (prot.res[i_res].n_conf <10000) {
                                n_fold = 36;
                            }
                            else if (prot.res[i_res].n_conf <100000) {
                                n_fold = 12;
                            }
                        }
                        */
                        for (i_fold=0; i_fold<n_fold; i_fold++) {
                            if (i_fold) {
                                if (!i_conf) break; /* do not add extra conf for backbone */
                                if (!handle_addconf) break; /* do not add extra conf if the flag is 0 */
                                
                                if (to_complete_atoms[0]->name[1] == 'H') {
                                    if (to_complete_atoms[1]->name[1] == 'H') {
                                        if (to_complete_atoms[2]->name[1] == 'H') break; /* do not add extra conf if added atoms are all protons (methyl) */
                                    }
                                }

                                ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                                if (ins == USERERR) return USERERR;
                                conf_p = &prot.res[i_res].conf[i_conf];
                                if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                                get_connect12_conf(i_res, res_p->n_conf-1, prot);
                                //printf("   Debugging! residue %s%4d,conformer %s to conformer %s\n", res_p->resName,res_p->resSeq,conf_p->confName,res_p->conf[res_p->n_conf-1].confName);
                            }
                            
                            torsion_angle = (env.PI + tors.gamma[0] + (i_fold)*2.*env.PI)/n_fold;
                            sp3_1known(atom_p->xyz,
                            known_atoms[0]->xyz,
                            back_atom_p->xyz,
                            &to_complete_atoms[0]->xyz,
                            NULL,
                            NULL,
                            bond_length,
                            bond_angle,
                            torsion_angle );
                            
                            //printf("   Debugging! i_fold %d, n_fold %f, residue %s%4d,history %s, pos %d\n",i_fold,tors.n_fold[0], res_p->resName,res_p->resSeq,conf_p->history,i_conf);
                            //printf("   Debugging2! residue%d %s%4d,conformer%d %s and last conformer %s\n", i_res, res_p->resName,res_p->resSeq,i_conf,conf_p->confName,res_p->conf[res_p->n_conf-1].confName);
                            to_complete_atoms[0]->on = 1;
                            n_added++;
                            get_connect12_conf(i_res, i_conf, prot);
                        }
                    }
                    else if (n_known == 0) {
                        int counter;
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        /* Define four corners of the box, atom would be placed on the conner. */
                        a = bond_length/sqrt(3.);
                        counter = 0;
                        v.x =  a; v.y =  a; v.z =  a; corners[0] = vector_vplusv(atom_p->xyz, v);
                        v.x = -a; v.y = -a; v.z =  a; corners[1] = vector_vplusv(atom_p->xyz, v);
                        v.x = -a; v.y =  a; v.z = -a; corners[2] = vector_vplusv(atom_p->xyz, v);
                        v.x =  a; v.y = -a; v.z = -a; corners[3] = vector_vplusv(atom_p->xyz, v);
                        for (i_corner = 0; i_corner < 4; i_corner++) {
                            
                            if (n_complete <=1) start = 3; else start = i_corner+1;
                            for (j_corner = start; j_corner < 4; j_corner++) {
                                
                                if (n_complete <=2) start = 3; else start = j_corner+1;
                                for (k_corner = start; k_corner < 4; k_corner++) {
                                    
                                    if (n_complete <=3) start = 3; else start = k_corner+1;
                                    for (l_corner = start; l_corner < 4; l_corner++) {
                                        
                                        to_complete_atoms[0]->xyz = corners[i_corner];
                                        to_complete_atoms[1]->xyz = corners[j_corner];
                                        to_complete_atoms[2]->xyz = corners[k_corner];
                                        to_complete_atoms[3]->xyz = corners[l_corner];
                                        
                                        to_complete_atoms[0]->on = 1;
                                        to_complete_atoms[1]->on = 1;
                                        to_complete_atoms[2]->on = 1;
                                        to_complete_atoms[3]->on = 1;
                                        n_added++;
                                        
                                        if (!handle_addconf) continue;
                                        
                                        ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                                        if (ins == USERERR) return USERERR;
                                        conf_p = &prot.res[i_res].conf[i_conf];
                                        if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                                        get_connect12_conf(i_res, i_conf, prot);
                                    }
                                }
                            }
                        }
                        
                        v.x = -a; v.y = -a; v.z = -a; corners[0] = vector_vplusv(atom_p->xyz, v);
                        v.x =  a; v.y =  a; v.z = -a; corners[1] = vector_vplusv(atom_p->xyz, v);
                        v.x =  a; v.y = -a; v.z =  a; corners[2] = vector_vplusv(atom_p->xyz, v);
                        v.x = -a; v.y =  a; v.z =  a; corners[3] = vector_vplusv(atom_p->xyz, v);
                        for (i_corner = 0; i_corner < 4; i_corner++) {
                            
                            if (n_complete <=1) start = 3; else start = i_corner+1;
                            for (j_corner = start; j_corner < 4; j_corner++) {
                                
                                if (n_complete <=2) start = 3; else start = j_corner+1;
                                for (k_corner = start; k_corner < 4; k_corner++) {
                                    
                                    if (n_complete <=3) start = 3; else start = k_corner+1;
                                    for (l_corner = start; l_corner < 4; l_corner++) {
                                        
                                        to_complete_atoms[0]->xyz = corners[i_corner];
                                        to_complete_atoms[1]->xyz = corners[j_corner];
                                        to_complete_atoms[2]->xyz = corners[k_corner];
                                        to_complete_atoms[3]->xyz = corners[l_corner];
                                        
                                        to_complete_atoms[0]->on = 1;
                                        to_complete_atoms[1]->on = 1;
                                        to_complete_atoms[2]->on = 1;
                                        to_complete_atoms[3]->on = 1;
                                        n_added++;
                                        
                                        if (!handle_addconf) continue;
                                        ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                                        if (ins == USERERR) return USERERR;
                                        conf_p = &prot.res[i_res].conf[i_conf];
                                        if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                                        get_connect12_conf(i_res, res_p->n_conf-1, prot);
                                    }
                                }
                            }
                        }
                        
                        if (handle_addconf) {
                            del_conf(res_p, res_p->n_conf-1);
                        }
                        conf_p = &prot.res[i_res].conf[i_conf];
                        continue;
                    }
                    else {
                        printf("   Error! place_missing(): number of known atoms can't be handled, sp3 with %i known atoms\n",n_known);
                        fatal++;
                    }
                }
                else if (!strcmp(orbital, "sp2")) {
                    if ( (n_known + n_complete) < 3 ) {
                        for (i_dummy = n_complete; i_dummy < (3-n_known); i_dummy++) {
                            to_complete_atoms[i_dummy] = &dummy_atom[i_dummy];
                        }
                    }
                    if (n_known == 2) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        /* get bond angle from tpl file, force the subroutine not to guess from orbital type */
                        bond_angle = get_bond_angle(conf_p, atom_p, known_atoms[0], to_complete_atoms[0], NULL);
                        //printf("%s %s %s %8.3f\n", atom_p->name, known_atoms[0]->name, to_complete_atoms[0]->name, bond_angle*180./env.PI);
                        if (bond_angle < 0) {
                            sp2_2known( atom_p->xyz,
                            known_atoms[0]->xyz,
                            known_atoms[1]->xyz,
                            &to_complete_atoms[0]->xyz,
                            bond_length);
                        }
                        else {
                            sp2_2known_with_angle( atom_p->xyz,
                            known_atoms[0]->xyz,
                            known_atoms[1]->xyz,
                            &to_complete_atoms[0]->xyz,
                            bond_length,
                            bond_angle );
                        }
                        
                        to_complete_atoms[0]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res, i_conf, prot);
                    }
                    else if (n_known == 1) {
                        int n_fold;
                        if ( param_get("TORSION",conf_p->confName, to_complete_atoms[0]->name, &tors) ) {
                            strncpy(resName, conf_p->confName, 3); resName[3] = '\0';
                            if ( param_get("TORSION",resName, to_complete_atoms[0]->name, &tors) ) {
                                memset(&tors,0,sizeof(TORS));
                                strcpy(tors.atom1, atom_p->name);
                                strcpy(tors.atom2, known_atoms[0]->name);
                                for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                                    if (!known_atoms[0]->connect12[i_connect]) break;
                                    if (!known_atoms[0]->connect12[i_connect]->on) continue;
                                    if (!strcmp(known_atoms[0]->connect12[i_connect]->name, atom_p->name)) continue;
                                    strcpy(tors.atom3, known_atoms[0]->connect12[i_connect]->name);
                                    break;
                                }
                                tors.V2[0]     = 0.;
                                tors.n_fold[0] = 2.;
                                tors.gamma[0]  = 180.*env.d2r;
                                tors.opt_hyd = 0;
                                param_sav("TORSION", conf_p->confName, to_complete_atoms[0]->name, &tors, sizeof(TORS));
                                debug_fp = fopen(env.debug_log, "a");
                                fprintf(debug_fp, "TORSION  %s %s %s %s %s  f  %9.1f %9.0f %9.2f\n",
                                conf_p->confName,
                                to_complete_atoms[0]->name,
                                tors.atom1,
                                tors.atom2,
                                tors.atom3,
                                tors.V2[0],
                                tors.n_fold[0],
                                tors.gamma[0]/env.d2r);
                                fclose(debug_fp);
                                error++;
                            }
                        }
                        if (strcmp(atom_p->name, tors.atom1)) {
                            printf("   Error! Cannot add atom \'%s\' onto \'%s\' because TORSION parameter of atom \'%s\' doesn't match the connected atoms in this protein, double check the TORSION parameter. (Res: %s %c%04d)\n",
                            to_complete_atoms[0]->name,atom_p->name,to_complete_atoms[0]->name,res_p->resName,res_p->chainID,res_p->resSeq);
                            continue;
                        }
                        if (strcmp(known_atoms[0]->name, tors.atom2)) {
                            printf("   Error! Cannot add atom \'%s\' onto \'%s\' because TORSION parameter of atom \'%s\' doesn't match the connected atoms in this protein, double check the TORSION parameter. (Res: %s %c%04d)\n",
                            to_complete_atoms[0]->name,atom_p->name,to_complete_atoms[0]->name,res_p->resName,res_p->chainID,res_p->resSeq);
                            continue;
                        }
                        back_atom_p = NULL;
                        for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                            if (!known_atoms[0]->connect12[i_connect]) break;
                            if (!known_atoms[0]->connect12[i_connect]->on) continue;
                            if (!strcmp(known_atoms[0]->connect12[i_connect]->name, tors.atom3)) {back_atom_p = known_atoms[0]->connect12[i_connect]; break;}
                        }
                        if (!back_atom_p) {
                            printf("   Error! Cannot add atom \'%s\' onto \'%s\' because TORSION parameter of atom \'%s\' doesn't match the connected atoms in this protein, double check the TORSION parameter. (Res: %s %c%04d)\n",
                            to_complete_atoms[0]->name,atom_p->name,to_complete_atoms[0]->name,res_p->resName,res_p->chainID,res_p->resSeq);
                            continue;
                        }
                        
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        bond_angle = get_bond_angle(conf_p, atom_p, known_atoms[0], to_complete_atoms[0], orbital);
                        //bond_angle = 120.;
                        //printf("%s %s %s %8.3f\n", atom_p->name, known_atoms[0]->name, to_complete_atoms[0]->name, bond_angle*180./env.PI);
                        n_fold = tors.n_fold[0];
                        /*
                        if (tors.n_fold[0] != 1 && handle_addconf ==2) {
                            printf("ires=%d, nconf=%d\n",i_res,prot.res[i_res].n_conf);
                            if (prot.res[i_res].n_conf <500) {
                                n_fold = 72;
                            }
                            else if (prot.res[i_res].n_conf <10000) {
                                n_fold = 36;
                            }
                            else if (prot.res[i_res].n_conf <100000) {
                                n_fold = 12;
                            }
                        }
                        */
                        for (i_fold=0; i_fold<n_fold; i_fold++) {
                            if (i_fold) {
                                if (!i_conf) break; /* do not add extra conf for backbone */
                                if (!handle_addconf) break; /* do not add extra conf if the flag is 0 */
                                
                                if (to_complete_atoms[0]->name[1] == 'H') {
                                    if (to_complete_atoms[1]->name[1] == 'H') break; /* do not add extra conf if added atoms are all protons */
                                }
                                
                                ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                                if (ins == USERERR) return USERERR;
                                conf_p = &prot.res[i_res].conf[i_conf];
                                if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                                get_connect12_conf(i_res, res_p->n_conf-1, prot);
                            }
                            torsion_angle = (env.PI + tors.gamma[0] + (i_fold)*2.*env.PI)/n_fold;
                            sp2_1known(atom_p->xyz,
                            known_atoms[0]->xyz,
                            back_atom_p->xyz,
                            &to_complete_atoms[0]->xyz,
                            bond_length,
                            bond_angle,
                            torsion_angle );
                            
                            to_complete_atoms[0]->on = 1;
                            n_added++;
                            get_connect12_conf(i_res,i_conf,prot);
                        }
                        
                        /* only 1 atom is added, rollback to add the second */
                        i_atom--;
                        continue;
                    }
                    else if (n_known == 0) {
                        printf("   Error! no rule to add atom \'%s\' onto \'%s\': (atom %s in conformer %s have no atom connected)\n",to_complete_atoms[0]->name, atom_p->name, atom_p->name,conf_p->confName);
                    }
                    else {
                        printf("   Error! place_missing(): number of connected atoms wrong, %i n_connected atoms\n",n_known);
                        fatal++;
                    }
                }
                else if (!strcmp(orbital, "sp2d")) {
                    if ( (n_known + n_complete) < 4 ) {
                        for (i_dummy = n_complete; i_dummy < (4-n_known); i_dummy++) {
                            to_complete_atoms[i_dummy] = &dummy_atom[i_dummy];
                        }
                    }
                    
                    if (n_known == 3) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        
                        sp2d_3known( atom_p->xyz,
                        known_atoms[0]->xyz,
                        known_atoms[1]->xyz,
                        known_atoms[2]->xyz,
                        &to_complete_atoms[0]->xyz,
                        bond_length );
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res,i_conf,prot);
                    }
                    else {
                        printf("   Error! place_missing(): number of known atoms can't be handled, sp2d with %i known atoms\n",n_known);
                        fatal++;
                    }
                }
                else if (!strcmp(orbital, "sp3d2")) {
                    if ( (n_known + n_complete) < 6 ) {
                        for (i_dummy = n_complete; i_dummy < (6-n_known); i_dummy++) {
                            to_complete_atoms[i_dummy] = &dummy_atom[i_dummy];
                        }
                    }
                    
                    if (n_known == 5) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        
                        sp3d2_5known( atom_p->xyz,
                        known_atoms[0]->xyz,
                        known_atoms[1]->xyz,
                        known_atoms[2]->xyz,
                        known_atoms[3]->xyz,
                        known_atoms[4]->xyz,
                        &to_complete_atoms[0]->xyz,
                        bond_length );
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res,i_conf,prot);
                    }
                    else {
                        printf("   Error! place_missing(): number of known atoms can't be handled, sp3d2 with %i known atoms\n",n_known);
                        fatal++;
                    }
                }
                else {
                    printf("   Error! place_missing(): Don't know yet how to deal with orbital type %s for atom %s in conformer %s\n",orbital,atom_p->name, conf_p->confName);
                }
                /* delete dulipcates probably causing memory trouble
                for (k_conf=prot.res[i_res].n_conf-1; k_conf>=1; k_conf--) {
                    int j_conf;
                    for (j_conf=1; j_conf<k_conf; j_conf++) {
                        if( strcmp(prot.res[i_res].conf[k_conf].confName, prot.res[i_res].conf[j_conf].confName)) continue;
                        if(!cmp_conf(prot.res[i_res].conf[k_conf], prot.res[i_res].conf[j_conf], 0.01)) {
                            del_conf(&prot.res[i_res], k_conf);
                            break;
                        }
                    }
                }
                */
                /* END of this atom */
            }
        }
        
        while (1) {
            int  counter, j_conf;
            for (i_conf=0; i_conf<res_p->n_conf; i_conf++) {
                if (!strncmp(res_p->conf[i_conf].history+6,"____",4)) break;
            }
            if (i_conf >= prot.res[i_res].n_conf) break;
            counter = 0;
            for (j_conf=0; j_conf<res_p->n_conf; j_conf++) {
                if (strncmp(res_p->conf[i_conf].history,res_p->conf[j_conf].history,6)) continue;
                res_p->conf[j_conf].history[6] = 'M';
                sprintf(sbuff,"%03d",counter);
                strncpy(res_p->conf[j_conf].history+7,sbuff,3);
                counter++;
            }
        }
    }
    
    /* Check for if there are still missing atoms after protonation. the way of treating NTR and CTR here is not good and standard */
    if (!n_added) {
        Missing = 0;
        for (kr=0; kr<prot.n_res; kr++) {
            for (kc=0; kc<prot.res[kr].n_conf; kc++) {
                for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
                    if (!prot.res[kr].conf[kc].atom[ka].on) {
                        sprintf(siatom, "%3d", ka);
                        if (param_get("ATOMNAME", prot.res[kr].conf[kc].confName, siatom, sbuff)) {
                            printf("   Missing ATOMNAME records for slot \"%d\" of conformer %s.\n",
                            ka, prot.res[kr].conf[kc].confName);
                            Missing++;
                        }
                        else {
                            if (!strcmp(prot.res[kr-1].resName, "NTR") || !strcmp(prot.res[kr-1].resName, "NTG")) {
                                if (!strncmp(sbuff+1,"HA",2)) sbuff[0] = ' ';
                                if (!strncmp(sbuff+1,"H ",2)) sbuff[0] = '1';
                                if (!strncmp(sbuff+1,"H\0",2)) sbuff[0] = '1';
                                if (!param_get("IATOM","NTRBK",sbuff,sbuff2)) continue;
                                if (!param_get("IATOM","NTR01",sbuff,sbuff2)) continue;
                                if (!param_get("IATOM","NTGBK",sbuff,sbuff2)) continue;
                                if (!param_get("IATOM","NTG01",sbuff,sbuff2)) continue;
                            }
                            if (!strcmp(prot.res[kr+1].resName, "CTR")) {
                                sbuff[0] = ' ';
                                if (!param_get("IATOM","CTRBK",sbuff,sbuff2)) continue;
                                if (!param_get("IATOM","CTR01",sbuff,sbuff2)) continue;
                            }
                            printf("   Missing atom %s at slot %4d of conf %s in \"%s %c %3d %02d\".\n",
                            sbuff, ka, prot.res[kr].conf[kc].confName, prot.res[kr].resName, prot.res[kr].chainID, prot.res[kr].resSeq,kc);
                            Missing++;
                        }
                    }
                }
            }
        }
    }
    if (error) {
       printf("   %d TORSION parameters were guessed and recorded in file %s.\n", error, env.debug_log);
       printf("   Modify and put them into a param file to have the full control.\n" );
    }
    return n_added;
}

int place_missing_res(PROT prot, int i_res, int handle_addconf) {
    int         i_conf, i_atom, ins;
    FILE        *debug_fp;
    CONNECT     connect;
    char        orbital[10],sbuffer[5],name[MAXCHAR_LINE];
    RES         *res_p;
    CONF        *conf_p;
    ATOM        *atom_p, *back_atom_p;
    ATOM        *known_atoms[MAX_CONNECTED], *to_complete_atoms[MAX_CONNECTED], dummy_atom[MAX_CONNECTED];
    int         n_known, n_complete;
    int         i_connect, i_dummy, i_fold, i_complete;
    int         i_corner, j_corner, k_corner, l_corner, start;
    VECTOR      corners[4], v;
    float       bond_length, bond_angle, torsion_angle, a;
    TORS        tors;
    int         fatal = 0;
    int error;
    int kc,ka;
    char sbuff[MAXCHAR_LINE],sbuff2[MAXCHAR_LINE], siatom[MAXCHAR_LINE];
    char resName[4];
    int  Missing, n_added=0;
    
    memset(dummy_atom,0,MAX_CONNECTED*sizeof(ATOM));
    
    res_p = &prot.res[i_res];
    for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
        conf_p = &prot.res[i_res].conf[i_conf];
        for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
            if (atom_p->on) continue;
            
            sprintf(sbuffer,"%d",i_atom);
            if( param_get("ATOMNAME", conf_p->confName, sbuffer, name) ) {
                strcpy(name, "    ");
                fatal++;
            }
            while (strlen(name)<4) strcat(name, " ");
            strncpy(atom_p->name, name, 4);
            atom_p->name[4] = '\0';
            if (atom_p->name[1] != 'H') { /* heavy atom */
                if(!param_get("CONNECT", conf_p->confName, atom_p->name, &connect) ) { /* find connectivity */
                    strip(orbital, connect.orbital);
                    if (strcmp(orbital,"ion")) { /* Not ion */
                        if (!strcmp(atom_p->name, " N  ") || !strcmp(atom_p->name, " CA ")) {
                            if (i_res>0) {
                                if (!strcmp(prot.res[i_res-1].resName, "NTR") || !strcmp(prot.res[i_res-1].resName, "NTG")) {
                                    continue;
                                }
                            }
                        }
                        if (!strcmp(atom_p->name, " C  ") || !strcmp(atom_p->name, " O  ")) {
                            if (i_res<prot.n_res-1) {
                                if (!strcmp(prot.res[i_res+1].resName, "CTR")){
                                    continue;
                                }
                            }
                        }
                        
                        //printf("    Warning! May add heavy atom %s onto residue \"%s %c%04d\"\n",atom_p->name, res_p->resName, res_p->chainID, res_p->resSeq);
                    }
                }
            }
        }
    }
    
    for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
        get_connect12_conf(i_res, i_conf, prot);
    }

    error = 0;
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            //int k_conf; /* used for checking duplicates later in the loop */
            conf_p = &prot.res[i_res].conf[i_conf];
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                if( param_get("CONNECT", conf_p->confName, atom_p->name, &connect) ) {
                    debug_fp = fopen(env.debug_log, "a");
                    fprintf(debug_fp, "   Error! place_missing(): Can't find CONNECT parameter of conformer \"%s\" atom \"%s\"\n", conf_p->confName, atom_p->name);
                    fclose(debug_fp);
                    continue;
                }
                //printf("orbital %s\n",connect.orbital);
                strip(orbital, connect.orbital);
                if (!strcmp(orbital,"ion")) continue;
                
                /* Collect known atoms and unknown atoms */
                memset(known_atoms,0,MAX_CONNECTED*sizeof(void *));
                memset(to_complete_atoms,0,MAX_CONNECTED*sizeof(void *));
                n_known = 0;
                n_complete = 0;
                for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                    if (!atom_p->connect12[i_connect]) break;
                    if (atom_p->connect12[i_connect]->on) {
                        n_known++;
                        known_atoms[n_known-1] = atom_p->connect12[i_connect];
                    }
                    else {
                        n_complete++;
                        to_complete_atoms[n_complete-1] = atom_p->connect12[i_connect];
                    }
                }
                
                if (!n_complete) continue;
                
                if (!strcmp(orbital, "sp3")) {
                    /* If total number connected atoms is less than 4 for sp3, use dummy atoms to complete 4 slots */
                    if ( (n_known + n_complete) < 4 ) {
                        for (i_dummy = n_complete; i_dummy < (4-n_known); i_dummy++) {
                            to_complete_atoms[i_dummy] = &dummy_atom[i_dummy];
                        }
                    }
                    
                    if (n_known == 3) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        
                        sp3_3known( atom_p->xyz,
                        known_atoms[0]->xyz,
                        known_atoms[1]->xyz,
                        known_atoms[2]->xyz,
                        &to_complete_atoms[0]->xyz,
                        bond_length );
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res,i_conf,prot);
                    }
                    else if (n_known == 2) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        //bond_angle = 109.;
                        bond_angle = get_bond_angle(conf_p, atom_p, known_atoms[0], to_complete_atoms[0], orbital);

                        //printf("   Debugging! Case sp3, n_known = 2\n");
                        
                        sp3_2known( atom_p->xyz,
                        known_atoms[0]->xyz,
                        known_atoms[1]->xyz,
                        &to_complete_atoms[0]->xyz,
                        &to_complete_atoms[1]->xyz,
                        bond_length,
                        bond_angle);
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res,i_conf,prot);
                        
                        if (!i_conf) continue; /* do not add extra conf for backbone */
                        //printf("   Debugging! Case sp3, n_known = 2, i_conf!=0\n");
                        if (!handle_addconf) continue; /* do not add extra conf if the flag is 0 */
                        if (to_complete_atoms[0]->name[1] == 'H') {
                            if (to_complete_atoms[1]->name[1] == 'H') continue; /* do not add extra conf if added atoms are all protons */
                        }
                        //printf("   Debugging! Case sp3, n_known = 2, handle_addconf !=0\n");
                        
                        ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                        if (ins == USERERR) return USERERR;
                        conf_p = &prot.res[i_res].conf[i_conf]; /* reassign conf_p because ins_conf changes the memory position of conf array */
                        if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                        get_connect12_conf(i_res, res_p->n_conf-1, prot);
                        
                        sp3_2known( atom_p->xyz,
                            known_atoms[0]->xyz,
                            known_atoms[1]->xyz,
                            &to_complete_atoms[1]->xyz,
                            &to_complete_atoms[0]->xyz,
                            bond_length,
                            bond_angle);
                        //printf("debug %s %s\n",prot.res[i_res].resName,to_complete_atoms[0]->name);
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res, i_conf, prot);
                    }
                    else if (n_known == 1) {
                        int n_fold;
                        //printf("   Debugging! residue %s%4d,conformer %s\n", res_p->resName,res_p->resSeq,conf_p->confName);
                        /* Get TORSION parameter, if not exit then make one */
                        if ( param_get("TORSION",conf_p->confName, to_complete_atoms[0]->name, &tors) ) {
                            strncpy(resName, conf_p->confName, 3); resName[3] = '\0';
                            if ( param_get("TORSION",resName, to_complete_atoms[0]->name, &tors) ) {
                                memset(&tors,0,sizeof(TORS));
                                strcpy(tors.atom1, atom_p->name);
                                strcpy(tors.atom2, known_atoms[0]->name);
                                for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                                    if (!known_atoms[0]->connect12[i_connect]) break;
                                    if (!known_atoms[0]->connect12[i_connect]->on) continue;
                                    if (!strcmp(known_atoms[0]->connect12[i_connect]->name, atom_p->name)) continue;
                                    strcpy(tors.atom3, known_atoms[0]->connect12[i_connect]->name);
                                    break;
                                }
                                tors.V2[0]     = 0.;
                                tors.n_fold[0] = 3.;
                                tors.gamma[0]  = 0.;
                                tors.opt_hyd = 0;
                                param_sav("TORSION", conf_p->confName, to_complete_atoms[0]->name, &tors, sizeof(TORS));
                                debug_fp = fopen(env.debug_log, "a");
                                fprintf(debug_fp, "TORSION  %s %s %s %s %s  f  %9.1f %9.0f %9.2f\n",
                                conf_p->confName,
                                to_complete_atoms[0]->name,
                                tors.atom1,
                                tors.atom2,
                                tors.atom3,
                                tors.V2[0],
                                tors.n_fold[0],
                                tors.gamma[0]/env.d2r);
                                fclose(debug_fp);
                                error++;
                            }
                        }
                        
                        if (strcmp(atom_p->name, tors.atom1)) {
                            printf("   Error! Atom %s is in torsion parameter of conformer %s atom %s, but connected atom %s was found\n",
                            tors.atom1,conf_p->confName,to_complete_atoms[0]->name,atom_p->name);
                            continue;
                        }
                        if (strcmp(known_atoms[0]->name, tors.atom2)) {
                            printf("   Error! Atom %s is in torsion parameter of conformer %s atom %s, but connected atom %s was found\n",
                            tors.atom2,conf_p->confName,to_complete_atoms[0]->name,known_atoms[0]->name);
                            continue;
                        }
                        back_atom_p = NULL;
                        for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                            if (!known_atoms[0]->connect12[i_connect]) break;
                            if (!known_atoms[0]->connect12[i_connect]->on) continue;
                            if (!strcmp(known_atoms[0]->connect12[i_connect]->name, tors.atom3)) {back_atom_p = known_atoms[0]->connect12[i_connect]; break;}
                        }
                        if (!back_atom_p) {
                            printf("   Error! place_missing(): cannot put hydrogen on atom %s in residue %s %d, because %s-%s- N/A is unexpected.\n",
                            atom_p->name, res_p->resName, res_p->resSeq, atom_p->name, known_atoms[0]->name);
                            continue;
                        }
                        
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        bond_angle  = get_bond_angle(conf_p, atom_p, known_atoms[0], to_complete_atoms[0], orbital);
                        
                        n_fold = tors.n_fold[0];
                        /*
                        if (tors.n_fold[0] != 1 && handle_addconf ==2) {
                            n_fold = N_ROT_REBUILD;
                        }
                        */
                        for (i_fold=0; i_fold<n_fold; i_fold++) {
                            if (i_fold) {
                                if (!i_conf) break; /* do not add extra conf for backbone */
                                if (!handle_addconf) break; /* do not add extra conf if the flag is 0 */
                                
                                if (to_complete_atoms[0]->name[1] == 'H') {
                                    if (to_complete_atoms[1]->name[1] == 'H') {
                                        if (to_complete_atoms[2]->name[1] == 'H') break; /* do not add extra conf if added atoms are all protons */
                                    }
                                }

                                ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                                if (ins == USERERR) return USERERR;
                                conf_p = &prot.res[i_res].conf[i_conf];
                                if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                                get_connect12_conf(i_res, res_p->n_conf-1, prot);
                                //printf("   Debugging! residue %s%4d,conformer %s to conformer %s\n", res_p->resName,res_p->resSeq,conf_p->confName,res_p->conf[res_p->n_conf-1].confName);
                            }
                            
                            torsion_angle = (env.PI + tors.gamma[0] + (i_fold)*2.*env.PI)/n_fold;
                            sp3_1known(atom_p->xyz,
                            known_atoms[0]->xyz,
                            back_atom_p->xyz,
                            &to_complete_atoms[0]->xyz,
                            NULL,
                            NULL,
                            bond_length,
                            bond_angle,
                            torsion_angle );
                            
                            //printf("   Debugging! i_fold %d, n_fold %f, residue %s%4d,history %s, pos %d\n",i_fold,tors.n_fold[0], res_p->resName,res_p->resSeq,conf_p->history,i_conf);
                            //printf("   Debugging2! residue%d %s%4d,conformer%d %s and last conformer %s\n", i_res, res_p->resName,res_p->resSeq,i_conf,conf_p->confName,res_p->conf[res_p->n_conf-1].confName);
                            to_complete_atoms[0]->on = 1;
                            n_added++;
                            get_connect12_conf(i_res, i_conf, prot);
                        }
                    }
                    else if (n_known == 0) {
                        int counter;
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        /* Define four corners of the box, atom would be placed on the conner. */
                        a = bond_length/sqrt(3.);
                        counter = 0;
                        v.x =  a; v.y =  a; v.z =  a; corners[0] = vector_vplusv(atom_p->xyz, v);
                        v.x = -a; v.y = -a; v.z =  a; corners[1] = vector_vplusv(atom_p->xyz, v);
                        v.x = -a; v.y =  a; v.z = -a; corners[2] = vector_vplusv(atom_p->xyz, v);
                        v.x =  a; v.y = -a; v.z = -a; corners[3] = vector_vplusv(atom_p->xyz, v);
                        for (i_corner = 0; i_corner < 4; i_corner++) {
                            
                            if (n_complete <=1) start = 3; else start = i_corner+1;
                            for (j_corner = start; j_corner < 4; j_corner++) {
                                
                                if (n_complete <=2) start = 3; else start = j_corner+1;
                                for (k_corner = start; k_corner < 4; k_corner++) {
                                    
                                    if (n_complete <=3) start = 3; else start = k_corner+1;
                                    for (l_corner = start; l_corner < 4; l_corner++) {
                                        
                                        to_complete_atoms[0]->xyz = corners[i_corner];
                                        to_complete_atoms[1]->xyz = corners[j_corner];
                                        to_complete_atoms[2]->xyz = corners[k_corner];
                                        to_complete_atoms[3]->xyz = corners[l_corner];
                                        
                                        to_complete_atoms[0]->on = 1;
                                        to_complete_atoms[1]->on = 1;
                                        to_complete_atoms[2]->on = 1;
                                        to_complete_atoms[3]->on = 1;
                                        n_added++;
                                        
                                        if (!handle_addconf) continue;
                                        
                                        ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                                        if (ins == USERERR) return USERERR;
                                        conf_p = &prot.res[i_res].conf[i_conf];
                                        if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                                        get_connect12_conf(i_res, i_conf, prot);
                                    }
                                }
                            }
                        }
                        
                        v.x = -a; v.y = -a; v.z = -a; corners[0] = vector_vplusv(atom_p->xyz, v);
                        v.x =  a; v.y =  a; v.z = -a; corners[1] = vector_vplusv(atom_p->xyz, v);
                        v.x =  a; v.y = -a; v.z =  a; corners[2] = vector_vplusv(atom_p->xyz, v);
                        v.x = -a; v.y =  a; v.z =  a; corners[3] = vector_vplusv(atom_p->xyz, v);
                        for (i_corner = 0; i_corner < 4; i_corner++) {
                            
                            if (n_complete <=1) start = 3; else start = i_corner+1;
                            for (j_corner = start; j_corner < 4; j_corner++) {
                                
                                if (n_complete <=2) start = 3; else start = j_corner+1;
                                for (k_corner = start; k_corner < 4; k_corner++) {
                                    
                                    if (n_complete <=3) start = 3; else start = k_corner+1;
                                    for (l_corner = start; l_corner < 4; l_corner++) {
                                        
                                        to_complete_atoms[0]->xyz = corners[i_corner];
                                        to_complete_atoms[1]->xyz = corners[j_corner];
                                        to_complete_atoms[2]->xyz = corners[k_corner];
                                        to_complete_atoms[3]->xyz = corners[l_corner];
                                        
                                        to_complete_atoms[0]->on = 1;
                                        to_complete_atoms[1]->on = 1;
                                        to_complete_atoms[2]->on = 1;
                                        to_complete_atoms[3]->on = 1;
                                        n_added++;
                                        
                                        if (!handle_addconf) continue;
                                        ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                                        if (ins == USERERR) return USERERR;
                                        conf_p = &prot.res[i_res].conf[i_conf];
                                        if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                                        get_connect12_conf(i_res, res_p->n_conf-1, prot);
                                    }
                                }
                            }
                        }
                        
                        if (handle_addconf) {
                            del_conf(res_p, res_p->n_conf-1);
                        }
                        conf_p = &prot.res[i_res].conf[i_conf];
                        continue;
                    }
                    else {
                        printf("   Error! place_missing(): number of known atoms can't be handled, sp3 with %i known atoms\n",n_known);
                        fatal++;
                    }
                }
                else if (!strcmp(orbital, "sp2")) {
                    if ( (n_known + n_complete) < 3 ) {
                        for (i_dummy = n_complete; i_dummy < (3-n_known); i_dummy++) {
                            to_complete_atoms[i_dummy] = &dummy_atom[i_dummy];
                        }
                    }
                    if (n_known == 2) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        /* get bond angle from tpl file, force the subroutine not to guess from orbital type */
                        bond_angle = get_bond_angle(conf_p, atom_p, known_atoms[0], to_complete_atoms[0], NULL);
                        //printf("%s %s %s %8.3f\n", atom_p->name, known_atoms[0]->name, to_complete_atoms[0]->name, bond_angle*180./env.PI);
                        if (bond_angle < 0) {
                            sp2_2known( atom_p->xyz,
                            known_atoms[0]->xyz,
                            known_atoms[1]->xyz,
                            &to_complete_atoms[0]->xyz,
                            bond_length);
                        }
                        else {
                            sp2_2known_with_angle( atom_p->xyz,
                            known_atoms[0]->xyz,
                            known_atoms[1]->xyz,
                            &to_complete_atoms[0]->xyz,
                            bond_length,
                            bond_angle );
                        }
                        
                        to_complete_atoms[0]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res, i_conf, prot);
                    }
                    else if (n_known == 1) {
                        int n_fold;
                        if ( param_get("TORSION",conf_p->confName, to_complete_atoms[0]->name, &tors) ) {
                            strncpy(resName, conf_p->confName, 3); resName[3] = '\0';
                            if ( param_get("TORSION",resName, to_complete_atoms[0]->name, &tors) ) {
                                memset(&tors,0,sizeof(TORS));
                                strcpy(tors.atom1, atom_p->name);
                                strcpy(tors.atom2, known_atoms[0]->name);
                                for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                                    if (!known_atoms[0]->connect12[i_connect]) break;
                                    if (!known_atoms[0]->connect12[i_connect]->on) continue;
                                    if (!strcmp(known_atoms[0]->connect12[i_connect]->name, atom_p->name)) continue;
                                    strcpy(tors.atom3, known_atoms[0]->connect12[i_connect]->name);
                                    break;
                                }
                                tors.V2[0]     = 0.;
                                tors.n_fold[0] = 2.;
                                tors.gamma[0]  = 180.*env.d2r;
                                tors.opt_hyd = 0;
                                param_sav("TORSION", conf_p->confName, to_complete_atoms[0]->name, &tors, sizeof(TORS));
                                debug_fp = fopen(env.debug_log, "a");
                                fprintf(debug_fp, "TORSION  %s %s %s %s %s  f  %9.1f %9.0f %9.2f\n",
                                conf_p->confName,
                                to_complete_atoms[0]->name,
                                tors.atom1,
                                tors.atom2,
                                tors.atom3,
                                tors.V2[0],
                                tors.n_fold[0],
                                tors.gamma[0]/env.d2r);
                                fclose(debug_fp);
                                error++;
                            }
                        }
                        if (strcmp(atom_p->name, tors.atom1)) {
                            printf("   Error! Cannot add atom \'%s\' onto \'%s\' because TORSION parameter of atom \'%s\' doesn't match the connected atoms in this protein, double check the TORSION parameter. (Res: %s %c%04d)\n",
                            to_complete_atoms[0]->name,atom_p->name,to_complete_atoms[0]->name,res_p->resName,res_p->chainID,res_p->resSeq);
                            continue;
                        }
                        if (strcmp(known_atoms[0]->name, tors.atom2)) {
                            printf("   Error! Cannot add atom \'%s\' onto \'%s\' because TORSION parameter of atom \'%s\' doesn't match the connected atoms in this protein, double check the TORSION parameter. (Res: %s %c%04d)\n",
                            to_complete_atoms[0]->name,atom_p->name,to_complete_atoms[0]->name,res_p->resName,res_p->chainID,res_p->resSeq);
                            continue;
                        }
                        back_atom_p = NULL;
                        for (i_connect=0; i_connect < MAX_CONNECTED; i_connect++) {
                            if (!known_atoms[0]->connect12[i_connect]) break;
                            if (!known_atoms[0]->connect12[i_connect]->on) continue;
                            if (!strcmp(known_atoms[0]->connect12[i_connect]->name, tors.atom3)) {back_atom_p = known_atoms[0]->connect12[i_connect]; break;}
                        }
                        if (!back_atom_p) {
                            printf("   Error! Cannot add atom \'%s\' onto \'%s\' because TORSION parameter of atom \'%s\' doesn't match the connected atoms in this protein, double check the TORSION parameter. (Res: %s %c%04d)\n",
                            to_complete_atoms[0]->name,atom_p->name,to_complete_atoms[0]->name,res_p->resName,res_p->chainID,res_p->resSeq);
                            continue;
                        }
                        
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        bond_angle = get_bond_angle(conf_p, atom_p, known_atoms[0], to_complete_atoms[0], orbital);
                        //bond_angle = 120.;
                        //printf("%s %s %s %8.3f\n", atom_p->name, known_atoms[0]->name, to_complete_atoms[0]->name, bond_angle*180./env.PI);
                        n_fold = tors.n_fold[0];
                        /*
                        if (tors.n_fold[0] != 1 && handle_addconf ==2) {
                            n_fold = N_ROT_REBUILD;
                        }
                        */
                        for (i_fold=0; i_fold<n_fold; i_fold++) {
                            if (i_fold) {
                                if (!i_conf) break; /* do not add extra conf for backbone */
                                if (!handle_addconf) break; /* do not add extra conf if the flag is 0 */
                                
                                if (to_complete_atoms[0]->name[1] == 'H') {
                                    if (to_complete_atoms[1]->name[1] == 'H') break; /* do not add extra conf if added atoms are all protons */
                                }
                                
                                ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                                if (ins == USERERR) return USERERR;
                                conf_p = &prot.res[i_res].conf[i_conf];
                                if (cpy_conf(&res_p->conf[res_p->n_conf-1], conf_p)) {printf("   Error! place_missing(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf_p->confName,res_p->resName, res_p->resSeq, res_p->n_conf-1); fatal++;}
                                get_connect12_conf(i_res, res_p->n_conf-1, prot);
                            }
                            torsion_angle = (env.PI + tors.gamma[0] + (i_fold)*2.*env.PI)/n_fold;
                            sp2_1known(atom_p->xyz,
                            known_atoms[0]->xyz,
                            back_atom_p->xyz,
                            &to_complete_atoms[0]->xyz,
                            bond_length,
                            bond_angle,
                            torsion_angle );
                            
                            to_complete_atoms[0]->on = 1;
                            n_added++;
                            get_connect12_conf(i_res,i_conf,prot);
                        }
                        
                        /* only 1 atom is added, rollback to add the second */
                        i_atom--;
                        continue;
                    }
                    else if (n_known == 0) {
                        printf("   Error! no rule to add atom \'%s\' onto \'%s\': (atom %s in conformer %s have no atom connected)\n",to_complete_atoms[0]->name, atom_p->name, atom_p->name,conf_p->confName);
                    }
                    else {
                        printf("   Error! place_missing(): number of connected atoms wrong, %i n_connected atoms\n",n_known);
                        fatal++;
                    }
                }
                else if (!strcmp(orbital, "sp2d")) {
                    if ( (n_known + n_complete) < 4 ) {
                        for (i_dummy = n_complete; i_dummy < (4-n_known); i_dummy++) {
                            to_complete_atoms[i_dummy] = &dummy_atom[i_dummy];
                        }
                    }
                    
                    if (n_known == 3) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        
                        sp2d_3known( atom_p->xyz,
                        known_atoms[0]->xyz,
                        known_atoms[1]->xyz,
                        known_atoms[2]->xyz,
                        &to_complete_atoms[0]->xyz,
                        bond_length );
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res,i_conf,prot);
                    }
                    else {
                        printf("   Error! place_missing(): number of known atoms can't be handled, sp2d with %i known atoms\n",n_known);
                        fatal++;
                    }
                }
                else if (!strcmp(orbital, "sp3d2")) {
                    if ( (n_known + n_complete) < 6 ) {
                        for (i_dummy = n_complete; i_dummy < (6-n_known); i_dummy++) {
                            to_complete_atoms[i_dummy] = &dummy_atom[i_dummy];
                        }
                    }
                    
                    if (n_known == 5) {
                        bond_length = get_bond_length(conf_p,atom_p,to_complete_atoms[0]);
                        
                        sp3d2_5known( atom_p->xyz,
                        known_atoms[0]->xyz,
                        known_atoms[1]->xyz,
                        known_atoms[2]->xyz,
                        known_atoms[3]->xyz,
                        known_atoms[4]->xyz,
                        &to_complete_atoms[0]->xyz,
                        bond_length );
                        
                        for (i_complete=0; i_complete<n_complete; i_complete++) to_complete_atoms[i_complete]->on = 1;
                        n_added++;
                        get_connect12_conf(i_res,i_conf,prot);
                    }
                    else {
                        printf("   Error! place_missing(): number of known atoms can't be handled, sp3d2 with %i known atoms\n",n_known);
                        fatal++;
                    }
                }
                else {
                    printf("   Error! place_missing(): Don't know yet how to deal with orbital type %s for atom %s in conformer %s\n",orbital,atom_p->name, conf_p->confName);
                }
                /* delete dulipcates probably causing memory trouble
                for (k_conf=prot.res[i_res].n_conf-1; k_conf>=1; k_conf--) {
                    int j_conf;
                    for (j_conf=1; j_conf<k_conf; j_conf++) {
                        if( strcmp(prot.res[i_res].conf[k_conf].confName, prot.res[i_res].conf[j_conf].confName)) continue;
                        if(!cmp_conf(prot.res[i_res].conf[k_conf], prot.res[i_res].conf[j_conf], 0.01)) {
                            del_conf(&prot.res[i_res], k_conf);
                            break;
                        }
                    }
                }
                */
                /* END of this atom */
            }
        }
        
        while (1) {
            int  counter, j_conf;
            for (i_conf=0; i_conf<res_p->n_conf; i_conf++) {
                if (!strncmp(res_p->conf[i_conf].history+6,"____",4)) break;
            }
            if (i_conf >= prot.res[i_res].n_conf) break;
            counter = 0;
            for (j_conf=0; j_conf<res_p->n_conf; j_conf++) {
                if (strncmp(res_p->conf[i_conf].history,res_p->conf[j_conf].history,6)) continue;
                res_p->conf[j_conf].history[6] = 'M';
                sprintf(sbuff,"%03d",counter);
                strncpy(res_p->conf[j_conf].history+7,sbuff,3);
                counter++;
            }
        }
    
    /* Check for if there are still missing atoms after protonation. the way of treating NTR and CTR here is not good and standard */
    if (!n_added) {
        Missing = 0;
            for (kc=0; kc<prot.res[i_res].n_conf; kc++) {
                for (ka=0; ka<prot.res[i_res].conf[kc].n_atom; ka++) {
                    if (!prot.res[i_res].conf[kc].atom[ka].on) {
                        sprintf(siatom, "%3d", ka);
                        if (param_get("ATOMNAME", prot.res[i_res].conf[kc].confName, siatom, sbuff)) {
                            printf("   Missing ATOMNAME records for slot \"%d\" of conformer %s.\n",
                            ka, prot.res[i_res].conf[kc].confName);
                            Missing++;
                        }
                        else {
                            if (!strcmp(prot.res[i_res-1].resName, "NTR") || !strcmp(prot.res[i_res-1].resName, "NTG")) {
                                if (!strncmp(sbuff+1,"HA",2)) sbuff[0] = ' ';
                                if (!strncmp(sbuff+1,"H ",2)) sbuff[0] = '1';
                                if (!strncmp(sbuff+1,"H\0",2)) sbuff[0] = '1';
                                if (!param_get("IATOM","NTRBK",sbuff,sbuff2)) continue;
                                if (!param_get("IATOM","NTR01",sbuff,sbuff2)) continue;
                                if (!param_get("IATOM","NTGBK",sbuff,sbuff2)) continue;
                                if (!param_get("IATOM","NTG01",sbuff,sbuff2)) continue;
                            }
                            if (!strcmp(prot.res[i_res+1].resName, "CTR")) {
                                sbuff[0] = ' ';
                                if (!param_get("IATOM","CTRBK",sbuff,sbuff2)) continue;
                                if (!param_get("IATOM","CTR01",sbuff,sbuff2)) continue;
                            }
                            printf("   Missing atom %s at slot %4d of conf %s in \"%s %c %3d %02d\".\n",
                            sbuff, ka, prot.res[i_res].conf[kc].confName, prot.res[i_res].resName, prot.res[i_res].chainID, prot.res[i_res].resSeq,kc);
                            Missing++;
                        }
                    }
                }
            }
    }
    if (error) {
       printf("   %d TORSION parameters were guessed and recorded in file %s.\n", error, env.debug_log);
       printf("   Modify and put them into a param file to have the full control.\n" );
    }
    return n_added;
}

void sp3_3known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR r3, VECTOR *r4, float bond_len_04) {

    /*
    .     r4
    .     |
    .     r0
    .    /|\
    .  r1 r2 r3
    .
    . r0 is sp3 type.
    . r0, r1, r2, r3's coordinates are known.
    . r4's coordinate is to be determined.
    */

    VECTOR n01;
    VECTOR n02;
    VECTOR n03;
    VECTOR n04;
    
    n01 = vector_normalize(vector_vminusv(r1, r0));
    n02 = vector_normalize(vector_vminusv(r2, r0));
    n03 = vector_normalize(vector_vminusv(r3, r0));
    n04 = vector_normalize(vector_neg(vector_sum3v(n01,n02,n03)));

    if (r4) *r4  = vector_vplusv(r0, vector_rescale(n04, bond_len_04));
}

void sp3_2known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, VECTOR *r4, float bond_len_03, float bond_angle_304) {
    
    /*
    .    r3 r4
    .     \/
    .     r0
    .    /|
    .  r1 r2
    .
    .
    . r0 is sp3 type.
    . r0, r1, r2's coordinates are known.
    . r3, r4's coordinates are to be determined.
    */

    float half_angle;

    VECTOR n01;
    VECTOR n02;
    VECTOR bisect304;	/* bisector of angle r3-r0-r4 */
    VECTOR norm102;	/* normal of plane r1-r0-r2 */

    n01 = vector_normalize(vector_vminusv(r1, r0));
    n02 = vector_normalize(vector_vminusv(r2, r0));

    bisect304 = vector_normalize(vector_neg(vector_vplusv(n01,n02)));

    norm102 = vector_normalize(vector_vxv(n01, n02));

    half_angle = bond_angle_304 / 2.;
    
    if (r3) *r3 = vector_vplusv(
                        r0,
                        vector_rescale(
                                vector_vplusv(
                                        vector_rescale(bisect304,cos(half_angle)),
                                        vector_rescale(norm102,sin(half_angle))),
                                bond_len_03));
    
    if (r4) *r4 = vector_vplusv(
                        r0,
                        vector_rescale(
                                vector_vminusv(
                                        vector_rescale(bisect304,cos(half_angle)),
                                        vector_rescale(norm102,sin(half_angle))),
                                bond_len_03));
    

    //if (r3) *r3 = r0[i] + bond_len_03 * (bisect304[i] * cos(half_angle) + norm102[i] * sin(half_angle));
    //if (r4) *r4 = r0[i] + bond_len_03 * (bisect304[i] * cos(half_angle) - norm102[i] * sin(half_angle));
}

void sp3_1known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, VECTOR *r4, VECTOR *r5, float bond_len_03, float bond_angle_103, float torsion_angle_3012) {
    
    /*
    .	r3 r4 r5
    .	  \|/
    .	   r0			  z (n10)
    .	   |			  |_ y (norm_0_1_n012)
    .	   r1			 /
    .	   \			x (norm012)
    .	    r2
    . r0 is sp3 type.
    . r0, r1, r2's coordinates are known.
    . r3, r4, r5's coordinates are to be determined.
    */


    VECTOR n10;
    VECTOR n12;
    VECTOR norm012;	/* normal of plane r0-r1-r2 */
    VECTOR norm_0_1_n012;	/* normal of plane r0-r1-norm012 */

    float theta;
    float phi;

    n10 = vector_normalize(vector_vminusv(r0, r1));
    n12 = vector_normalize(vector_vminusv(r2, r1));

    norm012 = vector_normalize(vector_vxv(n10,n12));
    norm_0_1_n012 = vector_normalize(vector_vxv(n10,norm012));
    
    theta = env.PI - bond_angle_103;
    if (r3) {
        phi = torsion_angle_3012 - env.PI/2.;
        *r3 = vector_vplusv(r0, vector_rescale(
                vector_sum3v(   vector_rescale(norm012,cos(phi)*sin(theta)),
                                vector_rescale(norm_0_1_n012,sin(phi)*sin(theta)),
                                vector_rescale(n10,cos(theta))),
                bond_len_03));
    }
    if (r4) {
        phi = torsion_angle_3012 + env.PI/6.;
        *r4 = vector_vplusv(r0, vector_rescale(
                vector_sum3v(   vector_rescale(norm012,cos(phi)*sin(theta)),
                                vector_rescale(norm_0_1_n012,sin(phi)*sin(theta)),
                                vector_rescale(n10,cos(theta))),
                bond_len_03));
    }
    if (r5) {
        phi = torsion_angle_3012 + 5.*env.PI/6.;
        *r5 = vector_vplusv(r0, vector_rescale(
                vector_sum3v(   vector_rescale(norm012,cos(phi)*sin(theta)),
                                vector_rescale(norm_0_1_n012,sin(phi)*sin(theta)),
                                vector_rescale(n10,cos(theta))),
                bond_len_03));
    }

    /* ^       ^                      ^                      ^              */
    /* r = r ( i*cos(phi)sin(theta) + j*sin(phi)sin(theta) + k*cos(theta) ) */
    /*
    if (r3) {
        for (i=0; i<3; i++) {
            r3[i] = r0[i] + bond_len_03 * (norm012[i]*cos(phi)*sin(theta)
            +                        norm_0_1_n012[i]*sin(phi)*sin(theta)
            +                                  n10[i]         *cos(theta) );
        }
    }

    if (r4) {
	phi = ( torsion_angle_3012 + 30 );
        for (i=0; i<3; i++) {
            r4[i] = r0[i] + bond_len_03 * (norm012[i]*cos(phi)*sin(theta) 
            +                        norm_0_1_n012[i]*sin(phi)*sin(theta) 
            +                                  n10[i]         *cos(theta) );
        }
    }

    if (r5) {
	phi = ( torsion_angle_3012 + 150 );
        for (i=0; i<3; i++) {
            r5[i] = r0[i] + bond_len_03 * (norm012[i]*cos(phi)*sin(theta)
            +                        norm_0_1_n012[i]*sin(phi)*sin(theta) 
            +                                  n10[i]         *cos(theta) );
        }
    }
    */
}

void sp2_2known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, float bond_len_03) {
    
    /*
    .     r3
    .     |
    .     r0
    .    / \
    .  r1   r2
    . 
    . r0 is sp3 type.
    . r0, r1, r2's coordinates are known.
    . r3's coordinate is to be determined.
    */
    
    VECTOR n01;
    VECTOR n02;
    VECTOR n03;

    n01 = vector_normalize(vector_vminusv(r1, r0));
    n02 = vector_normalize(vector_vminusv(r2, r0));
    n03 = vector_normalize(vector_neg(vector_vplusv(n01,n02)));

    if (r3) *r3  = vector_vplusv(r0, vector_rescale(n03, bond_len_03));
}

void sp2_2known_with_angle(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, float bond_len_03, float bond_angle_103) {
    
    /*
    .     r3
    .     |               y (norm_0_1_n012)
    .     r0               \
    .    / \               /
    .  r1   r2  	      x (n01)
    . 
    . r0 is sp3 type.
    . r0, r1, r2's coordinates are known.
    . r3's coordinate is to be determined.
    */
    
    VECTOR n01;
    VECTOR n02;
    VECTOR n03;
    VECTOR norm012;         /* normal of plane r0-r1-r2 */
    VECTOR norm_0_1_n012;	/* normal of plane r0-r1-norm012 */

    n01 = vector_normalize(vector_vminusv(r1, r0));
    n02 = vector_normalize(vector_vminusv(r2, r0));
    norm012  = vector_normalize(vector_vxv(n01,n02));
    norm_0_1_n012  = vector_normalize(vector_vxv(n01,norm012));
    n03 = vector_vplusv(vector_rescale(n01, cos(bond_angle_103)), vector_rescale(norm_0_1_n012, sin(bond_angle_103)));
    /* ^     ^                 ^              */
    /* n03 = i * cos(theta) +  j * sin(theta) */

    *r3  = vector_vplusv(r0, vector_rescale(n03, bond_len_03));
}

void sp2_1known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR *r3, float bond_len_03, float bond_angle_103, float torsion_angle_3012) {
    
    /*
    .	r3   r4
    .	  \ /
    .	   r0			  z (n10)
    .	   |			  |_ y (norm_0_1_n012)
    .	   r1			 /
    .	   \			x (norm012)
    .	    r2
    . r0 is sp3 type.
    . r0, r3, r4's coordinates are known.
    . r1, r2's coordinates are to be determined.
    
    */

    VECTOR n10;
    VECTOR n12;
    VECTOR norm012;	/* normal of plane r0-r1-r2 */
    VECTOR norm_0_1_n012;	/* normal of plane r0-r1-norm012 */

    float theta;
    float phi;

    n10 = vector_normalize(vector_vminusv(r0, r1));
    n12 = vector_normalize(vector_vminusv(r2, r1));

    norm012 = vector_normalize(vector_vxv(n10,n12));
    norm_0_1_n012 = vector_normalize(vector_vxv(n10,norm012));

    theta = env.PI - bond_angle_103;
    phi = torsion_angle_3012 - env.PI/2.;
    *r3 = vector_vplusv(r0, vector_rescale(
    vector_sum3v(   vector_rescale(norm012,cos(phi)*sin(theta)),
    vector_rescale(norm_0_1_n012,sin(phi)*sin(theta)),
    vector_rescale(n10,cos(theta))),
    bond_len_03));

    /* ^       ^                      ^                      ^              */
    /* r = r ( i*cos(phi)sin(theta) + j*sin(phi)sin(theta) + k*cos(theta) ) */
    /*
    for (i=0; i<3; i++) {
        r3[i] = r0[i] + bond_len_03 * (norm012[i]*cos(phi)*sin(theta)
        +                        norm_0_1_n012[i]*sin(phi)*sin(theta) 
        +                                  n10[i]         *cos(theta) );
    }

    if (r4) {
        phi = (torsion_angle_3012 + 90);
        for (i=0; i<3; i++) {
            r4[i] = r0[i] + bond_len_03 * (norm012[i]*cos(phi)*sin(theta)
            +                        norm_0_1_n012[i]*sin(phi)*sin(theta)
            +                                  n10[i]         *cos(theta) );
        }
    }
    */
}

void sp2d_3known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR r3, VECTOR *r4, float bond_len_04) {

    /*
    .     r4
    .     |
    .  r3-r0-r1
    .     |
    .     r2
    .
    . r0 is sp2d/dsp2 type.
    . r0, r1, r2, r3's coordinates are known.
    . r4's coordinate is to be determined.
    */

    VECTOR n01, n02, n03, n04, n13, swp;
    float product;
    int on_a_line;

    n01 = vector_normalize(vector_vminusv(r1, r0));
    n02 = vector_normalize(vector_vminusv(r2, r0));
    n03 = vector_normalize(vector_vminusv(r3, r0));

    /* find the two vectors on a line */
    on_a_line = 12;
    product = vdotv(n01, n02);
    if (vdotv(n01, n03) < product) {
        on_a_line = 13;
        product = vdotv(n01, n03);
    }
    if (vdotv(n02, n03) < product) {
        on_a_line = 23;
        product = vdotv(n02, n03);
    }
    /* make n01 and n03 on a line */
    if (on_a_line == 12) {
        swp = n02;
        n02 = n03;
        n03 = swp;
    }
    if (on_a_line == 23) {
        swp = n02;
        n02 = n01;
        n01 = swp;
    }
    
    /* consider if r1,r2,r3 are a little bended, try to guess the direction for r4 */
    n13 = vector_normalize(vector_vminusv(n03,n01)); /* average of n01 and n03 (if they are not exactly on a line */
    n02 = vector_normalize(vector_vminusv(n02,vector_rescale(n02, vdotv(n02,n13)))); /* n02 should be normal to n13 */
    n04 = vector_neg(n02);
    
    if (r4) *r4  = vector_vplusv(r0, vector_rescale(n04, bond_len_04));
}

void sp3d2_5known(VECTOR r0, VECTOR r1, VECTOR r2, VECTOR r3, VECTOR r4, VECTOR r5, VECTOR *r6, float bond_len_06) {

    /*
    .      r4  r3
    .       \ /
    .  r5-- r0 --r6
    .       / \
    .      r1  r2
    .
    . r0 is sp3d2 type.
    . r0, r1, r2, r3, r4, r5's coordinates are known.
    . r6's coordinate is to be determined.
    . In this figure, r0,r1,r2,r3,r4 are in a same plane (alpha)
    . r05, r06 are both normal to the plane alpha.
    */

    VECTOR n[6],n13,n24;
    float product, min_product,p13,p24;
    int i,j,perpend;

    n[1] = vector_normalize(vector_vminusv(r1, r0));
    n[2] = vector_normalize(vector_vminusv(r2, r0));
    n[3] = vector_normalize(vector_vminusv(r3, r0));
    n[4] = vector_normalize(vector_vminusv(r4, r0));
    n[5] = vector_normalize(vector_vminusv(r5, r0));

    /* find the four vectors on a plane (or find one vector that's perpendicular to all the others) */
    min_product = 1;
    for (i=1;i<6;i++) {
        product = 0;
        for (j=1;j<6;j++) {
            if (i==j) continue;
            /* remember the angle btw i and j that's lest close to 90 */
            if (fabs(vdotv(n[i],n[j])) > product) {
                product = fabs(vdotv(n[i],n[j]));
            }
        }
        if (product < min_product) {
            min_product = product;
            perpend = i;
        }
    }
    
    /* find r1-r3 and r2-r4 */
    p13 = 0;
    p24 = 0;
    for (i=1;i<6;i++) {
        if (i == perpend) continue;
        for (j=i+1;j<6;j++) {
            if (j == perpend) continue;
            
            if (fabs(vdotv(n[i],n[j])) > p13) {
                p24 = p13;
                n24 = n13;
                p13 = fabs(vdotv(n[i],n[j]));
                n13 = vector_normalize(vector_vminusv(n[i], n[j]));
            }
            else if (fabs(vdotv(n[i],n[j])) > p24) {
                p24 = fabs(vdotv(n[i],n[j]));
                n24 = vector_normalize(vector_vminusv(n[i], n[j]));
            }
        }
    }
    
    n[0] = vector_normalize(vector_vxv(n13,n24));
    if (vdotv(n[perpend],n[0]) > 0) {
        n[0] = vector_neg(n[0]);
    }

    if (r6) *r6  = vector_vplusv(r0, vector_rescale(n[0], bond_len_06));
}

VECTOR vector_rescale(VECTOR v, float c) {
    VECTOR z;
    z.x = v.x * c;
    z.y = v.y * c;
    z.z = v.z * c;
    return z;
}

VECTOR vector_sum3v(VECTOR v1, VECTOR v2, VECTOR v3) {
    VECTOR z;
    z.x = v1.x + v2.x + v3.x;
    z.y = v1.y + v2.y + v3.y;
    z.z = v1.z + v2.z + v3.z;
    return z;
}

VECTOR vector_neg(VECTOR v) {
    VECTOR z;
    z.x = -v.x;
    z.y = -v.y;
    z.z = -v.z;
    return z;
}

#define DEFAULT_RAD 0.75
float get_bond_length(CONF *conf_p, ATOM *atom1_p, ATOM *atom2_p) {
    float rad1=0;
    float rad2=0;
    char  elem1[5]="    ";
    char  elem2[5]="    ";
    char  bond_len[MAXCHAR_LINE], *sbuffer, resName[4];
    float len;
    int   exist;

    /* first search for specifically defined bond length */
    exist = 1;
    if ( param_get("BOND_LEN",conf_p->confName, atom1_p->name, bond_len) ) {
        strncpy(resName, conf_p->confName, 3); resName[3] = '\0';
        if ( param_get("BOND_LEN",resName, atom1_p->name, bond_len) ) {
            exist = 0;
        }
    }
    if (exist) {
        sbuffer = strstr(bond_len, atom2_p->name);
        if (sbuffer) {
            sbuffer[11] = '\0';
            len = atof(sbuffer + 5);
            return len;
        }
    }
    
    exist = 1;
    if ( param_get("BOND_LEN",conf_p->confName, atom2_p->name, bond_len) ) {
        strncpy(resName, conf_p->confName, 3); resName[3] = '\0';
        if ( param_get("BOND_LEN",resName, atom2_p->name, bond_len) ) {
            exist = 0;
        }
    }
    if (exist) {
        sbuffer = strstr(bond_len, atom1_p->name);
        if (sbuffer) {
            sbuffer[11] = '\0';
            len = atof(sbuffer + 5);
            return len;
        }
    }

    /* then use covalence radii to get bond length */
    elem1[0]=atom1_p->name[1];
    if (param_get("RADCOVAL", "SINGL", elem1, &rad1)) {
        printf("   Warning! get_bond_length(): Can't find covalent radius for %s, use default radius %f instead\n", atom1_p->name, DEFAULT_RAD);
        rad1 = DEFAULT_RAD;
        param_sav("RADCOVAL", "SINGL", elem1, &rad1, sizeof(float));
    }
    elem2[0]=atom2_p->name[1];
    if (param_get("RADCOVAL", "SINGL", elem2, &rad2)) {
        printf("   Warning! get_bond_length(): Can't find covalent radius for %s, use default radius %f instead\n", atom2_p->name, DEFAULT_RAD);
        rad2 = DEFAULT_RAD;
        param_sav("RADCOVAL", "SINGL", elem2, &rad2, sizeof(float));
    }
    return rad1+rad2;
}

float get_bond_angle(CONF *conf_p, ATOM *atom0_p, ATOM *atom1_p, ATOM *atom2_p, char *orbital) {
    char  bond_ang[MAXCHAR_LINE], *sbuff, *sbuffer1, *sbuffer2, resName[4];
    int   tpl_exist;
    
    tpl_exist = 1;
    if ( param_get("BOND_ANG",conf_p->confName, atom0_p->name, bond_ang) ) {
        strncpy(resName, conf_p->confName, 3); resName[3] = '\0';
        if ( param_get("BOND_ANG",resName, atom0_p->name, bond_ang) ) {
            tpl_exist = 0;
        }
    }
    
    /*
    .....................01234567890123456789012345678901234567890
    .0123456789012345678901234567890123456789012345678901234567890
    .                   |one block of def   |
    .#-------|-----|----|----|----|---------|----|----|---------|
    .BOND_ANG HEA   FE    N A  OH   109.0
    */
    if (tpl_exist) {
        sbuff = bond_ang;
        while (1) { /* keep searching */
            sbuffer1 = strstr(sbuff, atom1_p->name); /* look for atom1 */
            sbuffer2 = strstr(sbuff, atom2_p->name); /* look for atom2 */
            if (sbuffer1 && sbuffer2) {  /* both are found */
                if (fabs(strlen(sbuffer1)-strlen(sbuffer2)) <= 6) { /* within one block of angle definition */
                    /* get angle parameter */
                    if (strlen(sbuffer1)<strlen(sbuffer2)) {
                        sbuffer1[11] = '\0';
                        return env.d2r * atof(sbuffer1+5);
                    }
                    else {
                        sbuffer2[11] = '\0';
                        return env.d2r * atof(sbuffer2+5);
                    }
                }
                else { /* not in one block of definition, keep searching */
                    if (strlen(sbuffer1)<strlen(sbuffer2)) sbuff = sbuffer2+6;
                    else sbuff = sbuffer1+6;
                }
            }
            else break;
        }
    }
    
    if (!orbital) return -1;
    
    if (!strcmp(orbital, "sp2"))         return env.d2r * 120.;
    else if (!strcmp(orbital, "sp3"))    return env.d2r * 109.;
    else if (!strcmp(orbital, "sp3d2"))  return env.d2r * 90.;
    else if (!strcmp(orbital, "sp2d"))   return env.d2r * 90.;
    else {
        return env.d2r * 109;
    }
}

int ionization(PROT prot)
{  int kr, kc, ka, ic, ia, natoms, ins;
   STRINGS confs;
   ATOM *swap;

   for (kr=0; kr< prot.n_res; kr++) {
      if (param_get("CONFLIST", prot.res[kr].resName, "", &confs)) {
         printf("   FATAL: ionization(): No CONFLIST for residue %s\n", prot.res[kr].resName);
         return USERERR;
      }
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
          prot.res[kr].conf[kc].tmp_flag = 1; /* mark the existing conformers */
      }
      
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
          if (!prot.res[kr].conf[kc].n_atom) continue;
          if (!prot.res[kr].conf[kc].tmp_flag) continue; /* do not duplicate for conformers created here */
          
          /* each conformer will be multipied to the conformers in the conflist */
         for (ic=1; ic<confs.n; ic++) {
            if (strcmp(prot.res[kr].conf[kc].confName, confs.strings[ic])) {
               if (param_get("NATOM", confs.strings[ic], "", &natoms)) {
                  printf("   FATAL: ionization(): No NATOM for conformer %s\n", confs.strings[ic]);
                  return USERERR;
               }
               //ins = ins_conf(&prot.res[kr], prot.res[kr].n_conf, natoms);
               /* insert behind the parent conformer, instead of at the end, Yifan 01/25/07 */
               ins = ins_conf(&prot.res[kr], kc+1, natoms);
               if (ins == USERERR) return USERERR;

               /* inherite the parent conformer but the atoms */
               swap=prot.res[kr].conf[ins].atom; /* remember the position of new array created by ins_conf */
               memcpy(&prot.res[kr].conf[ins], &prot.res[kr].conf[kc], sizeof(CONF));
               prot.res[kr].conf[ins].tmp_flag = 0; /* new conformers are not marked */
               prot.res[kr].conf[ins].atom = swap; prot.res[kr].conf[ins].n_atom=natoms;

               /* update the confName and history */
               strcpy(prot.res[kr].conf[ins].confName, confs.strings[ic]); /* update confName */
               strncpy(prot.res[kr].conf[ins].history, confs.strings[ic]+3, 2); /* update history */

               for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
                   if (prot.res[kr].conf[kc].atom[ka].name[1] == 'H') {
                       char toggle[MAXCHAR_LINE];
                       if ( !param_get("DEL_HYDR", prot.res[kr].resName, "", toggle) ) { /* let mcce rebuilt H. (delete hydrogen is true) */
                           if (strchr(toggle,'t')) continue;
                           if (strchr(toggle,'T')) continue;
                       }
                       else { /* do not copy hydrogen, let mcce rebuild */
                           continue;
                       }
                   }
                       
                  ia = iatom(prot.res[kr].conf[ins].confName, prot.res[kr].conf[kc].atom[ka].name);
                  if (prot.res[kr].conf[kc].atom[ka].on && ia>=0) {
                     strcpy(prot.res[kr].conf[ins].atom[ia].confName, confs.strings[ic]); /* update confname */
                     strncpy(prot.res[kr].conf[ins].atom[ia].history, confs.strings[ic]+3, 2); /* update history */
                     prot.res[kr].conf[ins].atom[ia] = prot.res[kr].conf[kc].atom[ka];
                  }
               }
            }
         }
      }
   }

   return 0;
}

int rm_dupconf_hv(PROT prot)
{  int kr, kc, ic;
   int C = 0;

   for (kr=0; kr<prot.n_res; kr++) {
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
         for (ic=kc+1; ic<prot.res[kr].n_conf; ic++) {
            if(!cmp_conf_hv(prot.res[kr].conf[kc], prot.res[kr].conf[ic], env.prune_thr)) {
               //printf("debugging, %d,%d,%d\n",prot.res[kr].resSeq,kc,ic);
               del_conf(&prot.res[kr], ic);
               ic--;
               C++;
            }
         }
      }
   }

   return C;
}

int rm_dupconf(PROT prot, float prune_thr) {
    int i_res;
    int rm_counter=0;
    for (i_res=0; i_res<prot.n_res; i_res++) {
        rm_counter += rm_dupconf_res(prot, i_res, prune_thr);
    }
    return rm_counter;
}

int rm_dupconf_res(PROT prot, int i_res, float prune_thr) {
    int rm_counter=0;
    int i_conf;
    for (i_conf=prot.res[i_res].n_conf-1; i_conf>=1; i_conf--) {
        int j_conf;
        for (j_conf=1; j_conf<i_conf; j_conf++) {
            if( strcmp(prot.res[i_res].conf[i_conf].confName, prot.res[i_res].conf[j_conf].confName)) continue;
            if(!cmp_conf(prot.res[i_res].conf[i_conf], prot.res[i_res].conf[j_conf], prune_thr)) {
                del_conf(&prot.res[i_res], i_conf);
                rm_counter++;
                break;
            }
        }
    }
    return rm_counter;
}

int write_conflist(FILE *fp, PROT prot)
{  int kr, kc, ka;

   id_conf(prot);
   assign_vdw_param(prot);
   get_connect12(prot);
   get_vdw0(prot);
   get_vdw1(prot);
   for (kr=0; kr<prot.n_res; kr++) {
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
          prot.res[kr].conf[kc].E_torsion = torsion_conf(&prot.res[kr].conf[kc]);
          prot.res[kr].conf[kc].E_self = prot.res[kr].conf[kc].E_vdw0 + prot.res[kr].conf[kc].E_vdw1 + prot.res[kr].conf[kc].E_torsion;
          if (prot.res[kr].conf[kc].E_vdw0 > 999.) prot.res[kr].conf[kc].E_vdw0 = 999.;
          if (prot.res[kr].conf[kc].E_vdw1 > 999.) prot.res[kr].conf[kc].E_vdw1 = 999.;
          if (prot.res[kr].conf[kc].E_torsion > 999.) prot.res[kr].conf[kc].E_torsion = 999.;
          if (prot.res[kr].conf[kc].E_self > 999.) prot.res[kr].conf[kc].E_self = 999.;
      }
   }
   fprintf(fp, "CONFORMER         crg  vdw0  vdw1  tors  epol  self   occ  Hb     history    ASA\n");
   for (kr=0; kr<prot.n_res; kr++) {
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
         prot.res[kr].conf[kc].netcrg = 0.0;
         for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
            if (!prot.res[kr].conf[kc].atom[ka].on) continue;
            prot.res[kr].conf[kc].netcrg += prot.res[kr].conf[kc].atom[ka].crg;
         }

         if (fabs(prot.res[kr].conf[kc].netcrg - rintf(prot.res[kr].conf[kc].netcrg)) > 0.001) {
            printf("   WARNING: Conformer %s has non integer charge %6.3f\n", prot.res[kr].conf[kc].uniqID,
                                                                        prot.res[kr].conf[kc].netcrg);
            fflush(stdout);
         }

         fprintf(fp, "%s %6.3f%6.2f%6.2f%6.2f%6.2f%6.2f%6.3f  %2s  %10s %6.3f\n",prot.res[kr].conf[kc].uniqID,
                                                                    prot.res[kr].conf[kc].netcrg,
                                                                    prot.res[kr].conf[kc].E_vdw0,
                                                                    prot.res[kr].conf[kc].E_vdw1,
                                                                    prot.res[kr].conf[kc].E_torsion,
                                                                    prot.res[kr].conf[kc].E_epol,
                                                                    prot.res[kr].conf[kc].E_self,
                                                                    prot.res[kr].conf[kc].occ,
                                                                    prot.res[kr].conf[kc].history+8,
                                                                    prot.res[kr].conf[kc].history,
                                                                    prot.res[kr].conf[kc].sas_fraction);
      }
   }

   for (kr=0; kr<prot.n_res; kr++) {
         prot.res[kr].conf[0].netcrg = 0.0;
         for (ka=0; ka<prot.res[kr].conf[0].n_atom; ka++) {
            if (!prot.res[kr].conf[0].atom[ka].on) continue;
            prot.res[kr].conf[0].netcrg += prot.res[kr].conf[0].atom[ka].crg;
         }

         if (fabs(prot.res[kr].conf[0].netcrg - (int) prot.res[kr].conf[0].netcrg) > 0.001) {
            printf("   WARNING: Conformer %s has non integer charge %6.3f\n", prot.res[kr].conf[0].uniqID,
                                                                        prot.res[kr].conf[0].netcrg);
            fflush(stdout);
         }
   }

   return 0;
}

int get_demetri_out(void)
{ char line[MAXCHAR_LINE], num[2]=" ", atomname[6], conformerlable;
  int i;
  FILE *fp, *fp2;
  fp = fopen("step2_out.pdb","r+");
  fp2 = fopen("demetri_out.pdb","w");
  rewind(fp);
  rewind(fp2);

  while(fgets(line, sizeof(line), fp))
      {  if (strncmp("__", line+82, 2) || strncmp("__", line+84, 2)) {
            strncpy(line+7, "AB  ", 4);
         }
         else  {
            strncpy(line+7, "AA  ", 4);
         }
         /*deal with the atom name part in the next */
         strncpy(atomname, line+12, 5);
         atomname[5]=' ';

         if( !strncmp(line+17,"CTR",3) || !strncmp(line+17,"NTR",3))
            atomname[4] = 'T';
         if(atomname[0]<=57 && atomname[0]>=49)
         { atomname[5] = atomname[0];
            atomname[0] = ' ';
         }
         else if(atomname[0] != ' ')
         { for(i=4; i>=0 ; i--)
               atomname[i+1] = atomname[i];
            atomname[0] = ' ';
         }
         for(i=2; i<6; i++)
         { if(' ' == atomname[i-1] && atomname[i]!=' ')
            { atomname[i-1] = atomname[i];
               atomname[i] = ' ';
               i=i-2;
            }
         }
         strncpy(line+12, atomname,5);

/* Special Cases */
         strncpy(atomname, line+13, 6);
         if(!strncmp(atomname,"H ",4))
         { line[14] = 'N';
         }
         else if(!strncmp(atomname,"HXTTCT",6))
         { line[14] = 'T'; line[15] =' ' ; line[16] = ' ';
         }
         else if(!strncmp(atomname,"OXTTCT",6))
         { line[14] = 'T'; line[15] ='2' ; line[16] = ' ';
         }
         else if(!strncmp(atomname,"OT CT",6))
            line[15] = '1';
         else if(!strncmp(atomname,"HT1 NT",6))
            line[15] = '3';
         else if(!strncmp(atomname,"H1A PA",6) || !strncmp(atomname,"H2A PA",6))
         { line[14] = 'T'; line[15] =' ';
         }
         else if(!strncmp(atomname,"H1D PD",6) || !strncmp(atomname,"H2D PD",6))
         { line[14] = 'T'; line[15] =' ';
         }
         else if(!strncmp(atomname,"HE1 GL",6) || !strncmp(atomname,"HE2 GL",6))
         { line[14] = 'T'; line[15] =' ';
         }
         else if(!strncmp(atomname,"HD1 AS",6) || !strncmp(atomname,"HD2 AS",6))
         { line[14] = 'T'; line[15] =' ';
         }

         /*deal with the residue name part in the next*/
         if('+' == line[80])
            line[19] = '+';
         else if ('-' == line[80])
            line[19] = '-';

         /*deal with the conformor number part in the next*/
         for(i=23; ' '==line[i]; i++)
            line[i]='0';
         num[1] = line[29];
         num[0] = line[28];
         conformerlable = atoi(num);
         if(0 == conformerlable)
            conformerlable = 65;
         else if(conformerlable > 0 && conformerlable <= 26)
            conformerlable += 64;
         else if(conformerlable >26 && conformerlable <=52)
            conformerlable += 70;
         else
         { printf("\n   The number of conformers exceeds 52!\n");
            return USERERR;
         }
         line[26] = conformerlable;
         line[29] = ' ';
         line[28] = ' ';
         line[27] = ' ';

         for(i=80; i<90 ; i++)
            line[i] = ' ';

         fprintf(fp2,line);

      }

  rewind(fp);
  rewind(fp2);
  fclose(fp);
  fclose(fp2);
  return 0;
}

int load_listrot(char *fname, PROT prot)
{  FILE *fp;
   int i;
   char line[MAXCHAR_LINE];
   ATOM atom;

   if ((fp=fopen(fname, "r")) == NULL) {
      printf("   Read file %s error\n", fname);
      return USERERR;
   }

   memset(line, 0, sizeof(line));
   while (fgets(line, sizeof(line), fp) != NULL) {
      atom = pdbline2atom(line);

      for (i=0; i<prot.n_res; i++) {
         /* allow raw pdb comparison */
         if (atom.chainID == ' ') atom.chainID = '_';

         if (prot.res[i].iCode == atom.iCode && \
             prot.res[i].chainID == atom.chainID && \
             prot.res[i].resSeq == atom.resSeq && \
             !strncmp(prot.res[i].resName, atom.resName, 3)) {
            prot.res[i].do_rot = 1;
         }
      }

      memset(line, 0, sizeof(line));
   }

   return 0;
}

int prune_hdirected(PROT prot)
{  float d_low = 2.5;
   float d_high = 3.5;
   float vdw_limit = 20.0; /* kCal/mol */
   char  HaveON;
   int ires, iconf, iatom, jres, jconf, jatom;
   float dd, dd_low, dd_high, d_relax;
   char need_relax;
   float vdw_pair;
   int n_deleted=0;
   int i_nconf, j_nconf;

   dd_low = d_low*d_low;
   dd_high = d_high*d_high;

   for (ires=0; ires<prot.n_res; ires++) {
      prot.res[ires].i_res_prot = ires;
   }

   for (ires=0; ires<prot.n_res; ires++) {
      /* Do I have O or N?*/
      HaveON = 0;
      if (prot.res[ires].n_conf >= 2) {
         for (iatom = 0; iatom <prot.res[ires].conf[1].n_atom; iatom++) {
            if (prot.res[ires].conf[1].atom[iatom].on) {
               if (prot.res[ires].conf[1].atom[iatom].name[1] == 'O' || prot.res[ires].conf[1].atom[iatom].name[1] == 'N') {
                  HaveON = 1;
                  break;
               }
            }
         }
      }

      if (!HaveON) continue;  /* nex residue when no O or N atomes*/

      i_nconf = prot.res[ires].n_conf;
      for (iconf=1; iconf<i_nconf; iconf++) {
         if (prot.res[ires].conf[iconf].history[2] == 'H') break;

         for (jres=0; jres<prot.res[ires].n_ngh; jres++) {
            if (prot.res[ires].ngh[jres]->i_res_prot <= ires) continue;
            j_nconf = prot.res[ires].ngh[jres]->n_conf;
            
	    for (jconf=1; jconf<j_nconf; jconf++) {
               if (prot.res[ires].ngh[jres]->conf[jconf].history[2] == 'H') break;

               /* now, got the pair, check ON - ON distance */
               need_relax = 0; /* by default, no relaxation */
               for (iatom = 0; iatom <prot.res[ires].conf[iconf].n_atom; iatom++) {
                  if (need_relax) break;
                  if (!prot.res[ires].conf[iconf].atom[iatom].on) continue;
                  if (prot.res[ires].conf[iconf].atom[iatom].name[1] != 'N' && prot.res[ires].conf[iconf].atom[iatom].name[1] != 'O') continue;
                  for (jatom = 0; jatom <prot.res[ires].ngh[jres]->conf[jconf].n_atom; jatom++) {
                     if (!prot.res[ires].ngh[jres]->conf[jconf].atom[jatom].on) continue;
                     if (prot.res[ires].ngh[jres]->conf[jconf].atom[jatom].name[1] != 'N' && prot.res[ires].ngh[jres]->conf[jconf].atom[jatom].name[1] != 'O') continue;

                     dd = ddvv(prot.res[ires].conf[iconf].atom[iatom].xyz,prot.res[ires].ngh[jres]->conf[jconf].atom[jatom].xyz);

                     if (dd>dd_low && dd<dd_high) {
                        vdw_pair = vdw_conf_hv(ires, iconf, prot.res[ires].ngh[jres]->i_res_prot, jconf, prot);
                        if (vdw_pair<vdw_limit) {
                           need_relax = 1;
                           break;
                        }
                     }
                  }
               }
               if (need_relax) {
		  d_relax = relax_this_pair(prot, &prot.res[ires], iconf, prot.res[ires].ngh[jres], jconf);
                  /* DEBUG
                  printf("%9.6f %9.6f\n",sqrt(dd), d_relax);
                  fflush(stdout);
                  */
               }
            }
         }
      }

      /* for this polar/ionizable residue, delete non derived conformers
      for (iconf=2; iconf<prot.res[ires].n_conf; iconf++) {
         if (prot.res[ires].conf[iconf].history[2] == 'H') break;
         del_conf(&prot.res[ires], iconf);
         iconf--;
      }
      */

      for (iconf=2; iconf<prot.res[ires].n_conf; iconf++) {
         if (prot.res[ires].conf[iconf].history[2] == 'H') break;
      }

      if (prot.res[ires].n_conf - iconf >= 2) {

         /* delete conformers with high offset to meet conformer number limit */
         qsort((prot.res[ires].conf+iconf), prot.res[ires].n_conf-iconf, sizeof(CONF), cmp_Eself);
         /* DEBUG
         for (jconf=iconf; jconf<prot.res[ires].n_conf; jconf++) {
            printf("%6d, %8.4f %10s\n", jconf, prot.res[ires].conf[jconf].E_self, prot.res[ires].conf[jconf].history);
         }
      */
      }
      for (jconf=prot.res[ires].n_conf-1; jconf>=env.hdirlimt+iconf; jconf--) {
         del_conf(&prot.res[ires], jconf);
         n_deleted++;
      }
   }
   printf("   %d conformers are deleted to fit the limit (%d) of derived conformer number\n", n_deleted, env.hdirlimt);

   return 0;
}

float relax_this_pair(PROT prot, RES *res_a, int conf_a, RES *res_b, int conf_b)
{  int i, n_a, n_b, kc;
   int C;
   char C_str[5];
   ROTAMER rule;
   char sbuff[MAXCHAR_LINE];
   int iconf, jconf, iatom, jatom;
   float dd3, dd, dd3_min, dd_min;
   int i_min, j_min;
   int ic;

   /* remember the number of the current conformers */
   n_a = res_a->n_conf;
   n_b = res_b->n_conf;

   /* make a copy of passed in conformer */
   ins_conf(res_a, res_a->n_conf, res_a->conf[conf_a].n_atom);
   cpy_conf(&res_a->conf[res_a->n_conf-1], &res_a->conf[conf_a]);
   ins_conf(res_b, res_b->n_conf, res_b->conf[conf_b].n_atom);
   cpy_conf(&res_b->conf[res_b->n_conf-1], &res_b->conf[conf_b]);


   /* a side rotamer relaxation */
   C = 0;
   while (1) {
      sprintf(C_str, "%d", C);
      if (param_get("ROTAMER", res_a->resName, C_str, &rule)) break;
      sprintf(sbuff, " %s  ", rule.affected);
      sprintf(rule.affected, sbuff);
      kc = res_a->n_conf;
      for (i=n_a; i<kc; i++) {
         if (swing_conf_rot_rule(res_a, i, rule, 4.0)) {
            printf("      Error in H bond relaxation.\n");
         }
      }
      C++;
   }

   /* b side rotamer relaxation */
   C = 0;
   while (1) {
      sprintf(C_str, "%d", C);
      if (param_get("ROTAMER", res_b->resName, C_str, &rule)) break;
      sprintf(sbuff, " %s  ", rule.affected);
      sprintf(rule.affected, sbuff);
      kc = res_b->n_conf;
      for (i=n_b; i<kc; i++) {
         if (swing_conf_rot_rule(res_b, i, rule, 4.0)) {
            printf("      Error in H bond relaxation.\n");
         }
      }
      C++;
   }

   /* refine rotamer distance to 3.0 */
   i_min = n_a;
   j_min = n_b;
   dd3_min = 1000000.0;
   dd_min =  1000000.0;
   for (iconf=n_a; iconf<res_a->n_conf; iconf++) {
      for (jconf=n_b; jconf<res_b->n_conf; jconf++) {
         for (iatom=0; iatom<res_a->conf[iconf].n_atom; iatom++) {
            if (!res_a->conf[iconf].atom[iatom].on) continue;
            if (res_a->conf[iconf].atom[iatom].name[1] != 'N' && res_a->conf[iconf].atom[iatom].name[1] != 'O') continue;
            for (jatom=0; jatom<res_b->conf[jconf].n_atom; jatom++) {
               if (!res_b->conf[jconf].atom[jatom].on) continue;
               if (res_b->conf[jconf].atom[jatom].name[1] != 'N' && res_b->conf[jconf].atom[jatom].name[1] != 'O') continue;
               dd = ddvv(res_a->conf[iconf].atom[iatom].xyz, res_b->conf[jconf].atom[jatom].xyz);
               dd3 = fabs(dd-9.0);
               if (dd3<dd3_min) {
                  dd3_min = dd3;
                  dd_min  = dd;
                  strncpy(res_a->conf[i_min].history+2, "___", 3);
                  strncpy(res_b->conf[j_min].history+2, "___", 3);
                  strncpy(res_a->conf[iconf].history+2, "MIN", 3);
                  strncpy(res_b->conf[jconf].history+2, "MIN", 3);
                  i_min = iconf;
                  j_min = jconf;
                  /* DEBUG
                  printf ("%.6f, %d(%.3f, %.3f, %.3f) %d(%.3f, %.3f, %.3f)\n", sqrt(dd_min),
                     iconf,
                     res_a->conf[iconf].atom[iatom].xyz.x,
                     res_a->conf[iconf].atom[iatom].xyz.y,
                     res_a->conf[iconf].atom[iatom].xyz.z,
                     jconf,
                     res_b->conf[jconf].atom[jatom].xyz.x,
                     res_b->conf[jconf].atom[jatom].xyz.y,
                     res_b->conf[jconf].atom[jatom].xyz.z);
                  */
               }
            }
         }
      }
   }


   for (iconf=n_a; iconf<res_a->n_conf; iconf++) {
      if (strncmp(res_a->conf[iconf].history+2,"MIN",3)) {
         del_conf(res_a, iconf);
         iconf--;
      }
   }
   for (iconf=n_b; iconf<res_b->n_conf; iconf++) {
      if (strncmp(res_b->conf[iconf].history+2,"MIN",3)) {
         del_conf(res_b, iconf);
         iconf--;
      }
   }


   /* SECOND ROUND, repeat with swing angle 1 degree */
   C = 0;
   while (1) {
      sprintf(C_str, "%d", C);
      if (param_get("ROTAMER", res_a->resName, C_str, &rule)) break;
      sprintf(sbuff, " %s  ", rule.affected);
      sprintf(rule.affected, sbuff);
      kc = res_a->n_conf;
      for (i=n_a; i<kc; i++) {
         if (swing_conf_rot_rule(res_a, i, rule, 2.0)) {
            printf("      Error in H bond relaxation.\n");
         }
      }
      C++;
   }
   C = 0;
   while (1) {
      sprintf(C_str, "%d", C);
      if (param_get("ROTAMER", res_b->resName, C_str, &rule)) break;
      sprintf(sbuff, " %s  ", rule.affected);
      sprintf(rule.affected, sbuff);
      kc = res_b->n_conf;
      for (i=n_b; i<kc; i++) {
         if (swing_conf_rot_rule(res_b, i, rule, 2.0)) {
            printf("      Error in H bond relaxation.\n");
         }
      }
      C++;
   }
   i_min = n_a;
   j_min = n_b;
   dd3_min = 1000000.0;
   dd_min =  1000000.0;
   for (iconf=n_a; iconf<res_a->n_conf; iconf++) {
      for (jconf=n_b; jconf<res_b->n_conf; jconf++) {
         for (iatom=0; iatom<res_a->conf[iconf].n_atom; iatom++) {
            if (!res_a->conf[iconf].atom[iatom].on) continue;
            if (res_a->conf[iconf].atom[iatom].name[1] != 'N' && res_a->conf[iconf].atom[iatom].name[1] != 'O') continue;
            for (jatom=0; jatom<res_b->conf[jconf].n_atom; jatom++) {
               if (!res_b->conf[jconf].atom[jatom].on) continue;
               if (res_b->conf[jconf].atom[jatom].name[1] != 'N' && res_b->conf[jconf].atom[jatom].name[1] != 'O') continue;
               dd = ddvv(res_a->conf[iconf].atom[iatom].xyz, res_b->conf[jconf].atom[jatom].xyz);
               dd3 = fabs(dd-9.0);
               if (dd3<dd3_min) {
                  dd3_min = dd3;
                  dd_min  = dd;
                  strncpy(res_a->conf[i_min].history+2, "___", 3);
                  strncpy(res_b->conf[j_min].history+2, "___", 3);
                  strncpy(res_a->conf[iconf].history+2, "MIN", 3);
                  strncpy(res_b->conf[jconf].history+2, "MIN", 3);
                  i_min = iconf;
                  j_min = jconf;
               }
            }
         }
      }
   }
   for (iconf=n_a; iconf<res_a->n_conf; iconf++) {
      if (strncmp(res_a->conf[iconf].history+2,"MIN",3)) {
         del_conf(res_a, iconf);
         iconf--;
      }
   }
   for (iconf=n_b; iconf<res_b->n_conf; iconf++) {
      if (strncmp(res_b->conf[iconf].history+2,"MIN",3)) {
         del_conf(res_b, iconf);
         iconf--;
      }
   }


   /* THIRD ROUND, repeat with swing angle 1 degree */
   C = 0;
   while (1) {
      sprintf(C_str, "%d", C);
      if (param_get("ROTAMER", res_a->resName, C_str, &rule)) break;
      sprintf(sbuff, " %s  ", rule.affected);
      sprintf(rule.affected, sbuff);
      kc = res_a->n_conf;
      for (i=n_a; i<kc; i++) {
         if (swing_conf_rot_rule(res_a, i, rule, 1.0)) {
            printf("      Error in H bond relaxation.\n");
         }
      }
      C++;
   }
   C = 0;
   while (1) {
      sprintf(C_str, "%d", C);
      if (param_get("ROTAMER", res_b->resName, C_str, &rule)) break;
      sprintf(sbuff, " %s  ", rule.affected);
      sprintf(rule.affected, sbuff);
      kc = res_b->n_conf;
      for (i=n_b; i<kc; i++) {
         if (swing_conf_rot_rule(res_b, i, rule, 1.0)) {
            printf("      Error in H bond relaxation.\n");
         }
      }
      C++;
   }
   i_min = n_a;
   j_min = n_b;
   dd3_min = 1000000.0;
   dd_min =  1000000.0;
   for (iconf=n_a; iconf<res_a->n_conf; iconf++) {
      for (jconf=n_b; jconf<res_b->n_conf; jconf++) {
         for (iatom=0; iatom<res_a->conf[iconf].n_atom; iatom++) {
            if (!res_a->conf[iconf].atom[iatom].on) continue;
            if (res_a->conf[iconf].atom[iatom].name[1] != 'N' && res_a->conf[iconf].atom[iatom].name[1] != 'O') continue;
            for (jatom=0; jatom<res_b->conf[jconf].n_atom; jatom++) {
               if (!res_b->conf[jconf].atom[jatom].on) continue;
               if (res_b->conf[jconf].atom[jatom].name[1] != 'N' && res_b->conf[jconf].atom[jatom].name[1] != 'O') continue;
               dd = ddvv(res_a->conf[iconf].atom[iatom].xyz, res_b->conf[jconf].atom[jatom].xyz);
               dd3 = fabs(dd-9.0);
               if (dd3<dd3_min) {
                  dd3_min = dd3;
                  dd_min  = dd;
                  strncpy(res_a->conf[i_min].history+2, "___", 3);
                  strncpy(res_b->conf[j_min].history+2, "___", 3);
                  strncpy(res_a->conf[iconf].history+2, "MIN", 3);
                  strncpy(res_b->conf[jconf].history+2, "MIN", 3);
                  i_min = iconf;
                  j_min = jconf;
               }
            }
         }
      }
   }
   for (iconf=n_a; iconf<res_a->n_conf; iconf++) {
      if (strncmp(res_a->conf[iconf].history+2,"MIN",3)) {
         del_conf(res_a, iconf);
         iconf--;
      }
   }
   for (iconf=n_b; iconf<res_b->n_conf; iconf++) {
      if (strncmp(res_b->conf[iconf].history+2,"MIN",3)) {
         del_conf(res_b, iconf);
         iconf--;
      }
   }

   res_a->conf[n_a].history[2] = 'H';
   res_b->conf[n_b].history[2] = 'H';
   res_a->conf[n_a].E_self = dd3_min;
   res_b->conf[n_b].E_self = dd3_min;

   /* delete the identical confomers with higher E_self */
   for (ic=res_a->n_conf-2; ic>0; ic--) {
      if (res_a->conf[ic].history[2] != 'H') continue;
      if(!cmp_conf_hv(res_a->conf[res_a->n_conf-1], res_a->conf[ic], env.hdirdiff)) {
         if (res_a->conf[ic].E_self >= dd3_min) {
            del_conf(res_a, ic);
         }
         else {
            del_conf(res_a, res_a->n_conf-1);
            break;
         }
      }
   }

   for (ic=res_b->n_conf-2; ic>0; ic--) {
      if (res_b->conf[ic].history[2] != 'H') continue;
      if(!cmp_conf_hv(res_b->conf[res_b->n_conf-1], res_b->conf[ic], env.hdirdiff)) {
         if (res_b->conf[ic].E_self >= dd3_min) {
            del_conf(res_b, ic);
         }
         else {
            del_conf(res_b, res_b->n_conf-1);
            break;
         }
      }
   }

   return sqrt(dd_min);
}

int swing_conf_rot_rule(RES *res, int j, ROTAMER rule, float phi)
{  VECTOR v1, v2, v3;
   GEOM op;
   LINE axis;
   int ins;
   int k, l;
   ATOM atom1, atom2;
   char found;

   k=iatom(res->conf[j].confName, rule.atom2);

   if ((k=iatom(res->conf[j].confName, rule.atom2)) == -1) {
      printf("confname=\"%s\"\n",res->conf[j].confName);
      if ((k=iatom(res->conf[0].confName, rule.atom2)) == -1) {
         printf("   FATAL: place_rot_rule(): can't find atom \"%s\" in residue \"%s\"\n",
                 rule.atom2, res->resName);
         return USERERR;
      }
      else atom2 = res->conf[0].atom[k];
   }
   else atom2 = res->conf[j].atom[k];
   v2 = atom2.xyz;

   /* find the bond: atom1 is the first atom of the bond, it is one of the connected
    * atoms of atom2, but not necessary in this conformer or this residue */
   found = 0;
   if ((k=iatom(res->conf[j].confName, rule.atom1)) == -1) {
      l=0;
      while (atom2.connect12[l]!=NULL) {
         if (!strcmp(atom2.connect12[l]->name, rule.atom1)) {
            atom1 = *atom2.connect12[l];
            found = 1;
            break;
         }
         l++;
      }
   }
   else {
      atom1 = res->conf[j].atom[k];
      found = 1;
   }

   if (found)  v1 = atom1.xyz;
   else {
      printf("   FATAL: place_rot_rule(): can't find atom \"%s\" connected to \"%s\" in residue \"%s\"\n",
              rule.atom2, rule.atom1, res->resName);
      return USERERR;
   }

   axis = line_2v(v1, v2);

   /* swing left */
   geom_reset(&op);
   geom_roll(&op, -phi, axis);
   ins = ins_conf(res, res->n_conf, res->conf[j].n_atom);
   if (ins == USERERR) return USERERR;

   cpy_conf(&res->conf[ins], &res->conf[j]);
   strncpy(res->conf[ins].history+2, "HOSw", 4);
   for (k=0; k<res->conf[j].n_atom; k++) {
      if (strstr(rule.affected, res->conf[j].atom[k].name)) {
         v3 = res->conf[j].atom[k].xyz;
         geom_apply(op, &v3);
         res->conf[ins].atom[k].xyz = v3;
      }
   }

   /* swing right */
   geom_reset(&op);
   geom_roll(&op, phi, axis);
   ins = ins_conf(res, res->n_conf, res->conf[j].n_atom);
   if (ins == USERERR) return USERERR;

   cpy_conf(&res->conf[ins], &res->conf[j]);
   strncpy(res->conf[ins].history+2, "HOSw", 4);
   for (k=0; k<res->conf[j].n_atom; k++) {
      if (strstr(rule.affected, res->conf[j].atom[k].name)) {
         v3 = res->conf[j].atom[k].xyz;
         geom_apply(op, &v3);
         res->conf[ins].atom[k].xyz = v3;
      }
   }

   return 0;
}


int prune_pv(PROT prot, float c1, float c2, float c3)
{
    int n = 0; /* number of conformers deleted */
   float cutoff_geo = c1;
   float cutoff_ele = c2;
   float cutoff_vdw = c3;
   int ir, ic, ia, jc;
   int deleting;
   int i_res, j_res, i_conf, j_conf, i_atom, j_atom;
   long    idum;
   
   if (env.test_seed < 0) idum = time(NULL); //allows random numbers to be fixed for testing
   else idum = env.test_seed;
   
   int i;

    for (i=0;i<500;i++) {
         ran2(&idum);
    }

   id_conf(prot); /* is conf ID used here? - y*/
   
   /* calculate conformer netcrg */
   for (ir=0; ir<prot.n_res; ir++) {
      for (ic=1; ic<prot.res[ir].n_conf; ic++) {
         prot.res[ir].conf[ic].netcrg = 0.0;
         for (ia=0; ia<prot.res[ir].conf[ic].n_atom; ia++) {
            if (!prot.res[ir].conf[ic].atom[ia].on) continue;
            prot.res[ir].conf[ic].netcrg += prot.res[ir].conf[ic].atom[ia].crg;
         }
      }
   }
   
   /* set flags for all conformer */
   for (ir = 0; ir< prot.n_res; ir++) {
       for (ic = 1; ic< prot.res[ir].n_conf; ic++) {
           prot.res[ir].conf[ic].on =1;
       }
   }
   
   /* set up neighbor list */
   assign_rad(prot);
   assign_crg(prot);
   assign_vdw_param(prot);
   get_connect12(prot);
   get_vdw0(prot);
   get_vdw1(prot);
   for (ir=0; ir<prot.n_res; ir++) {
      for (ic=1; ic<prot.res[ir].n_conf; ic++) {
          prot.res[ir].conf[ic].E_torsion = torsion_conf(&prot.res[ir].conf[ic]);
          prot.res[ir].conf[ic].E_self = prot.res[ir].conf[ic].E_vdw0 + prot.res[ir].conf[ic].E_vdw1 + prot.res[ir].conf[ic].E_torsion;
      }
   }

   setup_vdw_fast(prot);
   for (i_res = 0; i_res < prot.n_res; i_res++) {
       prot.res[i_res].n_ngh = 0;
       prot.res[i_res].ngh = NULL;
       
       for (j_res = 0; j_res< prot.n_res; j_res++) {
           int found;
           if (i_res == j_res) continue;
           
           /* check if distance within the threshold */
           if (out_of_range(prot.res[i_res].r_min,prot.res[i_res].r_max,prot.res[j_res].r_min,prot.res[j_res].r_max,36.)) continue;  /* 8 away */
           found = 0;
           for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
               for (j_conf=0; j_conf<prot.res[j_res].n_conf; j_conf++) {
                   
                   for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                       if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                       for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
                           if (!prot.res[j_res].conf[j_conf].atom[j_atom].on) continue;
                           
                           if (ddvv(prot.res[i_res].conf[i_conf].atom[i_atom].xyz,prot.res[j_res].conf[j_conf].atom[j_atom].xyz) < 64.) { /* 8 away */
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
           
           prot.res[i_res].n_ngh++;
           prot.res[i_res].ngh = realloc(prot.res[i_res].ngh, prot.res[i_res].n_ngh * sizeof(RES *));
           prot.res[i_res].ngh[prot.res[i_res].n_ngh-1] = &prot.res[j_res];
       }
   }

   
   /* get pairwise vector, ele pairwise + vdw pairwise 
    * We may want to keep the first generation conformers, non-rotamers,
    * jmao Original conformer won't be pruned.
    */
    
   /* adding recursion here (originally outside of the subroutine, which makes self-energy calculations unnecessarily repeated) */
   deleting = 1;
   while (deleting) {
       deleting = 0;
   for (ir=0; ir<prot.n_res; ir++) {
      for (ic=1; ic<prot.res[ir].n_conf-1; ic++) {
         if (!prot.res[ir].conf[ic].on) continue;
         //if (!prot.res[ir].conf[ic].history[2] == 'E') continue; /* do not compare exposed conformers, they have too few interactions to distinguish */
         for (jc=ic+1; jc<prot.res[ir].n_conf; jc++) {
            if (!prot.res[ir].conf[jc].on) continue;
            //if (!prot.res[ir].conf[jc].history[2] == 'E') continue; /* do not compare exposed conformers, they have too few interactions to distinguish */
            
            //if (prot.res[ir].conf[jc].history[2] == 'O') continue;
            if (strcmp(prot.res[ir].conf[ic].confName, prot.res[ir].conf[jc].confName)) continue;
            if (fabs(prot.res[ir].conf[ic].E_self-prot.res[ir].conf[jc].E_self) > cutoff_vdw) continue;
            if (over_geo(prot.res[ir].conf[ic], prot.res[ir].conf[jc], cutoff_geo)) continue;
            if (over_ele(prot, ir, ic, jc, cutoff_ele)) continue;
            if (env.rot_mhd_prune == 1) { //execution path that allows for MHD's pruning function
            	if (over_vdw_mhd(prot, ir, ic, jc, cutoff_vdw)) continue;
            }
            else {
            	if (over_vdw(prot, ir, ic, jc, cutoff_vdw)) continue;
            }
            	
            //printf("Self energies of the kept and pruned conformers are %8.3f and %8.3f\n", prot.res[ir].conf[ic].E_self,prot.res[ir].conf[jc].E_self);
			//prot.res[ir].conf[jc].on = 0;
            
            deleting = 1;
            if (prot.res[ir].conf[ic].history[2] == 'O' && prot.res[ir].conf[jc].history[2] == 'O') {
			    if (ran2(&idum) < 0.5) prot.res[ir].conf[ic].on = 0;
				else prot.res[ir].conf[jc].on = 0;
                break;
			}
			else if (prot.res[ir].conf[ic].history[2] != 'O' && prot.res[ir].conf[jc].history[2] != 'O') {
			    if (ran2(&idum) < 0.5) prot.res[ir].conf[ic].on = 0;
				else prot.res[ir].conf[jc].on = 0;
                break;
			}
			else if (prot.res[ir].conf[ic].history[2] != 'O') {
                prot.res[ir].conf[ic].on = 0;
                break;
			}
            else {
                prot.res[ir].conf[jc].on = 0;
                break;
			}
         }
      }
   }
   }
   
   /* clean up neighbor list */
   for (i_res = 0; i_res < prot.n_res; i_res++) {
       prot.res[i_res].n_ngh = 0;
       free(prot.res[i_res].ngh);
   }
   
   /* delete conformers */
   for (ir=0; ir<prot.n_res; ir++) {
      for (ic=prot.res[ir].n_conf-1; ic>1; ic--) {
         if (!prot.res[ir].conf[ic].on) {
            del_conf(&prot.res[ir], ic);
            n++;
         }
      }
   }
   
   return n;
}

int over_geo(CONF conf1, CONF conf2, float cutoff)
{  int ia, ja;
   float dd_max = cutoff*cutoff;
   float dd;

   for (ia=0; ia<conf1.n_atom; ia++) {
      if (!conf1.atom[ia].on) continue;
      if (!conf2.atom[ia].on) continue;
      if ((ja=iatom(conf2.confName, conf1.atom[ia].name))<0) continue;
      if (!conf2.atom[ja].on) continue;
      dd = ddvv(conf1.atom[ia].xyz, conf2.atom[ja].xyz);
      if (dd>dd_max) return 1;
   }

   return 0;
}

int over_ele(PROT prot, int i_res, int i_conf, int j_conf, float cutoff)
{
    int k_res, k_conf, i_ngh;
    float Ei, Ej;
    
    for (i_ngh = 0; i_ngh < prot.res[i_res].n_ngh; i_ngh++) {
        k_res = prot.res[i_res].ngh[i_ngh]->i_res_prot;
        
        for (k_conf=0; k_conf < prot.res[k_res].n_conf; k_conf++) {
            if (!prot.res[k_res].conf[k_conf].on) continue;
            
            Ei = Ecoulomb_conf2conf(prot, i_res, i_conf, k_res, k_conf, env.epsilon_prot);
            Ej = Ecoulomb_conf2conf(prot, i_res, j_conf, k_res, k_conf, env.epsilon_prot);
            if (fabs(Ei-Ej)>cutoff) return 1;
        }
    }

   return 0;
}

int over_vdw(PROT prot, int i_res, int i_conf, int j_conf, float cutoff)
{
    int k_res, k_conf, i_ngh;
    float Ei, Ej;
    
    for (i_ngh = 0; i_ngh < prot.res[i_res].n_ngh; i_ngh++) {
        k_res = prot.res[i_res].ngh[i_ngh]->i_res_prot;
        
        for (k_conf=0; k_conf < prot.res[k_res].n_conf; k_conf++) {
            if (!prot.res[k_res].conf[k_conf].on) continue;
            
            Ei = Evdw_conf2conf(prot, i_res, i_conf, k_res, k_conf);
            Ej = Evdw_conf2conf(prot, i_res, j_conf, k_res, k_conf);
            if (!(Ei>10.0 && Ej>10.0) && fabs(Ei-Ej)>cutoff) return 1;
        }
    }
    
    return 0;
}

int label_exposed(PROT prot)
{
  	int kr, kc, i;
	float max_sas;

   //sas_native(prot);
   sas_ionizable(prot, env.radius_probe); /* changed to using exposure of terminal atoms of ionizable residues -Yifan */
   
   for (kr=0; kr<prot.n_res; kr++) {
       /*
       for (kc=1; kc<prot.res[kr].n_conf; kc++) {
           av_sas = 0.0;
           n_atom = 0;
           for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
               if (!prot.res[kr].conf[kc].atom[ka].on) continue;
               n_atom++;
               //av_sas += prot.res[kr].conf[kc].atom[ka].sas;
               av_sas += prot.res[kr].conf[kc].atom[ka].sas * 4 * 3.1415926 * prot.res[kr].conf[kc].atom[ka].rad * prot.res[kr].conf[kc].atom[ka].rad;
           }
           av_sas /= n_atom;
           
           prot.res[kr].conf[kc].E_rxn =av_sas;
       }
       */
       
       /* find the conformer with the maximum ASA */
       if (prot.res[kr].n_conf>1) {
           i=1;
           max_sas = prot.res[kr].conf[1].sas_fraction;
       }
       else i =0;
       for (kc=2; kc<prot.res[kr].n_conf; kc++) {
           if (prot.res[kr].conf[kc].sas_fraction > max_sas) {
               max_sas = prot.res[kr].conf[kc].sas_fraction;
               i=kc;
           }
           //printf("  %3s %c%04d%c%03d %8.3f %8.3f\n", prot.res[kr].resName,prot.res[kr].chainID,prot.res[kr].resSeq,prot.res[kr].iCode,kc,max_sas,prot.res[kr].conf[kc].sas_fraction);
       }
      
      /* 5% threshold is for the percentile exposed surface */
      //if (max_sas > 0.05 && i && prot.res[kr].conf[i].history[2] != 'O') {
          /* Here we found the non-native exposed side chain confomer
          ** We will relax this conformer to maximize the exposure.
          */
          /* printf("   %s%4d Before = %.3f", prot.res[kr].resName, prot.res[kr].resSeq, prot.res[kr].conf[i].sas_fraction);
          printf(" After = %.3f\n", max_confSAS(prot, kr, i)); */
          
          /* this function does not include self energy for optimization and uses a wrong ASA calculation, turned off by Yifan */
          //max_confSAS(prot, kr, i);
      //}

      /* mark conformers with big ASA -Yifan */
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
          if (prot.res[kr].conf[kc].history[2] == 'O') continue; /* skip original */
          if (prot.res[kr].conf[kc].history[2] == 'I') continue; /* skip original */
          if (prot.res[kr].conf[kc].history[2] == 'Y') continue; /* skip original */
          /* 20% threshold is for the percentile exposed surface
          - From Jin's paper 20% guarantees a <2pK units of desolvation energy -Yifan*/
          if (prot.res[kr].conf[kc].sas_fraction > 0.2) {
              prot.res[kr].conf[kc].history[2] = 'E';
              //printf("E  %3s %c%04d%c%03d %8.3f %8.3f\n", prot.res[kr].resName,prot.res[kr].chainID,prot.res[kr].resSeq,prot.res[kr].iCode,kc,max_sas,prot.res[kr].conf[kc].sas_fraction);
          }
      }
   }
   return 0;
}

float max_confSAS(PROT prot, int ir, int ic)
{  int nconf,i, kc, ka;
   int C;
   char C_str[5];
   ROTAMER rule;
   char sbuff[MAXCHAR_LINE];
   float av_sas, max_sas;
   int n_atom;

   for (i=0; i<prot.n_res; i++) {
      if (prot.res[i].n_conf > 1) prot.res[i].conf[1].on = 1;
      prot.res[i].conf[0].on = 1;
      for (kc=2; kc<prot.res[i].n_conf; kc++) prot.res[i].conf[kc].on = 0;
   }


   /* remember the number of the current conformers */
   nconf = prot.res[ir].n_conf;

   /* make a copy of passed in conformer */
   ins_conf(&prot.res[ir], nconf, prot.res[ir].conf[ic].n_atom);
   cpy_conf(&prot.res[ir].conf[nconf], &prot.res[ir].conf[ic]);

   /* conformer relaxation */
   C = 0;
   /* printf("      n_conf = %4d", prot.res[ir].n_conf);*/
   while (1) {
      sprintf(C_str, "%d", C);
      if (param_get("ROTAMER", prot.res[ir].resName, C_str, &rule)) break;
      sprintf(sbuff, " %s  ", rule.affected);
      sprintf(rule.affected, sbuff);
      kc = prot.res[ir].n_conf;
      for (i=nconf; i<kc; i++) {
         if (swing_conf_rot_rule(&prot.res[ir], i, rule, 5.0)) {
            printf("      Error in H bond relaxation.\n");
         }
      }
      C++;
   }

   /* compute SAS of conformers in this residue */
   surfw_res(prot, ir, env.radius_probe);
   for (kc = nconf; kc < prot.res[ir].n_conf; kc++) {
   //for (kc = 1; kc < prot.res[ir].n_conf; kc++) {
      av_sas = 0.0;
      n_atom = 0;
      for (ka=0; ka<prot.res[ir].conf[kc].n_atom; ka++) {
         if (!prot.res[ir].conf[kc].atom[ka].on) continue;
         n_atom++;
         av_sas += prot.res[ir].conf[kc].atom[ka].sas;
      }
      av_sas /= n_atom;
      prot.res[ir].conf[kc].E_rxn =av_sas;
   }

   i=ic;
   max_sas = prot.res[ir].conf[ic].E_rxn;
   for (kc=prot.res[ir].n_conf-1; kc >= nconf; kc--) {
      if (prot.res[ir].conf[kc].E_rxn > max_sas) {
         max_sas = prot.res[ir].conf[kc].E_rxn;
         i=kc;
      }
   }

   prot.res[ir].conf[i].history[2] = 'E';
   prot.res[ir].conf[ic].history[2] = 'E';

   for (kc=prot.res[ir].n_conf-1; kc >= nconf; kc--) {
      if (prot.res[ir].conf[kc].history[2] != 'E') del_conf(&prot.res[ir], kc);
   }
   //if (prot.res[ir].conf[ic].history[2] != 'E') del_conf(&prot.res[ir], ic);

   return max_sas;
}

int rebuild_sc(PROT prot) {
    int i_res;
    for (i_res=0; i_res<prot.n_res; ++i_res) {
        char copy_atoms[MAXCHAR_LINE];
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
            
            prot.res[i_res].conf[k_conf].history[2] = 'B'; /* Label */

            for (i_atom=0; i_atom<prot.res[i_res].conf[k_conf].n_atom; i_atom++) {
                if (!strstr(copy_atoms, prot.res[i_res].conf[k_conf].atom[i_atom].name)) {
                    /* if atom name does not match those defined in parameter, then do not keep */
                    prot.res[i_res].conf[k_conf].atom[i_atom].on = 0;
                }
            }
            
            /* rebuild */
            while(place_missing(prot,1) > 0);
            rm_dupconf(prot, 0.005);
            //printf("resi %4d, n_conf %4d\n",i_res,prot.res[i_res].n_conf);
        }
    }
        get_connect12(prot);
        assign_vdw_param(prot);
        get_vdw0_no_sas(prot);
        get_vdw1(prot);
        int kc;
        for (i_res=0; i_res<prot.n_res; ++i_res) {
            for (kc=1; kc<prot.res[i_res].n_conf; kc++) {
                prot.res[i_res].conf[kc].E_torsion = torsion_conf(&prot.res[i_res].conf[kc]);
                prot.res[i_res].conf[kc].E_self = prot.res[i_res].conf[kc].E_vdw0 + prot.res[i_res].conf[kc].E_vdw1 + prot.res[i_res].conf[kc].E_torsion;
            }
        }
        prune_by_vdw(prot, env.vdw_cutoff);
    
    return 0;
}

void del_non_common_h(PROT prot)
{
    int i_res, i_conf, i_atom;
    for (i_res=0; i_res<prot.n_res; ++i_res) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; ++i_conf) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; ++i_atom) {
                int common_h;
                ATOM * atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                if (prot.res[i_res].conf[i_conf].atom[i_atom].name[1] != 'H') continue;
                
                /* check if COMMON_H parameter has been saved */
                if (param_get("COMMON_H", prot.res[i_res].resName, atom_p->name, &common_h)) {
                    STRINGS conf_types;
                    int i_conf_type;
                    if (param_get("CONFLIST", prot.res[i_res].resName, "", &conf_types)) {
                        continue;
                    }
                    if (conf_types.n <=2) {
                        common_h = 1;
                        param_sav("COMMON_H", prot.res[i_res].resName, atom_p->name, &common_h, sizeof(int));
                        continue;
                    }
                    
                    common_h = 1;
                    for (i_conf_type=1; i_conf_type<conf_types.n; i_conf_type++) {
                        int buff;
                        if (param_get("IATOM", conf_types.strings[i_conf_type], atom_p->name, &buff)) {
                            common_h = 0;
                            break;
                        }
                    }
                    param_sav("COMMON_H", prot.res[i_res].resName, atom_p->name, &common_h, sizeof(int));
                }
                
                if (!common_h) {
                    /* delete if not common */
                    prot.res[i_res].conf[i_conf].atom[i_atom].on = 0;
                }
            }
        }
    }
}

int rand_conf_prune(PROT prot)
{
    int i_res, i_conf;
    
    if (env.test_seed < 0) idum = time(NULL); //allows random numbers to be fixed for testing
    else idum = env.test_seed;
    
    for (i_conf=0;i_conf<2000;i_conf++) ran2(&idum);
    
    for (i_res=0; i_res<prot.n_res; ++i_res) {
        int n_orig = 0;
        int n_hv_conf_limit = env.n_hv_conf_limit;
        
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; ++i_conf) {
            if (prot.res[i_res].conf[i_conf].history[2] == 'O') n_orig++;
            else if (prot.res[i_res].conf[i_conf].history[2] == 'X') n_orig++;
            else if (prot.res[i_res].conf[i_conf].history[2] == 'Y') n_orig++;
        }

        if (n_hv_conf_limit < n_orig) n_hv_conf_limit = n_orig;
        
        while (prot.res[i_res].n_conf-1 > n_hv_conf_limit) {
            i_conf = 1 + ran2(&idum) * (float) (prot.res[i_res].n_conf-1);
            if (prot.res[i_res].conf[i_conf].history[2] == 'O') continue;
            if (prot.res[i_res].conf[i_conf].history[2] == 'X') continue;
            if (prot.res[i_res].conf[i_conf].history[2] == 'Y') continue;
            
            del_conf(&prot.res[i_res], i_conf);
        }
    }
    return 0;
}

/* the following two functions are modified versions of the regular functions for the purpose
of facilitating MHD's pruning function */

int over_vdw_mhd(PROT prot, int i_res, int  i_conf, int j_conf, float cutoff)
{
    int k_res, k_conf, i_ngh, i_res_conf, k_res_conf;
    float Ei, Ej;

    
    for (i_ngh = 0; i_ngh < prot.res[i_res].n_ngh; i_ngh++) {
        k_res = prot.res[i_res].ngh[i_ngh]->i_res_prot;

        
	float min_vdw = 99999.;
	for (i_res_conf = 1; i_res_conf < prot.res[i_res].n_conf; i_res_conf++) {
	for (k_res_conf = 1; k_res_conf < prot.res[k_res].n_conf; k_res_conf++) {
	float pw_vdw = vdw_conf(i_res, i_res_conf, k_res, k_res_conf, prot);
	if (pw_vdw < min_vdw) {min_vdw = pw_vdw;}
			}	
		}	
	if (min_vdw < 0.) { min_vdw = 0.; }


        for (k_conf=0; k_conf < prot.res[k_res].n_conf; k_conf++) {
            if (!prot.res[k_res].conf[k_conf].on) continue;

            
            Ei = vdw_conf(i_res, i_conf, k_res, k_conf, prot);
            Ej = vdw_conf(i_res, j_conf, k_res, k_conf, prot);
            if (!(Ei>(min_vdw+10.0) && Ej>(min_vdw+10.0)) && fabs(Ei-Ej)>cutoff) return 1;
            //printf("Minghui find pruned conformers are %3s and %8d and %8d\n",prot.res[k_res].resName, prot.res[k_res].n_conf, k_conf);
        }
    }

    
    return 0;
}

int rand_conf_prune_mhd(PROT prot)
{
    int i_res, i_conf;
    
    if (env.test_seed < 0) idum = time(NULL); //allows random numbers to be fixed for testing
    else idum = env.test_seed;
    
    for (i_conf=0;i_conf<2000;i_conf++) ran2(&idum);
    
    for (i_res=0; i_res<prot.n_res; ++i_res) {
        int n_orig = 0;
        int n_hv_conf_limit = env.n_hv_conf_limit;
        
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; ++i_conf) {
            if (prot.res[i_res].conf[i_conf].history[2] == 'O') n_orig++;
            else if (prot.res[i_res].conf[i_conf].history[2] == 'X') n_orig++;
            else if (prot.res[i_res].conf[i_conf].history[2] == 'E') n_orig++;
        }

        if (n_hv_conf_limit < n_orig) n_hv_conf_limit = n_orig;
        
        while (prot.res[i_res].n_conf-1 > n_hv_conf_limit) {
            i_conf = 1 + ran2(&idum) * (float) (prot.res[i_res].n_conf-1);
            if (prot.res[i_res].conf[i_conf].history[2] == 'O') continue;
            if (prot.res[i_res].conf[i_conf].history[2] == 'X') continue;
            if (prot.res[i_res].conf[i_conf].history[2] == 'E') continue;
            
            del_conf(&prot.res[i_res], i_conf);
        }
    }
    return 0;
}

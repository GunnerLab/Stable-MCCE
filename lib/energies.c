#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "mcce.h"

/* precaculate energy terms:
 *       self vdw
 *       bkb  vdw
 *       torsion
 *       epol
 *       solvation
 *       pairwise vdw
 *       pairwise electristatic
 *
 * The results will be placed in one designated directory
 */


/* this is a global variable for bookkeeping the dielectric boundary
 * and the atoms linked to each each entry.
 */
typedef struct {
   int header;
/* jmao: IDL of gfortran needs a 4-byte int to store the length of a record */
   float x;
   float y;
   float z;
   float rad;
   float crg;
   int trailer;
} UNF;

#define NUNIQS 1000000

typedef struct {
   UNF  unf[NUNIQS];      /* this number limits the number of atoms */
   char head[NUNIQS][31]; /* head of the atom line */
   int n;
} ELE_BOUND;

ELE_BOUND ele_bound;
int del_runs;           /* depth of focusing */
int n_retry;
char del_folder[300];
char cur_folder[300];

int conf_energies(int kr, int kc, PROT prot);
int conf_rxn(int kr, int kc, PROT prot);
int define_boundary(PROT prot);
int find_atom(VECTOR v, int counter);
int delphi_depth();
void write_fort15();

int refresh_prot(PROT prot);
int add_dummies(PROT prot);
time_t time_start, nowA, nowB;

int energies()
{  FILE *fp;
   FILE *progress_fp;
   PROT prot;
   int i, j, k, counter, i_chosen;
   int n_conf;
   char sbuff[MAXCHAR_LINE];

   /* Load step 2 output pdb file */
   time_start = time(NULL);
   printf("   Load step 2 output file %s...\n", STEP2_OUT);
   if (!(fp=fopen(STEP2_OUT, "r"))) {
      printf("   FATAL: energies(): \"No step 2 output \"%s\".\n", STEP2_OUT);
      return USERERR;
   }
   prot = load_pdb(fp);
   if (prot.n_res == 0) {
      printf("   There are errors in pdb file, quiting ...\n");
      return USERERR;
   }
   printf("   Done\n\n"); fflush(stdout);
   fclose(fp);

   assign_vdw_param(prot);

   /* get unique name for each conformer */
   id_conf(prot);

   /* report net charge */
   if (env.reassign) {
       assign_rad(prot);
       assign_crg(prot);
   }
   counter = 0;
   printf("   Reporting non integer conformer charge ...\n");
   for (i=0; i<prot.n_res; i++) {
      for (j=0; j<prot.res[i].n_conf; j++) {
         prot.res[i].conf[j].netcrg = 0.0;
         for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
            if (!prot.res[i].conf[j].atom[k].on) continue;
            prot.res[i].conf[j].netcrg += prot.res[i].conf[j].atom[k].crg;
         }
         if (fabs(prot.res[i].conf[j].netcrg - rintf(prot.res[i].conf[j].netcrg)) > 0.001) {
            printf("   WARNING: Conformer %s has non integer charge %6.3f\n", prot.res[i].conf[j].uniqID,
                                                                        prot.res[i].conf[j].netcrg);
            counter ++;
            fflush(stdout);
         }
      }
   }
   if (counter) printf("   Make sure you intended to have these non integer charges.\n");
   else printf("   Passed, no non-integer charges\n");
   printf("   Done\n\n");
   fflush(stdout);

   printf("   Creating connectivity table...\n"); fflush(stdout);
   get_connect12(prot);
   printf("   Done\n\n"); fflush(stdout);

   printf("   Preparing DelPhi runs at epsilon %.2f ...\n", env.epsilon_prot); fflush(stdout);

   progress_fp = fopen(env.progress_log, "a");

   /* Propare run delphi in a different directory */
   if (env.pbe_folder[strlen(env.pbe_folder)-1] == '/') env.pbe_folder[strlen(env.pbe_folder)-1] = '\0';
   sprintf(del_folder, "%s/delphi_XXXXXX", env.pbe_folder);
   if (!mkdtemp(del_folder)) {
      printf("   Fatal: failed creating delphi tempoary folder \"%s\"\n", del_folder);
      return USERERR;
   }

   /* This piece of code will cause trouble in memprof, though it is 100% legal. To
    * do memory profiling, manully set cur_folder.
    */
   getcwd(cur_folder, sizeof(cur_folder));
   /*
   iostream = popen("pwd", "r");
   fgets(cur_folder, sizeof(cur_folder), iostream);
   if (cur_folder[strlen(cur_folder)-1] == '\n') cur_folder[strlen(cur_folder)-1] = '\0';
   pclose(iostream);
   */
   chdir(del_folder);


   /* define dielectric boundary */
   define_boundary(prot);
   if (ele_bound.n >= 100000) { /* atoms exceeding delphi limit - 100,000 lines */
      printf("   FATAL: %d unique atoms, exceeding 100,000 limit of delphi\n", ele_bound.n);
      return USERERR;
   }

   write_fort15();

   del_runs = delphi_depth();
   printf("      %d focusing runs required for this protein.\n", del_runs);

   /* DEBUG
   int k;
   for (i=0; i<prot.n_res; i++) {
      for (j=0; j<prot.res[i].n_conf; j++) {
         for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
            if (!prot.res[i].conf[j].atom[k].on) continue;
            printf ("%5d%8.3f%8.3f%8.3f%8.3f%8.3f\n", prot.res[i].conf[j].atom[k].serial,
                                                      ele_bound.unf[prot.res[i].conf[j].atom[k].serial].x,
                                                      ele_bound.unf[prot.res[i].conf[j].atom[k].serial].y,
                                                      ele_bound.unf[prot.res[i].conf[j].atom[k].serial].z,
                                                      ele_bound.unf[prot.res[i].conf[j].atom[k].serial].rad,
                                                      ele_bound.unf[prot.res[i].conf[j].atom[k].serial].crg);
         }
      }
   }
   */

   /* Compute self vdw, bkb vdw. torsion, epol, pairwise vdw and pairwise electrostatic energy
    * of each conformer and save to a file.
    */
   counter=0;
   for (i=0; i<prot.n_res; i++)
      counter += prot.res[i].n_conf-1;
   n_conf = counter;

   /* Set up index to head list */
   counter = 0;
   for (i=0; i<prot.n_res; i++) {
      for (j=0; j<prot.res[i].n_conf; j++) {
         prot.res[i].conf[j].iConf = counter;
         counter ++;
      }
   }
   printf("   Done\n\n");

   /* Correct start and end limit */
   if (env.pbe_start < 1) env.pbe_start = 1;
   if (env.pbe_end > n_conf) env.pbe_end = n_conf;

   printf("   Computing pairwise from conformer %d to %d of %d total conformers\n      see %s for progress...\n", env.pbe_start,env.pbe_end, n_conf, env.progress_log); fflush(stdout);
   counter = 0;
   for (i=0; i<prot.n_res; i++) {
      for (j=1; j<prot.res[i].n_conf; j++) {
         if (counter < env.pbe_start-1 || counter > env.pbe_end-1) {
            counter++;
            continue;
         }
         /* clean last run */
         for (k=1; k<=del_runs; k++) {
            sprintf(sbuff, "rxn%02d.log", k);
            remove(sbuff);
            sprintf(sbuff, "run%02d.frc", k);
            remove(sbuff);
         }
         remove("ARCDAT");

         counter++;
         nowA = time(NULL);
         fprintf(progress_fp, "   Doing pairwise %5d of %5d conformers.", counter, n_conf); fflush(progress_fp);

         if (conf_energies(i, j, prot)) {
            printf("   Fatal error reported by conf_energies(), QUIT.\n");
            fflush(stdout);
            return USERERR;
         }

         nowB = time(NULL);
         fprintf(progress_fp, "%5ld seconds\n", nowB-nowA);
      }
      /* turn back on the radii of conformers in this residue */
      for (j=1; j<prot.res[i].n_conf; j++) {
         for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
            if (!prot.res[i].conf[j].atom[k].on) continue;
            ele_bound.unf[prot.res[i].conf[j].atom[k].serial].rad = prot.res[i].conf[j].atom[k].rad;
         }
      }
   }
   printf("   Done\n\n"); fflush(stdout);

   /* reaction field energy */
   printf("   Computing RXN from conformer %d to %d of %d total conformers\n      see %s for progress...\n", env.pbe_start,env.pbe_end,n_conf, env.progress_log); fflush(stdout);

   /* set up single conformer dielectric boundary*/
   for (i=0; i<prot.n_res; i++) {

      if (prot.res[i].n_conf > 1) i_chosen = 1;
      else i_chosen = 0;
      for (j=1; j<prot.res[i].n_conf; j++) {
         if (fabs(prot.res[i].conf[j].netcrg)> 0.1) {
            i_chosen = j;
            break;
         }
      }

      for (j=1; j<prot.res[i].n_conf; j++) {
         for (k=0; k<prot.res[i].conf[j].n_atom; k++) { /* turn all side chains off */
            if (!prot.res[i].conf[j].atom[k].on) continue;
            ele_bound.unf[prot.res[i].conf[j].atom[k].serial].rad = 0.0;
         }
      }
      prot.res[i].i_bound = i_chosen;

      /* turn on one side chain */
      if (prot.res[i].n_conf>1) {
         for (k=0; k<prot.res[i].conf[i_chosen].n_atom; k++) {
            if (!prot.res[i].conf[i_chosen].atom[k].on) continue;
            ele_bound.unf[prot.res[i].conf[i_chosen].atom[k].serial].rad = prot.res[i].conf[i_chosen].atom[k].rad;
         }
      }
   }

   counter = 0;
   for (i=0; i<prot.n_res; i++) {
      for (j=1; j<prot.res[i].n_conf; j++) {
         if (counter < env.pbe_start-1 || counter > env.pbe_end-1) {
            counter++;
            continue;
         }
          /* clean records of the last run */
         for (k=1; k<=del_runs; k++) {
            sprintf(sbuff, "rxn%02d.log", k);
            remove(sbuff);
            sprintf(sbuff, "run%02d.frc", k);
            remove(sbuff);
         }
         remove("ARCDAT");

         counter++;
         nowA = time(NULL);
         fprintf(progress_fp, "   Doing rxn %5d of %5d conformers.", counter, n_conf); fflush(progress_fp);
         if (conf_rxn(i, j, prot)) {
            remove("ARCDAT");
            fprintf(progress_fp, "   Retry\n");
            n_retry++;
            fflush(stdout);
            j--;
         }
         nowB = time(NULL);
         fprintf(progress_fp, "%5ld seconds\n", nowB-nowA);
      }
      /* turn off all side chain of this residue but the one used for the dielectric boundary, prepare for the next residue */
      if (prot.res[i].n_conf>1) {
         for (k=0; k<prot.res[i].conf[prot.res[i].i_bound].n_atom; k++) { /* Restore the dielectric boundary */
            if (!prot.res[i].conf[prot.res[i].i_bound].atom[k].on) continue;
            ele_bound.unf[prot.res[i].conf[prot.res[i].i_bound].atom[k].serial].rad = prot.res[i].conf[prot.res[i].i_bound].atom[k].rad;
         }
      }
   }
   printf("   Done\n\n"); fflush(stdout);


   fclose(progress_fp);

   /* make dummy conformers */
   printf("   Add dummy conformers ... "); fflush(stdout);
   printf("%d dummy conformers were added.\n", add_dummies(prot));
   printf("   Done\n\n");

   /* write compressed pw table */
   printf("   Making pairwise energy matrices and conformer summary ... \n"); fflush(stdout);
   if (make_matrices(prot, cur_folder)) return USERERR;
   printf("   Done\n\n"); fflush(stdout);

   chdir(cur_folder);

   /* Restore directory */
   if (env.delphi_clean) {
      /* del_dir(del_folder); */
      sprintf(sbuff, "nohup rm -r %s > null", del_folder);
      system(sbuff);
   }


   if (!env.minimize_size) {
      fp = fopen(STEP3_OUT,"w");
      write_pdb(fp, prot);
      fclose(fp);
   }
   else {
      printf("   Skiping writing step 3 rotamer file %s.\n", STEP3_OUT);
   }


   nowB = time(NULL);
   printf("   Total time for step3 (energy calculation) is %ld seconds.\n\n", nowB - time_start);
   printf("   Output files (epsilon = %.2f):\n", env.epsilon_prot);
   printf("      %-16s: energy lookup table, use opp to decode the file\n", ENERGY_TABLE);
   printf("      %-16s: conformer summary\n", FN_CONFLIST3);
   printf("\n"); fflush(stdout);

   del_prot(&prot);
   return 0;
}

int conf_energies(int kr, int kc, PROT prot)
/* compute energies related to that conformer and write out to a file */
{  int i, j, k, counter, n_conf;
   char fname[MAXCHAR_LINE];
   float vdw0, vdw1,vdwt, torsion, ele0;
   float *pairwise_ele;
   FILE *fp, *fp2;
   float *potentials;      /* potentials of the submitted lines */
   char sbuff[MAXCHAR_LINE];
   float phi, fdummy;
   char notpassed, del_err;
   float weight;
   VECTOR center;


   n_retry = 0; /* reset delphi failure counter for this conformer */

   /* count conformers */
   counter=0;
   for (i=0; i<prot.n_res; i++) counter += prot.res[i].n_conf;
   n_conf = counter;

   if ((pairwise_ele = (float *) malloc(n_conf*sizeof(float))) == NULL) {
      printf("   FATAL: Memory error in conf_energies()\n");
      return USERERR;
   }

   if (!(potentials = (float *) calloc(sizeof(float), ele_bound.n))) {
      printf("   ERROR: memory error in conf_energies()\n");
      return USERERR;
   }

   if (!env.skip_ele) {
      /* eletrostatic potentials */
      /* set up the first run */
      /* fort.13 */
      for (i=1; i<prot.res[kr].n_conf; i++) { /* turn off this residue */
         for (j=0; j<prot.res[kr].conf[i].n_atom; j++) {
            if (!prot.res[kr].conf[i].atom[j].on) continue;
            ele_bound.unf[prot.res[kr].conf[i].atom[j].serial].rad = 0.0;
         }
      }
      for (i=0; i<prot.res[kr].conf[kc].n_atom; i++) { /* turn on this conformer */
         if (!prot.res[kr].conf[kc].atom[i].on) continue;
         ele_bound.unf[prot.res[kr].conf[kc].atom[i].serial].rad = prot.res[kr].conf[kc].atom[i].rad;
         ele_bound.unf[prot.res[kr].conf[kc].atom[i].serial].crg = prot.res[kr].conf[kc].atom[i].crg;
      }

      fp = fopen("fort.13", "w");
      fwrite(ele_bound.unf, sizeof(UNF), ele_bound.n, fp);
      fclose(fp);

      center.x = 0.0; center.y = 0.0; center.z = 0.0;  weight = 0.0;
      for (i=0; i<prot.res[kr].conf[kc].n_atom; i++) {
         if (!prot.res[kr].conf[kc].atom[i].on) continue;
         center.x += prot.res[kr].conf[kc].atom[i].xyz.x*fabs(prot.res[kr].conf[kc].atom[i].crg);
         center.y += prot.res[kr].conf[kc].atom[i].xyz.y*fabs(prot.res[kr].conf[kc].atom[i].crg);
         center.z += prot.res[kr].conf[kc].atom[i].xyz.z*fabs(prot.res[kr].conf[kc].atom[i].crg);
         weight += fabs(prot.res[kr].conf[kc].atom[i].crg);
      }
      center.x /= (weight+0.000001);
      center.y /= (weight+0.000001);
      center.z /= (weight+0.000001);


      fp = fopen("fort.27", "w");
      fprintf(fp, "ATOM  %5d  C   CEN  %04d    %8.3f%8.3f%8.3f\n", 1, 1,
                                                               center.x,
                                                               center.y,
                                                               center.z);
      fclose(fp);


      notpassed = 1;

      while (notpassed) {

         /* DEBUG formatted form of unformatted file sent to delphi
         printf("%s\n", prot.res[kr].conf[kc].uniqID);
         for (i=0; i<ele_bound.n; i++) {
            printf("ATOM  %5d  X   XXX  %04d    %8.3f%8.3f%8.3f%8.3f%8.3f\n", i, i,
                                                     ele_bound.unf[i].x,
                                                     ele_bound.unf[i].y,
                                                     ele_bound.unf[i].z,
                                                     ele_bound.unf[i].rad,
                                                     ele_bound.unf[i].crg);
         }
         printf("\n");
         system("cat fort.27");
         */

         fp = fopen("fort.10", "w");
         fprintf(fp, "gsize=%d\n", env.grids_delphi);
         fprintf(fp, "scale=%.2f\n", env.grids_per_ang/pow(2, del_runs-1));
         fprintf(fp, "in(unpdb,file=\"fort.13\")\n");
         fprintf(fp, "indi=%.1f\n", env.epsilon_prot);
         fprintf(fp, "exdi=%.1f\n", env.epsilon_solv);
         fprintf(fp, "ionrad=%.1f\n", env.ionrad);
         fprintf(fp, "salt=%.2f\n", env.salt);
         fprintf(fp, "bndcon=2\n");
         fprintf(fp, "center(777, 777, 0)\n");
         fprintf(fp, "out(frc,file=\"run01.frc\")\n");
         fprintf(fp, "out(phi,file=\"run01.phi\")\n");
         fprintf(fp, "site(a,c,p)\n");
         fprintf(fp, "energy(g,an,sol)\n");
         fclose(fp);


         sprintf(sbuff, "%s>delphi%02d.log 2>/dev/null", env.delphi_exe, 1);

         if (n_retry<5) {
	    system(sbuff);
	 }
         else {
            printf("   FATAL: too many failed delphi runs (%d), quitting...\n", n_retry);
            return USERERR;
         }

         for (i=1; i<del_runs; i++) {
            fp = fopen("fort.10", "w");
            fprintf(fp, "gsize=%d\n", env.grids_delphi);
            fprintf(fp, "scale=%.2f\n", env.grids_per_ang/pow(2,del_runs-1-i));
            fprintf(fp, "in(unpdb,file=\"fort.13\")\n");
            fprintf(fp, "in(phi,file=\"run%02d.phi\")\n", i);
            fprintf(fp, "indi=%.1f\n", env.epsilon_prot);
            fprintf(fp, "exdi=%.1f\n", env.epsilon_solv);
            fprintf(fp, "ionrad=%.1f\n", env.ionrad);
            fprintf(fp, "salt=%.2f\n", env.salt);
            fprintf(fp, "bndcon=3\n");
            fprintf(fp, "center(777, 777, 0)\n");
            fprintf(fp, "out(frc,file=\"run%02d.frc\")\n", i+1);
            fprintf(fp, "out(phi,file=\"run%02d.phi\")\n", i+1);
            fprintf(fp, "site(a,c,p)\n");
            fprintf(fp, "energy(g,an,sol)\n");
            fclose(fp);

            sprintf(sbuff, "%s>delphi%02d.log 2>/dev/null", env.delphi_exe, i+1);
            if (n_retry<5) {
               system(sbuff);
            }
            else {
               printf("   FATAL: too many failed delphi runs (%d), quitting...\n", n_retry);
               return USERERR;
            }
         }

         /* restore radii */
         for (i=1; i<prot.res[kr].n_conf; i++) { /* set radii of this conformer to be 0 */
            for (j=0; j<prot.res[kr].conf[i].n_atom; j++) {
               if (!prot.res[kr].conf[i].atom[j].on) continue;
               ele_bound.unf[prot.res[kr].conf[i].atom[j].serial].rad = prot.res[kr].conf[i].atom[j].rad;
            }
         }
         /* turn off the charge */
         for (i=0; i<prot.res[kr].conf[kc].n_atom; i++) {
            if (!prot.res[kr].conf[kc].atom[i].on) continue;
            ele_bound.unf[prot.res[kr].conf[kc].atom[i].serial].crg = 0.0;
         }

         /* Now collect and write the results */
         /* initialize the array with the first run */
         sprintf(fname, "run01.frc");
         if (!(fp = fopen(fname, "r"))) {
            printf("\n   WARNING: Delphi failed at focusing depth %d of %s, retry\n", 1, prot.res[kr].conf[kc].uniqID);
            n_retry++;
            remove("ARCDAT");
            del_err = 1;
            continue;
         }
         counter = 0;

         /* skip 12 lines */
         for (j=0; j<12; j++) fgets(sbuff, sizeof(sbuff), fp);

         while (fgets(sbuff, sizeof(sbuff), fp)) {
            if (strlen(sbuff)>39) {
               sscanf(sbuff+20, "%f %f", &phi, &fdummy);
               potentials[counter] = phi;
               counter++;
            }
         }
         fclose(fp);
         del_err = 0;

         /* use the first non 0 value for each atom */
         for (i=2; i<=del_runs; i++) {
            sprintf(fname, "run%02d.frc", i);
            if (!(fp = fopen(fname, "r"))) {
               printf("\n   WARNING: Delphi failed at focusing depth %d of %s, retry\n", i, prot.res[kr].conf[kc].uniqID);
               n_retry++;
               remove("ARCDAT");
               del_err = 1;
               break;
            }
            counter = 0;

            /* skip 12 lines */
            for (j=0; j<12; j++) fgets(sbuff, sizeof(sbuff), fp);

            while (fgets(sbuff, sizeof(sbuff), fp)) {
               if (strlen(sbuff)>39) {
                  sscanf(sbuff+20, "%f %f", &phi, &fdummy);
                  if (fabs(phi)>0.0001) potentials[counter] = phi;
                  counter++;
               }
            }
            fclose(fp);
            del_err = 0;
         }

         if (del_err) notpassed = 1;
         else notpassed = 0; /* so far so good */
      }
   }


   /* compute pairwise_ele, bkb moved to the bottom */
   counter = 0;
   for (i=0; i<prot.n_res; i++) {
      for (j=1; j<prot.res[i].n_conf; j++) {
         pairwise_ele[counter] = 0.0;
         for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
            if (!prot.res[i].conf[j].atom[k].on) continue;
            pairwise_ele[counter] += prot.res[i].conf[j].atom[k].crg
                                    *(potentials[prot.res[i].conf[j].atom[k].serial]);
         }
         pairwise_ele[counter] /= KCAL2KT;
         counter++;
      }
   }

   for (i=0; i<prot.n_res; i++) {
      pairwise_ele[counter] = 0.0;
      if (prot.res[i].n_conf == 0) {
         printf("   WARNING: no conformers in residue %s%c%04d\n", prot.res[i].resName,
                                                                   prot.res[i].chainID,
                                                                   prot.res[i].resSeq);
      }
      else {
         for (k=0; k<prot.res[i].conf[0].n_atom; k++) {
            if (!prot.res[i].conf[0].atom[k].on) continue;
            pairwise_ele[counter] += prot.res[i].conf[0].atom[k].crg
                               *(potentials[prot.res[i].conf[0].atom[k].serial]);
         }
         pairwise_ele[counter] /= KCAL2KT;

         counter++;
      }
   }

   free(potentials);

   /* get vdw0, intra conf vdw */
   vdw0 = 0.0;

   /* get vdw1, vdw to backbone */
   vdw1 = 0.0;

   /* torsion energy */
   torsion = torsion_conf(&prot.res[kr].conf[kc]);

   /* ele to the backbone */
   ele0 = 0.0;


   sprintf(fname, "%s.opp", prot.res[kr].conf[kc].uniqID);
   fp = fopen(fname, "w");
   counter = 0;
   for (i=0; i<prot.n_res; i++) {
      for (j=1; j<prot.res[i].n_conf; j++) {
         vdwt = vdw_conf(kr, kc, i, j, prot);
         if (vdwt > 999.0) vdwt = 999.0;
         fprintf(fp, "%05d%15s%10.3f%10.3f\n", counter,
                                            prot.res[i].conf[j].uniqID,
                                            pairwise_ele[counter],
                                            vdwt);
         if (kr==i && kc==j) {   /* I get the intra conf interaction */
            vdw0 = vdw_conf(kr, kc, i, j, prot);   /* vdw is counted */
                                                   /* ele is not counted */
         }
         counter ++;
      }
   }

   fprintf(fp, "\n");   /* separator of bkb */

   for (i=0; i<prot.n_res; i++) {
      if (prot.res[i].n_conf == 0) {
         printf("   WARNING: no conformers in residue %s%c%04d\n", prot.res[i].resName,
                                                                   prot.res[i].chainID,
                                                                   prot.res[i].resSeq);
      }
      else {
         vdwt = vdw_conf(kr, kc, i, 0, prot);
         fprintf(fp, "%05d%15s%10.3f%10.3f\n", counter,
                                            prot.res[i].conf[0].uniqID,
                                            pairwise_ele[counter],
                                            vdwt);
         vdw1 += vdwt;
         ele0 += pairwise_ele[counter];
         counter++;
      }
   }

   fprintf(fp, "\n");   /* separator of misc */

   fprintf(fp, "VDW_SELF  %10.3f\n", vdw0);
   fprintf(fp, "VDW_BKBN  %10.3f\n", vdw1);
   fprintf(fp, "ELE_BKBN  %10.3f\n", ele0);
   fprintf(fp, "TORSION   %10.3f\n", torsion);

   prot.res[kr].conf[kc].E_vdw0 = vdw0;
   prot.res[kr].conf[kc].E_vdw1 = vdw1;
   prot.res[kr].conf[kc].E_epol = ele0;
   prot.res[kr].conf[kc].E_tors = torsion;

   /* get rxn */
   fprintf(fp, "RXNMULTI  ");
   for (i=0; i<del_runs; i++) {
      sprintf(fname, "delphi%02d.log", i+1);
      fp2 = fopen(fname, "r");
      while (fgets(sbuff, sizeof(sbuff), fp2)) {
         if (strstr(sbuff, "corrected reaction field energy:")) {
            fprintf(fp, "%10.3f", atof(sbuff+34)/KCAL2KT);
            break;
         }
      }
      fclose(fp2);

   }
   fprintf(fp, "\n");



   fclose(fp);

   free(pairwise_ele);

   return 0;
}

int define_boundary(PROT prot)
{  int kr, kc, ka, i, counter;

   /* count unique atoms */
   ele_bound.n = 0;

   counter = 0;
   for (kr=0; kr<prot.n_res; kr++) {
      if (prot.res[kr].n_conf == 0) continue;
      for (ka=0; ka<prot.res[kr].conf[0].n_atom; ka++) {
         if (!prot.res[kr].conf[0].atom[ka].on) continue;
         ele_bound.unf[counter].header = 20;
         ele_bound.unf[counter].trailer = 20;
         ele_bound.unf[counter].x   = prot.res[kr].conf[0].atom[ka].xyz.x;
         ele_bound.unf[counter].y   = prot.res[kr].conf[0].atom[ka].xyz.y;
         ele_bound.unf[counter].z   = prot.res[kr].conf[0].atom[ka].xyz.z;
         ele_bound.unf[counter].rad = prot.res[kr].conf[0].atom[ka].rad;
         ele_bound.unf[counter].crg = 0.0;
         /* atom number may exceed 99999 */
         /*sprintf(ele_bound.head[counter], "ATOM  %5d %4s%c%3s %c%04d%c%03d", */
         sprintf(ele_bound.head[counter], "ATOM  %5d %4s%c%3s %c%04d%c%03d",
                            counter%100000, prot.res[kr].conf[0].atom[ka].name,
                            prot.res[kr].conf[0].altLoc,
                            prot.res[kr].resName,
                            prot.res[kr].chainID,
                            prot.res[kr].resSeq,
                            prot.res[kr].iCode,
                            0);

         prot.res[kr].conf[0].atom[ka].serial = counter;
         counter++;
      }

      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
         for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
            if (!prot.res[kr].conf[kc].atom[ka].on) continue;
            if ((i = find_atom(prot.res[kr].conf[kc].atom[ka].xyz, counter)) >= 0) {
               prot.res[kr].conf[kc].atom[ka].serial = i;
               if (ele_bound.unf[i].rad < prot.res[kr].conf[kc].atom[ka].rad)
                  ele_bound.unf[i].rad = prot.res[kr].conf[kc].atom[ka].rad;
            }
            else {
               ele_bound.unf[counter].header = 20;
               ele_bound.unf[counter].trailer = 20;
               ele_bound.unf[counter].x   = prot.res[kr].conf[kc].atom[ka].xyz.x;
               ele_bound.unf[counter].y   = prot.res[kr].conf[kc].atom[ka].xyz.y;
               ele_bound.unf[counter].z   = prot.res[kr].conf[kc].atom[ka].xyz.z;
               ele_bound.unf[counter].rad = prot.res[kr].conf[kc].atom[ka].rad;
               sprintf(ele_bound.head[counter], "ATOM  %5d %4s%c%3s %c%04d%c%03d",
                            counter%100000, prot.res[kr].conf[kc].atom[ka].name,
                            prot.res[kr].conf[kc].altLoc,
                            prot.res[kr].resName,
                            prot.res[kr].chainID,
                            prot.res[kr].resSeq,
                            prot.res[kr].iCode,
                            prot.res[kr].conf[kc].iConf);

               prot.res[kr].conf[kc].atom[ka].serial = counter;
               counter++;
            }
         }
      }
   }
   ele_bound.n = counter;
   return 0;
}

int find_atom(VECTOR v, int counter)
{  int i;

   for (i=0; i<counter; i++) {
      if (fabs(v.x-ele_bound.unf[i].x) > 0.01) continue;
      if (fabs(v.y-ele_bound.unf[i].y) > 0.01) continue;
      if (fabs(v.z-ele_bound.unf[i].z) > 0.01) continue;
      return i;
   }

   return -1;
}

int delphi_depth()
{  float x_min, x_max;
   float y_min, y_max;
   float z_min, z_max;
   float dx, dy, dz, dm;
   int i;
   float scale;
   int depth;

   x_min = x_max = ele_bound.unf[0].x;
   y_min = y_max = ele_bound.unf[0].y;
   z_min = z_max = ele_bound.unf[0].z;
   for (i=1; i<ele_bound.n; i++) {
      if (x_min > ele_bound.unf[i].x) x_min = ele_bound.unf[i].x;
      if (x_max < ele_bound.unf[i].x) x_max = ele_bound.unf[i].x;
      if (y_min > ele_bound.unf[i].y) y_min = ele_bound.unf[i].y;
      if (y_max < ele_bound.unf[i].y) y_max = ele_bound.unf[i].y;
      if (z_min > ele_bound.unf[i].z) z_min = ele_bound.unf[i].z;
      if (z_max < ele_bound.unf[i].z) z_max = ele_bound.unf[i].z;
   }

   dx = x_max - x_min;
   dy = y_max - y_min;
   dz = z_max - z_min;

   dm = (dx > dy) ? dx : dy;
   dm = (dm > dz) ? dm : dz;
   dm += env.radius_probe * 2.0 + 3.4;

   scale = env.grids_per_ang/(env.grids_delphi/(2.0*dm));

   if (scale <= 1.0) depth = 1;   /* first resolution good enough */
   else depth = (int) (log(scale)/log(2.0)) + 2;

   return depth;
}

void write_fort15()
{  int i;
   FILE *fp;

   fp = fopen("fort.15", "w");
   for (i=0; i<ele_bound.n; i++) {
      fprintf(fp, "%30s%8.3f%8.3f%8.3f\n", ele_bound.head[i],
                                                  ele_bound.unf[i].x,
                                                  ele_bound.unf[i].y,
                                                  ele_bound.unf[i].z);
   }

   fclose(fp);
   return;
}

int conf_rxn(int kr, int kc, PROT prot)
{  int i, ic, j, k, counter;
   char fname[MAXCHAR_LINE];
   float *potentials = NULL;
   FILE *fp, *fp2;
   char sbuff[MAXCHAR_LINE];
   char rxnmulti[MAXCHAR_LINE];
   float rxn[100]; /*rxn at focusing runs*/
   float rxn_min;
   VECTOR center;
   float weight;
   char uniqID[MAXCHAR_LINE];
   float k_single_multi;
   char notpassed, del_err;
   float phi, fdummy;
   float vdw0, vdw1,vdwt, torsion, ele0;
   float reliability;

   n_retry = 0; /* reset delphi failure counter for this conformer */

   if (!(potentials = (float *) calloc(sizeof(float), ele_bound.n))) {
      printf("   ERROR: memory error in conf_energies()\n");
      return USERERR;
   }

   if (!env.skip_ele) {
      /* rxn */
      /* set up the first run */
      /* fort.13 */
      /* turn off all sidechain atoms and turn on this conformer */
      for (ic=1; ic<prot.res[kr].n_conf; ic++) {
         for (i=0; i<prot.res[kr].conf[ic].n_atom; i++) {
            if (!prot.res[kr].conf[ic].atom[i].on) continue;
            ele_bound.unf[prot.res[kr].conf[ic].atom[i].serial].rad = 0.0;
            ele_bound.unf[prot.res[kr].conf[ic].atom[i].serial].crg = 0.0;
         }
      }


      for (i=0; i<prot.res[kr].conf[kc].n_atom; i++) {
         if (!prot.res[kr].conf[kc].atom[i].on) continue;
         ele_bound.unf[prot.res[kr].conf[kc].atom[i].serial].rad = prot.res[kr].conf[kc].atom[i].rad;
         ele_bound.unf[prot.res[kr].conf[kc].atom[i].serial].crg = prot.res[kr].conf[kc].atom[i].crg;
      }


      fp = fopen("fort.13", "w");
      fwrite(ele_bound.unf, sizeof(UNF), ele_bound.n, fp);
      fclose(fp);


      center.x = 0.0; center.y = 0.0; center.z = 0.0;  weight = 0.0;
      for (i=0; i<prot.res[kr].conf[kc].n_atom; i++) {
         if (!prot.res[kr].conf[kc].atom[i].on) continue;
         center.x += prot.res[kr].conf[kc].atom[i].xyz.x*fabs(prot.res[kr].conf[kc].atom[i].crg);
         center.y += prot.res[kr].conf[kc].atom[i].xyz.y*fabs(prot.res[kr].conf[kc].atom[i].crg);
         center.z += prot.res[kr].conf[kc].atom[i].xyz.z*fabs(prot.res[kr].conf[kc].atom[i].crg);
         weight += fabs(prot.res[kr].conf[kc].atom[i].crg);
      }
      center.x /= (weight+0.000001);
      center.y /= (weight+0.000001);
      center.z /= (weight+0.000001);

      fp = fopen("fort.27", "w");
      fprintf(fp, "ATOM  %5d  C   CEN  %04d    %8.3f%8.3f%8.3f\n", 1, 1,
                                                               center.x,
                                                               center.y,
                                                               center.z);
      fclose(fp);


      fp = fopen("fort.10", "w");
      fprintf(fp, "gsize=%d\n", env.grids_delphi);
      fprintf(fp, "scale=%.2f\n", env.grids_per_ang/pow(2, del_runs-1));
      fprintf(fp, "in(unpdb,file=\"fort.13\")\n");
      fprintf(fp, "indi=%.1f\n", env.epsilon_prot);
      fprintf(fp, "exdi=%.1f\n", env.epsilon_solv);
      fprintf(fp, "ionrad=%.1f\n", env.ionrad);
      fprintf(fp, "salt=%.2f\n", env.salt);
      fprintf(fp, "bndcon=2\n");
      fprintf(fp, "center(777, 777, 0)\n");
      fprintf(fp, "out(frc,file=\"run01.frc\")\n");
      fprintf(fp, "out(phi,file=\"run01.phi\")\n");
      fprintf(fp, "site(a,c,p)\n");
      fprintf(fp, "energy(g,an,sol)\n");
      fclose(fp);

      /* DEBUG formatted form of unformatted file sent to delphi
       printf("%s\n", prot.res[kr].conf[kc].uniqID);
       for (i=0; i<ele_bound.n; i++) {
          printf("ATOM  %5d  X   XXX  %04d    %8.3f%8.3f%8.3f%8.3f%8.3f\n", i, i,
                                                  ele_bound.unf[i].x,
                                                  ele_bound.unf[i].y,
                                                  ele_bound.unf[i].z,
                                                  ele_bound.unf[i].rad,
                                                  ele_bound.unf[i].crg);
       }
       printf("\n");
      */

      sprintf(sbuff, "%s>rxn%02d.log", env.delphi_exe, 1);
      if (n_retry<5) system(sbuff);
      else {
         printf("   FATAL: too many failed delphi runs (%d), quitting...\n", n_retry);
         return USERERR;
      }


      for (i=1; i<del_runs; i++) {
         fp = fopen("fort.10", "w");
         fprintf(fp, "gsize=%d\n", env.grids_delphi);
         fprintf(fp, "scale=%.2f\n", env.grids_per_ang/pow(2,del_runs-1-i));
         fprintf(fp, "in(unpdb,file=\"fort.13\")\n");
         fprintf(fp, "in(phi,file=\"run%02d.phi\")\n", i);
         fprintf(fp, "indi=%.1f\n", env.epsilon_prot);
         fprintf(fp, "exdi=%.1f\n", env.epsilon_solv);
         fprintf(fp, "ionrad=%.1f\n", env.ionrad);
         fprintf(fp, "salt=%.2f\n", env.salt);
         fprintf(fp, "bndcon=3\n");
         fprintf(fp, "center(777, 777, 0)\n");
         fprintf(fp, "out(frc,file=\"run%02d.frc\")\n", i+1);
         fprintf(fp, "out(phi,file=\"run%02d.phi\")\n", i+1);
         fprintf(fp, "site(a,c,p)\n");
         fprintf(fp, "energy(g,an,sol)\n");
         fclose(fp);
         sprintf(sbuff, "%s>rxn%02d.log", env.delphi_exe, i+1);
         if (n_retry<5) system(sbuff);
         else {
            printf("   FATAL: too many failed delphi runs (%d), quitting...\n", n_retry);
            return USERERR;
         }
      }

      /* turn off this conformer */
      for (i=0; i<prot.res[kr].conf[kc].n_atom; i++) {
         if (!prot.res[kr].conf[kc].atom[i].on) continue;
         ele_bound.unf[prot.res[kr].conf[kc].atom[i].serial].rad = 0.0;
         ele_bound.unf[prot.res[kr].conf[kc].atom[i].serial].crg = 0.0;
      }

      /* Now collect and write the results */
      /* initialize the array with the first run */
      sprintf(fname, "run01.frc");
      if (!(fp = fopen(fname, "r"))) {
         printf("\n   WARNING: Delphi failed at focusing depth %d of %s, retry\n", 1, prot.res[kr].conf[kc].uniqID);
         n_retry++;
         remove("ARCDAT");
         del_err = 1;
         return (USERERR);
      }

      /* skip 12 lines */
      for (j=0; j<12; j++) fgets(sbuff, sizeof(sbuff), fp);
      counter = 0;
      while (fgets(sbuff, sizeof(sbuff), fp)) {
         if (strlen(sbuff)>39) {
            sscanf(sbuff+20, "%f %f", &phi, &fdummy);
            potentials[counter] = phi;
            counter++;
         }
      }
      fclose(fp);
      del_err = 0;

      /* use the first non 0 value for each atom */
      for (i=2; i<=del_runs; i++) {
         sprintf(fname, "run%02d.frc", i);
         if (!(fp = fopen(fname, "r"))) {
            printf("\n   WARNING: Delphi failed at focusing depth %d of %s, retry\n", i, prot.res[kr].conf[kc].uniqID);
            n_retry++;
            remove("ARCDAT");
            del_err = 1;
            break;
         }
         counter = 0;

         /* skip 12 lines */
         for (j=0; j<12; j++) fgets(sbuff, sizeof(sbuff), fp);

         while (fgets(sbuff, sizeof(sbuff), fp)) {
            if (strlen(sbuff)>39) {
               sscanf(sbuff+20, "%f %f", &phi, &fdummy);
               if (fabs(phi)>0.0001) potentials[counter] = phi;
               counter++;
            }
         }
         fclose(fp);
         del_err = 0;
      }

      if (del_err) notpassed = 1;
      else notpassed = 0; /* so far so good */
   }


   for (i=0; i<prot.n_res; i++) {
      /* compute pairwise_ele to the dielectric boundary conformer */
      prot.res[i].pw_bound = 0.0;
      for (k=0; k<prot.res[i].conf[prot.res[i].i_bound].n_atom; k++) {
         if (!prot.res[i].conf[prot.res[i].i_bound].atom[k].on) continue;
         prot.res[i].pw_bound += prot.res[i].conf[prot.res[i].i_bound].atom[k].crg
                                 *(potentials[prot.res[i].conf[prot.res[i].i_bound].atom[k].serial])/KCAL2KT;
      }
      /* compute pairwise_ele to backbone at single conformation condition */
      prot.res[i].conf[0].tmp_pw_epol = 0.0;
      for (k=0; k<prot.res[i].conf[0].n_atom; k++) {
         if (!prot.res[i].conf[0].atom[k].on) continue;
         prot.res[i].conf[0].tmp_pw_epol += prot.res[i].conf[0].atom[k].crg
                                 *(potentials[prot.res[i].conf[0].atom[k].serial])/KCAL2KT;
      }
   }

   free(potentials);

   /* get vdw0, intra conf vdw */
   vdw0 = 0.0;
   vdw0 = vdw_conf(kr, kc, kr, kc, prot);

   /* get vdw1, vdw to backbone */
   vdw1 = 0.0;

   /* torsion energy */
   torsion = torsion_conf(&prot.res[kr].conf[kc]);

   /* ele to the backbone */
   ele0 = 0.0;

   sprintf(fname, "%s.opp", prot.res[kr].conf[kc].uniqID);

   /* read the existing opp file */
   fp = fopen(fname, "r");

   /* read */
   /* side chain */
   for (i=0; i<prot.n_res; i++) {
      for (j=1; j<prot.res[i].n_conf; j++) {
         fgets(sbuff, sizeof(sbuff), fp);
         sscanf(sbuff, "%d %s %f %f", &prot.res[i].conf[j].iConf, \
                                      uniqID, \
                                      &prot.res[i].conf[j].tmp_pw_ele, \
                                      &prot.res[i].conf[j].tmp_pw_vdw);
         if (strcmp(uniqID, prot.res[i].conf[j].uniqID)) {
             printf("   ERROR: Mismatch %s and %s in %s.opp\n", prot.res[i].conf[j].uniqID, uniqID, prot.res[kr].conf[kc].uniqID);
             return USERERR;
         }

      }
   }
   fgets(sbuff, sizeof(sbuff), fp);

   for (i=0; i<prot.n_res; i++) {
      fgets(sbuff, sizeof(sbuff), fp);
      sscanf(sbuff, "%d %s %f %f", &counter, \
                                   uniqID, \
                                   &prot.res[i].conf[0].tmp_pw_ele, \
                                   &prot.res[i].conf[0].tmp_pw_vdw);
   }

   rxnmulti[0] = '\n';
   while (fgets(sbuff, sizeof(sbuff), fp)) {
      if (!strncmp("RXNMULTI", sbuff, 8)) {
         strcpy(rxnmulti, sbuff);
         break;
      }
   }

   fclose(fp);


   /* write */
   fp = fopen(fname, "w");
   for (i=0; i<prot.n_res; i++) {
      if (i==kr)
         k_single_multi = 0.0;
      else if (fabs(prot.res[i].conf[prot.res[i].i_bound].tmp_pw_ele)> 0.1)
         k_single_multi = prot.res[i].pw_bound/prot.res[i].conf[prot.res[i].i_bound].tmp_pw_ele;
      else
         k_single_multi = 1.0;

      if (k_single_multi > 1.0) k_single_multi = 1.0;

      prot.res[i].single_multi_bound = k_single_multi;



      for (j=1; j<prot.res[i].n_conf; j++) {
         /* correction on this scaling factor, 0 to 1 with 1 is reliable.
          * (pw-ref)/ref        reliability  k_single_multi = 10.00  5.00   1.00   0.50   0.10
          *  0.01                99%                           9.77  4.92   1.00   0.50   0.10
          *  0.1                 90%                           8.03  4.29   1.00   0.53   0.10
          *  1.0                 37%                           2.33  1.81   1.00   0.77   0.43
          *  2.0                 14%                           1.37  1.24   1.00   0.91   0.73
          *  5.0                  1%                           1.02  1.01   1.00   1.00   0.98
          * 10.0                  0%                           1.00  1.00   1.00   1.00   1.00

         reliability = exp(-fabs((prot.res[i].conf[prot.res[i].i_bound].tmp_pw_ele-prot.res[i].conf[j].tmp_pw_ele))\
                                 /(fabs(prot.res[i].conf[prot.res[i].i_bound].tmp_pw_ele)+0.01));
         fprintf(fp, "%05d%15s%10.3f%10.3f%10.3f", prot.res[i].conf[j].iConf,
                                            prot.res[i].conf[j].uniqID,
                                            prot.res[i].conf[j].tmp_pw_ele * \
                                            (k_single_multi*reliability+(1-reliability)),
                                            prot.res[i].conf[j].tmp_pw_vdw,
                                            prot.res[i].conf[j].tmp_pw_ele);
         */
         /* No additional correction on k */
         fprintf(fp, "%05d%15s%10.3f%10.3f%10.3f", prot.res[i].conf[j].iConf,
                                            prot.res[i].conf[j].uniqID,
                                            prot.res[i].conf[j].tmp_pw_ele * k_single_multi,
                                            prot.res[i].conf[j].tmp_pw_vdw,
                                            prot.res[i].conf[j].tmp_pw_ele);

         if (prot.res[i].conf[prot.res[i].i_bound].tmp_pw_vdw > 50.0) fprintf(fp, "?");
         if (j==prot.res[i].i_bound) fprintf(fp, "*\n");
         else fprintf(fp, "\n");
      }
   }


   fprintf(fp, "\n");   /* separator of bkb */

   for (i=0; i<prot.n_res; i++) {
      if (prot.res[i].n_conf == 0) {
         printf("   WARNING: no conformers in residue %s%c%04d\n", prot.res[i].resName,
                                                                   prot.res[i].chainID,
                                                                   prot.res[i].resSeq);
      }
      else {
         vdwt = vdw_conf(kr, kc, i, 0, prot);
         fprintf(fp, "%05d%15s%10.3f%10.3f%10.3f\n", i,
                                            prot.res[i].conf[0].uniqID,
                                            prot.res[i].conf[0].tmp_pw_epol, /* backbone recalculated */
                                            prot.res[i].conf[0].tmp_pw_vdw,
                                            prot.res[i].conf[0].tmp_pw_ele);
         vdw1 += vdwt;
         ele0 += prot.res[i].conf[0].tmp_pw_epol;
      }
   }

   fprintf(fp, "\n");   /* separator of misc */

   fprintf(fp, "VDW_SELF  %10.3f\n", vdw0);
   fprintf(fp, "VDW_BKBN  %10.3f\n", vdw1);
   fprintf(fp, "ELE_BKBN  %10.3f\n", ele0);
   fprintf(fp, "TORSION   %10.3f\n", torsion);
   fprintf(fp, "%s", rxnmulti);

   prot.res[kr].conf[kc].E_vdw0 = vdw0;
   prot.res[kr].conf[kc].E_vdw1 = vdw1;
   prot.res[kr].conf[kc].E_epol = ele0;
   prot.res[kr].conf[kc].E_tors = torsion;

   /* get rxn at single conformation environment */
   fprintf(fp, "RXNSINGLE ");
   for (i=0; i<del_runs; i++) {
      if (!env.skip_ele) {
         sprintf(fname, "rxn%02d.log", i+1);
         if (!(fp2 = fopen(fname, "r"))) {
            printf("\n   WARNING: Delphi(rxn) failed at focusing %d of %s, retry\n", i, prot.res[kr].conf[kc].uniqID);
            n_retry++;
            return USERERR;
         }

         while (fgets(sbuff, sizeof(sbuff), fp2)) {
            if (strstr(sbuff, "corrected reaction field energy:")) {
               /* printf("\n%s", sbuff);
               printf("%s", prot.res[kr].conf[kc].uniqID);
               */
               rxn[i] = atof(sbuff+34)/KCAL2KT;
               if (strstr(sbuff, "NaN")) rxn[i] = 999.99;
               fprintf(fp, "%10.3f", rxn[i]);
               break;
            }
         }
         fclose(fp2);
      }
      else {
         rxn[i] = 0.0;
         fprintf(fp, "%10.3f", rxn[i]);
      }
   }
   fprintf(fp, "\n");


   rxn_min = 0.0;
   if (del_runs < 3) i=0;
   else i = del_runs-3;
   for (; i<del_runs; i++) {
      if (rxn_min > rxn[i]) rxn_min = rxn[i];
   }
   prot.res[kr].conf[kc].E_rxn = rxn_min; /* report the most negative in the last 3 runs*/

   fclose(fp);
   return 0;
}

int add_dummies(PROT prot)
{  int n = 0;
   int i, kr, kc;
   STRINGS confs;
   char ins;
   int natom;

   for (kr=0; kr<prot.n_res; kr++) {
      if (param_get("CONFLIST", prot.res[kr].resName, "", &confs)) {
         printf("   WARNING: add_dummies() says \"No CONFLIST for %s\"\n", prot.res[kr].resName);
      }
      if (prot.res[kr].iCode == '\0' || prot.res[kr].iCode == ' ') {
         ins = '_';
      }
      else ins = prot.res[kr].iCode;
      for (i=0; i<confs.n; i++) {
         if (param_get("NATOM", confs.strings[i], "", &natom)) {
            printf("   WARNING: No NATOM for %s\n", confs.strings[i]);
            break;
         }
         if (natom == 0 && strncmp(confs.strings[i]+3, "BK", 2)) { /* dummy conformer */
            n++;
            kc = ins_conf(&prot.res[kr], prot.res[kr].n_conf, 0);
            memset(&prot.res[kr].conf[kc], 0, sizeof(CONF));
            strncpy(prot.res[kr].conf[kc].confName, confs.strings[i], 5);
            strncpy(prot.res[kr].conf[kc].history, confs.strings[i]+3, 2);

            sprintf(prot.res[kr].conf[kc].uniqID, "%3s%c%c%c%04d%c%03d", prot.res[kr].resName,
                                                                  prot.res[kr].conf[kc].history[0],
                                                                  prot.res[kr].conf[kc].history[1],
                                                                  prot.res[kr].chainID,
                                                                  prot.res[kr].resSeq,
                                                                  ins, kc);
            prot.res[kr].conf[kc].uniqID[sizeof(prot.res[kr].conf[kc].uniqID)-1] = '\0';

         }
      }
   }

   return n;
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

/* Standalone program of prining based on pairwise interaction vector
 * test comb/1lse/ on hestia with full CPU yfsong's version: user    46m5.200s
 */
int prune_pv(PROT prot, float c1, float c2, float c3);
int over_geo(CONF conf1, CONF conf2, float cutoff);
int over_ele(PROT prot, int ir, int ic, int jc, float cutoff);
int over_vdw(PROT prot, int ir, int ic, int jc, float cutoff);

int main(int argc, char *argv[])
{  FILE *pdb_fp, *fp;
   PROT prot;
   float c1, c2, c3;
   int N;

   db_open();
   if (get_env()) {
      printf("   This program needs run.prm in working directory.\n");
      return USERERR;
   }

   if (argc < 4) {
      c1 = 0.5;
      c2 = 1.0;
      c3 = 8.0;
   }
   else {
      c1 = atof(argv[1]);
      c2 = atof(argv[2]);
      c3 = atof(argv[3]);
   }

   /* load parameters */
   if ((fp=fopen(env.new_tpl, "r"))) {
      fclose(fp);
      load_param(env.new_tpl);
      printf("%s loaded.\n",env.new_tpl);
   }
   if (strlen(env.param)) {
      if (load_all_param(env.param)) {
         printf("   Can't load parameter files in %s\n",env.param);
         return USERERR;
      }
   }

   /* load pdb */
   if ((pdb_fp=fopen(STEP2_OUT, "r"))) {
      prot = load_pdb(pdb_fp);
      fclose(pdb_fp);
      if (prot.n_res == 0) {
         printf("   Fail to load pdb file: %s\n",STEP2_OUT);
         return USERERR;
      }
   }
   else {
      printf("   Specified PDB file \"%s\" was not found\n",STEP2_OUT);
      return USERERR;
   }

   id_conf(prot);
   assign_vdw_param(prot);
   N=prune_pv(prot, c1, c2, c3);

   printf("geo_cutoff = %.3f, ele_cutoff = %.3f, vdw_cutoff = %.3f, %d conf deleted\n", c1, c2, c3, N);
   write_pdb(stdout, prot);

   db_close();
   return 0;
}


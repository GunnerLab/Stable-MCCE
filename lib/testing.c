#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mcce.h"

time_t nowA, nowB, nowStart, nowEnd;
long    idum;

int testing()
{
    FILE *fp;
    PROT prot;
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

   initial_relaxation(prot);
   
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
   printf("\n"); fflush(stdout);

   del_prot(&prot);

   return 0;
}

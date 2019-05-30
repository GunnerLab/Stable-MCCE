#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int load_occtable(PROT prot);

int main()
{
    FILE *fp;
    PROT prot;
    int kr, kc;
    
   /* Do step 0, initialization */
   db_open();
   if (init()) {
      db_close();
      printf("Help message: double check file \"run.prm\" in current directory.\n");
   }
   
   if (!(fp=fopen("step2_out.pdb", "r"))) {
      printf("   FATAL: No input \"step2_out.pdb\".\n");
      return USERERR;
   }
   prot = load_pdb(fp);
   if (prot.n_res == 0) {
      printf("   There are errors in pdb file, quiting ...\n");
      return USERERR;
   }
   fclose(fp);
   
   id_conf(prot);
   assign_rad(prot);
   
   fp = fopen("acc.cnf","w");
   sas_native(prot);
   load_occtable(prot);
   for (kr=0; kr<prot.n_res; kr++) {
       float sas_min = 99999;
       float sas_max = -1;
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
          if (prot.res[kr].conf[kc].sas > sas_max) sas_max = prot.res[kr].conf[kc].sas;
          if (prot.res[kr].conf[kc].sas < sas_min) sas_min = prot.res[kr].conf[kc].sas;
      }
      
      if (sas_max - sas_min < 25.) continue;
      for (kc=1; kc<prot.res[kr].n_conf; kc++) {
          fprintf(fp, "%s %8.2f %8.3f %8.3f %8.3f %8.3f",
              prot.res[kr].conf[kc].uniqID,
              prot.res[kr].conf[kc].occ,
              prot.res[kr].conf[kc].sas,
              prot.res[kr].conf[kc].sas_fraction,
              sas_min,
              sas_max);
          if (sas_max - sas_min > 1e-3) fprintf(fp, " %8.3f\n", (prot.res[kr].conf[kc].sas - sas_min) / (sas_max - sas_min) );
          else fprintf(fp, " %8.3f\n", 0.);
      }
   }
   fclose(fp);
   db_close();
   return 0;
}

int load_occtable(PROT prot)
{
    FILE *fp;
    char sbuff[160], confID[15];
    int n_titra, i_titra;
    int kr, kc, found;
    fp = fopen("fort.38", "r");
    if (!fp) {
        printf("   FATAL: Can't open file fort.38\n");
        return USERERR;
    }
    fgets(sbuff, sizeof(sbuff), fp); /* the first line */
    //if (strstr(sbuff, "ph")) titration_type = 1; /*  pH titration */
    n_titra = (int) ((strlen(sbuff) -14.) / 6.);
    for (kr=0; kr<prot.n_res; kr++) {
        for (kc=1; kc<prot.res[kr].n_conf; kc++) {
            prot.res[kr].conf[kc].occ = 0.;
        }
    }
    
    while(fgets(sbuff, sizeof(sbuff), fp)) {
        strncpy(confID, sbuff, 14); confID[14] = '\0';
   
        found = 0;
        for (kr=0; kr<prot.n_res; kr++) {
            for (kc=1; kc<prot.res[kr].n_conf; kc++) {
                if (!strcmp(prot.res[kr].conf[kc].uniqID, confID)) { /* found conf */
                    found = 1;
                    prot.res[kr].conf[kc].occ_table = malloc(n_titra * sizeof(float));
                    for (i_titra = 0; i_titra<n_titra; i_titra++) {
                        char temp_str[10];
                        strncpy(temp_str, sbuff+14+6*i_titra, 6); temp_str[6] = '\0';
                        prot.res[kr].conf[kc].occ_table[i_titra] = atof(strtok(temp_str, " "));
                        prot.res[kr].conf[kc].occ += prot.res[kr].conf[kc].occ_table[i_titra];
                    }
                    break;
                }
            }
            if (found) break;
        }
    }
    fclose(fp);
    return 0;
}

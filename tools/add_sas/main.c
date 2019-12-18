



#include <stdio.h>
#include <string.h>
#include "mcce.h"

extern PROT monte2_load_conflist(char *fname);

int main()
{
    FILE *fp;
    PROT prot, prot_head3;
    int ic, i_res, i_conf;
    
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
   get_connect12(prot);
   setup_vdw_fast(prot);
   assign_vdw_param(prot);
   get_vdw0(prot);
   for (i_res=0; i_res<prot.n_res; i_res++) {
       for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
           prot.res[i_res].conf[i_conf].E_torsion = torsion_conf(&prot.res[i_res].conf[i_conf]);
       }
   }
   
   prot_head3 = monte2_load_conflist("head3.lst");
   
   for (ic = 0; ic < prot_head3.nc; ic++) {
       
       for (i_res = 0; i_res< prot.n_res; i_res++) {
           for (i_conf = 1; i_conf < prot.res[i_res].n_conf; i_conf++) {
               
               if (!strcmp(prot_head3.conf[ic]->uniqID, prot.res[i_res].conf[i_conf].uniqID)) {
                   prot_head3.conf[ic]->E_vdw0 = prot.res[i_res].conf[i_conf].E_vdw0;
                   prot_head3.conf[ic]->E_tors = prot.res[i_res].conf[i_conf].E_torsion;
               }
           }
       }
   }
   
   fp = fopen("head3.lst", "w");
   fprintf(fp, "iConf CONFORMER     FL  occ    crg   Em0  pKa0 ne nH    vdw0    vdw1    tors    epol   dsolv   extra    history\n");
   for (ic=0; ic<prot_head3.nc; ic++) {
         fprintf(fp, "%05d %s %c %4.2f %6.3f %5.0f %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %11s\n",
                                                                    ic+1,
                                                                    prot_head3.conf[ic]->uniqID,
                                                                    'f', 0.00,
                                                                    prot_head3.conf[ic]->netcrg,
                                                                    prot_head3.conf[ic]->Em,
                                                                    prot_head3.conf[ic]->pKa,
                                                                    prot_head3.conf[ic]->e,
                                                                    prot_head3.conf[ic]->H,
                                                                    prot_head3.conf[ic]->E_vdw0,
                                                                    prot_head3.conf[ic]->E_vdw1,
                                                                    prot_head3.conf[ic]->E_tors,
                                                                    prot_head3.conf[ic]->E_epol,
                                                                    prot_head3.conf[ic]->E_dsolv,
                                                                    prot_head3.conf[ic]->E_extra,
                                                                    prot_head3.conf[ic]->history);
   }
   
   db_close();
   return 0;
}


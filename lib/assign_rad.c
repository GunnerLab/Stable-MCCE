#include <stdio.h>
#include "mcce.h"

int assign_rad(PROT prot)
{
    int i, j, k;
    float r;
    FILE *debug_fp;
    int err=0;
    
    for (i=0; i<prot.n_res; i++) {
        for (j=0; j<prot.res[i].n_conf; j++) {
            for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
                if (prot.res[i].conf[j].atom[k].on) {
                    if (param_get("RADIUS", prot.res[i].conf[j].confName, prot.res[i].conf[j].atom[k].name, &r)) { /* first try on conf level */
                        if (param_get("RADIUS", prot.res[i].resName, prot.res[i].conf[j].atom[k].name, &r)) {
                            debug_fp = fopen(env.debug_log,"a");
                            fprintf(debug_fp, "RADIUS   %s   %s   %4.2f\n", prot.res[i].resName, prot.res[i].conf[j].atom[k].name, env.default_radius);
                            fclose(debug_fp);
                            err=1;
                            prot.res[i].conf[j].atom[k].rad = env.default_radius;
                            param_sav("RADIUS", prot.res[i].resName, prot.res[i].conf[j].atom[k].name, &env.default_radius, sizeof(float));
                        }
                        else prot.res[i].conf[j].atom[k].rad = r;
                    }
                    else prot.res[i].conf[j].atom[k].rad = r;
                }
            }
        }
    }
    if (err) printf("   Warning! assign_rad():      missing parameter(s), default value is used and saved in %s.\n", env.debug_log); 
    return 0;
}

int assign_vdw_param(PROT prot)
{  int i, j, k;
   float val;
   FILE *debug_fp;
   int err=0;

   for (i=0; i<prot.n_res; i++) {
       for (j=0; j<prot.res[i].n_conf; j++) {
           for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
               if (prot.res[i].conf[j].atom[k].on) {
                   if (param_get("VDW_RAD", prot.res[i].conf[j].confName, prot.res[i].conf[j].atom[k].name, &val)) {
                       /* amber94 values*/
                       if (prot.res[i].conf[j].atom[k].name[1] == 'C') { /* carbon of a methyl group */
                           prot.res[i].conf[j].atom[k].vdw_rad = 1.908;
                       }
                       else if (prot.res[i].conf[j].atom[k].name[1] == 'O') { /* oxygen of a carbonyl group */
                           prot.res[i].conf[j].atom[k].vdw_rad = 1.6612;
                       }
                       else if (prot.res[i].conf[j].atom[k].name[1] == 'H') { /* proton of a methyl group */
                           prot.res[i].conf[j].atom[k].vdw_rad = 1.487;
                       }
                       else if (prot.res[i].conf[j].atom[k].name[1] == 'N') { /* amide */
                           prot.res[i].conf[j].atom[k].vdw_rad = 1.824;
                       }
                       else if (prot.res[i].conf[j].atom[k].name[1] == 'S') { /* cys/met */
                           prot.res[i].conf[j].atom[k].vdw_rad = 2.0;
                       }
                       else {                                                   /* using carbon value */
                           prot.res[i].conf[j].atom[k].vdw_rad = 1.908;
                       }
                       
                       err=1;
                       param_sav("VDW_RAD", prot.res[i].conf[j].confName, prot.res[i].conf[j].atom[k].name, &prot.res[i].conf[j].atom[k].vdw_rad, sizeof(float));
                       debug_fp = fopen(env.debug_log,"a");
                       fprintf(debug_fp, "VDW_RAD  %s %s   %4.2f\n",
                       prot.res[i].conf[j].confName, prot.res[i].conf[j].atom[k].name, prot.res[i].conf[j].atom[k].vdw_rad);
                       fclose(debug_fp);
                   }
                   else prot.res[i].conf[j].atom[k].vdw_rad = val;

                   if (param_get("VDW_EPS", prot.res[i].conf[j].confName, prot.res[i].conf[j].atom[k].name, &val)) {
                       if (prot.res[i].conf[j].atom[k].name[1] == 'C') {
                           prot.res[i].conf[j].atom[k].vdw_eps = 0.1094;
                       }
                       else if (prot.res[i].conf[j].atom[k].name[1] == 'O') {
                           prot.res[i].conf[j].atom[k].vdw_eps = 0.21;
                       }
                       else if (prot.res[i].conf[j].atom[k].name[1] == 'H') {
                           prot.res[i].conf[j].atom[k].vdw_eps = 0.0157;
                       }
                       else if (prot.res[i].conf[j].atom[k].name[1] == 'N') {
                           prot.res[i].conf[j].atom[k].vdw_eps = 0.17;
                       }
                       else if (prot.res[i].conf[j].atom[k].name[1] == 'S') {
                           prot.res[i].conf[j].atom[k].vdw_eps = 0.25;
                       }
                       else {
                           prot.res[i].conf[j].atom[k].vdw_eps = 0.1094;
                       }
                       
                       err=1;
                       debug_fp = fopen(env.debug_log,"a");
                       fprintf(debug_fp, "VDW_EPS  %s %s   %4.2f\n", prot.res[i].conf[j].confName, prot.res[i].conf[j].atom[k].name, prot.res[i].conf[j].atom[k].vdw_eps);
                       fclose(debug_fp);
                       param_sav("VDW_EPS", prot.res[i].conf[j].confName, prot.res[i].conf[j].atom[k].name, &prot.res[i].conf[j].atom[k].vdw_eps, sizeof(float));
                   }
                   else prot.res[i].conf[j].atom[k].vdw_eps = val;
                   
               }
           }
      }
   }
   if (err) printf("   Warning! assign_vdw_param(): missing parameter(s), default value is used and saved in %s.\n", env.debug_log); 
   return 0;
}


#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

float vdw_conf(int ires, int iconf, int jres, int jconf, PROT prot)
{  float  e = 0.0;
   int    iatom, jatom;
   ATOM   *atomi, *atomj;
   ATOM   *connect123[MAX_CONNECTED2];
   ATOM   *connect14[MAX_CONNECTED3];
   int    connect123_res[MAX_CONNECTED2];
   int    connect14_res[MAX_CONNECTED3];
   int    n_connect123, n_connect14;
   int    c1, c2, c3;

   for (iatom=0; iatom<prot.res[ires].conf[iconf].n_atom; iatom++) {
      atomi=&prot.res[ires].conf[iconf].atom[iatom];
      if (!atomi->on) continue;
      n_connect123 = 0;
      n_connect14  = 0;
      memset(connect123, 0,MAX_CONNECTED2*sizeof(ATOM *));
      memset(connect14,  0,MAX_CONNECTED3*sizeof(ATOM *));

      for (c1=0; c1<MAX_CONNECTED; c1++) {
         if (!atomi->connect12[c1]) break;
         n_connect123++;
         connect123[n_connect123-1] = atomi->connect12[c1];
         connect123_res[n_connect123-1] = atomi->connect12_res[c1];

         for (c2= 0; c2<MAX_CONNECTED; c2++) {
            if (!atomi->connect12[c1]->connect12[c2]) break;
            n_connect123++;
            connect123[n_connect123-1] = atomi->connect12[c1]->connect12[c2];
            connect13_res[n_connect13-1] = atomi->connect12[c1]->connect12_res[c2];

            for (c3 = 0; c3 < MAX_CONNECTED; c3++) {
               if (!atomi->connect12[c1]->connect12[c2]->connect12[c3]) break;
               n_connect14++;
               connect14[n_connect14-1] = atomi->connect12[c1]->connect12[c2]->connect12[c3];
               connect14_res[n_connect14-1] = atomi->connect12[c1]->connect12[c2]->connect12_res[c3];
            }
         }
      }

      for (jatom=0; jatom<prot.res[jres].conf[jconf].n_atom; jatom++) {
         atomj = &prot.res[jres].conf[jconf].atom[jatom];
         if (!atomj->on) continue;

         for (c1 = 0; c1<n_connect123; c1++) {
            if (jres == connect123_res[c1] && !strcmp(atomj->name,connect123[c1]->name)) break;
         }
         if (c1 < n_connect13) continue;

         for (c1= 0; c1 < n_connect14; c1++) {
            if (jres == connect14_res[c1] && !strcmp(atomj->name,connect14[iconnect]->name)) break;
         }
         if (iconnect < n_connect14) {
            e += vdw(*iatom_p, *jatom_p) * env.factor_14lj;
            continue;
         }
         e += vdw(*iatom_p, *jatom_p);
      }
   }

    return e;
}

#define VDW_CUTOFF_NEAR  1
#define VDW_ELIMIT_NEAR  999
#define VDW_CUTOFF_FAR   10
#define VDW_ELIMIT_FAR   0

/* This function uses AMBER C12 and C6 parameters to calculate
* Van der Waals potential. The parameters are read from paramter
* file amber.vdw and were originally cited from:
* http://solar.uits.indiana.edu/ad_html/Using_AutoDock_305.a.shtml#29613
*/

float cutoff_near2 = VDW_CUTOFF_NEAR * VDW_CUTOFF_NEAR;
float cutoff_far2  = VDW_CUTOFF_FAR  * VDW_CUTOFF_FAR;

float vdw(ATOM atom1, ATOM atom2)
{
    float       e = 0;
    float       d2, d6, d12;      /* distance with power 2, 6, and 12 */
    float       C12, C6;
    char        pair[4];
    VECTOR      v1, v2;
    FILE        *debug_fp;

    if (!atom1.on) return e;
    if (!atom2.on) return e;

    v1 = atom1.xyz;
    v2 = atom2.xyz;

    d2 = ddvv(v1, v2);

    if (d2 > cutoff_far2)
        e = VDW_ELIMIT_FAR;             /* Cutoff */
    else if (d2 < cutoff_near2)
        e = VDW_ELIMIT_NEAR;            /* Cutoff */
    else
    {                                   /* 12-6 L-J potential */
        pair[0] = atom1.name[1];        /* element name */
        pair[1] = '-';
        pair[2] = atom2.name[1];        /* element name */
        pair[3] = '\0';


        if (param_get("VDWAMBER", "C12", pair, &C12)) {
            param_get("VDWAMBER", "C12", "X-X", &C12);
            debug_fp = fopen(env.debug_log,"a");
            fprintf(debug_fp, "   Warning! vdw_amber(): No C12 vdw param for pair %s, default = %15.3f.\n", pair, C12);
            fclose(debug_fp);
        }

        if(param_get("VDWAMBER", "C6", pair, &C6)) {
            param_get("VDWAMBER", "C6", "X-X", &C6);
            debug_fp = fopen(env.debug_log,"a");
            fprintf(debug_fp, "   Warning! vdw_amber(): No C6  vdw param for pair %s, default = %15.3f.\n", pair, C6);
            fclose(debug_fp);
        }

        d6 = d2*d2*d2;
        d12 = d6*d6;
        e = C12/d12 - C6/d6;
    }

    return e;
}

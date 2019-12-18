#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

#define VDW_CUTOFF_NEAR  1
#define VDW_ELIMIT_NEAR  999
#define VDW_CUTOFF_FAR   10
#define VDW_ELIMIT_FAR   0

static float cutoff_near2 = VDW_CUTOFF_NEAR * VDW_CUTOFF_NEAR;
static float cutoff_far2  = VDW_CUTOFF_FAR  * VDW_CUTOFF_FAR;

float a2a_coulomb(ATOM atom1, ATOM atom2)
/* return unit is Kcal/mol */
{  float dx, dy, dz, d;
   dx = atom1.xyz.x - atom2.xyz.x;
   dy = atom1.xyz.y - atom2.xyz.y;
   dz = atom1.xyz.z - atom2.xyz.z;

   d = sqrt(dx*dx+dy*dy+dz*dz);
   if (d<0.5) return 999.0;
   else return 331.5*atom1.crg*atom2.crg/d;
}

float CoulombBySAS(ATOM atom1, ATOM atom2)
/* return unit is Kcal/mol */
{  float dx, dy, dz, d;
   dx = atom1.xyz.x - atom2.xyz.x;
   dy = atom1.xyz.y - atom2.xyz.y;
   dz = atom1.xyz.z - atom2.xyz.z;

   d = sqrt(dx*dx+dy*dy+dz*dz);
   if (d<0.5) return 999.0;
   else return 331.5*atom1.crg*atom2.crg/d * exp(-d/7.8) * (1.25-atom1.sas)*(1.25-atom2.sas)*0.44;
}

float RxnBySAS(CONF conf)
{  float dsolv, rxn0;
   int ia, n_atom;
   float av_sas, av_crg, rxn;

   n_atom = 0;
   av_sas = 0.0;
   av_crg = 0.0;

   for (ia=0; ia<conf.n_atom; ia++) {
      if (!conf.atom[ia].on) continue;
      n_atom++;
      av_sas += conf.atom[ia].sas;
      av_crg += fabs(conf.atom[ia].crg);
      //factor += conf.atom[ia].sas*fabs(conf.atom[ia].crg);
   }

   av_sas /= n_atom;
   av_crg /= n_atom;

   dsolv = (-400*av_sas+64)/env.epsilon_prot*fabs(conf.netcrg);

   if (param_get("RXN", conf.confName, "", &rxn0)) {
      printf("   WARNING: No RXN entry for %s, set to 0.\n", conf.confName);
      fflush(stdout);
      rxn0 = 0.0;
      param_sav("RXN", conf.confName, "", &rxn0, sizeof(float));
   }

   //return dsolv;
   if (n_atom<6) n_atom=6;

   rxn =(dsolv+rxn0)*0.8*log(n_atom-2.5);

   if ((rxn - rxn0)<0) return rxn0;
   else return rxn;
}

float Ecoulomb_conf2conf(PROT prot, int ir, int ic, int kr, int kc, float epsilon)
{  int ia, ka;
   float E = 0;
   
   for (ia=0; ia<prot.res[ir].conf[ic].n_atom; ia++) {
      if (!prot.res[ir].conf[ic].atom[ia].on) continue;
      if (fabs(prot.res[ir].conf[ic].atom[ia].crg) < 0.01) continue;
      for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
         if (!prot.res[kr].conf[kc].atom[ka].on) continue;
         if (fabs(prot.res[kr].conf[kc].atom[ka].crg) < 0.01) continue;
         E+=a2a_coulomb(prot.res[ir].conf[ic].atom[ia], prot.res[kr].conf[kc].atom[ka]);
      }
   }

   return E/epsilon;
}

float Evdw_conf2conf(PROT prot, int ir, int ic, int kr, int kc)
{  int ia, ka;
   float E = 0;
   
   for (ia=0; ia<prot.res[ir].conf[ic].n_atom; ia++) {
      if (!prot.res[ir].conf[ic].atom[ia].on) continue;
      for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
         if (!prot.res[kr].conf[kc].atom[ka].on) continue;
         if (fabs(prot.res[ir].conf[ic].atom[ia].xyz.x-prot.res[kr].conf[kc].atom[ka].xyz.x)>10.0) continue;
         E+=vdw(prot.res[ir].conf[ic].atom[ia], prot.res[kr].conf[kc].atom[ka]);
      }
   }
   
   return E;
}


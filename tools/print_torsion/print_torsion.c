#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mcce.h"
extern float torsion_angle(VECTOR v0, VECTOR v1, VECTOR v2, VECTOR v3);
extern int water_orient(VECTOR v0, VECTOR v1, VECTOR v2, float *theta, float *phi, float *psi);
float PI,d2r;

int main(int argc, char *argv[]) {
    PROT prot;
    FILE *fp;
    int kr,kc,ka,warn,nc;
    int i_res,i_conf,i_atom;
    ATOM *atom_p[4];
    float angle,theta,phi,psi;
    int tag;
    PI = 4.*atan(1.);
    d2r = PI/180.;        
	if (argc < 2) {
		printf("%s pdb_file\n",argv[0]);
		return USERERR;
	}

    db_open();
    if (init()) {
        db_close();
        printf("Help message: mcce requires a file \"run.prm\" in current directory.\n");
        return USERERR;
    }

   if (!(fp=fopen(argv[1], "r"))) {
       db_close();
      printf("   No file %s\n", argv[1]);
      return USERERR;
   }
   prot=load_pdb(fp);
   fclose(fp);

//   assign_crg(prot);
//   assign_rad(prot);
   get_connect12(prot);
   
   id_conf(prot);
   
   for (i_res=0; i_res<prot.n_res; i_res++) {
       for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
           printf("%s ",prot.res[i_res].conf[i_conf].uniqID);
           if (strcmp(prot.res[i_res].resName,"HOH")) {
               for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
				TORS tors;
				if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                   if (prot.res[i_res].conf[i_conf].atom[i_atom].name[1] == 'H') continue;
                   if (torsion_atoms(&prot.res[i_res].conf[i_conf],i_atom,&atom_p[0],&atom_p[1],&atom_p[2],&atom_p[3],&tors,1)) {
                       //printf("%s",prot.res[i_res].conf[i_conf].atom[i_atom].name);
                       if (torsion_atoms(&prot.res[i_res].conf[i_conf],i_atom,&atom_p[0],&atom_p[1],&atom_p[2],&atom_p[3],&tors, 0)) {
                           continue;
                       }
                   }
                   
                   angle = torsion_angle(atom_p[0]->xyz,atom_p[1]->xyz,atom_p[2]->xyz,atom_p[3]->xyz);
                   //printf("%f,%s-%s-%s-%s\n",angle,atom_p[0]->name,atom_p[1]->name,atom_p[2]->name,atom_p[3]->name);
                   tag = (int) ((angle/d2r)+0.5); if (tag == 360) tag = 0;
                   printf("%4s %03d|",atom_p[0]->name,tag);
               }
               for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
				   TORS tors;
                   if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                   if (prot.res[i_res].conf[i_conf].atom[i_atom].name[1] != 'H') continue;
                   if (torsion_atoms(&prot.res[i_res].conf[i_conf],i_atom,&atom_p[0],&atom_p[1],&atom_p[2],&atom_p[3],&tors,1)) continue;
                   angle = torsion_angle(atom_p[0]->xyz,atom_p[1]->xyz,atom_p[2]->xyz,atom_p[3]->xyz);
                   tag = (int) ((angle/d2r) + 0.5); if (tag == 360) tag = 0;
                   printf("%4s %03d|",atom_p[0]->name,tag);
               }
           }
           else {
               if (prot.res[i_res].conf[i_conf].n_atom>=3) {
                   water_orient(prot.res[i_res].conf[i_conf].atom[0].xyz,
                   prot.res[i_res].conf[i_conf].atom[1].xyz,
                   prot.res[i_res].conf[i_conf].atom[2].xyz,
                   &theta, &phi, &psi);
                   tag = (int) ((theta/d2r) + 0.5); if (tag == 360) tag = 0;
                   printf("thet %03d|",tag);
                   tag = (int) ((phi/d2r) + 0.5); if (tag == 360) tag = 0;
                   printf("phi  %03d|",tag);
                   tag = (int) ((psi/d2r) + 0.5); if (tag == 360) tag = 0;
                   printf("psi  %03d|",tag);
               }
           }
           printf("\n");
       }
   }

   db_close();

   return 0;
}


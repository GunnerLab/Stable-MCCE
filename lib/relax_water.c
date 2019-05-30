#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"
#define  RELAX_THR2 env.water_relax_thr * env.water_relax_thr

int relax_water(PROT prot) {
    int i_res,i_conf,i_atom,j_res,j_conf,j_atom,k_res,k_conf,k_atom;
    int n_conf, add, counter;
    FILE *debug_fp;

    for (i_res=0; i_res<prot.n_res; i_res++) {
        RES *ires_p = &prot.res[i_res];
        if (strcmp(prot.res[i_res].resName, "HOH")) continue;
	n_conf = prot.res[i_res].n_conf;
        counter = 0;
        for (i_conf=1; i_conf<n_conf; i_conf++) {
            CONF *iconf_p = &prot.res[i_res].conf[i_conf];
	    if (!iconf_p->n_atom) continue;
            
            for (j_res=0; j_res<prot.n_res; j_res++) {
                if (!strcmp(prot.res[j_res].resName, "HOH")) continue;
                for (j_conf=0; j_conf<prot.res[j_res].n_conf; j_conf++) {
                    CONF *jconf_p = &prot.res[j_res].conf[j_conf];
                    ATOM *iatom_p;
                    
                    CONF conf;
                    conf.n_atom = iconf_p->n_atom;
                    conf.atom  = malloc(conf.n_atom*sizeof(ATOM));
                    cpy_conf(&conf,iconf_p);
                    for (i_atom=0;i_atom<conf.n_atom;i_atom++) {
                        if (!conf.atom[i_atom].on) continue;
                        if (conf.atom[i_atom].name[1] == 'O') break;
                    }
                    if (i_atom<conf.n_atom) iatom_p = &conf.atom[i_atom];
                    else break;
                    
                    add = 0;
                    for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
                        ATOM *jatom_p = &jconf_p->atom[j_atom];
			if (!jatom_p->on) continue;
                        if (jatom_p->name[1] == 'H') continue;
                        if (ddvv(iatom_p->xyz,jatom_p->xyz) > RELAX_THR2) continue;
                        add = 1;

                        for (i_atom=0;i_atom<conf.n_atom;i_atom++) {
                            if (!conf.atom[i_atom].on) continue;
                            
                            conf.atom[i_atom].xyz = 
                            vector_vplusv(jatom_p->xyz,
                            vector_rescale(vector_normalize(
                            vector_vminusv(conf.atom[i_atom].xyz, jatom_p->xyz)), env.water_relax_thr));
                        }
                        
                        /*
                        ins_conf(ires_p, ires_p->n_conf, iconf_p->n_atom);
                        iconf_p = &prot.res[i_res].conf[i_conf];
                        if (cpy_conf(&ires_p->conf[ires_p->n_conf-1], &ires_p->conf[i_conf])) {printf("   Error! relax_water(): couldn't copy the conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",ires_p->conf[i_conf].confName,ires_p->resName, ires_p->resSeq, ires_p->n_conf-1); return USERERR;}

                        
                        for (i_conf2=1; i_conf2<ires_p->n_conf-1; i_conf2++) {
                            if (!cmp_conf(ires_p->conf[i_conf2], ires_p->conf[ires_p->n_conf-1],0.05)) break;
                        }
                        if ( i_conf2 < ires_p->n_conf-1 ) {
                            del_conf(ires_p, ires_p->n_conf-1);
                            iconf_p = &prot.res[i_res].conf[i_conf];
                            continue;
                        }
                        */
                        //printf("   Debugging! residue %s%4d nconf=%d, residue %s%4d, distance %8.3f\n", prot.res[i_res].resName,prot.res[i_res].resSeq,ires_p->n_conf, prot.res[j_res].resName,prot.res[j_res].resSeq,dvv(iatom_p->xyz,jatom_p->xyz));
                    }
                    
                    if (add) {
                    for (j_atom=0; j_atom<prot.res[j_res].conf[j_conf].n_atom; j_atom++) {
                        ATOM *jatom_p = &jconf_p->atom[j_atom];
			if (!jatom_p->on) continue;
                        if (jatom_p->name[1] == 'H') continue;
                        if (ddvv(iatom_p->xyz,jatom_p->xyz) > RELAX_THR2) continue;
                        add = 0;
                    }
                    }
                    
                    if (add) {
                    for (k_conf=1; k_conf<ires_p->n_conf; k_conf++) {
                        if (!cmp_conf(ires_p->conf[k_conf], conf, 0.05)) break;
                    }
                    if ( k_conf<ires_p->n_conf ) add = 0;
                    }
                    
                    if (add) {
                    for (k_res=0; k_res<prot.n_res; k_res++) {
                        k_conf = 0;
                        for (k_atom=0; k_atom<prot.res[k_res].conf[k_conf].n_atom; k_atom++) {
                            ATOM *katom_p = &prot.res[k_res].conf[k_conf].atom[k_atom];
                            if (!katom_p->on) continue;
                            if (katom_p->name[1] == 'H') continue;
                            if (ddvv(iatom_p->xyz,katom_p->xyz) > RELAX_THR2) continue;
                            add = 0;
                        }
                    }
                    }

                    if (add) {
                        counter++;
                        ins_conf(ires_p, ires_p->n_conf, conf.n_atom);
                        iconf_p = &prot.res[i_res].conf[i_conf];
                        if (cpy_conf(&ires_p->conf[ires_p->n_conf-1], &conf)) {printf("   Error! relax_water(): couldn't copy a new conformer \"%s\" in residue %s %d, to new position k_conf = %d\n",conf.confName,ires_p->resName, ires_p->resSeq, ires_p->n_conf-1); return USERERR;}
                        debug_fp = fopen(env.debug_log,"a");
                        fprintf(debug_fp,"add conformer to HOH %c%04d to relax against %s %c%04d\n",ires_p->chainID,ires_p->resSeq,prot.res[j_res].resName,prot.res[j_res].chainID,prot.res[j_res].resSeq);
                        fclose(debug_fp);
                    }
                    
                    free(conf.atom);
                }
            }
            
        }
    }
    return 0;
}

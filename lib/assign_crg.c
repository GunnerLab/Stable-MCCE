#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int assign_crg(PROT prot) {
    int k_res,k_conf,k_atom;
    CONF *conf_p;
    ATOM *atom_p;
    FILE *debug_fp;
    int  err=0;
    for (k_res=0; k_res<prot.n_res; k_res++) {
        for (k_conf=0; k_conf<prot.res[k_res].n_conf; k_conf++) {
            conf_p = &prot.res[k_res].conf[k_conf];
            for (k_atom=0; k_atom<prot.res[k_res].conf[k_conf].n_atom; k_atom++) {
                atom_p = &conf_p->atom[k_atom];
                if (!atom_p->on) continue;
                
                if ( param_get( "CHARGE", conf_p->confName, atom_p->name, &atom_p->crg) )
                {
                    debug_fp = fopen(env.debug_log,"a");
                    if (!debug_fp) debug_fp = stdout;
                    fprintf(debug_fp, "CHARGE   %s %s     0.000\n", conf_p->confName, atom_p->name);
                    fclose(debug_fp);
                    err = 1;
                    atom_p->crg = 0;
                    param_sav("CHARGE", conf_p->confName, atom_p->name, &atom_p->crg, sizeof(float));
                }
                //printf("   Debugging! residue %s%4d,conformer %s atom %s crg %f\n", prot.res[k_res].resName,prot.res[k_res].resSeq,conf_p->confName,atom_p->name,atom_p->crg);
            }
        }
    }
    if (err) printf("   Warning! assign_crg():      missing parameter(s), default value is used and saved in %s.\n", env.debug_log);
    return 0;
}


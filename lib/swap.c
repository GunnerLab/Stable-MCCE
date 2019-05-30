#include <math.h>
#include "mcce.h"

int rot_swap(PROT prot) {
    int i_res, i_conf, n_conf, ins, i_atom, j_atom, i_swap, counter;
    SWAP_RULE swap_rule;
    char sbuff[MAXCHAR_LINE];
    VECTOR swap_xyz;
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        RES *res_p = &prot.res[i_res];
        
        counter = 0;
        while (1) {
            /* go over each ROT_SWAP parameter associated with this residue */
            sprintf(sbuff, "%d", counter);
            if (param_get("ROT_SWAP", res_p->resName, sbuff, &swap_rule)) break;
            counter++;
            
            /* for each conformer, make a new conformer and apply swapping rule */
            n_conf = prot.res[i_res].n_conf;
            for (i_conf=1; i_conf<n_conf; i_conf++) {
                CONF *conf_p = &prot.res[i_res].conf[i_conf];

                /* make new conformer */
                ins = ins_conf(res_p, res_p->n_conf, conf_p->n_atom);
                conf_p = &prot.res[i_res].conf[i_conf];
                if (ins == USERERR) {
                    printf("   FATAL: rot_swap(): Can't add conformer. Maybe out of memory!\n");
                    return USERERR;
                }
                cpy_conf(&res_p->conf[ins], conf_p);
                //strncpy(res_p->conf[ins].history+2, "SP", 2);
                
                /* apply swapping rule */
                for (i_swap=0; i_swap<swap_rule.n_swap; i_swap++) {
                    if (param_get("IATOM",res_p->conf[ins].confName,swap_rule.swap_atom1[i_swap],&i_atom)) continue;
                    if (param_get("IATOM",res_p->conf[ins].confName,swap_rule.swap_atom2[i_swap],&j_atom)) continue;
                    
                    swap_xyz = res_p->conf[ins].atom[i_atom].xyz;
                    res_p->conf[ins].atom[i_atom].xyz = res_p->conf[ins].atom[j_atom].xyz;
                    res_p->conf[ins].atom[j_atom].xyz = swap_xyz;
                }
            }
        }
    }
    return 0;
}



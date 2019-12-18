#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

int main(int argc, char *argv[]) {
    int i_arg, i_res,i_conf, j_res,j_conf;
    FILE *fp;
    char jconf_uniqID[20];
    PROT prot;
    
    memset(jconf_uniqID,   0,20*sizeof(char));

    if (argc < 2) {
        printf("   Usage: atom_e confID [-p param_dir] [-vdw confID]\n");
        printf("   -p param_dir: define the parameter folder if there is no run.prm\n");
        printf("   -vdw confID:  define the second conformer to calculate the pairwise vdw\n");
        return -1;
    }
    
    db_open(); /* initialize gdbm database */
    if (get_env()) { /* load run.prm */
        printf("   No run.prm, use default settings\n");
    }
    
    /* set up command line parameters */
    for (i_arg=2; i_arg<argc; i_arg++) {
        if (!strcmp(argv[i_arg],"-p")) {
            if (argc > i_arg+1) {
                strcpy(env.param, argv[i_arg+1]);
            }
        }
        if (!strcmp(argv[i_arg],"-vdw")) {
            if (argc > i_arg+1) {
                strcpy(jconf_uniqID, argv[i_arg+1]);
            }
        }
        
    }
    
    /* load new.tpl */
    if ((fp=fopen(env.new_tpl, "r"))) {
        fclose(fp);
        load_param(env.new_tpl);
        printf("%s loaded.\n",env.new_tpl);
    }
    
    /* load parameter folder */
    if (strlen(env.param)) {
        if (load_all_param(env.param)) {
            printf("   Can't load parameter files in %s\n",env.param);
            return USERERR;
        }
    }
    
    /* load pdb file */
    if (!(fp=fopen(STEP2_OUT, "r"))) {
        printf("   FATAL: pdb file \"%s\".\n", STEP2_OUT);
        return USERERR;
    }
    prot = load_pdb(fp);
    fclose(fp);
    
    /* set up conformer ID */
    id_conf(prot);
    
    assign_vdw_param(prot);
    get_connect12(prot);
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            float E_torsion;
            if (strcmp(prot.res[i_res].conf[i_conf].uniqID, argv[1]) != 0) continue;
            
            setup_vdw_fast(prot);
            setup_connect_res(prot,i_res);
            
            E_torsion = torsion_conf_print(&prot.res[i_res].conf[i_conf]);
            printf("Torsion: %s total, e=%8.3f\n",prot.res[i_res].conf[i_conf].uniqID, E_torsion);

            /* pairwise */
            if (strlen(jconf_uniqID)) {
                for (j_res=0; j_res<prot.n_res; j_res++) {
                    for (j_conf=0; j_conf<prot.res[j_res].n_conf; j_conf++) {
                        float E_vdw;
                        if (strcmp(prot.res[j_res].conf[j_conf].uniqID, jconf_uniqID) != 0) continue;
                        printf("\n");
                        E_vdw = vdw_conf_fast_print(i_res,i_conf,j_res,j_conf,prot);
                        printf("VDW total btw \"     %s\" and \"     %s\": %8.3f\n",
                            prot.res[i_res].conf[i_conf].uniqID,
                            prot.res[j_res].conf[j_conf].uniqID,
                            E_vdw);
                    }
                }
                
            }
        }
    }
    db_close();
    
    return 0;
}

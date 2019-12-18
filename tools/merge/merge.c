#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "mcce.h"
extern int rm_comment(char *target, char *str);

int main(int argc, char *argv[]) {
    int i_arg, i, i_res,i_conf,k_res,k_conf,n_atom,n_prot;
    FILE *fp, *pdb_fp;
    char outpdb[256], list[256], sbuff[256], line[256];
    PROT prot, merged_prot;
    
    memset(outpdb, 0,256*sizeof(char));
    memset(list,   0,256*sizeof(char));
    memset(sbuff,  0,256*sizeof(char));
    memset(line,   0,256*sizeof(char));

    printf("Usage: merge_pdb -i input_list -o output_pdb -p param_dir\n");

    db_open();
    if (get_env()) {
        printf("   No run.prm, use default settings\n");
    }
    
    for (i_arg=1; i_arg<argc; i_arg++) {
        if (!strcmp(argv[i_arg],"-p")) {
            if (argc > i_arg+1) {
                strcpy(env.param, argv[i_arg+1]);
            }
        }
        if (!strcmp(argv[i_arg],"-i")) {
            if (argc > i_arg+1) {
                strcpy(list, argv[i_arg+1]);
            }
        }
        if (!strcmp(argv[i_arg],"-o")) {
            if (argc > i_arg+1) {
                strcpy(outpdb, argv[i_arg+1]);
            }
        }
    }
    if ((fp=fopen(env.new_tpl, "r"))) {
        fclose(fp);
        load_param(env.new_tpl);
        printf("%s loaded.\n",env.new_tpl);
    }
    if (strlen(env.param)) {
        if (load_all_param(env.param)) {
            printf("   Can't load parameter files in %s\n",env.param);
            return USERERR;
        }
    }
    if (!(fp=fopen(list, "r"))) {
        printf("   Can't load input pdb list: %s\n",list);
        db_close();
        exit(-1);
    }

    n_prot=0;
    memset(&merged_prot,0,sizeof(PROT));
    while (fgets(sbuff, sizeof(line), fp)) {
        rm_comment(line, sbuff);
        for(i=strlen(line)-1;i>=0;i--) {
            if (line[i]==' ') line[i]='\0';
            else if (line[i]=='\t') line[i]='\0';
            else break;
        }
        if (!strlen(line)) continue;

        if ((pdb_fp=fopen(line, "r"))) {
            prot = load_pdb(pdb_fp);
            fclose(pdb_fp);
            if (prot.n_res == 0) {
                printf("   Fail to load pdb file %s\n",line);
                continue;
            }
        }
        else {
            printf("   Can't load pdb file: %s\n",line);
            continue;
        }

        n_prot++;
        for (i_res = 0; i_res < prot.n_res; i_res++) {
            /* find matched residue */
            for (k_res = 0; k_res < merged_prot.n_res; k_res++) {
                if (!strcmp(prot.res[i_res].resName, merged_prot.res[k_res].resName) &&
                    prot.res[i_res].chainID == merged_prot.res[k_res].chainID &&
                prot.res[i_res].resSeq  == merged_prot.res[k_res].resSeq  &&
                prot.res[i_res].iCode   == merged_prot.res[k_res].iCode  )
                {
                    break;
                }
            }
            if (k_res >= merged_prot.n_res) {
                k_res = ins_res(&merged_prot, merged_prot.n_res);
                strcpy(merged_prot.res[k_res].resName, prot.res[i_res].resName);
                merged_prot.res[k_res].chainID = prot.res[i_res].chainID;
                merged_prot.res[k_res].resSeq  = prot.res[i_res].resSeq;
                merged_prot.res[k_res].iCode   = prot.res[i_res].iCode;

                /* Insert a backbone conformer */
                n_atom = prot.res[i_res].conf[0].n_atom;
                ins_conf(&merged_prot.res[k_res], 0, n_atom);
                cpy_conf(&merged_prot.res[k_res].conf[0], &prot.res[i_res].conf[0]);

                /*
                sprintf(sbuff, "F%03X", n_prot);
                strncpy(merged_prot.res[k_res].conf[merged_prot.res[k_res].n_conf-1].history+2, sbuff, 4);
                strncpy(merged_prot.res[k_res].conf[merged_prot.res[k_res].n_conf-1].history+6, "C000", 4);
                */
            }

            /* paste side chains */
            for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
                /* Find matched conformer */
                for (k_conf = 1; k_conf < merged_prot.res[k_res].n_conf; k_conf++) {
                    /* same conformer */
                    if (!strcmp(prot.res[i_res].conf[i_conf].confName, merged_prot.res[k_res].conf[k_conf].confName)) {
                        if (!cmp_conf(prot.res[i_res].conf[i_conf], merged_prot.res[k_res].conf[k_conf], 0.001)) break;
                    }
                }
                if (k_conf < merged_prot.res[k_res].n_conf) {
                    if (strncmp("O000",prot.res[i_res].conf[i_conf].history+2,4)) {
                        merged_prot.res[k_res].conf[k_conf].history[2] = 'C';
                    }
                    continue;
                }
                if (merged_prot.res[k_res].n_conf>=1000) continue;
                ins_conf(&merged_prot.res[k_res], merged_prot.res[k_res].n_conf, prot.res[i_res].conf[i_conf].n_atom);
                cpy_conf(&merged_prot.res[k_res].conf[merged_prot.res[k_res].n_conf-1], &prot.res[i_res].conf[i_conf]);

                if (strncmp("O000",prot.res[i_res].conf[i_conf].history+2,4)) {
                sprintf(sbuff, "F%03X", n_prot);
                strncpy(merged_prot.res[k_res].conf[merged_prot.res[k_res].n_conf-1].history+2, sbuff, 4);
                sprintf(sbuff, "C%s", prot.res[i_res].conf[i_conf].confID);
                strncpy(merged_prot.res[k_res].conf[merged_prot.res[k_res].n_conf-1].history+6, sbuff, 4);
                }
            }
        }
        del_prot(&prot);
    }
    fclose(fp);

    sort_res(merged_prot);
    sort_conf(merged_prot);

    if (!(fp=fopen(outpdb, "w"))) {
        fp = stdout;
    }
    write_pdb(fp, merged_prot);
    if (fp != stdout) fclose(fp);
    db_close();
    return 0;
}

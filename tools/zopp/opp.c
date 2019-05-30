#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "mcce.h"

RES      conflist;

int load_conflist_dir(char *dir);

int main(int argc, char *argv[])
{
    EMATRIX ematrix, ematrix2;
    int i, j, counter,args;
    struct stat buf;
    int is_dir;
    char option;
    char param[256];
    char match;
    char fname[256];
    char sbuff[256];
    char tempID[256];
    int int_dummy;
    float temp_ele,temp_vdw, temp_crt, temp_ori;
    char  temp_mark[4];
    FILE *fp;
    int verbose;
    char **arguments;
    
    db_open();
    if (get_env()) {
        printf("   This program needs run.prm in working directory.\n");
        db_close();
        return USERERR;
    }
    
    if (argc<3) {
		printf("Zopp V2.01 \n\
            This program manipilates energy lookup table in working directory.\n\
            \n\
            Usage:\n\
            zopp <options> [-x dir | -a <dir1...dirN> | -d conf | -u dir]\n\
            \n\
            Options:\n\
            -v         Verbose\n\
            Energy matrix operations:\n\
            -x:        extract the energy lookup table to the directory named dir\n\
            -a:        append \"energies.opp\" from dir1...dirN to the current one in working directory\n\
            -u:        update energies.opp and head3.lst with a new calculation\n\
            -d conf:   display pairwise interaction of conf (ID as in head3.lst)\n\
            -c:        Create energies.opp from opp files and head3.lst in a directory\n");
            return 0;
    }
    else {
        verbose = 0;
        for (i=0; i<argc; i++) {
			if (argv[i][0] == '-') {
				if (argv[i][1] == 'v') {
					verbose = 1;
					//printf("---- Verbose-Mode ----\n");
				}
				else {
					option = argv[i][1];
					if (option == 'a' || option == 'u') {
						//printf(" Append mode\n");
						args = argc-i-1;
						arguments = malloc(args * sizeof(char *));
						printf(" Read in %d arguments\n",args);
						for (j=0; j<args; j++) {
							//printf("Dir-argument %d is %s\n",j , argv[i+j+1]);
							arguments[j] = malloc ((strlen(argv[j+i+1])+1) * sizeof(char));
							strcpy(arguments[j], argv[j+i+1]);
						}
					}
					else {
                        args = 1;
                        j=0;
                        //printf("Forced %d argument\n",args);
                        arguments = malloc(args * sizeof(char *));
                        arguments[j] = malloc ((strlen(argv[j+i+1])+1) * sizeof(char));
                        strcpy(arguments[j], argv[j+i+1]);
					}
				}
			}
        }
    }
    
    
    
    switch (option) {
    case 'x':
        ematrix.n = 0;
        if (load_energies(&ematrix, ".", verbose) < 0) {
            printf("   Error loading energy table\n");
            return USERERR;
        }
        if ( verbose == 1 ) {
			printf("creating directory %s\n", arguments[0]); fflush(stdout);
        }/* create a directory */
        if (stat(arguments[0], &buf)) { /* no such file or dir */
            if (mkdir(arguments[0], 0755)) {
                printf("   FATAL: Failed creating directory %s, no write permission.\n", arguments[0]);
                return USERERR;
            }
        }
        else { /* there is a file or dir named as this */
            is_dir = S_ISDIR(buf.st_mode);
            if (!is_dir) { /* it is not a directory */
                printf("   FATAL: %s exists but is not a directory.\n", arguments[0]);
                return USERERR;
            }
        }
        if (extract_matrix(&ematrix, arguments[0], verbose)) {
            printf("   Error extracting energy table\n");
            return USERERR;
        }
        free_ematrix(&ematrix);
        break;
        
    case 'a':
        ematrix.n = 0;
        if ( verbose == 1 ) {
			printf(" Loading first matrix (energies.opp in the current dir)\n"); fflush(stdout);
        }
		if (load_energies(&ematrix, ".", verbose) < 0) {
            printf("   Error loading energy table\n");
            return USERERR;
        }
		for (i=0; i<args; i++) {
            if ( verbose == 1 ) {
                printf(" Loading argument/matrix %d:%s\n",i, arguments[i]); fflush(stdout);
            }
            /* load the second one */
            if (load_energies(&ematrix, arguments[i], verbose) < 0) {
                printf("   Error loading energy table at append mode\n");
                return USERERR;
            }
            
		}
        
        /* write the mergied */
        write_energies(&ematrix, ".", verbose);
        
        free_ematrix(&ematrix);
        break;
    case 'u':
        ematrix.n = 0;
        if (load_energies(&ematrix, ".", verbose ) < 0) {
            printf("   Error loading energy table\n");
            return USERERR;
        }
        /* load the second one */
        ematrix2.n = 0;
        if (load_energies(&ematrix2, argv[2], verbose) < 0) {
            printf("   Error loading energy table\n");
            return USERERR;
        }
        /* update, will use the conformers in matrix 2 */
        for (i=0; i<ematrix2.n; i++) {
            if (ematrix2.conf[i].on == 't') { /* always the new calculation */
                memcpy(&ematrix.conf[i], &ematrix2.conf[i], sizeof(CONF_HEAD));
                memcpy(ematrix.pw[i], ematrix2.pw[i], ematrix.n*sizeof(PAIRWISE));
            }
            else { /* get a number from the old matrix */
                //pass;
            }
        }
        /* write the updated */
        write_energies(&ematrix2, ".", verbose);
        free_ematrix(&ematrix);
        free_ematrix(&ematrix2);
        break;
    case 'd':
        ematrix.n = 0;
        if (load_energies(&ematrix, ".", verbose) < 0) {
            printf("   Error loading energy table\n");
            return USERERR;
        }
        
        
        match = '\0';
        for (i=0; i<ematrix.n; i++) {
            if (!strcmp(ematrix.conf[i].uniqID, param)) {
                counter = 0;
                for (j=0; j<ematrix.n; j++) {
                    counter++;
                    printf("%05d %s %8.3f%8.3f%8.3f%8.3f %s\n", counter, ematrix.conf[j].uniqID, ematrix.pw[i][j].ele, ematrix.pw[i][j].vdw, ematrix.pw[i][j].crt, ematrix.pw[i][j].ori, ematrix.pw[i][j].mark);
                }
                match = 't';
            }
        }
        if (match != 't') printf("conformer %s not found in energy table\n", param);
        free_ematrix(&ematrix);
        
        break;
    case 'c': /* creat energies from a directory of text opp files */
        /* read in head3.lst */
        load_conflist_dir(argv[2]);
        
        /* make an empty matrix */
        if (!(ematrix.conf = (CONF_HEAD *) calloc(conflist.n_conf, sizeof(CONF_HEAD)))) {
            printf("   Memory error in A load_energies\n");
            return 0; /* none loaded and memory cleared */
        }
        if (!(ematrix.pw = (PAIRWISE **) calloc(conflist.n_conf, sizeof(PAIRWISE *)))) {
            printf("   Memory error in B load_energies\n");
            return 0; /* none loaded and memory cleared */
        }
        for (i=0; i<conflist.n_conf; i++) {
            if (!(ematrix.pw[i] = (PAIRWISE *) calloc(conflist.n_conf, sizeof(PAIRWISE)))) {
                printf("   Memory error in C load_energies\n");
                return 0; /* none loaded and memory cleared */
            }
        }
        ematrix.n = conflist.n_conf;
        
        /* load in self energy */
        for (i=0; i<conflist.n_conf; i++) {
            strncpy(ematrix.conf[i].uniqID, conflist.conf[i].uniqID, 14); ematrix.conf[i].uniqID[14] = '\0';
            ematrix.conf[i].netcrg = conflist.conf[i].netcrg;
            ematrix.conf[i].Em = conflist.conf[i].Em;
            ematrix.conf[i].pKa = conflist.conf[i].pKa;
            ematrix.conf[i].e = conflist.conf[i].e;
            ematrix.conf[i].H = conflist.conf[i].H;
            ematrix.conf[i].E_vdw0 = conflist.conf[i].E_vdw0;
            ematrix.conf[i].E_vdw1 = conflist.conf[i].E_vdw1;
            ematrix.conf[i].E_epol = conflist.conf[i].E_epol;
            ematrix.conf[i].E_rxn0 = conflist.conf[i].E_rxn0;
            ematrix.conf[i].E_rxn = conflist.conf[i].E_rxn;
            ematrix.conf[i].E_tors = conflist.conf[i].E_tors;
            ematrix.conf[i].E_dsolv = conflist.conf[i].E_dsolv;
            ematrix.conf[i].E_extra = conflist.conf[i].E_extra;
            strncpy(ematrix.conf[i].history, conflist.conf[i].history, 10);  ematrix.conf[i].history[10]= '\0';
        }
        
        /* load opp files */
        for (i=0; i<conflist.n_conf; i++) {
            sprintf(fname, "%s/%s.opp", argv[2], ematrix.conf[i].uniqID);
            if (!(fp=fopen(fname, "r"))) { /* not exist, assume no interactions  */
                printf("   opp file %s does not exist, assumed no interactins\n", fname);
                ematrix.conf[i].on = 'f';
            }
            else {
                ematrix.conf[i].on = 't';
                counter = 0;
                while(fgets(sbuff, sizeof(sbuff), fp)) {
                    sscanf(sbuff, "%d %s %f %f %f %f %s", &int_dummy, tempID, &temp_ele, &temp_vdw, &temp_crt, &temp_ori, temp_mark);
                    while (strcmp(tempID, ematrix.conf[counter].uniqID)) { /* not match, proceed to the next one that matches */
                        counter ++;
                        if (counter >= ematrix.n) { /* reached the end but no match in head3.lst was found */
                            printf("   Conformer %s in file %s was not matched by head3.lst\n", tempID, fname);
                            return USERERR;
                        }
                    }
                    if (counter < ematrix.n) { /* a match was found */
                        ematrix.pw[i][counter].ele = temp_ele;
                        ematrix.pw[i][counter].vdw = temp_vdw;
                        ematrix.pw[i][counter].crt = temp_crt;
                        ematrix.pw[i][counter].ori = temp_ori;
                        strcpy(ematrix.pw[i][counter].mark, temp_mark);
                    }
                    else { /* this line doesn't have match in ematrix */
                        printf("      The line \"%s\" in file %s does not have a match in %s/%s\n", sbuff, fname, argv[2],FN_CONFLIST3);
                        return USERERR;
                    }
                    
                    counter ++;
                }
                fclose(fp);
            }
        }
        
        
        /* write the mergied */
        write_energies(&ematrix, ".", verbose);
        
        free_ematrix(&ematrix);
        break;
        
        
        
    default:
        printf("unrecognized option -%c", option);
        return 0;
    }
    
    
    db_close();
    return 0;
}

int load_conflist_dir(char *dir)
{
    FILE *fp;
    char sbuff[MAXCHAR_LINE];
    char stemp[MAXCHAR_LINE];
    char fname[256];
    CONF conf_temp;
    int iconf;
    int counter;
    
    conflist.n_conf = 0;
    conflist.conf   = NULL;
    
    
    sprintf(fname, "%s/%s", dir, FN_CONFLIST3);
    if (!(fp=fopen(fname, "r"))) {
        printf("   FATAL: Can't open file %s\n", FN_CONFLIST3);
        return USERERR;
    }
    fgets(sbuff, sizeof(sbuff), fp); /* skip the first line */
    counter = 0;
    while(fgets(sbuff, sizeof(sbuff), fp)) {
        /* load this line to a conf template */
        if (strlen(sbuff) < 20) continue;
        sscanf(sbuff, "%d %s %c %f %f %f %f %d %d %f %f %f %f %f %f %s", &conf_temp.iConf,
            conf_temp.uniqID,
            &conf_temp.on,
            &conf_temp.occ,
            &conf_temp.netcrg,
            &conf_temp.Em,
            &conf_temp.pKa,
            &conf_temp.e,
            &conf_temp.H,
            &conf_temp.E_vdw0,
            &conf_temp.E_vdw1,
            &conf_temp.E_tors,
            &conf_temp.E_epol,
            &conf_temp.E_dsolv,
            &conf_temp.E_extra,
            conf_temp.history);
        
        strncpy(conf_temp.resName, conf_temp.uniqID, 3); conf_temp.resName[3] = '\0';
        strncpy(conf_temp.confName, conf_temp.uniqID, 5); conf_temp.confName[5] = '\0';
        conf_temp.chainID = conf_temp.uniqID[5];
        strncpy(stemp, conf_temp.uniqID+6, 4); stemp[4] = '\0';
        conf_temp.resSeq = atoi(stemp);
        conf_temp.iCode = conf_temp.uniqID[10];
        conf_temp.n_atom = 0;
        if (conf_temp.on == 't' || conf_temp.on == 'T') conf_temp.on = 't';
        else conf_temp.on = 'f';
        conf_temp.iConf = counter;
        /* creating conflist */
        iconf = ins_conf(&conflist, conflist.n_conf, 0);
        cpy_conf(&conflist.conf[iconf], &conf_temp);
        counter++;
    }
    
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

#define BOND_THR   2.2
#define BOND_THR_H 1.8

int create_param(FILE *pdb_fp, int k_line) {
    CONF conf;
    int cnt, ichar, i_atom, j_atom;
    char sbuffer[MAXCHAR_LINE], line[MAXCHAR_LINE];
    FILE *param_fp;
    
    /* Collect a list of atoms from the same residue. */
    rewind(pdb_fp);
    cnt = -1;
    while (cnt < k_line && fgets(sbuffer, MAXCHAR_LINE, pdb_fp) ) {
        cnt++;
    }
    if (cnt != k_line) {
        printf(" Error! create_param(): Error in reaching defined line number. line counter = %d, defined line number = %d\n",cnt,k_line);
        return USERERR;
    }
    if (strncmp("ATOM  ", sbuffer, 6) && strncmp("HETATM", sbuffer, 6)) {
        printf(" Error! create_param(): input line is not an atom line. \n");
        return USERERR;
    }
    
    strcpy(line, sbuffer);
    
    memset(&conf,0,sizeof(CONF));
    rewind(pdb_fp);
    while ( fgets(sbuffer, MAXCHAR_LINE, pdb_fp) ) {
        if (strncmp("ATOM  ", sbuffer, 6) && strncmp("HETATM", sbuffer, 6)) continue;
        if (!strncmp(sbuffer+16, line+16, 4) &&
            !strncmp(sbuffer+21, line+21, 8) &&
        !strncmp(sbuffer+80, line+80, 2) ) {
            conf.n_atom++;
            conf.atom = realloc(conf.atom, conf.n_atom*sizeof(ATOM));
            memset(&conf.atom[conf.n_atom-1],0,sizeof(ATOM));
            conf.atom[conf.n_atom-1] = pdbline2atom(sbuffer);
        }
    }
    
    for (ichar=0;ichar<3;ichar++) {
        if (conf.atom[0].resName[ichar] == ' ') {
            printf(" Error! create_param(): A space found in residue name %s\n",conf.atom[0].resName);
            return USERERR;
        }
    }
    
    param_fp = fopen(env.new_tpl,"a");
    /* CONFLSIT */
    strcpy(sbuffer, "CONFLIST ");
    strcat(sbuffer, conf.atom[0].resName);
    strcat(sbuffer, "        ");
    
    strcat(sbuffer, conf.atom[0].resName);
    strcat(sbuffer, "BK ");
    
    strcat(sbuffer, conf.atom[0].resName);
    strcat(sbuffer, "01 ");
    
    fprintf(param_fp,"%s\n\n",sbuffer);
    
    /* NATOM */
    strcpy(sbuffer, "NATOM    ");
    strcat(sbuffer, conf.atom[0].resName);
    strcat(sbuffer, "BK      0");
    fprintf(param_fp,"%s\n",sbuffer);
    
    strcpy(sbuffer, "NATOM    ");
    strcat(sbuffer, conf.atom[0].resName);
    strcat(sbuffer, "01      ");
    sprintf(sbuffer+20, "%d", conf.n_atom);
    fprintf(param_fp,"%s\n\n",sbuffer);
    
    /* IATOM */
    for (i_atom=0; i_atom<conf.n_atom; i_atom++) {
        strcpy(sbuffer, "IATOM    ");
        strcat(sbuffer, conf.atom[0].resName);
        strcat(sbuffer, "01 ");
        strcat(sbuffer, conf.atom[i_atom].name);
        sprintf(sbuffer+19, " %4d", i_atom);
        fprintf(param_fp,"%s\n",sbuffer);
    }
    fprintf(param_fp,"\n");
    
    /* ATOMNAME */
    for (i_atom=0; i_atom<conf.n_atom; i_atom++) {
        strcpy(sbuffer, "ATOMNAME ");
        strcat(sbuffer, conf.atom[0].resName);
        sprintf(sbuffer+12, "01 %4d", i_atom);
        strcat(sbuffer, " ");
        strcat(sbuffer, conf.atom[i_atom].name);
        fprintf(param_fp,"%s\n",sbuffer);
    }
    fprintf(param_fp,"\n");
    
    /* CONNECT */
    for (i_atom=0; i_atom<conf.n_atom; i_atom++) {
        strcpy(sbuffer, "CONNECT  ");
        strcat(sbuffer, conf.atom[0].resName);
        strcat(sbuffer, "01 ");
        strcat(sbuffer, conf.atom[i_atom].name);
        strcat(sbuffer, " ion      ");
        
        for (j_atom=0; j_atom<conf.n_atom; j_atom++) {
            if (j_atom == i_atom) continue;
            if (conf.atom[i_atom].name[1] == 'H' || conf.atom[j_atom].name[1] == 'H') {
                if (dvv(conf.atom[i_atom].xyz, conf.atom[j_atom].xyz) > BOND_THR_H) continue;
            }
            else {
                if (dvv(conf.atom[i_atom].xyz, conf.atom[j_atom].xyz) > BOND_THR) continue;
            }
            strcat(sbuffer, "  0   ");
            strcat(sbuffer, conf.atom[j_atom].name);
        }
        fprintf(param_fp,"%s\n",sbuffer);
    }
    fprintf(param_fp,"\n");
    
    fclose(param_fp);
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

PROT load_pdb(FILE *fp)
{   PROT  prot;
    ATOM  atom;
    char  line[MAXCHAR_LINE];   /* line buffer */
    char  confName[6];
    int   k_res;
    int   k_conf;
    int   n_atom;
    int   k_atom;
    char  Fatal = 0;

    memset(&prot, 0, sizeof(PROT));

    /*
    int c;
    while ((c=fgetc(fp)) != EOF) fputc(c, stdout);
    fputc(EOF, stdout);

    rewind(fp);
    */

    while ( memset(line,0,MAXCHAR_LINE*sizeof(char)),
    fgets(line, MAXCHAR_LINE, fp)) {
        if (strncmp("ATOM  ", line, 6) && strncmp("HETATM", line, 6)) continue;
        if (line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0'; /* this is to make sure conformer history not end with a new line */
        atom = pdbline2atom(line);

        /* search for the residue: each unique combination of
        residue name, chain ID, residue number and insertion code
        defines one residue. */
        for (k_res = prot.n_res - 1; k_res >= 0; k_res--) {
            if (!strcmp(atom.resName, prot.res[k_res].resName) &&
                atom.chainID == prot.res[k_res].chainID &&
                atom.resSeq  == prot.res[k_res].resSeq  &&
                atom.iCode   == prot.res[k_res].iCode  )
            {
                break;
            }
        }
        /* If couldn't find the residue, add a new one */
        if (k_res == -1) {
            k_res = ins_res(&prot, prot.n_res);
            strcpy(prot.res[k_res].resName, atom.resName);
            prot.res[k_res].chainID = atom.chainID;
            prot.res[k_res].resSeq  = atom.resSeq;
            prot.res[k_res].iCode   = atom.iCode;
            if (strncmp("HETATM", line, 6)) prot.res[k_res].groupID   = 0;
            else prot.res[k_res].groupID = 1;


            /* Insert a backbone conformer */
            strcpy(confName, atom.resName);
            strcat(confName, "BK");
            if (param_get("NATOM", confName, "", &n_atom)) {
                printf("   WARNING: load_pdb(): \"No NATOM for conformer %s, set to 0.\"\n",confName);
                n_atom=0;
                continue;
            }
            ins_conf(&prot.res[k_res], 0, n_atom);
            strcpy(prot.res[k_res].conf[0].confName, confName);
            strcpy(prot.res[k_res].conf[0].confID, "000");
            if (!strncmp(atom.history,"BK",2)) strcpy(prot.res[k_res].conf[0].history, atom.history);
            else strcpy(prot.res[k_res].conf[0].history,"BK________");
            prot.res[k_res].conf[0].altLoc = ' ';
        }

        /* search for the conformer */
        if (!strcmp(atom.confName+3, "BK"))
            k_conf = 0;
        else {
            for (k_conf = prot.res[k_res].n_conf - 1; k_conf >= 0; k_conf--) {
                if (!strcmp(atom.confID, prot.res[k_res].conf[k_conf].confID)) {
                    break;
                }
            }
        }
        /* If couldn't find the conformer, add a new one */
        if (k_conf == -1) {
            if (param_get("NATOM", atom.confName, "", &n_atom)) {
                printf("   WARNING: load_pdb(): \"No NATOM for conformer %s, set to 0.\"\n",atom.confName);
                continue;
            }
            k_conf = ins_conf(&prot.res[k_res], prot.res[k_res].n_conf, n_atom);
            strcpy(prot.res[k_res].conf[k_conf].confName, atom.confName);
            strcpy(prot.res[k_res].conf[k_conf].confID, atom.confID);
            prot.res[k_res].conf[k_conf].altLoc = atom.altLoc;
            strcpy(prot.res[k_res].conf[k_conf].history, atom.history);
        }

        if (param_get("IATOM", atom.confName, atom.name, &k_atom)) {
            printf("   FATAL: load_pdb(): \"No IATOM of conformer \"%s\" atom \"%s\".\n", atom.confName, atom.name);
            Fatal = 1;
            continue;
        }

        if (k_atom >= prot.res[k_res].conf[k_conf].n_atom) {
            printf("   FATAL: load_pdb(): Confliction in  IATOM of conformer \"%s\" atom \"%s\" and NATOM of conformer \"%s\".\n", atom.confName, atom.name, atom.confName);
            Fatal = 1;
            continue;
        }


        if (prot.res[k_res].conf[k_conf].atom[k_atom].on) {
            if (k_conf == 0) {  /* backbone atoms, ignore and continue */
                printf("   Warning: load_pdb(): \"Duplicate backbone atom ignored, \"%3s %3s %c %3d\".\n",atom.name, atom.resName, atom.chainID, atom.resSeq);
            }
            else {
                printf("   FATAL: load_pdb(): \"Atom already exists, \"%3s %3s %c %3d\".\n",
                atom.name, atom.resName, atom.chainID, atom.resSeq);
                Fatal = 1;
            }
            continue;
        }

        prot.res[k_res].conf[k_conf].atom[k_atom] = atom;
        prot.res[k_res].conf[k_conf].atom[k_atom].on = 1;

        if ( !strlen(prot.res[k_res].conf[k_conf].history) )
            strcpy(prot.res[k_res].conf[k_conf].history, atom.history);
    }

    if (Fatal) {
       del_prot(&prot);
       printf("   STOP: Fix fatal errors before proceeding.\n");
       printf("   Try delete local new.tpl, then check the loaded tpl files\n");
    }

    return prot;
}


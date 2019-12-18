#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

PROT load_pdb_no_param(FILE *fp)
{   PROT  prot;
    ATOM  atom;
    char  line[MAXCHAR_LINE];   /* line buffer */
    char  confName[6];
    int   k_res;
    int   k_conf;
    int   k_atom;
    memset(&prot, 0, sizeof(PROT));


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

            /* Insert a backbone conformer */
            strcpy(confName, atom.resName);
            strcat(confName, "BK");
            ins_conf(&prot.res[k_res], 0, 0);
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
            k_conf = ins_conf(&prot.res[k_res], prot.res[k_res].n_conf, 0);
            strcpy(prot.res[k_res].conf[k_conf].confName, atom.confName);
            strcpy(prot.res[k_res].conf[k_conf].confID, atom.confID);
            prot.res[k_res].conf[k_conf].altLoc = atom.altLoc;
            strcpy(prot.res[k_res].conf[k_conf].history, atom.history);
        }

        /* Insert atom */
        /*--- increase the atom list size by 1 */
        if (!(prot.res[k_res].conf[k_conf].atom = (ATOM *) realloc(prot.res[k_res].conf[k_conf].atom, (prot.res[k_res].conf[k_conf].n_atom + 1) * sizeof(ATOM)))) {
            printf("ins_atom(): Fails resizing memory.\n");
            del_prot(&prot);
            return prot;
        }
        prot.res[k_res].conf[k_conf].n_atom++;
        k_atom = prot.res[k_res].conf[k_conf].n_atom-1;
        
        prot.res[k_res].conf[k_conf].atom[k_atom] = atom;
        prot.res[k_res].conf[k_conf].atom[k_atom].on = 1;

        if ( !strlen(prot.res[k_res].conf[k_conf].history) )
            strcpy(prot.res[k_res].conf[k_conf].history, atom.history);
    }
    return prot;
}

int atom2pdbline(char *line, ATOM atom)
{
    
    sprintf(line, "ATOM        %4s%c%3s %c%04d%c%3s%8.3f%8.3f%8.3f %7.3f      %6.3f      %-10s",
            atom.name,
            atom.altLoc,
            atom.resName,
            atom.chainID,
            atom.resSeq,
            atom.iCode,
            atom.confID,
            atom.xyz.x,
            atom.xyz.y,
            atom.xyz.z,
            atom.rad,
            atom.crg,
            atom.history);
    
            return 0;
    
}


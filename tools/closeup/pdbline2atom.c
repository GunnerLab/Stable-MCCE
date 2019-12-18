/*
 * NAME
 *       pdbline2atom - convert a PDB formatted string to atom structure
 *
 * SYNOPSIS
 *       #include <mcce.h>
 *       ...
 *       ATOM pdbline2atom(char *line);
 *
 * DESCRIPTION
 *       The pdbline2atom function convert a string, for example a line in PDB file,
 *       into ATOM structure. And return the converted ATOM structure.
 *
 * SEE ALSO
 *       load_pdb
 *
 * EXAMPLE
 *       #include <stdio.h>
 *       #include <stdlib.h>
 *       #include <string.h>
 *       #include "mcce.h"
 *
 *       int main() {
 *           char line[] = "ATOM     10  N   ASP     2 BK  -11.857  10.111   6.057  1.00  0.00           N  BK"
 *           ATOM atom;
 *
 *           atom = pdbline2atom(line);
 *       }
 *
 * AUTHOR
 *       Yifan Song, 06/17/2003, Edited by Junjun Mao 08/28/2003
 */

#include <stdlib.h>
#include <string.h>
#include "mcce.h"

ATOM pdbline2atom(char *line)
{   ATOM atom;
    char sbuffer[MAXCHAR_LINE];

    memset(&atom, 0, sizeof(ATOM));

    strncpy(sbuffer, (line+6), 5);
    sbuffer[5] = '\0';
    atom.serial = atoi(sbuffer);
    
    strncpy(atom.name, (line+12), 4);
    atom.name[4] = '\0';

    atom.altLoc = line[16];
    
    strncpy(atom.confName,  (line+17), 3);
    strncpy(atom.confName+3,(line+80), 2);
    atom.confName[5] = '\0';

    strncpy(atom.resName, (line+17), 3);
    atom.resName[3] = '\0';

    atom.chainID = line[21];

    strncpy(sbuffer, (line+22), 4);
    sbuffer[4] = '\0';
    atom.resSeq = atoi(sbuffer);

    atom.iCode = line[26];
    if (atom.iCode == ' ') atom.iCode = '_';

    strncpy(atom.confID, (line+27), 3);
    atom.confID[3] = '\0';

    strncpy(sbuffer, (line+30), 8);
    sbuffer[8] = '\0';
    atom.xyz.x = atof(sbuffer);

    strncpy(sbuffer, (line+38), 8); 
    sbuffer[8] = '\0'; 
    atom.xyz.y = atof(sbuffer);
    
    strncpy(sbuffer, (line+46), 8); 
    sbuffer[8] = '\0'; 
    atom.xyz.z = atof(sbuffer);
    
    strncpy(sbuffer, (line+55), 7); 
    sbuffer[7] = '\0';
    atom.rad = atof(sbuffer);

    strncpy(sbuffer, (line+68), 6);
    sbuffer[6] = '\0';
    atom.crg = atof(sbuffer);
    
    strncpy(atom.history, (line+80), 10);
    atom.history[10] = '\0';

    return atom;
}

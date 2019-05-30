#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mcce.h"

int main(int argc, char *argv[]) {
    FILE *fp;
    char line[MAXCHAR_LINE], dist[MAXCHAR_LINE];
    VECTOR origin;
    ATOM atom;
    float threshold;
    
    if (argc < 6) {
        printf("   Usage: closeup pdbfile x y z threshold\n");
        return 0;
    }
    
    fp = fopen(argv[1],"r");
    origin.x = atof(argv[2]);
    origin.y = atof(argv[3]);
    origin.z = atof(argv[4]);
    threshold = atof(argv[5]);

    while ( memset(line,0,MAXCHAR_LINE*sizeof(char)),
    fgets(line, MAXCHAR_LINE, fp)) {
        if (strncmp("ATOM  ", line, 6) && strncmp("HETATM", line, 6)) continue;

        atom = pdbline2atom(line);
        if (dvv(atom.xyz, origin) < threshold) {
            sprintf(dist, "%5.2f",dvv(atom.xyz, origin));
            strncpy(line+6,dist,strlen(dist));
            printf("%s",line);
        }
    }
    fclose(fp);
    return 0;
}


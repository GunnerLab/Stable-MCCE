#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mcce.h"

extern int rm_comment(char *target, char *str);

int check_tpl(char *fname) {
    FILE *fp;                              /* handler of the loaded parameter file */
    char sbuff[MAXCHAR_LINE];              /* string buffer */
    char line[MAXCHAR_LINE];               /* line buffer of a parameter line */
    char key1[LEN_KEY1+1];                 /* key1 of the parameter entry */
    char key2[LEN_KEY2+1];                 /* key2 of the parameter entry */
    char key3[LEN_KEY3+1];                 /* key3 of the parameter entry */
    int i;

    if ((fp=fopen(fname, "r")) == NULL) {
        printf("\ncheck_tpl(): failed oppening file \"%s\".\n", fname);
        return USERERR;
    }
    
    while (fgets(sbuff, MAXCHAR_LINE*sizeof(char), fp)) {
        /* strip off the comment of the line */
        memset(line,0,MAXCHAR_LINE*sizeof(char));
        rm_comment(line, sbuff);
        if (strncmp(line, "ATOMNAME ",9)) {
            for(i=strlen(line)-1;i>=0;i--) {
                if (line[i]==' ') line[i]='\0';
                else if (line[i]=='\t') line[i]='\0';
                else break;
            }
        }
        
        /* if the line is shorter than the total length of 3 keys,
        * then it has no value, we proceed to the next line */
        if (strlen(line) < (LEN_KEY1+LEN_KEY2+LEN_KEY3)) continue;
        for(i=strlen(line)-1;i>=0;i--) {
            if (line[i]=='\t') {
                printf("\n   Tab charater is not recognized by MCCE\n");
                printf("   Please check line \"%s\" in file \"%s\" and replace tab with spaces\n",line,fname);
                return USERERR;
            }
        }
        
        /* else we create 3 keys */
        strncpy(key1,line,LEN_KEY1); key1[LEN_KEY1] = '\0';
        strncpy(key2,line+LEN_KEY1,LEN_KEY2); key2[LEN_KEY2] = '\0';
        strncpy(key3,line+LEN_KEY1+LEN_KEY2,LEN_KEY3); key3[LEN_KEY3] = '\0';
        
        if (param_exist(key1, key2, key3)) {
            printf("\n   In file %s, parameter \"%s\" is already loaded somewhere else.\n",fname,line);
            printf("   Try delete this entry and run MCCE again\n");
            return 0;
        }
    }
    return 0;
}

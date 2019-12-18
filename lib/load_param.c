/*******************************************************************************
 * NAME
 *        load_param - load a parameter file into parameter database and convert 
 *                     the parameter value to appropiate data type.
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        int load_param(char *fname)
 *
 * DESCRIPTION
 *        The load_param() function loads a parameter file to parameter database.
 *        Each entry is a pair of key and value, which can be retrieved later by
 *        the key.
 *
 *        The key is composed of 3 key strings called key1, key2  and key3.   In
 *        general key1 is the record name,  key2 is the conformer name,  and key3
 *        is the atom name.  The leading and ending spaces of 3 keys are stripped
 *        off at the time the parameter entry is saved.
 *
 *        The value has  been either  converted  to the the  appropriate types if 
 *        listed below or left as a string and you can process it later.  Here is
 *        a table of the converted values:
 *
 *        record      data type
 *        CONNECT     CONNECT
 *        CONFLIST    STRINGS
 *        IATOM       int
 *        PROTON      int
 *        ELECTRON    int
 *        PKA         float
 *        RXN         float
 *        CHARGE      float
 *        RADIUS      float
 *        all other   STRING
 *
 *        A db_open()  should be called separately before this function.  This is
 *        because multiple parameter files are loaded,  but db_open() can be only
 *        called once.
 *
 * RETURN VALUE
 *        Return integer 0 on success or -1 on failure.
 *
 * SEE ALSO
 *        db_open, db_close, param_sav, param_get, iatom
 *
 * EXAMPLE
 *        #include <mcce.h>
 *
 *        int main()
 *        {  db_open();
 *           load_param("asp.txt");
 *           load_param("asn.txt");
 *           db_close();
 *           return 0;
 *        }
 *
 * DEPENDENCY
 *        strip_comment, strip_spc, param_sav
 *
 * AUTHOR
 *        Junjun Mao, 06/06/2003
 *******************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mcce.h"


int rm_comment(char *target, char *str)
/* remove comment */
{  char *pt;

   strcpy(target, str);
   if ((pt=strstr(target, "REMARK")) ||
       (pt=strchr(target, '#'))      ||
       (pt=strchr(target, '\n'))       )  *pt = '\0';

   return 0;
}


int load_param(char *fname)
{  FILE *fp;                              /* handler of the loaded parameter file */
   char sbuff[MAXCHAR_LINE];              /* string buffer */
   char line[MAXCHAR_LINE];               /* line buffer of a parameter line */
   char *ptr;                             /* position pointer of the above line */
   char key1[LEN_KEY1+1];                 /* key1 of the parameter entry */
   char key2[LEN_KEY2+1];                 /* key2 of the parameter entry */
   char key3[LEN_KEY3+1];                 /* key3 of the parameter entry */
   int Counter;                           /* atom counter */
   int i;

   if (check_tpl(fname)) {
       printf("FATAL: There is an error in tpl file %s\n", fname);
       return USERERR;
   }
   
   /* open the parameter file, return -1 on failure */
   if ((fp=fopen(fname, "r")) == NULL) {
      printf("load_param(): failed oppening file \"%s\".\n", fname);
      return USERERR;
   }

   /* process each valid line in the parameter file */
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

      /* else we create 3 keys */
      strncpy(key1,line,LEN_KEY1); key1[LEN_KEY1] = '\0';
      strncpy(key2,line+LEN_KEY1,LEN_KEY2); key2[LEN_KEY2] = '\0';
      strncpy(key3,line+LEN_KEY1+LEN_KEY2,LEN_KEY3); key3[LEN_KEY3] = '\0';

      /* Then we convert the value string to appropiate type */
      if (!strcmp(key1, "CONNECT  ")) { /* convert to CONNECT structure */
         CONNECT connect;              /* temporary CONNECT record */
         CONNECTED_ATOM connected;     /* temporary CONNECTED_ATOM record */

         /* get the orbital name from column LEN_KEY1+LEN_KEY2+LEN_KEY3, 10 characters */
         strncpy(sbuff, line+LEN_KEY1+LEN_KEY2+LEN_KEY3, 10); sbuff[10] = '\0';
         strip(connect.orbital, sbuff);

         /* set the line pointer to the first atom */
         if (strlen(line) <= LEN_KEY1+LEN_KEY2+LEN_KEY3 + 10)
            ptr = line+strlen(line);                    /* no atom, step to the end */
         else
            ptr = line+LEN_KEY1+LEN_KEY2+LEN_KEY3 + 10; /* 1st atom, step to the 1st atom */

         /* set atom counter */
         Counter = 0;

         /* start looping over of atoms */
         while (strlen(ptr)) {
            strncpy(sbuff, ptr, 5); sbuff[5] = '\0';  /* cut the first 5 char */
            connected.res_offset = atoi(sbuff);       /* convert to integer */
            if (strstr(sbuff, "LIG")) connected.ligand = 1;/* ligand flag is on */
            else connected.ligand = 0;                     /* ligand flag is off */

            strncpy(connected.name, ptr+5, 4); /* cut the next 4 char as the atom name */
            connected.name[4] = '\0';
            while (strlen(connected.name)<4)   /* If atom name is too short, complete with space */
                strcat(connected.name," ");

            connect.atom[Counter] = connected; /* save this connected atom */
            Counter ++;                        /* Counter steps forward by 1 */
            if (strlen(ptr) > 10) ptr += 10;   /* string pointer steps forward by 1 atom */
            else ptr += strlen(ptr);           /* string pointer steps forward to the end */
         }

         /* store the value */
         connect.n = Counter;
         param_sav(key1, key2, key3, &connect, sizeof(CONNECT));

      }
      else if (!strcmp(key1, "CONFLIST ")) { /* convert to STRINGS */
         STRINGS value;
         strip(sbuff, line+LEN_KEY1+LEN_KEY2+LEN_KEY3);
         value.n = strlen(sbuff)/6+1;  /* 5 char for each conformer name, 1 separation space, # of conformers */

         /* allocate string array from the number of conf names */
         value.strings = (char **) malloc(value.n*sizeof(char *));

         /* loop over all conformer names in the value string */
         for (Counter=0; Counter<value.n; Counter++) {
            /* allocate the string for 5 char in each conf entry */
            value.strings[Counter] = (char *) malloc(6*sizeof(char));
            /* copy the conformer name */
            strncpy(value.strings[Counter], sbuff+Counter*6, 5);
            value.strings[Counter][5] = '\0';
         }
         /* save this STRINGS structure */
         param_sav(key1, key2, key3, &value, sizeof(value));
      }
      /* Edited for transfer hb.tpl to the parameter database --Cai 20170812 */
      else if (!strcmp(key1, "HDONOR   ") ||
               !strcmp(key1, "HACCEPT  ")
              ) { /* convert to STRINGS */
         STRINGS value;
         strip(sbuff, line+LEN_KEY1+LEN_KEY2+LEN_KEY3);
         value.n = strlen(sbuff)/5+1;  /* 4 char for each atom name, 1 separation space, # of donor/acceptor atom (hydrogen) */

         /* allocate string array from the number of atom names */
         value.strings = (char **) malloc(value.n*sizeof(char *));

         /* loop over all donor/acceptor atom names in the value string */
         for (Counter=0; Counter<value.n; Counter++) {
            /* allocate the string for 4 char in each h_donor/h_acceptor entry */
            value.strings[Counter] = (char *) malloc(5*sizeof(char));
            /* copy the donor/acceptor atom name */
            strncpy(value.strings[Counter], sbuff+Counter*5, 4);
            value.strings[Counter][4] = '\0';
         }
         /* save this STRINGS structure */
         param_sav(key1, key2, key3, &value, sizeof(value));
      }

      else if (!strcmp(key1, "PROTON   ") ||
               !strcmp(key1, "IATOM    ") ||
               !strcmp(key1, "NATOM    ") ||
               !strcmp(key1, "ELECTRON ")
              ) {    /* if falls into one of them, convert string to float */
         int value;
         value = atoi(line+LEN_KEY1+LEN_KEY2+LEN_KEY3);
         param_sav(key1, key2, key3, &value, sizeof(int));
      }
      else if (!strcmp(key1, "PKA      ") ||
               !strcmp(key1, "CHARGE   ") ||
               !strcmp(key1, "RADIUS   ") ||
               !strcmp(key1, "RXN      ") ||
               !strcmp(key1, "VDWAMBER ") ||
               !strcmp(key1, "VDW_RAD  ") ||
               !strcmp(key1, "VDW_EPS  ") ||
               !strcmp(key1, "RADCOVAL ") ||
               !strcmp(key1, "EM       ") ||
               !strcmp(key1, "EXTRA    ") ||
               !strcmp(key1, "IPECE_SC ") ||
               !strcmp(key1, "SCALING  ")
              ) {    /* if falls into one of them, convert string to float */
         float value;
         value = atof(line+LEN_KEY1+LEN_KEY2+LEN_KEY3);
         param_sav(key1, key2, key3, &value, sizeof(float));
      }
      else if (!strcmp(key1, "CAL_VDW  ") ||
               !strcmp(key1, "RELAX    ")
              ) {
          int value;
          if (strchr(line+LEN_KEY1+LEN_KEY2+LEN_KEY3, 't') ||
              strchr(line+LEN_KEY1+LEN_KEY2+LEN_KEY3, 'T') ) value = 1;
          else value = 0;
         param_sav(key1, key2, key3, &value, sizeof(int));
      }
      else if (!strcmp(key1, "ROTAMER  ") ||
               !strcmp(key1, "HYDROXYL ")) { /* convert to ROTAMER */
         ROTAMER value;
         strncpy(value.atom1, line+LEN_KEY1+LEN_KEY2+LEN_KEY3, 4); value.atom1[4] = '\0';
         strncpy(value.atom2, line+LEN_KEY1+LEN_KEY2+LEN_KEY3+5, 4); value.atom2[4] = '\0';
         strip(value.affected, line+LEN_KEY1+LEN_KEY2+LEN_KEY3+10);
         /* save this STRINGS structure */
         param_sav(key1, key2, key3, &value, sizeof(value));
      }
      else if (!strcmp(key1, "ROT_SWAP ")) {
         SWAP_RULE value;
         int i_swap;
         char *sbuffer;
         sbuffer = line+LEN_KEY1+LEN_KEY2+LEN_KEY3;
         value.n_swap = 0;
         while (strlen(sbuffer) > 5) {
             value.n_swap++;
             i_swap = value.n_swap - 1;
             strncpy(value.swap_atom1[i_swap], sbuffer,   4);
             value.swap_atom1[i_swap][4] = '\0';
             strncpy(value.swap_atom2[i_swap], sbuffer+5, 4);
             value.swap_atom2[i_swap][4] = '\0';
             
             if (strlen(sbuffer)>=10) sbuffer+=10;
             else sbuffer+=(strlen(sbuffer));
         }
         param_sav(key1, key2, key3, &value, sizeof(value));
      }
      else if (!strcmp(key1, "TORSION  ")) {
          TORS value;
          char  sbuffer[10];
          int i_term;
          strncpy(value.atom1, line+LEN_KEY1+LEN_KEY2+LEN_KEY3,    4); value.atom1[4] = '\0';
          strncpy(value.atom2, line+LEN_KEY1+LEN_KEY2+LEN_KEY3+5,  4); value.atom2[4] = '\0';
          strncpy(value.atom3, line+LEN_KEY1+LEN_KEY2+LEN_KEY3+10, 4); value.atom3[4] = '\0';
          strncpy(sbuffer,     line+LEN_KEY1+LEN_KEY2+LEN_KEY3+15, 4);     sbuffer[4] = '\0';
          if (strchr(sbuffer, 't') ||
              strchr(sbuffer, 'T') ) value.opt_hyd = 1;
          else value.opt_hyd = 0;
          for(i=strlen(line)-1;i>=0;i--) {
              if (line[i]==' ') line[i]='\0';
              else if (line[i]=='\t') line[i]='\0';
              else break;
          }
          value.n_term = ((float)(strlen(line)-(LEN_KEY1+LEN_KEY2+LEN_KEY3+20)))/30.+1;
          if (value.n_term < 1) value.n_term = 1;
          for (i_term=0;i_term<value.n_term; i_term++) {
              strncpy(sbuffer,     line+LEN_KEY1+LEN_KEY2+LEN_KEY3+20+i_term*30, 9);     sbuffer[9] = '\0'; value.V2[i_term]    = atof(sbuffer);
              strncpy(sbuffer,     line+LEN_KEY1+LEN_KEY2+LEN_KEY3+30+i_term*30, 9);     sbuffer[9] = '\0'; value.n_fold[i_term]  = atof(sbuffer);
              strncpy(sbuffer,     line+LEN_KEY1+LEN_KEY2+LEN_KEY3+40+i_term*30, 9);     sbuffer[9] = '\0'; value.gamma[i_term] = env.d2r * atof(sbuffer);
          }
         /* save */
         param_sav(key1, key2, key3, &value, sizeof(TORS));
      }
      else {         /* others, convert STRING structure */
         strcpy(sbuff, line+LEN_KEY1+LEN_KEY2+LEN_KEY3);
         param_sav(key1, key2, key3, sbuff, sizeof(sbuff));
         /* save with the string terminating character */
      }
   }


   fclose(fp);
   return 0;
}

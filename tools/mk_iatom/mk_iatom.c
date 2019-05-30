/*
NAME
        mk_iatom - make indexed atom list

SYNOPSIS
        #include <stdio.h>
        #include <stdlib.h>
        #include <string.h>
        #define MAX 256
        #define MAX_T_CONF 16

        int mk_iatom(char *par_fname, FILE *par_iatom)

DESCRIPTION
        mk_iatom - is passed with a parameter file, "*.tpl", and a file stream on which
        all the necessary information will be written for other processes.
        "mk_iatom" first opens and scans the passed parameter file and copies necessary portions
        to the string. As it scans through the file it remembers conformer names and counts
        the numbers of atoms for each conformer as well in order to make iatom with proper
        indexing and to list conformer types.

RETURN VALUE
        returns "0" if everything goes smoothly, "-1" otherwise.

SEE ALSO
        iatom, param_get

EXAMPLE
        #include <stdio.h>
        #include <string.h>
        #include <mcce.h>

        int main(int argc, char *argv[]){
           FILE *fpOut;
           char new_name[256] = "";
           strcpy(new_name, argv[1]);
           strcat(new_name, ".iatom");
           fpOut = fopen(new_name, "w");

           if(argc < 2){
              printf("mk_iatom param_file\n");
              return 0;
           }

           mk_iatom(argv[1], fpOut);
           fclose(fpOut);
           return 0;
        }
AUTHOR
        Jinrang Kim, 06/20/2003
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX 256
#define MAX_T_CONF 16

int mk_iatom(char *par_fname, FILE *par_iatom);

int main(int argc, char *argv[]){
   FILE *fpOut;
   FILE *temp_fp;
   char new_name[256] = "";
   char line[256];
   
   strcpy(new_name, argv[2]);
   temp_fp = tmpfile();

   if(argc < 3){
      printf("mk_iatom param_file out_file\n");
      return 0;
   }

   mk_iatom(argv[1], temp_fp);
   
   rewind(temp_fp);
   fpOut = fopen(new_name, "w");
   while ( fgets(line, 256, temp_fp) )
       fprintf(fpOut,"%s",line);
   fclose(temp_fp);
   fclose(fpOut);
   return 0;
}

/*This function is passed a parameter file, deletes existing "IATOM" list,
creates new "CONFLIST" of the residue, and recreates "IATOM" lists with index
numbers starting from 0 for each conformer. */

int mk_iatom(char *par_fname, FILE *par_iatom){
    FILE *iptr, *tmpstr;        /*pointers to hold passed file, output file*/
    char key1[10] = "CONFLIST ";
    char key2[10] = "IATOM    ";
    char key3[10] = "CONNECT  ";
    char key4[10] = "NATOM    ";
    char key5[10] = "ATOMNAME ";
    char atom_index[12];
    char confname[6];
    char atomname[5];
    char line[MAX+1];                  /*using a char array to hold lines of parameter file*/
    char chk_name[6];                /*using a char array to check entries and record conformers*/
    char buffer[MAX_T_CONF*6];         /*a char array of buffer is used to hold all the conformer names for 
comparison*/
    int counters[MAX_T_CONF];
    int count = 0;
    int i;
    
    /*get the filename and open the file for reading*/
    if((iptr = fopen(par_fname, "r")) == NULL){
        printf("mk_iatom(): Cannot open file to read \"%s\"\n", par_fname);
        return -1;
    }

    /*opens a temporary file*/
    /*if end of the file not reached, reads a line in and checks if it starts with "CONFLIST" or "IATOM*/
    /*if not, copies it to the tempfile for further processing*/
    /*if a line starts with "CONNECT", keeps the conformer name and compares it with ones in the buffer string*/
    /*if it's a new, then concatanates it to the buffer string and increasess the count by one*/
    /*closes the file and appends "END" to the end of the tempfile*/
    tmpstr = tmpfile();
    buffer[0] = '\0';
    while(fgets(line, MAX, iptr)){
       if(((strncmp(line, key1, 9)) != 0) && ((strncmp(line, key2, 9)) != 0) && ((strncmp(line, key4, 9)) != 0) && ((strncmp(line, key5, 9)) != 0)){
          fprintf(tmpstr, line);
          if(!strncmp(line, key3, 9)){
             strncpy(chk_name, line+9, 5); 
             chk_name[5] = '\0';
             if(!strstr(buffer, chk_name)){   /*compares new conformers against existing names*/
                strcat(buffer, chk_name);         /*if it's new, attach it at the end of the buffer array*/
                count++;
             }
          }
       }
    }
    fclose(iptr);
    fprintf(tmpstr, "END");

    /*reads the buffer string and writes out "CONFLIST"*/
    chk_name[3] = '\0'; 
    fprintf(par_iatom, "%-9s%3s        ", key1, chk_name);
    /*checks for backbone conformer, creates a dummy one if there is none*/
    strcat(chk_name, "BK"); chk_name[5] = '\0'; 
    if (strstr(buffer, chk_name)){
        for (i=0; i<count; i++) {
           strncpy(chk_name, buffer+i*5, 5); 
           chk_name[5] = '\0';
           fprintf(par_iatom, "%-6s", chk_name);
        }
    }   
    else{
       fprintf(par_iatom, "%-6s", chk_name); /*the dummy backbone conformer*/
       for(i=0; i<count; i++){
           strncpy(chk_name, buffer+i*5, 5);
           chk_name[5] = '\0';
           fprintf(par_iatom, "%-6s", chk_name);
       }
    }
    
    fprintf(par_iatom, "\n");

    /*set counters for indices to be 0*/
    for(i=0; i<count; i++){
       counters[i] = 0;
    }
    
    /*loop over "CONNECT" and count the number of atom of each conformer*/
    rewind(tmpstr);
    while(fgets(line, sizeof(line), tmpstr)) {
       if (!strncmp("END", line, 3)) break;
       if (!strncmp(line, key3, 9)) {
          strncpy(chk_name, line+9, 5); chk_name[5] = '\0';
          for (i=0; i<count; i++)
          if (!strncmp(buffer+i*5, chk_name, 5)) break;
          counters[i]++;
       }
    }
    fprintf(par_iatom, "\n");
    chk_name[3] = '\0';
    /*checks for backbone conformer, creates a dummy one with "NATOM = 0" if there is none*/
    strcat(chk_name, "BK"); chk_name[5] = '\0';
    if (strstr(buffer, chk_name)){
        for(i=0; i<count; i++){
           strncpy(chk_name, buffer+i*5, 5);
           chk_name[5] = '\0';
           fprintf(par_iatom, "%-9s%-11s%d\n", key4, chk_name, counters[i]);
        }
    }
    else{
        fprintf(par_iatom, "%-9s%-11s%d\n", key4, chk_name, 0); /*the dummy backbone conformer with zero atom*/
        for(i=0; i<count; i++){
            strncpy(chk_name, buffer+i*5, 5);
            chk_name[5] = '\0';
            fprintf(par_iatom, "%-9s%-11s%d\n", key4, chk_name, counters[i]);
        }
    }
    
    fprintf(par_iatom, "\n");
    
    /*set the counters to be "0" again for atom indices*/
    for(i=0; i<count; i++){
       counters[i] = 0;
    }
    
    /*loop over "CONNECT" and print out line by line with atom index number*/
    rewind(tmpstr);
    while(fgets(line, sizeof(line), tmpstr)) {
       if (!strncmp("END", line, 3)) break;
       if (!strncmp(line, key3, 9)) {
          strncpy(chk_name, line+9, 5); chk_name[5] = '\0';
          for (i=0; i<count; i++){
              if (!strncmp(buffer+i*5, chk_name, 5)) break;
          }
          strncpy(atom_index, line + 9, 11); 
          atom_index[11] = '\0';
          fprintf(par_iatom, "IATOM    %-11s%d\n", atom_index, counters[i]);
          counters[i]++;
       }
    }
    fprintf(par_iatom, "\n");

    for(i=0; i<count; i++){
       counters[i] = 0;
    }
    rewind(tmpstr);
    while(fgets(line, sizeof(line), tmpstr)) {
       if (!strncmp("END", line, 3)) break;
       if (!strncmp(line, key3, 9)) {
          strncpy(chk_name, line+9, 5); chk_name[5] = '\0';
          for (i=0; i<count; i++){
              if (!strncmp(buffer+i*5, chk_name, 5)) break;
          }
          strncpy(atom_index, line + 9, 11); 
          atom_index[11] = '\0';
          strncpy(confname, atom_index, 5);
          confname[5] = '\0';
          strncpy(atomname, atom_index+6, 4);
          atomname[4] = '\0';
          fprintf(par_iatom, "ATOMNAME %s %4d %s\n", confname, counters[i], atomname);
          counters[i]++;
       }
    }

    /*rewinds tempfile and copies lines to the new file*/
    /*fprintf(par_iatom, "\n");*/
    rewind(tmpstr);
    while(fgets(line, sizeof(line), tmpstr)) {
       if (!strncmp("END", line, 3)){ 
           break;
       }
       fprintf(par_iatom, line);
    }

    fclose(tmpstr);

    return 0;
}


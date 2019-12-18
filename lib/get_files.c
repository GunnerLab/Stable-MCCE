/*
 * NAME
 *    	 get_files - list of file names in a directory
 *    	  
 * SYNOPSIS
 *    	 #include <mcce.h>
 *    	 
 *    	 STRINGS get_files(char *dir);
 *    	 
 * DESCRIPTION
 *       The get_files() function returns a STRINGS  object with total number of
 *       files, and a list of file names.  
 *
 *       The STRINGS structure is declared as follows:
 *             typedef struct {
 *                int  n;     	           - total number of strings 
 *                char **strings;          - the list of strings
 *             } STRINGS;
 *
 *       You should be aware that the field 'strings' is dynamically allocated, 
 *       and it is your responsibily to deallocate if necessary.  The function 
 *       free_strings() is designed to free a STRING object, and reset n to 0.
 *
 *       This version of get_files() has been tested on Mac OS X and Linux.    
 *       The version for irix is get_files_irix(), which is different in:
 *         1) it includes different header files.
 *         2) it determines whether a filename is a directory based on *stat, 
 *            rather than a *dirent structure.
 *
 * RETURN VALUE
 *    	 The get_files() function returns a STRINGS variable.
 * 
 * SEE ALSO
 *       free_strings
 *
 * EXAMPLE
 *    	 #include <stdio.h>
 *    	 #include <mcce.h>
 *       ...
 *    	 char    *mydir = "funalgo"; 
 *       int     i;
 *    	 STRINGS files;
 *  
 *    	 files = get_files(mydir);
 *       for (i = 0; i < files.n; i++) {
 *          printf("%d, file : %s\n", i, files.strings[i]);
 *       }
 *       free_strings(files);
 * 
 * AUTHOR
 *    	 Yanjun Wang, 06/09/2003
 */


#include <stdio.h>
#include <stdlib.h>
#include <sys/dir.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include "mcce.h"


STRINGS get_files(char *dir_name)
{
   STRINGS file_names;          /* the value to be returned  */
   DIR     *dir;      	       /* DIR returned by opendir() */
   struct  dirent *dir_entry;  /* returned by readdir(dir)  */
   FILE    *temp;     	       /* a stream to temporarily hold filenames */
   char    s[256];   	       /* a buffer to hold one filename */
   int     i;

   temp = tmpfile();
   file_names.n  = 0;

   /*
    * read the directory, excluding current directory (.), parent directory (..),
    * and sub-directories.
    *
    * if failed reading the directory:
    *    return to the caller immediately
    * else
    *    save number of files in 'file_names.n'
    *    save all filenames to   'temp'
    *    malloc space to hold     file_names.strings
    */
   if ((dir = opendir(dir_name)) == NULL) {
      printf("   FATAL: get_files(): \"failed reading directory %s.\"\n", dir_name);
      file_names.strings = NULL;
   }
   else {
      while ((dir_entry  = readdir(dir)) != NULL) {
         if ( strcmp(dir_entry->d_name, ".")  != 0 &&
              strcmp(dir_entry->d_name, "..") != 0 &&
              dir_entry->d_type != DT_DIR) {
            file_names.n ++;
            fprintf(temp, "%s\n", dir_entry->d_name);
         }
      }
      rewind(temp);

      /* request physical memory to put filenames */
      file_names.strings = (char **) malloc (file_names.n * sizeof (char *));

      /* copy filenames from the stream to the STRINGS structure */
      for (i = 0; i < file_names.n; i++) {
         fgets(s, sizeof(s), temp);
         *strchr(s, '\n') = '\0';
         file_names.strings[i] = (char *) malloc ((1 + strlen(s)) * sizeof (char));
         strcpy(file_names.strings[i], s);
      }

      fclose(temp);		 /* close the stream */
      closedir(dir);		 /* close the directory */
   }

   return file_names;
}

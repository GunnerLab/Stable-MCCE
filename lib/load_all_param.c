/*******************************************************************************
 * NAME
 *        load_all_param - load all files to parameter database from a directory
 *
 * SYNOPSIS
 *        #include <mcce.h>
 *
 *        int load_all_param(char *dirname)
 *
 * DESCRIPTION
 *        The  load_param() function loads all parameter files  in a directory to
 *        the parameter database.  It is built up on  two functions:  get_files()
 *        and load_param().  A separate  call of  function db_open()  to open the
 *        parameter database is required before using this function.
 *
 * RETURN VALUE
 *        Return integer 0 on success or -1 on failure.
 *
 * SEE ALSO
 *        load_param, db_open, db_close
 *
 * EXAMPLE
 *        #include <mcce.h>
 *
 *        int main()
 *        {  db_open();
 *           load_all_param("param_dir");
 *           db_close();
 *           return 0;
 *        }
 *
 * DEPENDENCY
 *        get_files, load_param
 *
 * AUTHOR
 *        Junjun Mao, 06/11/2003
 *******************************************************************************/

#include <stdio.h>
#include <string.h>
#include "mcce.h"

int load_all_param(char *dirname)
{  int i;
   STRINGS files;
   char fullname[MAXCHAR_LINE];

   /* Open and read the parameter file names in a directory */
   files = get_files(dirname);
   if(files.n == 0) {
      printf("   FATAL: load_all_param(): \"Invalid or empty directory.\"\n");
      return USERERR;
   }

   /* load each file in this directory */
   for (i=0; i<files.n; i++) {
      /* create the full file name with path */
      //if (!strstr(files.strings[i], ".tpl")) continue;
      if (strlen(files.strings[i])<=4) continue;
      if (strcmp(files.strings[i]+(strlen(files.strings[i])-4),".tpl")) continue;
      sprintf(fullname, "%s/%s", dirname, files.strings[i]);

      /* load the file */
      /* printf("   loading \"%s\"\n", fullname);*/
      if (load_param(fullname)) {
         printf("   FATAL: load_all_param(): Failed loading file \"%s\".\n", fullname);
         return USERERR;
      }
   }

   free_strings(&files);
   return 0;
}

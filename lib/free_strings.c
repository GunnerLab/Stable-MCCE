/*
 * NAME
 *   free_strings - deallocate the dynamic memory of a STRINGS object
 *
 * SYNOPSIS
 *   #include <mcce.h>
 *
 *   void free_strings(STRINGS *s);
 *
 * DESCRIPTION
 *   The free_strings() function frees memory dynamically allocated to the
 *        'strings' field of a STRINGS object.
 *
 * SEE ALSO
 *   get_files
 *
 * EXAMPLE
 *   #include <mcce.h>
 *   ...
 *        char    *dir = "funalgo";
 *   STRINGS files;
 *
 *   files = get_files(dir);
 *   ...
 *   free_strings(&files);
 *
 * AUTHOR
 *   Yanjun Wang, 06/09/2003
 */


#include <stdlib.h>
#include "mcce.h"

void free_strings(STRINGS *s)
{  int i;

   for (i = 0; i < s->n; i++) free(s->strings[i]);
   if (s->n > 0) {
   free(s -> strings);
   s -> n = 0;
   }
   return;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include "mcce.h"
#define KEY_LEN 21
#define VALUE_LEN 360

/* This pointer is a static variable that points to parameter database
 * and pairwise.
 */
void *param_root = NULL;
void *pw_root = NULL;

typedef struct
{
    char key[KEY_LEN];
    char value[VALUE_LEN];
    int size;
} RECORD;


int dbcmp (const void *r1, const void *r2)
{
    return strcmp (((RECORD *) r1)->key, ((RECORD *) r2)->key);
}

int param_sav (char *key1, char *key2, char *key3, void *value, int s)
{
    RECORD *data, *val;
    char key[MAXCHAR_LINE];
    char sbuff[MAXCHAR_LINE];

    /* convert 3 key strings to one key, leading and ending spaces stripped */
    strip (key, key1);
    strip (sbuff, key2);
    strcat (key, sbuff);
    strip (sbuff, key3);
    strcat (key, sbuff);

    if (strlen (key) > KEY_LEN)
    {
        printf ("\n   param_sav(): key %s can not exceed %d characters\n", key,
                KEY_LEN);
        return EXIT_FAILURE;
    }
    if (s > VALUE_LEN)
    {
        printf
        ("\n   param_sav(): value (len=%d) can not exceed %d characters\n", s,
         VALUE_LEN);
        return EXIT_FAILURE;
    }

    /* construct key, value pair */
    data = malloc (sizeof (RECORD));
    strcpy (data->key, key);
    memcpy (data->value, value, s);
    data->value[s]='\0';
    data->size = s;
    val = tsearch (data, &param_root, dbcmp);
    return EXIT_SUCCESS;
}


int param_get(char *key1, char *key2, char *key3, void *value)
{
    RECORD *data;
    void *node;
    char key[MAXCHAR_LINE];
    char sbuff[MAXCHAR_LINE];

    /* convert 3 key strings to one key, leading and ending spaces stripped */
    strip (key, key1);
    strip (sbuff, key2);
    strcat (key, sbuff);
    strip (sbuff, key3);
    strcat (key, sbuff);
    if (strlen (key) > KEY_LEN)
    {
        printf ("\n   param_get(): key %s can not exceed %d characters\n", key, KEY_LEN);
        return -1; /* failure */
    }
    data = malloc (sizeof (RECORD));
    strcpy (data->key, key);
    node = (void *) tfind ((void *) data, &param_root, dbcmp);
    free (data);
    if (!node)
    {
       return -1; /* failure */
    }
    else /* success */
    {
       memcpy(value, (*(RECORD **)node)->value, (*(RECORD **)node)->size);
       return 0;
    }
}

int iatom(char *conf_name, char *atom_name)
/* iatom is a special case of param_get.
 * While param_get returns a general pointer, iatom decipers the pointer to an integer.
 */
{  int val, i_atom;
   if (param_get("IATOM", conf_name, atom_name, &val) < 0) {
      return -1; /* failure */
  }
  else {
      return val; /* an index number based from 0 */
  }
}


int param_exist(char *key1, char *key2, char *key3)
{  char val[VALUE_LEN];

   if (param_get(key1, key2, key3, &val) < 0) {
      return 0; /* failure */
   }
   else {
      return 1; /* success, yes exists */
   }
}

/* release database memory */
void free_param() {
   tdestroy(param_root, free);
   return;
}





#include <string.h>
#include <stdlib.h>
#include "mcce.h"

int ins_conf(RES *res, int ins, int n_atom)
{
    if (ins > res->n_conf)
    {
        printf("ins_conf(): off range insertion.\n");
        return USERERR;
    }

    /*--- increase the conformer list size by 1 */
    if (!(res->conf = (CONF *) realloc(res->conf, (res->n_conf + 1) * sizeof(CONF))))
    {
        printf("ins_conf(): Fails resizing memory\n");
        return USERERR;
    }
    res->n_conf++;

    /* move contents after position "ins" forward by 1 */
    memmove(res->conf+ins+1, res->conf+ins, (res->n_conf - ins - 1) * sizeof(CONF));

    /* initialize the inserted conformer */
    memset(&res->conf[ins], 0, sizeof(CONF));
    res->conf[ins].n_atom = n_atom;
    if (n_atom) {
        if (!(res->conf[ins].atom = (ATOM *) malloc(n_atom*sizeof(ATOM))))
        {
            printf("ins_conf(): Fails allocating memory for atoms.\n");
            return USERERR;
        }
        memset( res->conf[ins].atom, 0, (n_atom * sizeof(ATOM)) );
    }

    return ins;
}

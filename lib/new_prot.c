#include <string.h>
#include "mcce.h"

PROT new_prot()
{  PROT prot;
   memset(&prot, 0, sizeof(PROT));
   return prot;
}

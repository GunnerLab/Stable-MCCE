/* call to the C function realloc is necessary as memory is adjusted at
   run time 			S. Sridharan. Sep 95 */
#include <stdlib.h>

#ifdef IRIX
#define memalloc memalloc_
#elif defined (LINUX)
#define memalloc memalloc_
#elif defined (CRAY)
#define memalloc MEMALLOC
#elif defined (PC)
#define memalloc (_stdcall MEMALLOC)
#endif

void *memalloc( void **ptr,int *entry_size, int *new_size )
{       void *newptr;
/*      printf("\n Entry_size: %d, new_size: %d\n",*entry_size, *new_size);*/
	if(!*new_size)
    {
       if(*ptr) free(*ptr);
       return NULL;
	}

	if(!*ptr) {
                  newptr=calloc((size_t)(*new_size),(size_t)(*entry_size));
                  }
	else      {
            	// jmao realloc gives segmentation fault
                //newptr=realloc(*ptr, (size_t)((*new_size)*(*entry_size)));
	             free(*ptr);
	             newptr=calloc((size_t)(*new_size),(size_t)(*entry_size));
                 }





	if (newptr == 0 && *new_size != 0)
        {
#ifdef PC
          perror("memalloc");
#else
          perror("memalloc");
#endif
	  exit(EXIT_FAILURE);
	}

	return newptr;
}

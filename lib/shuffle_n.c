#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void shuffle_n(int *array, int n)
{
   int i;
   int x;
   int temp;

   for (i=0; i<n; i++) array[i] = i;
   for (i=n-1; i>1; i--) {           /* n-th thru 2nd element */
      x = rand() / (RAND_MAX/(i+1)+1);      /* 0 <= num <= i, C faq, 13.16 */
      temp = array[i];                  /* exchange array[i] and array[x] */
      array[i] = array[x];
      array[x] = temp;
   }
   return;
}

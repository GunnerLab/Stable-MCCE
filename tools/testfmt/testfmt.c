#include <stdio.h>

int main(int argc, char *argv[])
{  FILE *fp;
   struct {
      float x;
      float y;
      float z;
      float c;
      float r;
   } data;

   if (argc < 2) {
      printf("readout unpdb\n");
      return 0;
   }

   fp = fopen(argv[1], "r");

   while (fread(&data, sizeof(data), 1, fp)) {
      printf("%8.3f%8.3f%8.3f%8.3f%8.3f\n", data.x, data.y, data.z, data.c, data.r);
   }

   fclose(fp);

   return 0;
}

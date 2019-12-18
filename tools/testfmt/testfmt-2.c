#include <stdio.h>

int main(int argc, char *argv[])
{  FILE *fp;
   struct {
      int header;
      int header2;
      float x;
      float y;
      float z;
      float c;
      float r;
      int trailer;
      int trailer2;
   } data;

   if (argc < 2) {
      printf("readout unpdb\n");
      return 0;
   }

   fp = fopen(argv[1], "r");

   while (fread(&data, sizeof(data), 1, fp)) {
      printf("%8d %8d %8.3f%8.3f%8.3f%8.3f%8.3f %8d %8d\n", data.header,data.header2,data.x, data.y, data.z, data.c, data.r, data.trailer,data.trailer2);
   }

   fclose(fp);

   return 0;
}

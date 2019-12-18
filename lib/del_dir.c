#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

int del_dir(char *dir_name)
{  struct dirent *dentry;
   DIR *dp;
   
   if (!(dp = opendir(dir_name))) return -1;

   while ((dentry = readdir(dp))) {
      /*printf("removing ==%s==\n", dentry->d_name);*/
      if (!strcmp(dentry->d_name, ".")) continue;
      if (!strcmp(dentry->d_name, "..")) continue;
      remove(dentry->d_name);
   }
   rmdir(dir_name);

   return 0;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mcce.h"
int write_ms_gold(FILE *stream, PROT prot)
{  int i, j, count,flag_write;
   char res_name[3];

   //resName that will not output in ms_gold file
   char * not_like[14]={"BKB", "CTR", "NTR", "MEM", "PHE", "LEU", "ILE", "VAL", "ALA", "PRO", "GLY", "SEC", "CYS", "MET"};
   count=sizeof(not_like)/sizeof(not_like[0]); 	
   for (i=0; i<prot.n_res; i++) {
	strncpy(res_name,prot.res[i].resName,3);
	flag_write=1;
	for (j=0; j < count; j++) {
//		printf("res_name: %s, not_like: %s\n", res_name, not_like[j]);
// if (!( (res_name=="BKB") || (res_name=="CTR") || (res_name=="NTR") || (res_name=="MEM") || (res_name=="PHE") || (res_name=="LEU") || (res_name=="ILE") || (res_name=="VAL") || (res_name=="ALA") || (res_name=="PRO") || (res_name=="GLY") || (res_name=="SEC") || (res_name=="CYS") || (res_name=="MET") ))
		if ((strcmp(res_name, not_like[j])==0))
			{
				flag_write=0;break;
			}
   	}
	if (flag_write){
		fprintf(stream, "%3s%c%04d\n",
                                prot.res[i].resName,
                                prot.res[i].chainID,
                                prot.res[i].resSeq);
	}
   }
   return 0;
}



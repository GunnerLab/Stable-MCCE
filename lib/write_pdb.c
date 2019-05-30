#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mcce.h"

int write_pdb(FILE *stream, PROT prot)
{  int i, j, k, iConf, c;
   c = 0;
   for (i=0; i<prot.n_res; i++) {
      iConf =0;
      for (j=0; j<prot.res[i].n_conf; j++) {
         for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
            if (!(prot.res[i].conf[j].atom[k].on)) continue;
            if (c<99999) c++;
            fprintf(stream, "ATOM  %5d %4s%c%3s %c%04d%c%03d%8.3f%8.3f%8.3f %7.3f      %6.3f      %-11s\n",
                            c, prot.res[i].conf[j].atom[k].name,
                            prot.res[i].conf[j].altLoc,
                            prot.res[i].resName,
                            prot.res[i].chainID,
                            prot.res[i].resSeq,
                            prot.res[i].iCode,
                            iConf,
                            prot.res[i].conf[j].atom[k].xyz.x,
                            prot.res[i].conf[j].atom[k].xyz.y,
                            prot.res[i].conf[j].atom[k].xyz.z,
                            prot.res[i].conf[j].atom[k].rad,
                            prot.res[i].conf[j].atom[k].crg,
                            prot.res[i].conf[j].history);
         }
         iConf++;
      }
   }
   return 0;
}

int write_full_header(FILE *stream, int i, int j, int *c, int *iConf, PROT prot) {
	int k;
	float x,y,z;//use temporary float variables since the ATOM structure has xyz as "DOUBLE" which are not necessary
	//write out all the atoms
        for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
               	 if (!(prot.res[i].conf[j].atom[k].on)) {
			continue;
		 }
                 if ((*c)<99999) (*c)++;
                 fwrite(c,sizeof(int),1,stream);
                 fwrite(prot.res[i].conf[j].atom[k].name, sizeof(char),4,stream);
                 fwrite(&prot.res[i].conf[j].altLoc,sizeof(char),1,stream);
                 fwrite(prot.res[i].resName,sizeof(char),3,stream);
                 fwrite(&prot.res[i].chainID,sizeof(char), 1, stream);
                 fwrite(&prot.res[i].resSeq,sizeof(int),1,stream);
                 fwrite(&prot.res[i].iCode,sizeof(char),1,stream);
                 fwrite(iConf,sizeof(int),1,stream);
		 x = prot.res[i].conf[j].atom[k].xyz.x; y = prot.res[i].conf[j].atom[k].xyz.y; z = prot.res[i].conf[j].atom[k].xyz.z;
                 fwrite(&x,sizeof(float),1,stream); fwrite(&y,sizeof(float),1,stream); fwrite(&z,sizeof(float),1,stream);
                 fwrite(&prot.res[i].conf[j].atom[k].rad,sizeof(float),1,stream);
                 fwrite(&prot.res[i].conf[j].atom[k].crg,sizeof(float),1,stream);
                 fwrite(prot.res[i].conf[j].history,sizeof(char),11,stream);
	}
	(*iConf)++;
	return 0;
}

int write_coordinates(FILE* stream, int i, int j, int *iConf, PROT prot) {
	int k;
	float x,y,z;
	//for all conformers of the same residue, the properties that change from one conformer to another are:
	//-xyz coordinates
	//-iConf index value; this value is currently only supported to be in range 0-999 (which is a problem). It needs to be a 0-10000000!
        for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
                 if (!(prot.res[i].conf[j].atom[k].on)) continue;
		 fwrite(iConf,sizeof(int),1,stream);
		 x = prot.res[i].conf[j].atom[k].xyz.x; y = prot.res[i].conf[j].atom[k].xyz.y; z = prot.res[i].conf[j].atom[k].xyz.z;
                 fwrite(&x,sizeof(float),1,stream); fwrite(&y,sizeof(float),1,stream); fwrite(&z,sizeof(float),1,stream);
		 fwrite(prot.res[i].conf[j].history,sizeof(char),11,stream);
        }
	(*iConf)++;
	return 0;
}

int number_atoms(int i, int j, PROT prot) {
	//compute the total number of atoms of the backbone to be written out
        int nb_atoms = 0;
	int k;
        for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
                if (!(prot.res[i].conf[j].atom[k].on)) continue;
	        nb_atoms++;
        }
	return nb_atoms;
}

int write_pdb_binary(FILE *stream, PROT prot) {
	int i, j, iConf, c, nb_atoms;
   	c = 0;
	nb_atoms = 0;

	//write out number of residues to read back in as the write int value of the binary pdb file
	fwrite(&prot.n_res,sizeof(int),1,stream);
	printf("Nb of residues written:%i\n",prot.n_res);

	//for every residue, 1st write out the total number of conformers
	//for conformers other than the backbone, write out the total number of atoms that each conformer has.
	for (i=0; i<prot.n_res; i++) {
		iConf = 0;
		fwrite(&prot.res[i].n_conf,sizeof(int),1,stream);
		//fwrite(&prot.res[i].original_index,sizeof(int),1,stream);
		//printf("Nb of conformers res[%i]:%i\n",i,prot.res[i].n_conf);
		//includes 1 back-bone and 1 conformer
		//in this case, we have to write out full headers for both conformers.
		//conformer index '0' is backbone, conformer index '1' is 1st and only conformer
		if(prot.res[i].n_conf<=2) {
        		//get the total number of atoms of the backbone to be written out
        		nb_atoms = number_atoms(i,0,prot);
			//printf("Nb of atoms res[%i] conf[%i]: %i\n",i,0,nb_atoms);
			fwrite(&nb_atoms,sizeof(int),1,stream);
			write_full_header(stream,i,0,&c,&iConf,prot);
			if(prot.res[i].n_conf==2) {
				nb_atoms = number_atoms(i,1,prot);
				//printf("Nb of atoms res[%i] conf[%i]: %i\n",i,1,nb_atoms);
				fwrite(&nb_atoms,sizeof(int),1,stream);
				write_full_header(stream,i,1,&c,&iConf,prot);
			}
		}
		//includes 1 back-bone, 1st conformer and other conformers with the same number of atoms as the 1st conformer (non back-bone one)
		else {
			nb_atoms = number_atoms(i,0,prot);
			fwrite(&nb_atoms,sizeof(int),1,stream);
			write_full_header(stream,i,0,&c,&iConf,prot);
			//printf("Nb of atoms res[%i] conf[%i]: %i\n",i,0,nb_atoms);

                        nb_atoms = number_atoms(i,1,prot);
                        fwrite(&nb_atoms,sizeof(int),1,stream);
                        write_full_header(stream,i,1,&c,&iConf,prot);
			//printf("Nb of atoms res[%i] conf[%i]: %i N_ATOM:%i\n",i,1,nb_atoms,prot.res[i].conf[1].n_atom);
			//all following conformers, just print 3 coordinates of each atom - retain the same nomenclature as the 1st conformer (index '1').
			for(j=2; j<prot.res[i].n_conf; j++) {
				//printf("Res[%i] conf[%i]_N_ATOM:%i\n",i,j,prot.res[i].conf[j].n_atom);
				write_coordinates(stream,i,j,&iConf,prot);
			}
		}
	}
	return 0;
}

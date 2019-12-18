#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

PROT load_pdb_no_param(FILE *fp)
{   PROT  prot;
    ATOM  atom;
    char  line[MAXCHAR_LINE];   /* line buffer */
    char  confName[6];
    int   k_res;
    int   k_conf;
    int   k_atom;

    memset(&prot, 0, sizeof(PROT));

    while ( memset(line,0,MAXCHAR_LINE*sizeof(char)),
    fgets(line, MAXCHAR_LINE, fp)) {
        if (strncmp("ATOM  ", line, 6) && strncmp("HETATM", line, 6)) continue;
        if (line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0'; /* this is to make sure conformer history not end with a new line */
        atom = pdbline2atom(line);

        /* search for the residue: each unique combination of
        residue name, chain ID, residue number and insertion code
        defines one residue. */
        for (k_res = prot.n_res - 1; k_res >= 0; k_res--) {
            if (!strcmp(atom.resName, prot.res[k_res].resName) &&
                atom.chainID == prot.res[k_res].chainID &&
                atom.resSeq  == prot.res[k_res].resSeq  &&
                atom.iCode   == prot.res[k_res].iCode  )
            {
                break;
            }
        }
        
        /* If couldn't find the residue, add a new one */
        if (k_res == -1) {
            k_res = ins_res(&prot, prot.n_res);
            strcpy(prot.res[k_res].resName, atom.resName);
            prot.res[k_res].chainID = atom.chainID;
            prot.res[k_res].resSeq  = atom.resSeq;
            prot.res[k_res].iCode   = atom.iCode;

            /* Insert a backbone conformer */
            strcpy(confName, atom.resName);
            strcat(confName, "BK");
            ins_conf(&prot.res[k_res], 0, 0);
            strcpy(prot.res[k_res].conf[0].confName, confName);
            strcpy(prot.res[k_res].conf[0].confID, "0000");
            if (!strncmp(atom.history,"BK",2)) strcpy(prot.res[k_res].conf[0].history, atom.history);
            else strcpy(prot.res[k_res].conf[0].history,"BK________");
            prot.res[k_res].conf[0].altLoc = ' ';
        }

        /* search for the conformer */
        if (!strcmp(atom.confName+3, "BK"))
            k_conf = 0;
        else {
            for (k_conf = prot.res[k_res].n_conf - 1; k_conf >= 0; k_conf--) {
                if (!strcmp(atom.confID, prot.res[k_res].conf[k_conf].confID)) {
                    break;
                }
            }
        }
        /* If couldn't find the conformer, add a new one */
        if (k_conf == -1) {
            k_conf = ins_conf(&prot.res[k_res], prot.res[k_res].n_conf, 0);
            strcpy(prot.res[k_res].conf[k_conf].confName, atom.confName);
            strcpy(prot.res[k_res].conf[k_conf].confID, atom.confID);
            prot.res[k_res].conf[k_conf].altLoc = atom.altLoc;
            strcpy(prot.res[k_res].conf[k_conf].history, atom.history);
        }

        /* Insert atom */
        /*--- increase the atom list size by 1 */
        if (!(prot.res[k_res].conf[k_conf].atom = (ATOM *) realloc(prot.res[k_res].conf[k_conf].atom, (prot.res[k_res].conf[k_conf].n_atom + 1) * sizeof(ATOM)))) {
            printf("ins_atom(): Fails resizing memory.\n");
            del_prot(&prot);
            return prot;
        }
        prot.res[k_res].conf[k_conf].n_atom++;
        k_atom = prot.res[k_res].conf[k_conf].n_atom-1;
        
        prot.res[k_res].conf[k_conf].atom[k_atom] = atom;
        prot.res[k_res].conf[k_conf].atom[k_atom].on = 1;

        if ( !strlen(prot.res[k_res].conf[k_conf].history) )
            strcpy(prot.res[k_res].conf[k_conf].history, atom.history);
    }
    return prot;
}

int check_param(ATOM atom, PROT *prot) {
        int   k_res;
        int   k_conf;
        int   k_atom;
	int   n_atom;
        char  confName[6];
        /* search for the residue: each unique combination of
 	*  residue name, chain ID, residue number and insertion code
	*  defines one residue. */
	//printf("ATOM RES: %s\n",atom.resName);
        for (k_res = prot->n_res - 1; k_res >= 0; k_res--) {
            if (!strcmp(atom.resName, prot->res[k_res].resName) &&
                atom.chainID == prot->res[k_res].chainID &&
                atom.resSeq  == prot->res[k_res].resSeq  &&
                atom.iCode   == prot->res[k_res].iCode  )
            {
                break;
            }
        }
        /* If couldn't find the residue, add a new one */
        if (k_res == -1) {
            k_res = ins_res(prot, prot->n_res);
            strcpy(prot->res[k_res].resName, atom.resName);
            prot->res[k_res].chainID = atom.chainID;
            prot->res[k_res].resSeq  = atom.resSeq;
            prot->res[k_res].iCode   = atom.iCode;
            //if (strncmp("HETATM", line, 6)) prot->res[k_res].groupID   = 0;
            //else prot->res[k_res].groupID = 1;
            //assuming that we only write/read ATOM formatted atoms, not HETATM
            prot->res[k_res].groupID = 1; 

            /* Insert a backbone conformer */
            strcpy(confName, atom.resName);
            strcat(confName, "BK");
            if (param_get("NATOM", confName, "", &n_atom)) {
                printf("   WARNING: load_pdb(): \"No NATOM for conformer %s, set to 0.\"\n",confName);
                n_atom=0;
		return -1;
            }
            ins_conf(&prot->res[k_res], 0, n_atom);
            strcpy(prot->res[k_res].conf[0].confName, confName);
            strcpy(prot->res[k_res].conf[0].confID, "0000");
            if (!strncmp(atom.history,"BK",2)) strcpy(prot->res[k_res].conf[0].history, atom.history);
            else strcpy(prot->res[k_res].conf[0].history,"BK________");
            prot->res[k_res].conf[0].altLoc = ' ';
        }
        /* search for the conformer */
        if (!strcmp(atom.confName+3, "BK"))
            k_conf = 0;
        else {
            for (k_conf = prot->res[k_res].n_conf - 1; k_conf >= 0; k_conf--) {
                if (!strcmp(atom.confID, prot->res[k_res].conf[k_conf].confID)) {
                    break;
                }
            }
        }
        /* If couldn't find the conformer, add a new one */
        if (k_conf == -1) {
	    //printf("Param get for: %s\n",atom.confName);
            if (param_get("NATOM", atom.confName, "", &n_atom)) {
                printf("   WARNING: load_pdb(): \"No NATOM for conformer %s, set to 0.\"\n",atom.confName);
		return -1;
            }
            k_conf = ins_conf(&prot->res[k_res], prot->res[k_res].n_conf, n_atom);
            strcpy(prot->res[k_res].conf[k_conf].confName, atom.confName);
            strcpy(prot->res[k_res].conf[k_conf].confID, atom.confID);
            prot->res[k_res].conf[k_conf].altLoc = atom.altLoc;
            strcpy(prot->res[k_res].conf[k_conf].history, atom.history);
        }

        if (param_get("IATOM", atom.confName, atom.name, &k_atom)) {
            printf("   FATAL: load_pdb(): \"No IATOM of conformer \"%s\" atom \"%s\".\n", atom.confName, atom.name);
            printf("Now exiting...\n");
	    exit(0);
        }

        if (k_atom >= prot->res[k_res].conf[k_conf].n_atom) {
            printf("   FATAL: load_pdb(): Confliction in  IATOM of conformer \"%s\" atom \"%s\" and NATOM of conformer \"%s\".\n", atom.confName, atom.name, atom.confName);
	    printf("Now exiting...\n");
	    exit(0);
        }

        if (prot->res[k_res].conf[k_conf].atom[k_atom].on) {
            if (k_conf == 0) {  /* backbone atoms, ignore and continue */
                printf("   Warning: load_pdb(): \"Duplicate backbone atom ignored, \"%3s %3s %c %3d\".\n",atom.name, atom.resName, atom.chainID, atom.resSeq);
            }
            else {
                printf("   FATAL: load_pdb(): \"Atom already exists, \"%3s %3s %c %3d\".\n",
                atom.name, atom.resName, atom.chainID, atom.resSeq);
            }
            return -1;
        }

        prot->res[k_res].conf[k_conf].atom[k_atom] = atom;
        prot->res[k_res].conf[k_conf].atom[k_atom].on = 1;

        if ( !strlen(prot->res[k_res].conf[k_conf].history) ) strcpy(prot->res[k_res].conf[k_conf].history, atom.history);
	return k_res;
};


/*wherever 'check_param()' is called, replace it with this one if you wish to NOT load param*/
void check(ATOM atom, PROT *prot) {
	int   k_res;
        int   k_conf;
        int   k_atom;
	char  confName[6];

        /*  search for the residue: each unique combination of
 	**  residue name, chain ID, residue number and insertion code
 	**  defines one residue. */
        for (k_res = prot->n_res - 1; k_res >= 0; k_res--) {
            if (!strcmp(atom.resName, prot->res[k_res].resName) &&
                atom.chainID == prot->res[k_res].chainID &&
                atom.resSeq  == prot->res[k_res].resSeq  &&
                atom.iCode   == prot->res[k_res].iCode  )
            {
                break;
            }
        }

        /* If couldn't find the residue, add a new one */
        if (k_res == -1) {
            k_res = ins_res(prot, prot->n_res);
            strcpy(prot->res[k_res].resName, atom.resName);
            prot->res[k_res].chainID = atom.chainID;
            prot->res[k_res].resSeq  = atom.resSeq;
            prot->res[k_res].iCode   = atom.iCode;

            /* Insert a backbone conformer */
            strcpy(confName, atom.resName);
            strcat(confName, "BK");
            ins_conf(&prot->res[k_res], 0, 0);
            strcpy(prot->res[k_res].conf[0].confName, confName);
            strcpy(prot->res[k_res].conf[0].confID, "0000");
            if (!strncmp(atom.history,"BK",2)) strcpy(prot->res[k_res].conf[0].history, atom.history);
            else strcpy(prot->res[k_res].conf[0].history,"BK________");
            prot->res[k_res].conf[0].altLoc = ' ';
        }

        /* search for the conformer */
        if (!strcmp(atom.confName+3, "BK"))
            k_conf = 0;
        else {
            for (k_conf = prot->res[k_res].n_conf - 1; k_conf >= 0; k_conf--) {
                if (!strcmp(atom.confID, prot->res[k_res].conf[k_conf].confID)) {
                    break;
                }
            }
        }
        /* If couldn't find the conformer, add a new one */
        if (k_conf == -1) {
            k_conf = ins_conf(&prot->res[k_res], prot->res[k_res].n_conf, 0);
            strcpy(prot->res[k_res].conf[k_conf].confName, atom.confName);
            strcpy(prot->res[k_res].conf[k_conf].confID, atom.confID);
            prot->res[k_res].conf[k_conf].altLoc = atom.altLoc;
            strcpy(prot->res[k_res].conf[k_conf].history, atom.history);
        }

        /* Insert atom */
        /*--- increase the atom list size by 1 */
        if (!(prot->res[k_res].conf[k_conf].atom = (ATOM *) realloc(prot->res[k_res].conf[k_conf].atom, (prot->res[k_res].conf[k_conf].n_atom + 1) * sizeof(ATOM)))) {
            printf("ins_atom(): Fails resizing memory. Now exiting...\n");
            del_prot(prot);
            exit(0);
        }
        prot->res[k_res].conf[k_conf].n_atom++;
        k_atom = prot->res[k_res].conf[k_conf].n_atom-1;

        prot->res[k_res].conf[k_conf].atom[k_atom] = atom;
        prot->res[k_res].conf[k_conf].atom[k_atom].on = 1;

        if ( !strlen(prot->res[k_res].conf[k_conf].history) ) strcpy(prot->res[k_res].conf[k_conf].history, atom.history);
}

ATOM read_full_header(FILE* fp) {
	ATOM  tmp_atom;
	float x,y,z;
	memset(&tmp_atom, 0, sizeof(ATOM));//reset temporary storage for 1 atom to null
	fread(&tmp_atom.serial,sizeof(int),1,fp); fread(tmp_atom.name, sizeof(char),4,fp);
        fread(&tmp_atom.altLoc,sizeof(char),1,fp); fread(tmp_atom.resName,sizeof(char),3,fp);
        fread(&tmp_atom.chainID,sizeof(char), 1,fp); fread(&tmp_atom.resSeq,sizeof(int),1,fp);
        fread(&tmp_atom.iCode,sizeof(char),1,fp); fread(&tmp_atom.iConf,sizeof(int),1,fp);
        fread(&x,sizeof(float),1,fp); tmp_atom.xyz.x = x;
        fread(&y,sizeof(float),1,fp); tmp_atom.xyz.y = y;
        fread(&z,sizeof(float),1,fp); tmp_atom.xyz.z = z;
        fread(&tmp_atom.rad,sizeof(float),1,fp); fread(&tmp_atom.crg,sizeof(float),1,fp);
        fread(tmp_atom.history,sizeof(char),11,fp);
        strcpy(tmp_atom.confName,tmp_atom.resName); strncpy(tmp_atom.confName+3,tmp_atom.history,2); tmp_atom.confName[5] = '\0';
	sprintf(tmp_atom.confID,"%i",tmp_atom.iConf);//length 8
	return tmp_atom;
}

void read_coordinates(FILE *fp, int *iatm, int count, int nb_atoms, PROT *prot) {
        int k;
        float x,y,z;
	int iConf;
	ATOM atom;
        for (k=0; k<nb_atoms; k++) {
                fread(&iConf,sizeof(int),1,fp);
                fread(&x,sizeof(float),1,fp); fread(&y,sizeof(float),1,fp); fread(&z,sizeof(float),1,fp);
		atom.serial = 0;
		strcpy(atom.name,prot->res[count].conf[1].atom[iatm[k]].name);
		atom.altLoc = prot->res[count].conf[1].atom[iatm[k]].altLoc;
		strcpy(atom.resName,prot->res[count].conf[1].atom[iatm[k]].resName);
		atom.chainID = prot->res[count].conf[1].atom[iatm[k]].chainID;
		atom.resSeq = prot->res[count].conf[1].atom[iatm[k]].resSeq;
		atom.iCode = prot->res[count].conf[1].atom[iatm[k]].iCode;

		atom.iConf = iConf;
		atom.xyz.x = x; atom.xyz.y = y; atom.xyz.z = z;
		atom.rad = prot->res[count].conf[1].atom[iatm[k]].rad;
		atom.crg = prot->res[count].conf[1].atom[iatm[k]].crg;
		//strcpy(atom.history,prot->res[count].conf[1].atom[iatm[k]].history);
		fread(atom.history, sizeof(char), 11, fp);
		//sprintf(atom.history,"%s","01R000M000");//fixed 03-31/2009 Pascal
		//printf("Atom History read:%s\n",atom.history);	
	
		strcpy(atom.confName,atom.resName); strncpy(atom.confName+3,atom.history,2); atom.confName[5] = '\0';
		sprintf(atom.confID,"%i",iConf);//length 8
		
		check_param(atom,prot);
        }
};

PROT load_pdb_binary_no_param(FILE *fp) {
    PROT  prot; memset(&prot, 0, sizeof(PROT));
    ATOM  atom;
    int nb_atoms;
    int nb_residues;
    int nb_conformers;
    int i,j,k;
    int i_atom;
    //start reading binary PDB file
    fread(&nb_residues,sizeof(int),1,fp);
    //printf("%i RESIDUES to read in from binary PDB File. Please wait...\n",nb_residues);
    int count = 0;
    int *iatm;
    int k_RES;
    int original_index;
    for(count=0; count<nb_residues; count++) {
	    fread(&nb_conformers,sizeof(int),1,fp);/*read in the number of conformers for this residue*/
	    //fread(&original_index, sizeof(int),1,fp);
	    //printf("Number of conformers for res#%i:%i\n",count,nb_conformers);
	    if(nb_conformers<=2) {
		fread(&nb_atoms,sizeof(int),1,fp);/*read in the number of atoms for the backbone conformer*/
		//printf("Nb of atoms:%i\n",nb_atoms);
		for(i=0; i<nb_atoms; i++) {
        		atom = read_full_header(fp);
			k_RES = check_param(atom,&prot);
		}
		if(nb_conformers==2) {
			fread(&nb_atoms,sizeof(int),1,fp);/*read in the number of atoms for the 1st conformer*/
			//printf("Nb of atoms:%i\n",nb_atoms);
   		        for(i=0; i<nb_atoms; i++) {
                        	atom = read_full_header(fp);
				k_RES = check_param(atom,&prot);
                	}
		}
	    }
	    else {
            	fread(&nb_atoms,sizeof(int),1,fp);/*read in the number of atoms for the backbone conformer*/
		for(i=0; i<nb_atoms; i++) {
                        atom = read_full_header(fp);
			k_RES = check_param(atom,&prot);
                }
		fread(&nb_atoms,sizeof(int),1,fp);/*read in the number of atoms for the 1st conformer*/
				
		iatm = malloc(sizeof(int)*nb_atoms);//required to keep track of the indices the atoms of the conformer belong to
                for(i=0; i<nb_atoms; i++) {
                        atom = read_full_header(fp);
			iatm[i] = iatom(atom.confName, atom.name);
			k_RES = check_param(atom,&prot);
                }

                for(j=2; j<nb_conformers; j++) {
			read_coordinates(fp,iatm,k_RES,nb_atoms,&prot);/*FIXED to handle different residue locations than 'count', by using the k_RES variable*/
                }
		free(iatm);
	    }
	    //prot.res[k_RES].original_index = original_index;
    }

    //FIX SERIAL numbering of all atoms, starting from the 1st atom of the 1st conformer...
    int serial=1;
    for(i=0; i<prot.n_res; i++) {
    	for(j=0; j<prot.res[i].n_conf; j++) {
		for(k=0; k<prot.res[i].conf[j].n_atom; k++) {
			prot.res[i].conf[j].atom[k].serial = serial;
			serial++;
		}
	}
    }
    return prot;
}

int atom2pdbline(char *line, ATOM atom) {
    sprintf(line, "ATOM        %4s%c%3s %c%04d%c%8s%8.3f%8.3f%8.3f %7.3f      %6.3f      %-11s\n",
            atom.name,
            atom.altLoc,
            atom.resName,
            atom.chainID,
            atom.resSeq,
            atom.iCode,
            atom.confID,
            atom.xyz.x,
            atom.xyz.y,
            atom.xyz.z,
            atom.rad,
            atom.crg,
            atom.history);
            return 0;
};


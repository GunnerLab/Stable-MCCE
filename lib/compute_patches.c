#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "mcce.h"

typedef struct {
	char flagged;
} PFLAG;

typedef struct {
	float mdist;
	int cube;
	int res_id;//residue index in protein structure "PROT"
} DIST;

typedef struct {
	int res;
	int flag;
} RES_LIST;

typedef struct {
	float cube_size;
        float x_min, y_min, z_min;
        float x_max, y_max, z_max;
        int nb_bonds;
        int nb_res;
        RES_LIST *reslist;
}CUBE;

int compareCUBE(const void *a, const void *b) {
	CUBE *aa = (CUBE*)a;
	CUBE *bb = (CUBE*)b;
	return ( aa->nb_res - bb->nb_res );
}

//used for qsort ascending order
int compare (const void * a, const void * b) {
  DIST *aa = (DIST*)a;
  DIST *bb = (DIST*)b;
  return (int)( (aa->mdist*10000.f) - (bb->mdist*10000.f) ) ;
}

//we need to find the average minimum distance from prot.res[first_pick].conf[1] to all other non-flagged residues' conformer
//in this computation, we only take into account the atoms of the conformer[1] that are turned on
void build_distance_array(PROT prot, int first_pick, DIST *avg_min_dist, int nb_conformers, float eps) {
	/*we need to find the average minimum distance from prot.res[first_pick].conf[1] to all other non-flagged residues' conformer
	in this computation, we only take into account the atoms of the conformer[1] that are turned on
	return the minimal distance within eps distance from all atoms
	otherwise, return the closest furtest away residue*/
	int i,j,k;
	float tmp;
	ATOM *atm1,*atm2;
	int counter;
	float min_pairwise;
	for(i=0; i<nb_conformers; i++) {
		counter=0;
		min_pairwise = FLT_MAX;
		for(j=0; j<prot.res[avg_min_dist[i].res_id].conf[1].n_atom; j++) {
			atm1 = &prot.res[avg_min_dist[i].res_id].conf[1].atom[j];
			if(!atm1->on) continue;
			for(k=0; k<prot.res[first_pick].conf[1].n_atom; k++) {
				atm2 = &prot.res[first_pick].conf[1].atom[k];
				if(!atm2->on) continue;
				tmp = (atm1->xyz.x-atm2->xyz.x) * (atm1->xyz.x-atm2->xyz.x) + (atm1->xyz.y-atm2->xyz.y) * (atm1->xyz.y-atm2->xyz.y) + (atm1->xyz.z-atm2->xyz.z) * (atm1->xyz.z-atm2->xyz.z);
				tmp = sqrt(tmp);
				if(tmp < min_pairwise) min_pairwise = tmp;
				if(tmp < eps) {
					avg_min_dist[i].mdist += tmp;
					counter++;
				}
			}
		}
		if(counter!=0) avg_min_dist[i].mdist /= counter;
		else {avg_min_dist[i].mdist = min_pairwise;}
	}
	qsort(avg_min_dist,nb_conformers,sizeof(DIST),compare);
};

void build_distance_array_mode3(PROT prot, float Xc, float Yc, float Zc, DIST *avg_min_dist, int nb_conformers) {
        int i,j;
        float tmp;
        ATOM *atm1;
        float min;
        for(i=0; i<nb_conformers; i++) {
                min = FLT_MAX;
                for(j=0; j<prot.res[avg_min_dist[i].res_id].conf[1].n_atom; j++) {
                        atm1 = &prot.res[avg_min_dist[i].res_id].conf[1].atom[j];
                        if(!atm1->on) continue;
                        tmp = (atm1->xyz.x-Xc) * (atm1->xyz.x-Xc) + (atm1->xyz.y-Yc) * (atm1->xyz.y-Yc) + (atm1->xyz.z-Zc) * (atm1->xyz.z-Zc);
                        tmp = sqrt(tmp);
			if(tmp<min) {min=tmp;}
                }
                avg_min_dist[i].mdist = min;
        }
        qsort(avg_min_dist,nb_conformers,sizeof(DIST),compare);
};

float find_distance(PROT prot,int res_id1, int res_id2, float eps) {
	//return the minimal distance within eps distance from all atoms
	//otherwise, return the closest furtest away residue
	int j,k;
	ATOM *atm1, *atm2;
	int counter = 0;
	float tmp;
	float min_dist=0.f;
	float min_pairwise = FLT_MAX;
        for(j=0; j<prot.res[res_id1].conf[1].n_atom; j++) {
                atm1 = &prot.res[res_id1].conf[1].atom[j];
                if(!atm1->on) continue;
                for(k=0; k<prot.res[res_id2].conf[1].n_atom; k++) {
                	atm2 = &prot.res[res_id2].conf[1].atom[k];
                        if(!atm2->on) continue;
                        tmp = (atm1->xyz.x-atm2->xyz.x) * (atm1->xyz.x-atm2->xyz.x) + (atm1->xyz.y-atm2->xyz.y) * (atm1->xyz.y-atm2->xyz.y) + (atm1->xyz.z-atm2->xyz.z) * (atm1->xyz.z-atm2->xyz.z);
                        tmp = sqrt(tmp);
			if (tmp<min_pairwise) min_pairwise = tmp;
			if(tmp<=eps) {min_dist += tmp; counter++;}
                }
        }
	if(counter!=0) {return min_dist /= counter;}
        else {return min_pairwise;}
};

float find_fartest_distance(PROT prot, float Xc, float Yc, float Zc, int res_id) {
	int i,j,k;
	ATOM *atm1;
	float max_dist,tmp;
	max_dist = FLT_MIN;
	for(i=0; i<prot.res[res_id].conf[1].n_atom; i++) {
		atm1 = &prot.res[res_id].conf[1].atom[i];
		if(!atm1->on) continue;
		tmp = (atm1->xyz.x-Xc) * (atm1->xyz.x-Xc) + (atm1->xyz.y-Yc) * (atm1->xyz.y-Yc) + (atm1->xyz.z-Zc) * (atm1->xyz.z-Zc);
		tmp = sqrt(tmp);
		if(tmp>max_dist) {max_dist=tmp;}
	}
	return max_dist;
};

int write_one_patch_as_pdb(char *fname, PROT prot, PRES *patch, char p_number) {
	FILE *stream;
	stream = fopen(fname,"wt");
	int i, j, k, iConf, c,ID;
   	c = 0;
   	for (i=0; i<patch->res_count; i++) {
      		iConf = 0;
		ID = patch->res_ids[i];
      		for (j=0; j<prot.res[ID].n_conf; j++) {
         		for (k=0; k<prot.res[ID].conf[j].n_atom; k++) {
            			if (!(prot.res[ID].conf[j].atom[k].on)) continue;
            			if (c<99999) c++;
            			fprintf(stream, "ATOM  %5d %4s%c%3s %c%04d%c%03d%8.3f%8.3f%8.3f %7.3f      %6.3f      %-11s\n",
                            	c, prot.res[ID].conf[j].atom[k].name,
                            	prot.res[ID].conf[j].altLoc,
                            	prot.res[ID].resName,
                            	p_number,/*prot.res[ID].chainID,*/
                            	prot.res[ID].resSeq,
                            	prot.res[ID].iCode,
                            	iConf,
                            	prot.res[ID].conf[j].atom[k].xyz.x,
                            	prot.res[ID].conf[j].atom[k].xyz.y,
                            	prot.res[ID].conf[j].atom[k].xyz.z,
                            	prot.res[ID].conf[j].atom[k].rad,
                            	prot.res[ID].conf[j].atom[k].crg,
                            	prot.res[ID].conf[j].history);
         		}
         		iConf++;
      		}
   	}
	fclose(stream);
   	return 0;
}

int write_full_patch_as_pdb(FILE *stream, PROT prot, PATCH *patches) {
        int i, j, k,m, iConf, c,ID;
	char p_number;
        c = 0;
	for(m=0; m<patches->count; m++) {
		p_number = m+65;
        	for (i=0; i<patches->pat[m].res_count; i++) {
                	iConf = 0;
                	ID = patches->pat[m].res_ids[i];
                	for (j=0; j<prot.res[ID].n_conf; j++) {
				for (k=0; k<prot.res[ID].conf[j].n_atom; k++) {
                                	if (!(prot.res[ID].conf[j].atom[k].on)) continue;
                                	if (c<99999) c++;
                                	fprintf(stream, "ATOM  %5d %4s%c%3s %c%04d%c%03d%8.3f%8.3f%8.3f %7.3f      %6.3f      %-11s\n",
                                	c, prot.res[ID].conf[j].atom[k].name,
                                	prot.res[ID].conf[j].altLoc,
                               		prot.res[ID].resName,
                                	p_number,/*prot.res[ID].chainID,*/
                                	prot.res[ID].resSeq,
                                	prot.res[ID].iCode,
                                	iConf,
                                	prot.res[ID].conf[j].atom[k].xyz.x,
                                	prot.res[ID].conf[j].atom[k].xyz.y,
                                	prot.res[ID].conf[j].atom[k].xyz.z,
                                	prot.res[ID].conf[j].atom[k].rad,
                                	prot.res[ID].conf[j].atom[k].crg,
                                	prot.res[ID].conf[j].history);
                        	}
                        	iConf++;
                	}
        	}
		fprintf(stream,"TER\n");
	}
        return 0;
}

float full_atom_distance(PROT prot, int res1, CUBE *cube) {
	int i,j,k,m;
	int res2;
	ATOM *atom1, *atom2;
	float distance = 0.0f;
	float tmp;
	float new_distance = FLT_MAX;
	for(i=0; i<cube->nb_res; i++) {
		if(cube->reslist[i].flag) continue;
		res2 = cube->reslist[i].res;
		distance = 0.0f;
		for(j=0; j<prot.res[res2].conf[1].n_atom; j++) {
			atom2 = &prot.res[res2].conf[1].atom[j];
			if(!atom2->on) continue;
			for(k=0; k<prot.res[res1].conf[1].n_atom; k++) {
				atom1 = &prot.res[res1].conf[1].atom[k];
				if(!atom1->on) continue;
				tmp = (atom1->xyz.x-atom2->xyz.x) * (atom1->xyz.x-atom2->xyz.x) + (atom1->xyz.y-atom2->xyz.y) * (atom1->xyz.y-atom2->xyz.y) + (atom1->xyz.z-atom2->xyz.z) * (atom1->xyz.z-atom2->xyz.z);
				tmp = sqrt(tmp);
				distance += tmp;
			}
		}
		if(distance < new_distance) {
			new_distance = distance;
		}
	}
	return new_distance;
};


int atom_atom_full_distance(PROT prot, int res1, int res2, float param) {
	int i,j;
	ATOM *atm1, *atm2;
	float tmp;
	for(i=0; i<prot.res[res1].conf[1].n_atom; i++) {
		atm1 = &prot.res[res1].conf[1].atom[i];
		if(!atm1->on) continue;
		for(j=0; j<prot.res[res2].conf[1].n_atom; j++) {
			atm2 = &prot.res[res2].conf[1].atom[i];
			if(!atm2->on) continue;
			tmp = (atm1->xyz.x-atm2->xyz.x) * (atm1->xyz.x-atm2->xyz.x) + (atm1->xyz.y-atm2->xyz.y) * (atm1->xyz.y-atm2->xyz.y) + (atm1->xyz.z-atm2->xyz.z) * (atm1->xyz.z-atm2->xyz.z);
			tmp = sqrt(tmp);
			if (tmp < param) return 1;
		}
	}
	return 0;
};


void compute_boundary1(PROT prot,PATCH *patches) {
	int i,j,k;
	PFLAG *initial = malloc(sizeof(PFLAG)*prot.n_res);
	PFLAG *list = malloc(sizeof(PFLAG)*prot.n_res);
        for(i=0; i<prot.n_res; i++) {
                if(prot.res[i].n_conf>=2) {initial[i].flagged = 0;}
		else {initial[i].flagged = 1;}
        }
	
	for(i=0; i<patches->count; i++) {
		/*for each patch, build list of residues that does not include the ones already in the patch*/
		memmove(list,initial,sizeof(PFLAG)*prot.n_res);/*copy initial state into list to be processed by this conformer*/
		for(j=0; j<patches->pat[i].res_count; j++) {list[patches->pat[i].res_ids[j]].flagged = 1;}/*ignore the ones that are already in this patch, so flag them*/
		patches->pat[i].b1_count = 0;
		patches->pat[i].boundary1 = calloc(0,sizeof(int));
				
		for(j=0; j<prot.n_res; j++) {
			if(list[j].flagged) continue;/*ignore the ones that we just flagged*/
			for(k=0; k<patches->pat[i].res_count; k++) {
				/*if one atom was found to be within the b1_param distance then add the residue to boundarylist and flag this 'j' residue*/
				if (atom_atom_full_distance(prot, j, patches->pat[i].res_ids[k], patches->b1_param)) {
					patches->pat[i].boundary1 = realloc(patches->pat[i].boundary1,sizeof(int)*(patches->pat[i].b1_count+1));
					patches->pat[i].boundary1[patches->pat[i].b1_count] = j;
					patches->pat[i].b1_count++;
					break;
				}
			}
		}
		//printf("Patch[%i]: %i RES in boundary1\n",i,patches->pat[i].b1_count);
	}
	free(initial); free(list);
} 

void compute_boundary2(PROT prot, PATCH *patches) {
        int i,j,k;
        PFLAG *initial = malloc(sizeof(PFLAG)*prot.n_res);
        PFLAG *list = malloc(sizeof(PFLAG)*prot.n_res);
        for(i=0; i<prot.n_res; i++) {
                if(prot.res[i].n_conf>=2) {initial[i].flagged = 0;}
                else {initial[i].flagged = 1;}
        }

        for(i=0; i<patches->count; i++) {
                /*for each patch, build list of residues that does not include the ones already in the patch*/
                memmove(list,initial,sizeof(PFLAG)*prot.n_res);/*copy initial state into list to be processed by this conformer*/
                /*ignore the ones that are already in this patch, so flag them*/
		for(j=0; j<patches->pat[i].res_count; j++) {
			list[patches->pat[i].res_ids[j]].flagged = 1;
		}
		/*ignore the ones that are also in the 1st boundary in this patch, so flag them*/
		for(j=0; j<patches->pat[i].b1_count; j++) {
			list[patches->pat[i].boundary1[j]].flagged = 1;
		}
		
                patches->pat[i].b2_count = 0;
                patches->pat[i].boundary2 = calloc(0,sizeof(int));

                for(j=0; j<prot.n_res; j++) {
                        if(list[j].flagged) continue;/*ignore the ones that we just flagged*/
                        for(k=0; k<patches->pat[i].res_count; k++) {
                                /*if one atom was found to be within the b1_param distance then add the residue to boundarylist and flag this 'j' residue*/
                                if (atom_atom_full_distance(prot, j, patches->pat[i].res_ids[k], patches->b1_param + patches->b2_param)) {
                                        patches->pat[i].boundary2 = realloc(patches->pat[i].boundary2,sizeof(int)*(patches->pat[i].b2_count+1));
                                        patches->pat[i].boundary2[patches->pat[i].b2_count] = j;
                                        patches->pat[i].b2_count++;
                                        break;
                                }
                        }
                }
		//printf("Patch[%i]: %i RES in boundary2\n",i,patches->pat[i].b2_count);
        }
	free(initial); free(list);
}


// mode=0: overlapping patches
// mode=1: non-overlapping patches
// mode=2: proper grid patching
// mode=3: new spherical focus mode for isolation of mutant residues etc...
// patch results will vary; need to experiment to find which one produces best results
// 'N' = number of maximum residues to include in a patch
// 'eps' = epsilon value. You want at least 1 residue already in the patch within this distance value of the next one to be inserted in the patch
PATCH compute_patches(PROT prot, int mode, int N, float eps, float b1_param, float b2_param, int resSeq_center, int iso_res) {
	printf("Now Computing patch subtrees...\n");
	int i,j,k,m;
	int number_conformers=0;
	int number_flagged_conformers=0;
        int first_pick;
        int to_insert_flag;
        int inserted;
        int closest_index;
        float closest_dist;
        float tmp;
	PFLAG *flags = malloc(sizeof(PFLAG)*prot.n_res);
	PATCH patches;	patches.count=0; patches.pat = calloc(0,sizeof(PRES));
	DIST *avg_min_dist;
	patches.b1_param = b1_param;
        patches.b2_param = b2_param;
	//count the total number of conformers we have; we check for equality of 2; the backbone (1) and the first initial state conformer (2)
	//if we have at least 2, then we know for sure we have 1 initial state conformer for this residue.
	//We need to know of the total number of conformers so we can keep track of the number of patch-flagged conformers
	//We don't include a residue in computation if it doesn't have an initial state conformer and flag it as already been patched
	for(i=0; i<prot.n_res; i++) {
		if(prot.res[i].n_conf>=2) {number_conformers++; flags[i].flagged = 0;}
		else {flags[i].flagged = 1;}
	}
	printf("Total number of conformers to patch in subunits: %i\n",number_conformers);

	/*pick the first non flagged residue as the one to start the patching from*/
        for(i=0; i<prot.n_res; i++) {
        	if(!flags[i].flagged) {first_pick = i; flags[first_pick].flagged = 1; break;}
        }
	//first_pick = 66;
	//flags[66].flagged = 1;

	if(mode==0) {
		/*while we have more conformers to patch up together*/
		while(number_flagged_conformers!=number_conformers) {
			avg_min_dist = malloc((number_conformers-1)*sizeof(DIST));
			/*we now skip the first one and assign the index of the protein structure of the avg_min_dist array structure*/
			/*so there is a 1-1 mapping*/
			j=0;
			for(i=0; i<prot.n_res; i++) {
				if( (i!=first_pick) && (prot.res[i].n_conf>=2) ) {
					avg_min_dist[j].res_id = i;
					avg_min_dist[j].mdist = 0.f;
					j++;
				}
			}
			build_distance_array(prot,first_pick,avg_min_dist,number_conformers-1,eps);

			/*insert the 'first_pick' conformer inside the patch*/
			patches.pat = realloc(patches.pat,sizeof(PRES)*(patches.count+1));
			patches.pat[patches.count].res_count = 1; patches.pat[patches.count].res_ids = calloc(1,sizeof(int));
			patches.pat[patches.count].res_ids[0] = first_pick;/*just need to save res_id from original protein structure*/
			flags[first_pick].flagged = 1; number_flagged_conformers++;

			inserted = 1;
			for(i=0; i<number_conformers-1; i++) {
				to_insert_flag = 0;
				for(j=0; j<patches.pat[patches.count].res_count; j++) {
					if( find_distance(prot,patches.pat[patches.count].res_ids[j],avg_min_dist[i].res_id,eps) <= eps && (flags[avg_min_dist[i].res_id].flagged!=1) ) {
						to_insert_flag = 1; break;
					}
				}
				/*OK we're good, found at least 1 within eps distance*/
				/*insert 'i' into patch, set flag as inserted*/
				if(to_insert_flag) {
					patches.pat[patches.count].res_ids = realloc(patches.pat[patches.count].res_ids,sizeof(int)*(patches.pat[patches.count].res_count+1));
					patches.pat[patches.count].res_ids[patches.pat[patches.count].res_count] = avg_min_dist[i].res_id;
					patches.pat[patches.count].res_count++;
					if(!flags[avg_min_dist[i].res_id].flagged) number_flagged_conformers++;
					flags[avg_min_dist[i].res_id].flagged = 1;
					inserted++;
					if(inserted==N) break;/*maximum number of residues for core of a patch*/
				}
			}
		
			/*now look for the closest non-flagged residue to this patch and use this residue as the first_pick for the next patch.*/
			/*if there are no non-flagged residues, you're done.*/
			if (number_flagged_conformers!=number_conformers) {
				/*printf("Nb Conformers: %i Nb Flagged Conformers: %i\n",number_conformers,number_flagged_conformers);*/
				closest_dist = FLT_MAX; closest_index = 0;
               			for(i=0; i<prot.n_res; i++) {
              				if(!flags[i].flagged) {
					for(j=0; j<patches.pat[patches.count].res_count; j++) {
							tmp = find_distance(prot,patches.pat[patches.count].res_ids[j],i,eps);
							if( tmp < closest_dist ) {
								closest_dist = tmp; closest_index = i; k = patches.pat[patches.count].res_ids[j];
							}
						}
					}
               			}
				first_pick = closest_index;
			}
			/*
			printf("Patch[%i]: ",patches.count); printf("\tInserted %i residues (by id): ",inserted);
			for(i=0; i<patches.pat[patches.count].res_count; i++) {
				printf("%i ",patches.pat[patches.count].res_ids[i]);
			}
			printf("\n");
			printf("\t\t\tNext closest residue: %4i from residue %4i in patch. Their distance: %8.3f\n",first_pick,k,closest_dist);
			*/
			free(avg_min_dist);
			patches.count++;
			/*exit(0);*/
		}
	}
	else if(mode==1) {
		/*while we have more conformers to patch up together*/
                while(number_flagged_conformers!=number_conformers) {
                        avg_min_dist = malloc((number_conformers-number_flagged_conformers-1)*sizeof(DIST));
                        /*we now skip the first one and assign the index of the protein structure of the avg_min_dist array structure*/
                        /*so there is a 1-1 mapping*/
                        j=0;
                        for(i=0; i<prot.n_res; i++) {
				if(flags[i].flagged) continue;
                                if( (i!=first_pick) && (prot.res[i].n_conf>=2) ) {
                                        avg_min_dist[j].res_id = i;
                                        avg_min_dist[j].mdist = 0.f;
                                        j++;
                                }
                        }
                        build_distance_array(prot,first_pick,avg_min_dist,number_conformers-number_flagged_conformers-1,eps);

                        /*insert the 'first_pick' conformer inside the patch*/
                        patches.pat = realloc(patches.pat,sizeof(PRES)*(patches.count+1));
                        patches.pat[patches.count].res_count = 1; patches.pat[patches.count].res_ids = calloc(1,sizeof(int));
                        patches.pat[patches.count].res_ids[0] = first_pick;/*just need to save res_id from original protein structure*/
                        flags[first_pick].flagged = 1;

                        inserted = 1;
                        for(i=0; i<number_conformers-number_flagged_conformers-1; i++) {
                                to_insert_flag = 0;
                                for(j=0; j<patches.pat[patches.count].res_count; j++) {
                                        if(find_distance(prot,patches.pat[patches.count].res_ids[j],avg_min_dist[i].res_id,eps) <= eps) {
                                                to_insert_flag = 1; break;
                                        }
                                }
                                /*OK we're good, found at least 1 within eps distance*/
                                /*insert 'i' into patch, set flag as inserted*/
                                if(to_insert_flag) {
                                        patches.pat[patches.count].res_ids = realloc(patches.pat[patches.count].res_ids,sizeof(int)*(patches.pat[patches.count].res_count+1));
                                        patches.pat[patches.count].res_ids[patches.pat[patches.count].res_count] = avg_min_dist[i].res_id;
                                        patches.pat[patches.count].res_count++;
                                        flags[avg_min_dist[i].res_id].flagged = 1;
                                        inserted++;
                                        if(inserted==N) break;/*maximum number of residues for a patch*/
                                }
                        }
			
			number_flagged_conformers += inserted;//update the number of flagged conformers based on the number of insertions done			

                        /*now look for the closest non-flagged residue to this patch and use this residue as the first_pick for the next patch.*/
                        /*if there are no non-flagged residues, you're done.*/
                        if (number_flagged_conformers!=number_conformers) {
                                /*printf("Nb Conformers: %i Nb Flagged Conformers: %i\n",number_conformers,number_flagged_conformers);*/
                                closest_dist = FLT_MAX; closest_index = 0;
                                for(i=0; i<prot.n_res; i++) {
                                        if(!flags[i].flagged) {
                                        for(j=0; j<patches.pat[patches.count].res_count; j++) {
                                                        tmp = find_distance(prot,patches.pat[patches.count].res_ids[j],i,eps);
                                                        if( tmp < closest_dist ) {
                                                                closest_dist = tmp; closest_index = i; k = patches.pat[patches.count].res_ids[j];
                                                        }
                                                }
                                        }
                                }
                                first_pick = closest_index;
                        }
			/*
                        printf("Patch[%i]: ",patches.count); printf("\tInserted {%3i} residues (by id): ",inserted);
                        for(i=0; i<patches.pat[patches.count].res_count; i++) {
                                printf("%i ",patches.pat[patches.count].res_ids[i]);
                        }
                        printf("\n");
                        printf("\t\t\tNext closest residue: %4i from residue %4i in patch. Their distance: %8.3f\n",first_pick,k,closest_dist);
			printf("Number Conformers: %i Number Flagged Conformers: %i\n",number_conformers,number_flagged_conformers);
			*/
                        free(avg_min_dist);
                        patches.count++;
                        /*exit(0);*/
		}
	}
	else if(mode==2) {

		int merge_cutoff = env.patch_merge_cutoff;
		double cube_size = env.patch_cube_size;//Angstrom
		double x_min, y_min, z_min, x_max, y_max, z_max;
		int x_size, y_size, z_size;
		double tmp;
		CUBE ***cubes;
		x_min = y_min = z_min = FLT_MAX;
		x_max = y_max = z_max = FLT_MIN;
		double x,y,z;
		for(i=0; i<prot.n_res; i++) {
			for(j=0; j<prot.res[i].n_conf; j++) {
				for(k=0; k<prot.res[i].conf[j].n_atom; k++) {
					if(!prot.res[i].conf[j].atom[k].on) continue;
					if(prot.res[i].conf[j].atom[k].xyz.x < x_min) {x_min = prot.res[i].conf[j].atom[k].xyz.x;}
					if(prot.res[i].conf[j].atom[k].xyz.x > x_max) {x_max = prot.res[i].conf[j].atom[k].xyz.x;}
                                               
					if(prot.res[i].conf[j].atom[k].xyz.y < y_min) {y_min = prot.res[i].conf[j].atom[k].xyz.y;}
                                       	if(prot.res[i].conf[j].atom[k].xyz.y > y_max) {y_max = prot.res[i].conf[j].atom[k].xyz.y;}
                                       	if(prot.res[i].conf[j].atom[k].xyz.z < z_min) {z_min = prot.res[i].conf[j].atom[k].xyz.z;}
                                       	if(prot.res[i].conf[j].atom[k].xyz.z > z_max) {z_max = prot.res[i].conf[j].atom[k].xyz.z;}
				}
			}
		}	
		x_min -= 0.5f; y_min -= 0.5f; z_min -= 0.5f; x_max += 0.5f; y_max += 0.5f; z_max += 0.5f;
		
		//need to re-size x_max and x_min
		x_size = (int)ceil((x_max-x_min)/cube_size);
		y_size = (int)ceil((y_max-y_min)/cube_size);
		z_size = (int)ceil((z_max-z_min)/cube_size);	

                printf("Grid corners:\n");
                printf("\t 1st Corner x:%f y:%f z:%f\n",x_min, y_min, z_min); printf("\t 2nd Corner x:%f y:%f z:%f\n",x_min, y_max, z_min);
                printf("\t 3th Corner x:%f y:%f z:%f\n",x_min, y_min, z_max); printf("\t 4th Corner x:%f y:%f z:%f\n",x_min, y_max, z_max);
                printf("\t 5th Corner x:%f y:%f z:%f\n",x_max, y_min, z_min); printf("\t 6th Corner x:%f y:%f z:%f\n",x_max, y_max, z_min);
                printf("\t 7th Corner x:%f y:%f z:%f\n",x_max, y_min, z_max); printf("\t 8th Corner x:%f y:%f z:%f\n",x_max, y_max, z_max);
	
		printf("Number of cubes in X: %i in Y: %i in Z: %i Total: %i\n",x_size, y_size, z_size, x_size*y_size*z_size);
		
		cubes = calloc(x_size, sizeof(CUBE**));
		if(cubes!=NULL) {
			for(i=0; i<x_size; i++) {
				cubes[i] = calloc(y_size,sizeof(CUBE*));
				if(cubes[i]!=NULL) {
					for(j=0; j<y_size; j++) {
						cubes[i][j] = calloc(z_size,sizeof(CUBE));
						if(cubes[i][j]==NULL) {
							printf("Memory allocation error in CUBE z_size[%i][%i].\n",i,j); exit(0);
						}
					}
				}
				else {
					printf("Memory allocation error in CUBE y_size[%i].\n",i); exit(0);
				}
			}
		}
		else {
			printf("Memory allocation error in CUBE x_size.\n"); exit(0);
		}
		
		FILE *grid_pdb;
		grid_pdb = fopen("patch_grid.pdb","wt"); 
		int count=0;
		int cs[8];
		x = x_min;
		for(i=0; i<x_size; i++) {
			y = y_min;
			for(j=0; j<y_size; j++) {
				z = z_min;
				for(k=0; k<z_size; k++) {
					cubes[i][j][k].cube_size = cube_size;
					cubes[i][j][k].x_min = x;
					cubes[i][j][k].y_min = y;
					cubes[i][j][k].z_min = z;
					cubes[i][j][k].x_max = x+cube_size;
					cubes[i][j][k].y_max = y+cube_size;
					cubes[i][j][k].x_max = z+cube_size;
					cubes[i][j][k].nb_res = 0;
					cubes[i][j][k].nb_bonds = 0;
					
					cs[0] = ++count; cs[1] = ++count; cs[2] = ++count; cs[3] = ++count; cs[4] = ++count; cs[5] = ++count; cs[6] = ++count; cs[7] = ++count;
					
					fprintf(grid_pdb,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)cs[0],"GRD","COR",(long int)1,x,          y,	                z+cube_size,1.0,0.0);
					fprintf(grid_pdb,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)cs[1],"GRD","COR",(long int)1,x,          y+cube_size,        z+cube_size,1.0,0.0);
					fprintf(grid_pdb,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)cs[2],"GRD","COR",(long int)1,x+cube_size,y,                  z+cube_size,1.0,0.0);
					fprintf(grid_pdb,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)cs[3],"GRD","COR",(long int)1,x+cube_size,y+cube_size,        z+cube_size,1.0,0.0);
					fprintf(grid_pdb,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)cs[4],"GRD","COR",(long int)1,x,          y,                  z,          1.0,0.0);
					fprintf(grid_pdb,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)cs[5],"GRD","COR",(long int)1,x,          y+cube_size,        z,          1.0,0.0);
					fprintf(grid_pdb,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)cs[6],"GRD","COR",(long int)1,x+cube_size,y,                  z,          1.0,0.0);
					fprintf(grid_pdb,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)cs[7],"GRD","COR",(long int)1,x+cube_size,y+cube_size,        z,          1.0,0.0);

					//printf("Min:%8.3f %8.3f %8.3f Max:%8.3f %8.3f %8.3f\n",x,y,z, x+cube_size,y+cube_size,z+cube_size);
					z += cube_size;
				}
				y += cube_size;
			}
			x += cube_size;
		}
		count=0;
		for(i=0; i<x_size; i++) {
			for(j=0; j<y_size; j++) {
				for(k=0; k<z_size; k++) {
					cs[0] = ++count; cs[1] = ++count; cs[2] = ++count; cs[3] = ++count; cs[4] = ++count; cs[5] = ++count; cs[6] = ++count; cs[7] = ++count;
					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[4],cs[5]);
					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[4],cs[6]);
					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[5],cs[7]);
					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[6],cs[7]);

					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[0],cs[1]);
                                        fprintf(grid_pdb, "CONECT%5i%5i\n",cs[0],cs[2]);
                                        fprintf(grid_pdb, "CONECT%5i%5i\n",cs[1],cs[3]);
                                        fprintf(grid_pdb, "CONECT%5i%5i\n",cs[2],cs[3]);
					
					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[0],cs[4]);
					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[1],cs[5]);
					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[2],cs[6]);
					fprintf(grid_pdb, "CONECT%5i%5i\n",cs[3],cs[7]);
				}
			}
		}
		fprintf(grid_pdb, "END\n");
		fclose(grid_pdb);

		int tot_atoms, cube_count, inserted, largest_idx;
		int x_idx, y_idx, z_idx;
		char str[8];
		float largest;
		typedef struct {
			int x_idx, y_idx, z_idx;
			float weight;
			char str[8];
		}CUBE_ATM;
		/*go through all residues of the protein which have >1 conformers and assign them to the cube they belong to*/
		for(i=0; i<prot.n_res; i++) {//for all residues
			if(prot.res[i].n_conf>1) {
				tot_atoms = cube_count = 0;
				CUBE_ATM *atms = calloc(0,sizeof(CUBE_ATM));
				for(j=0; j<prot.res[i].conf[1].n_atom; j++) {//all heavy atoms in this initial conformer
					if(!prot.res[i].conf[1].atom[j].on) continue;
					inserted = 0;
					x_idx = ((int)((prot.res[i].conf[1].atom[j].xyz.x-x_min) / (cube_size)));
					y_idx = ((int)((prot.res[i].conf[1].atom[j].xyz.y-y_min) / (cube_size)));
					z_idx = ((int)((prot.res[i].conf[1].atom[j].xyz.z-z_min) / (cube_size)));
					sprintf(str,"%02i%02i%02i",x_idx, y_idx, z_idx);
					for(k=0; k<cube_count; k++) {
						if( strcmp(str, atms[k].str)==0 ) {
							atms[k].weight += 1.0f; inserted=1; break;
						}
					}
					if(!inserted) {
						atms = realloc(atms,sizeof(CUBE_ATM)*(cube_count+1));
						atms[cube_count].weight = 1.0f;
						atms[cube_count].x_idx = x_idx; atms[cube_count].y_idx = y_idx; atms[cube_count].z_idx = z_idx;
						sprintf(atms[cube_count].str,"%02i%02i%02i",x_idx, y_idx, z_idx);
						cube_count++;
					}
					tot_atoms += 1;
				}
				largest = FLT_MIN; largest_idx = 0;
				for(j=0; j<cube_count; j++) {atms[j].weight /= tot_atoms;}
				for(j=0; j<cube_count; j++) {
					if(atms[j].weight > largest) {largest = atms[j].weight; largest_idx = j;}
				}
				x_idx = atms[largest_idx].x_idx; y_idx = atms[largest_idx].y_idx; z_idx = atms[largest_idx].z_idx;
				cubes[x_idx][y_idx][z_idx].reslist = realloc(cubes[x_idx][y_idx][z_idx].reslist, sizeof(RES_LIST)*(cubes[x_idx][y_idx][z_idx].nb_res+1));
				cubes[x_idx][y_idx][z_idx].reslist[cubes[x_idx][y_idx][z_idx].nb_res].res = i;
				cubes[x_idx][y_idx][z_idx].nb_res++;
				free(atms);
			}//if-statement
		}//for-loop
		count=0;
		int m, p;
		
		
		for(i=0; i<x_size; i++) {
                        for(j=0; j<y_size; j++) {
                                for(k=0; k<z_size; k++) {
                                	if(cubes[i][j][k].nb_res>0) {
                                		FILE *cube_out;
						char out_str[64];
						sprintf(out_str,"./pdb_patches/patch_%02i.pdb",count);
						cube_out = fopen(out_str,"wt");
						for(m=0; m<cubes[i][j][k].nb_res; m++) {
							int ID = cubes[i][j][k].reslist[m].res;
							for(p=0; p<prot.res[ID].conf[1].n_atom; p++) {
								fprintf(cube_out, "ATOM  %5d %4s%c%3s %c%04d%c%03d%8.3f%8.3f%8.3f %7.3f      %6.3f      %-11s\n",
                                        			0, prot.res[ID].conf[1].atom[p].name,
                                        			prot.res[ID].conf[1].altLoc,
                                        			prot.res[ID].resName,
                                        			prot.res[ID].chainID,
                                        			prot.res[ID].resSeq,
                                        			prot.res[ID].iCode,
                                        			0,
                                        			prot.res[ID].conf[1].atom[p].xyz.x,
                                        			prot.res[ID].conf[1].atom[p].xyz.y,
                                        			prot.res[ID].conf[1].atom[p].xyz.z,
                                        			prot.res[ID].conf[1].atom[p].rad,
                                        			prot.res[ID].conf[1].atom[p].crg,
                                        			prot.res[ID].conf[1].history);
							}
						}
						fclose(cube_out);
						count++;
					}
				}
                        }
                }
		

		CUBE *merge_cubes = calloc(count, sizeof(CUBE));
		count=0;
		for(i=0; i<x_size; i++) {
                        for(j=0; j<y_size; j++) {
                                for(k=0; k<z_size; k++) {
                                        if(cubes[i][j][k].nb_res>0) {
						memmove(&merge_cubes[count], &cubes[i][j][k], sizeof(CUBE));
						merge_cubes[count].reslist = calloc(merge_cubes[count].nb_res,sizeof(RES_LIST));
						for(m=0; m<cubes[i][j][k].nb_res; m++) {
							merge_cubes[count].reslist[m].res = cubes[i][j][k].reslist[m].res;
							merge_cubes[count].reslist[m].flag = 0;
						}
						count++;
					}
				}
			}
		}
		qsort(merge_cubes, count, sizeof(CUBE), compareCUBE);//sort cubes in ascending order of nb_res

		//all cubes with less than 6 resides, put in nearest cube to which the atom - atom full sum distance is closest
		int c_count;
		int trial;
		int smallest_idx;
		int smallest;
		int nb_active_atoms = 0;
		float atom_dist_cutoff = env.patch_merge_atm_cutoff;
		merge:
		for(i=0; i<count; i++) {
			if(merge_cubes[i].nb_res < merge_cutoff) {
				for(j=0; j<merge_cubes[i].nb_res; j++) {//assign each residue to the cube for which its distance is minimal
					DIST *cdist = calloc(count-1,sizeof(DIST));
					c_count=0;
					for(k=0; k<count; k++) {
						if(k==i) continue;//skip the cube we're in right now (since we want to redistribute its current residues)
						cdist[c_count].mdist = full_atom_distance(prot, merge_cubes[i].reslist[j].res, &merge_cubes[k]);//returns the minimum distance residue<->residue
						cdist[c_count].cube = k;
						c_count++;
					}
					qsort(cdist, c_count, sizeof(DIST), compare);//sort the distances in ascending order
					
					nb_active_atoms = 0;
					for(k=0; k<prot.res[merge_cubes[i].reslist[j].res].conf[1].n_atom; k++) {
						if(!prot.res[merge_cubes[i].reslist[j].res].conf[1].atom[k].on) continue;
						nb_active_atoms++;
					}
					int *list = calloc(c_count,sizeof(int));
					trial = 1; list[0] = cdist[0].cube;
					for(k=trial; k<c_count; k++) {
						//printf("Cube[%i] - Cube[0] = %f - %f = %f < %f ?\n",k, cdist[k].mdist, cdist[0].mdist, cdist[k].mdist - cdist[0].mdist, nb_active_atoms*(atom_dist_cutoff));
						if( (cdist[k].mdist - cdist[0].mdist) < nb_active_atoms*(atom_dist_cutoff) ) {
							list[trial] = cdist[k].cube; trial++;
						} 
					}
					smallest = 30000000; smallest_idx = 0;
					for(k=0; k<trial; k++) {
						if(merge_cubes[list[k]].nb_res < smallest) {
							smallest = merge_cubes[list[k]].nb_res; smallest_idx = k;
						}
					}
					merge_cubes[list[smallest_idx]].reslist = realloc(merge_cubes[list[smallest_idx]].reslist,sizeof(RES_LIST)* (merge_cubes[list[smallest_idx]].nb_res+1));
					merge_cubes[list[smallest_idx]].reslist[merge_cubes[list[smallest_idx]].nb_res].res =  merge_cubes[i].reslist[j].res;
					merge_cubes[list[smallest_idx]].reslist[merge_cubes[list[smallest_idx]].nb_res].flag = 1;
					merge_cubes[list[smallest_idx]].nb_res++;
					free(cdist);
					free(list);
				}
				//remove this cube
				CUBE *tmp_merge = calloc(count-1, sizeof(CUBE));
				trial=0;
				for(j=0; j<count; j++) {
					if (j==i) continue;
					memmove(&tmp_merge[trial],&merge_cubes[j], sizeof(CUBE));
					tmp_merge[trial].reslist = calloc(merge_cubes[j].nb_res, sizeof(RES_LIST));
					for(k=0; k<merge_cubes[j].nb_res; k++) {
						tmp_merge[trial].reslist[k].res = merge_cubes[j].reslist[k].res;
						tmp_merge[trial].reslist[k].flag = 0;
					}
					trial++;
				}
				free(merge_cubes);
				merge_cubes = calloc(count-1,sizeof(CUBE));
				count = count-1;
				for(j=0; j<count; j++) {
				        memmove(&merge_cubes[j], &tmp_merge[j], sizeof(CUBE));
                                        merge_cubes[j].reslist = calloc(merge_cubes[j].nb_res, sizeof(RES_LIST));
                                        for(k=0; k<merge_cubes[j].nb_res; k++) {
                                                merge_cubes[j].reslist[k].res = tmp_merge[j].reslist[k].res;
                                                merge_cubes[j].reslist[k].flag = 0;
                                        }
				}
				free(tmp_merge);
				qsort(merge_cubes, count, sizeof(CUBE), compareCUBE);
				goto merge;
			}
			else {break;}
		}
		
		/*
 		printf("Number of cubes: %i\n",count);
		for(i=0; i<count; i++) {
			printf("Cube %i - %i residues\n",i+1, merge_cubes[i].nb_res);
		}
		*/
		//finally create the patches for generation of boundary1 and boundary2
		patches.count = count;
	 	patches.pat = realloc(patches.pat,sizeof(PRES)*(count));
		for(i=0; i<count; i++) {
			patches.pat[i].res_count = merge_cubes[i].nb_res;
			patches.pat[i].res_ids = calloc(merge_cubes[i].nb_res, sizeof(int));
			for(j=0; j<merge_cubes[i].nb_res; j++) {
				patches.pat[i].res_ids[j] = merge_cubes[i].reslist[j].res;
			}
		}
	 	for(i=0; i<count; i++) {free(merge_cubes[i].reslist);}
		free(merge_cubes);

	}//outter if-statement
	else if(mode==3) {
		float Xo,Yo,Zo;//center of mass
		Xo = Yo = Zo = 0.0f;
		int n_atm;
		for(k=0; k<prot.res[iso_res].conf[1].n_atom; k++) {
			if(!prot.res[iso_res].conf[1].atom[k].on) continue;
			n_atm++;
			Xo += prot.res[iso_res].conf[1].atom[k].xyz.x;
			Yo += prot.res[iso_res].conf[1].atom[k].xyz.y;
			Zo += prot.res[iso_res].conf[1].atom[k].xyz.z;
		}
		Xo /= n_atm; Yo /= n_atm; Zo /= n_atm;//center from which to include residues
		
		//next, find the initial probe radius from the center of mass to include residues
		float probe_radius = find_fartest_distance(prot, Xo, Yo, Zo, iso_res);
		
		if(probe_radius<20.0f) {probe_radius =20.0f;}		

		printf("Using probe radius of %f from center %f %f %f\n",probe_radius,Xo,Yo,Zo);

		//build list of all residue distances from center
		avg_min_dist = malloc(number_conformers*sizeof(DIST));
                j=0;
		for(i=0; i<prot.n_res; i++) {
                	if( prot.res[i].n_conf>=2 ) {
                        	avg_min_dist[j].res_id = i;
                        	avg_min_dist[j].mdist = 0.f;
                        	j++;
			}
		}
		build_distance_array_mode3(prot, Xo, Yo, Zo, avg_min_dist, number_conformers);
	
		//calculate which residues are to be included in the core of the patch
		int nb_patch_res=0;
		for(i=0; i<number_conformers; i++) {
			if( avg_min_dist[i].mdist <= probe_radius ) continue;
			else {
				nb_patch_res=i; break;
			}
		}
		printf("%i residues to include in patch:\n",nb_patch_res);
		
		patches.pat = realloc(patches.pat,sizeof(PRES)*1);
		patches.pat[patches.count].res_ids = calloc(nb_patch_res,sizeof(int));
		patches.pat[patches.count].res_count = nb_patch_res;
		for(i=0; i<nb_patch_res; i++) {
			printf("\tresSeq#%i (dist:%f)\n",prot.res[avg_min_dist[i].res_id].resSeq,avg_min_dist[i].mdist);
			patches.pat[patches.count].res_ids[i] = avg_min_dist[i].res_id;
		}
		patches.count=1;

	}//outter if-statement

	compute_boundary1(prot,&patches);
	compute_boundary2(prot,&patches);

	FILE *stream;
	stream = fopen("./pdb_patches/pdb_patch.pdb","wt");
	write_full_patch_as_pdb(stream, prot, &patches);
	fclose(stream);

	return patches;
};

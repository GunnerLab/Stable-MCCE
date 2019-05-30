/* Pascal Comte, Brock University, St. Catharines, Ontario, Canada - 2010
 * Master's Thesis Project: Bio-Inspired optimization and sampling technique for sidechain packing in MCCE
 *
 * Read thesis chapter 4 before attempting to modify code
 * One needs a clear understanding of genetic algorithm & simulated annealing theories and genetic operators,
 * Boltzmann & exponential laws. We used C random number generator although we have
 * a setup for GSL random number.  There are a few papers in the literature that
 * show using a poor 'C random number generator' performs better for genetic algorithms
 * than using a more elaborate one.  Therefore, change the randomization with GREAT CARE.
 *
 * IMPORTANT:
 * This Algorithm for sidechain packing is MEANT to execute on a system with multiple CORES
 * and operates using openMP.  We suggest a system with at least a QUAD-CORE CPU.
 * Original research was compiled using intel compiler v10.0+ on an iCore7 (4 dual SMP=8 cores) machine with 12GB RAM
 *
 * Works that use this algorithm need to cite paper reference: []
 *
 * IMPORTANT:
 * CTR AND NTR TERMINALS ARE NOT USED IN THIS ALGORITHM. PLEASE CHECK RUN.PRM AND TURN OFF LABELING OF TERMINALS
 */

#include "mcce.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <float.h>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>

//#include <plplot.h>

#include "ga_engine.h"

#define VDW_CUTOFF_NEAR  1
#define VDW_ELIMIT_NEAR  999
#define VDW_CUTOFF_FAR   10
#define VDW_ELIMIT_FAR   0

const static float ccutoff_near2 = VDW_CUTOFF_NEAR * VDW_CUTOFF_NEAR;
const static float ccutoff_far2  = VDW_CUTOFF_FAR  * VDW_CUTOFF_FAR;

extern int initial_relaxation(PROT prot);
extern int get_resnbrs(PROT prot, float dlimit);
extern int relax_water(PROT prot);
extern int rm_dupconf_hv(PROT prot);   

/*MCCE initialization for GA Execution
 *The functions called here are the same as in rotamers() of step2
 *The only difference here is that we don't want to re-generate all the rotamers
 *We have a list of already known rotamers
 *The protein structure 'prot' is then ready to be accessed by the GA
 */

GA_STRUCTURE ga;

void update_VDW_RADEPS_UA(PROT *prot) {
        
	typedef struct {
                char name[4]; float vdw_rad; float vdw_eps;
        }TYPE;

        typedef struct {
                char name[4]; char type[4]; int typeI;
        }UATOM;
        typedef struct {
                char name[4]; int nb_atoms; UATOM *atoms;
        }UARES;

        int nb_types = 3;
        FILE *input, *file;
        int nb_res = 0;
        int nb_atoms = 0;
        int i,j,k,L,m;
        char buffer[64];
	char atom_name[4];
	char atom_type[4];
	char res_name[4];

        UARES *uas;
        TYPE *types = calloc(nb_types,sizeof(TYPE));
        sprintf(types[0].name,"C1"); types[0].vdw_rad = 1.9580; types[0].vdw_eps = 0.0994;
        sprintf(types[1].name,"C2"); types[1].vdw_rad = 2.0580; types[1].vdw_eps = 0.1094;
        sprintf(types[2].name,"C3"); types[2].vdw_rad = 2.0580; types[2].vdw_eps = 0.1494;

        input = fopen("./mcce_data/UA/List.txt","rt"); fscanf(input,"%i\n",&nb_res);
        uas = calloc(nb_res, sizeof(UARES));
        chdir("./mcce_data/UA");

        for(i=0; i<nb_res; i++) {
                fscanf(input, "%s\n", buffer);
                memmove(uas[i].name, buffer, 3*sizeof(char));
                printf("%i - %s", i+1, uas[i].name);
                file = fopen(buffer,"rt");
                fscanf(file,"%i\n",&nb_atoms);
                uas[i].atoms = calloc(0,sizeof(UATOM));
                printf(" %i atoms...\n",nb_atoms);
                uas[i].nb_atoms = 0;
                for(j=0; j<nb_atoms; j++) {
                        fscanf(file,"%s %s %s\n", res_name, atom_name, atom_type);
                        for(k=0; k<nb_types; k++) {
                                if ( strncmp(atom_type, types[k].name, 3) == 0) {
                                        uas[i].atoms = realloc(uas[i].atoms, sizeof(UATOM)*(uas[i].nb_atoms+1));
                                        memmove(uas[i].atoms[uas[i].nb_atoms].name, atom_name, sizeof(char)*4);
                                        memmove(uas[i].atoms[uas[i].nb_atoms].type, atom_type, sizeof(char)*4);
                                        uas[i].atoms[uas[i].nb_atoms].typeI = k;
					printf("\t %i - %s %s\n", uas[i].nb_atoms, uas[i].atoms[uas[i].nb_atoms].name, uas[i].atoms[uas[i].nb_atoms].type);
                                        uas[i].nb_atoms++;
                                        break;
                                }
                        }
                }
                fclose(file);
        }
        fclose(input);

	for(i=0; i<nb_res; i++) {
		for(j=0; j<prot->n_res; j++) {
			if( strncmp(uas[i].name, prot->res[j].resName, 3) == 0 ) {
				for(k=1; k<prot->res[j].n_conf; k++) {/*skipping backbone*/
					for(L=0; L<prot->res[j].conf[k].n_atom; L++) {
						if(!prot->res[j].conf[k].atom[L].on) {continue;}
						for(m=0; m<uas[i].nb_atoms; m++) {
							if( strncmp(&uas[i].atoms[m].name[0], &prot->res[j].conf[k].atom[L].name[1], strlen(uas[i].atoms[m].name)) == 0 ) {
								prot->res[j].conf[k].atom[L].vdw_rad = types[uas[i].atoms[m].typeI].vdw_rad;
								prot->res[j].conf[k].atom[L].vdw_eps = types[uas[i].atoms[m].typeI].vdw_eps;
								break;
							}
						}
					}
				}
			}
		}
	}

        chdir("../../");
};

extern int atom2pdbline(char *line, ATOM atom);

int compareINT(const void *A, const void *B) {
	return *(int*)A - *(int*)B;
}

int compareSOBJ(const void *A, const void *B) {
        //return ( (int)( ((long double)(*(SOBJ*)A).fitness)*10000.0f - ((long double)(*(SOBJ*)B).fitness)*10000.0f) );
	if ( (*(SOBJ*)A).fitness - (*(SOBJ*)B).fitness > 0.0f ) return 1;
        else if ( fabs((*(SOBJ*)A).fitness - (*(SOBJ*)B).fitness) < 0.000000000000001 ) return 0;
	else return -1;
}

int compareFLT(const void *A, const void *B) {
	//return ( (int)( ((long double)(*(CHROMOSOME*)A).fitness)*10000.0f - ((long double)(*(CHROMOSOME*)B).fitness)*10000.0f) );
	if ( (*(CHROMOSOME*)A).fitness - (*(CHROMOSOME*)B).fitness > 0.0f ) return 1;
	else if ( fabs((*(CHROMOSOME*)A).fitness - (*(CHROMOSOME*)B).fitness) < 0.000000000000001 ) return 0;
	else return -1;
}

void GA_SETUP(PROT *prot, GA_STRUCTURE gs, int patch_nb) {

   ga.dE = gs.dE;
   ga.pop_size = gs.pop_size;
   ga.mutation_rate = gs.mutation_rate;
   ga.migration_rate = gs.migration_rate;
   ga.xover_rate = gs.xover_rate;
   ga.elitism = gs.elitism;
   ga.seed = gs.seed;
   ga.k_select = gs.k_select;
   ga.generations = gs.generations;
   ga.phase = gs.phase;
   ga.shift = gs.shift;
   ga.sphere_focus_resid = gs.sphere_focus_resid;
   ga.sphere_focus_probe_radius = gs.sphere_focus_probe_radius;
   ga.dist_center = gs.dist_center;
   ga.pop_bucket = gs.pop_bucket;
   ga.dist_center_eps = gs.dist_center_eps;
   ga.residue_check = gs.residue_check;
   ga.rand_cut_points = gs.rand_cut_points;
   ga.occupancy = gs.occupancy;

   int i, j, k, L, m, n, kr, kc;
   int n_threads;
   #pragma omp parallel
   {
	n_threads = omp_get_num_threads();
   }
   printf("GA: Using %i threads to run Sidechain Optimization Algorithm\n",n_threads);

   /* load CONFLIST1 */
   load_headlst((*prot));
   get_resnbrs((*prot), 5.0);

   get_connect12((*prot));
   printf("GA: (MCCE) Deleting H atoms...%d H atoms were deleted.\n", delete_h((*prot))); fflush(stdout);
   //assign_vdw_param((*prot));
   printf("GA: (MCCE) Assigning radii.\n"); fflush(stdout); assign_rad((*prot));
   //printf("   Estimating Solvent Accessible Surface (SAS).\n"); fflush(stdout); surfw((*prot), 1.4);
   //write_headlst((*prot));
   printf("GA: (MCCE) Done.\n\n");
  
   struct timeb start,end;

   printf("GA: (MCCE) Creating connectivity table...\n"); fflush(stdout);
   get_connect12((*prot));
   printf("GA: (MCCE) Computing self VDW potential. It may take a while...\n"); fflush(stdout);
   assign_vdw_param((*prot));
   get_vdw0_no_sas((*prot));
   get_vdw1((*prot));
   for (kr=0; kr<prot->n_res; kr++) {
      for (kc=1; kc<prot->res[kr].n_conf; kc++) {
	  prot->res[kr].conf[kc].E_torsion = torsion_conf(&prot->res[kr].conf[kc]);
          prot->res[kr].conf[kc].E_self = prot->res[kr].conf[kc].E_vdw0 + prot->res[kr].conf[kc].E_vdw1 + prot->res[kr].conf[kc].E_torsion;
      }
   }

   printf("GA: Determining which residues are part of the GA\n");
   /*establish which residues are part of the GA run (rotamers etc)*/

   /*keep track of the residue index number from the protein structure in 'ga.list'*/
   /*keep track of the number of conformers for a residue in 'ga.values'*/
   /*there is a 1-to-1 mapping between 'ga.values' and 'ga.list'*/
   /*example: ga.list[0] points to a residue int index value for identifying the residue in 'prot.res'*/
   /*and ga.values[0] points to the number of conformers for this residue*/
   ga.list = calloc(0,sizeof(int));
   ga.values = calloc(0,sizeof(int));
   ga.N = 0;
   int nrg_index = 0;
   int tot_conf = 0;
   for(i=0; i<prot->n_res; i++) {
	if( (prot->res[i].n_conf>1) && (strncmp(prot->res[i].resName,"CTR",3)!=0) && (strncmp(prot->res[i].resName,"NTR",3)!=0)) {
		ga.list = realloc(ga.list,sizeof(int)*(ga.N+1));
		ga.values = realloc(ga.values,sizeof(int)*(ga.N+1));
		ga.values[ga.N] = prot->res[i].n_conf;
		ga.list[ga.N] = i; ga.N++;
		for(j=0; j<prot->res[i].n_conf; j++) {
	 		prot->res[i].conf[j].on_atm_idx = calloc(0,sizeof(int));
			prot->res[i].conf[j].on_atoms = 0;
		}
	}
   }
  
   printf("GA: Building contiguous indexing array of turned 'on' atoms\n"); 
   /*build contiguous index of atoms so we don't have to check whether the atoms are turned on (major speedup in energy fitness function)*/
   /*can save memory if saved at the residue level but saved at the conformer level for now*/
   /*we completely ignore NTR and CTR terminals*/
   for(i=0; i<prot->n_res; i++) {
	/*need at least 2 conformers*/
	if( (prot->res[i].n_conf<2) || (strncmp(prot->res[i].resName,"CTR",3)==0) || (strncmp(prot->res[i].resName,"NTR",3)==0)) continue;
   	tot_conf += prot->res[i].n_conf-1;//not to backbone
	for(j=0; j<prot->res[i].n_conf; j++) {
		if(j>0) {prot->res[i].conf[j].nrg_idx = nrg_index; nrg_index++;}//not to backbone
		for(k=0; k<prot->res[i].conf[j].n_atom; k++) {
			if(prot->res[i].conf[j].atom[k].on) {
				prot->res[i].conf[j].on_atm_idx = realloc(prot->res[i].conf[j].on_atm_idx,sizeof(int)*(prot->res[i].conf[j].on_atoms+1));
				prot->res[i].conf[j].on_atm_idx[prot->res[i].conf[j].on_atoms] = k;
				prot->res[i].conf[j].on_atoms++;
			}
		}
	}
   }
  
   long phy_tot = sysconf(_SC_PHYS_PAGES);
   long psize = sysconf(_SC_PAGE_SIZE);
   long phy_av = sysconf(_SC_AVPHYS_PAGES);

   printf("GA: Total Conformers (without backbones):%i\n",tot_conf);
   printf("GA: Total free RAM memory: [%f] / [%f]\n",phy_av*psize/(float)1048576, phy_tot*psize/(float)1048576);

   int conf_count = 0;
   int s_total=0;

   for(i=0; i<prot->n_res; i++) {
   	if((prot->res[i].n_conf<2) || (strncmp(prot->res[i].resName,"CTR",3)==0) || (strncmp(prot->res[i].resName,"NTR",3)==0)){continue;}
        	for(j=1; j<prot->res[i].n_conf; j++) {//remember we don't compute against backbone
                	conf_count = 0;
                	for(k=i+1; k<prot->n_res; k++) {//triangular style
                        	if( (prot->res[k].n_conf<2) || (strncmp(prot->res[k].resName,"CTR",3)==0) || (strncmp(prot->res[k].resName,"NTR",3)==0)){continue;}
                                for(m=1; m<prot->res[k].n_conf; m++) {
                                        conf_count++;
                                }
                        }
			s_total += conf_count;
                }
   }
  
   printf("GA: Memory required for pre-calculation of pair-wise conformer VDW interactions: [%f] Megs.\n", (4*s_total)/(float)1048576);
	
   if((4*s_total)/(float)1048576 < phy_av*psize/(float)1048576) {
	fprintf(stdout,"GA: Pre-calculating pairwise conformer VDW interactions please wait... "); fflush(stdout);
   	double sum = 0.0f;
   	int i_idx, j_idx, tmp_idx;
	struct timeb start, end;
   	float SMIN, EPS_MIN, SIG_D2, SIG_D6, SIG_D12, dx,dy,dz,d2;
   	ATOM *atm0, *atm1;

	ga.tot_conf = tot_conf;
   	ga.conf_nrg = (float**)calloc(tot_conf,sizeof(float*));
	
	/* square array - waste of memory */
  	/*
	for(i=0; i<tot_conf; i++) {
        	ga.conf_nrg[i] = (float*)calloc(tot_conf,sizeof(float));
		for(j=0; j<tot_conf; j++) {ga.conf_nrg[i][j] = 0.0f;}
   	}
	*/
	
	/*this memory allocation algorithm properly allocates only required memory*/
	/*instead of full NxN matrix. It forms a ragged-triangular matrix layout*/
	/*and is more efficient than the NxN.  We first need to compute the total*/
	/*number of conformers for a specific pair-wise calculation and then allocate*/
	/*memory for this amount.  Then, we fill the matrix by re-mapping the 'j' index*/
	/*into a correct [0..max-1] sequence.*/
	conf_count = 0;
	s_total=0;//has a different meaning here
        for(i=0; i<prot->n_res; i++) {
                if((prot->res[i].n_conf<2) || (strncmp(prot->res[i].resName,"CTR",3)==0) || (strncmp(prot->res[i].resName,"NTR",3)==0)){continue;}
		for(j=1; j<prot->res[i].n_conf; j++) {//ignore backbone
			conf_count = 0;
			for(k=i+1; k<prot->n_res; k++) {//triangular style
				if( (prot->res[k].n_conf<2) || (strncmp(prot->res[k].resName,"CTR",3)==0) || (strncmp(prot->res[k].resName,"NTR",3)==0)){continue;}
				for(m=1; m<prot->res[k].n_conf; m++) {//ignore backbone
					conf_count++;
				}
			}
			if(conf_count>0) {ga.conf_nrg[s_total] = calloc(conf_count,sizeof(float));}
			else {ga.conf_nrg[s_total] = NULL;}
			s_total++;
		}
	}

	ftime(&start);
	/*parallelize this loop*/
	/*important speedup achieved using this parallel for loop*/
	/*scales almost linearly with # of cores*/
	for(i=0; i<prot->n_res; i++) {
   		if((prot->res[i].n_conf<2) || (strncmp(prot->res[i].resName,"CTR",3)==0) || (strncmp(prot->res[i].resName,"NTR",3)==0)){continue;}
		#pragma omp parallel for shared(i) private(k,m,j_idx,i_idx,sum,atm0,atm1,dx,dy,dz,d2,SMIN,EPS_MIN,SIG_D2,SIG_D6,SIG_D12)
		for(j=1; j<prot->res[i].n_conf; j++) {//this conformer pairwise against all other conformers of other residues, avoid backbone
			i_idx = prot->res[i].conf[j].nrg_idx;
			for(k=i+1; k<prot->n_res; k++) {//triangular matrix style
				if( (prot->res[k].n_conf<2) || (strncmp(prot->res[k].resName,"CTR",3)==0) || (strncmp(prot->res[k].resName,"NTR",3)==0)){continue;}
				for(m=1; m<prot->res[k].n_conf; m++) {//against this residue's conformer 'm', avoid backbone
					j_idx = tot_conf - prot->res[k].conf[m].nrg_idx - 1;//need to re-map the j index properly from [0..max-1]
					sum = 0.0f;
					for(L=0; L<prot->res[i].conf[j].on_atoms; L++) {
						atm0 = &prot->res[i].conf[j].atom[prot->res[i].conf[j].on_atm_idx[L]];
						for(n=0; n<prot->res[k].conf[m].on_atoms; n++) {
							atm1 = &prot->res[k].conf[m].atom[prot->res[k].conf[m].on_atm_idx[n]];
							dx = atm1->xyz.x- atm0->xyz.x; dy = atm1->xyz.y- atm0->xyz.y; dz = atm1->xyz.z- atm0->xyz.z;
							d2 = dx*dx+dy*dy+dz*dz;
							/*Truncated LJ12-6 potential*/
							if (__builtin_expect(d2 > ccutoff_far2,0)) {sum += VDW_ELIMIT_FAR; continue;}
							if (__builtin_expect(d2 < ccutoff_near2,0)) {sum += VDW_ELIMIT_NEAR; continue;}
							SMIN = atm0->vdw_rad + atm1->vdw_rad;
							EPS_MIN = sqrt(atm0->vdw_eps*atm1->vdw_eps);
							SIG_D2 = SMIN*SMIN / d2; SIG_D6 = SIG_D2*SIG_D2*SIG_D2; SIG_D12 = SIG_D6*SIG_D6;
							sum += EPS_MIN*(SIG_D12 - 2.*SIG_D6);
						}
	
					}
					ga.conf_nrg[i_idx][j_idx] = sum;
				}
			}
		}
   	}
	ftime(&end);
	printf("%6.2f seconds to pre-calculate interactions.\n",(1000.0*(end.time-start.time)+(end.millitm-start.millitm))/1000);
  	//update_VDW_RADEPS_UA(prot);/*United Atom Hydrogen for C1 C2 and C3 heavy atoms.*/
	RUN_GA(prot, patch_nb, 1);
	return;
   }
   else {
	printf("GA: Not enough memory! Calculating pairwise conformer VDW on-the-fly (slower performance)\n");
   	//update_VDW_RADEPS_UA(prot);/*United Atom Hydrogen for C1 C2 and C3 heavy atoms.*/
   	RUN_GA(prot, patch_nb, 0);
	return;
   }
};

/*1 to M+1*/
__inline int random_number_1_M(int M) {
	return (int)((double)rand()/ ( (double)RAND_MAX + (double)(1)) * M)+1;
}

/*0 to M*/
__inline int random_number_0_M(int M) {
        return (int)((double)rand()/ ( (double)RAND_MAX + (double)(1)) * M);
}

/*0.0 to < 1.0*/
__inline float random_number_0_1(void) {
        return (double)rand() / ( (double)RAND_MAX + (double)(1));
}

void initialize_populations(CHROMOSOME *pop, CHROMOSOME *new_pop) {
	int i,j;
        for(i=0; i<ga.pop_size; i++) {
                pop[i].genes = calloc(ga.N,sizeof(int));
		pop[i].fitness = 0.0f;
		pop[i].flag = 0;
		pop[i].mated = 0;
		pop[i].copied = 0;
		pop[i].pw_vdw = malloc(ga.N*sizeof(VDWCONF));
	
		new_pop[i].genes = calloc(ga.N,sizeof(int));
		new_pop[i].fitness = 0.0f;
		new_pop[i].flag = 0;
		new_pop[i].mated = 0;
		new_pop[i].copied = 0;
		new_pop[i].pw_vdw = malloc(ga.N*sizeof(VDWCONF));

                for(j=0; j<ga.N; j++) {
                        pop[i].genes[j] = random_number_1_M(ga.values[j]-1);
			pop[i].pw_vdw[j].vdw = 0.0;
			pop[i].pw_vdw[j].min = FLT_MAX;
		}
        }
};

void evaluate_population2(PROT prot, CHROMOSOME *pop) {
	int i,j,k,m,i_idx,j_idx,tmp_idx,p;
	#pragma omp parallel for private(i, j, k, m, i_idx, j_idx, tmp_idx)
	for(p=0; p<ga.pop_size; p++) {
		if(pop[p].flag==1) {continue;}/*skip the ones which we already know the fitness*/
		pop[p].fitness = 0;
		for(i=0; i<ga.N; i++) {/*for all residues in chromosome*/
			k = pop[p].genes[i];/*conformer*/
			pop[p].fitness += prot.res[ga.list[i]].conf[k].E_self;
			pop[p].pw_vdw[i].vdw = prot.res[ga.list[i]].conf[k].E_self;
			for(j=i+1; j<ga.N; j++) {/*for all other residues in chromosome*/
				m = pop[p].genes[j];/*conformer*/
				//avoid conformers from the same residue - this can happen with MCCE structure
				if(prot.res[ga.list[i]].resSeq == prot.res[ga.list[j]].resSeq){continue;}
				i_idx = prot.res[ga.list[i]].conf[k].nrg_idx;
				j_idx = ga.tot_conf - prot.res[ga.list[j]].conf[m].nrg_idx - 1;//re-mapping to [0..max-1]

				pop[p].fitness += ga.conf_nrg[i_idx][j_idx];
				pop[p].pw_vdw[i].vdw += ga.conf_nrg[i_idx][j_idx];
				pop[p].pw_vdw[j].vdw += ga.conf_nrg[i_idx][j_idx];
			}
		}
	}
	return;
}

/*openMP parallelized for multi-core CPUs*/
/*algorithm without pre-calculated pairwise VDW interactions*/
/*very slow compared to the other version*/
void evaluate_population(PROT prot, CHROMOSOME *pop) {
	float SMIN, EPS_MIN, SIG_D2, SIG_D6, SIG_D12, dx,dy,dz,d2, vdw;
	int i,j,k,L,m,n,p,r;
	ATOM *atm0, *atm1;

	#pragma omp parallel for private(atm0, atm1, SMIN, EPS_MIN, SIG_D2, SIG_D6, SIG_D12, i, j, k, L, m, n, dx, dy, dz, d2, vdw)
	for(p=0; p<ga.pop_size; p++) {
		if(pop[p].flag==1) {continue;}/*skip the ones which we already know the fitness*/
		pop[p].fitness = 0;
		for(i=0; i<ga.N; i++) {/*for all residues in chromosome*/
			k = pop[p].genes[i];/*conformer*/
			pop[p].fitness += prot.res[ga.list[i]].conf[k].E_self;
			pop[p].pw_vdw[i].vdw += prot.res[ga.list[i]].conf[k].E_self;
			for(j=i+1; j<ga.N; j++) {/*for all other residues in chromosome*/
				m = pop[p].genes[j];/*conformer*/
                        	if(prot.res[ga.list[i]].resSeq == prot.res[ga.list[j]].resSeq) continue;
				for(L=0; L<prot.res[ga.list[i]].conf[k].on_atoms; L++) {/*for all turned on atoms of conformer 'k'*/
					atm0 = &prot.res[ga.list[i]].conf[k].atom[prot.res[ga.list[i]].conf[k].on_atm_idx[L]];
					for(n=0; n<prot.res[ga.list[j]].conf[m].on_atoms; n++) {/*for all turned on atoms of conformer 'm'*/
						atm1 = &prot.res[ga.list[j]].conf[m].atom[prot.res[ga.list[j]].conf[m].on_atm_idx[n]];
						dx = atm1->xyz.x- atm0->xyz.x; dy = atm1->xyz.y- atm0->xyz.y; dz = atm1->xyz.z- atm0->xyz.z;	
						d2 = dx*dx+dy*dy+dz*dz;
						/*Truncated LJ12-6 potential*/
						if (__builtin_expect(d2 > ccutoff_far2,0)) {
							pop[p].fitness += VDW_ELIMIT_FAR;
							pop[p].pw_vdw[i].vdw += VDW_ELIMIT_FAR;
							pop[p].pw_vdw[j].vdw += VDW_ELIMIT_FAR;
							continue;
						}
    						if (__builtin_expect(d2 < ccutoff_near2,0)) {
							pop[p].fitness += VDW_ELIMIT_NEAR;
							pop[p].pw_vdw[i].vdw += VDW_ELIMIT_NEAR;
							pop[p].pw_vdw[j].vdw += VDW_ELIMIT_NEAR;
							continue;
						}
						SMIN = atm0->vdw_rad + atm1->vdw_rad;
						EPS_MIN = sqrt(atm0->vdw_eps*atm1->vdw_eps);
						SIG_D2 = SMIN*SMIN / d2; SIG_D6 = SIG_D2*SIG_D2*SIG_D2; SIG_D12 = SIG_D6*SIG_D6;
						vdw = EPS_MIN*(SIG_D12 - 2.*SIG_D6);
						pop[p].fitness += vdw;
						pop[p].pw_vdw[i].vdw += vdw;
						pop[p].pw_vdw[j].vdw += vdw;
					}
                		}
        		}
   		}
	}
	return;
};

void evaluate_ind(PROT prot, CHROMOSOME *ind) {
        float SMIN, EPS_MIN, SIG_D2, SIG_D6, SIG_D12, dx,dy,dz,d2,vdw;
        int i,j,k,L,m,n,p;
        ATOM *atm0, *atm1;
        ind->fitness = 0;
	//ind->flag = 1;
        for(i=0; i<ga.N; i++) {/*for all residues in chromosome*/
        	k = ind->genes[i];/*conformer*/
                ind->fitness += prot.res[ga.list[i]].conf[k].E_self;
		ind->pw_vdw[i].vdw += prot.res[ga.list[i]].conf[k].E_self;
		for(j=i+1; j<ga.N; j++) {/*for all other residues in chromosome*/
                	m = ind->genes[j];/*conformer*/
                        if(prot.res[ga.list[i]].resSeq == prot.res[ga.list[j]].resSeq) continue;
			for(L=0; L<prot.res[ga.list[i]].conf[k].n_atom; L++) {/*for all atoms of conformer 'k'*/
				atm0 = &prot.res[ga.list[i]].conf[k].atom[L];
                                if(!atm0->on) continue;
				for(n=0; n<prot.res[ga.list[j]].conf[m].n_atom; n++) {/*for all atoms of conformer 'm'*/
                                	atm1 = &prot.res[ga.list[j]].conf[m].atom[n];
                                        if(!atm1->on) continue;
					dx = atm1->xyz.x- atm0->xyz.x; dy = atm1->xyz.y- atm0->xyz.y; dz = atm1->xyz.z- atm0->xyz.z;
                                        d2 = dx*dx+dy*dy+dz*dz;
					/*Truncated LJ12-6 potential*/
                                        if (d2 > ccutoff_far2) {
						vdw = VDW_ELIMIT_FAR;
						ind->fitness += vdw;
						ind->pw_vdw[i].vdw += vdw;
						ind->pw_vdw[j].vdw += vdw;
					} /* Cutoff */
                                        else if (d2 < ccutoff_near2) {
						vdw = VDW_ELIMIT_NEAR;
						ind->fitness += vdw;
						ind->pw_vdw[i].vdw += vdw;
						ind->pw_vdw[j].vdw += vdw;
					} /* Cutoff */
                                        else {
                                        	SMIN = atm0->vdw_rad + atm1->vdw_rad;
                                                EPS_MIN = sqrt(atm1->vdw_eps*atm0->vdw_eps);//ga.eps_atm[atm0->index][atm1->index];
                                                SIG_D2 = SMIN*SMIN / d2; SIG_D6 = SIG_D2*SIG_D2*SIG_D2; SIG_D12 = SIG_D6*SIG_D6;
						vdw = EPS_MIN*(SIG_D12 - 2.*SIG_D6);
						ind->fitness += vdw;
						ind->pw_vdw[i].vdw += vdw;
						ind->pw_vdw[j].vdw += vdw;
					}
                                }
                        }
                }
        }
};

__inline void mutate_population(CHROMOSOME *new_pop) {
	/*don't mutate genes that have values[i] == 2; it means that only the initial state conformer is present for this residue*/
	int i, gene, new_value;
	/*don't mutate the elite ones*/
	int start = 1;//(int)(ga.pop_size*ga.elitism);
	for(i=start; i<ga.pop_size; i++) {
		if(random_number_0_1() < ga.mutation_rate) {
			do {
				gene = random_number_0_M(ga.N);
			} while(ga.values[gene]==2);//check to ensure mutation selects a residue that can actually be mutated
			do {
				new_value = random_number_1_M(ga.values[gene]-1);
			} while(new_value==new_pop[i].genes[gene]);//check to ensure mutation generates a different value than the current one
			new_pop[i].genes[gene] = new_value;
			new_pop[i].flag = 0; new_pop[i].fitness = 0.0f;
		}
	}
};

/*selection pressure used for driving convergence of GA algorithm*/
__inline int Boltzmann_Pressure_Converge(CHROMOSOME *pop, float dE, float E_min) {
	int parent, parent2, parent3;
	float roll;
	parent = random_number_0_M(ga.pop_size);
	do {
		parent2 = random_number_0_M(ga.pop_size);
		
		roll = (-dE)*log(random_number_0_1());
		if( pop[parent].fitness > pop[parent2].fitness ) parent = parent2;
		
		if( fabs(pop[parent].fitness-E_min) < roll ) {
			return parent;
		}
	} while(1);
};

/*selection pressure used for driving evolutionary sampling AFTER convergence of GA*/
/*pay careful attention; this is NOT the same as the above one*/
/*Here, we randomly select individuals falling below the fixed bar and apply*/
/*selective pressure only above the fixed bar*/
__inline int Boltzmann_Pressure_Converge2(CHROMOSOME *pop, float dE, float E_min) {
        int parent, parent2, parent3;
        float roll;
        do {
                parent2 = random_number_0_M(ga.pop_size);

                roll = (-dE)*log(random_number_0_1());
                
                if(fabs(pop[parent2].fitness-E_min) < ga.dist_center) return parent2;
                else if( fabs(pop[parent2].fitness-E_min) < roll ) {
                	return parent2;
                }
        } while(1);
};

/*this is the one used in evolutionary sampling AFTER convergence of GA*/
/*pay attention to the signs of inequality for divergence*/
__inline int Boltzmann_Pressure_Diverge2(CHROMOSOME *pop, float dE, float E_min) {
        int parent, parent2;
        float roll;
        do {
                parent2 = random_number_0_M(ga.pop_size);

                roll = (-dE)*log(random_number_0_1());

                if(fabs(pop[parent2].fitness-E_min) > ga.dist_center) return parent2;
                else if( fabs(pop[parent2].fitness-E_min) > roll ) {
               		return parent2;
                }
        } while(1);
};

void elitism(int *new_pop_count, CHROMOSOME *pop, CHROMOSOME *new_pop) {
	int i, IDX;
	int nb_elite = 1;//hard-coded as THE best
	//int nb_elite = (int)(ga.elitism*ga.pop_size);
	int count=0;
	
	SOBJ *tmp_pop = malloc(sizeof(SOBJ)*ga.pop_size);
	for(i=0; i<ga.pop_size; i++) {
		tmp_pop[i].index = i;
		tmp_pop[i].fitness = pop[i].fitness;
	}
	qsort(tmp_pop, ga.pop_size, sizeof(SOBJ), compareSOBJ);
	
	i=0;
	count=0;
	IDX = 0;
	while(i<nb_elite) {
		if(i>0) {
			if ( genetic_difference(pop[tmp_pop[IDX].index].genes, pop[tmp_pop[count].index].genes) != 0  ) {
				memmove(new_pop[*new_pop_count].genes,pop[tmp_pop[count].index].genes,sizeof(int)*ga.N);
				memmove(new_pop[*new_pop_count].pw_vdw,pop[tmp_pop[count].index].pw_vdw,sizeof(VDWCONF)*ga.N);
				new_pop[*new_pop_count].flag = 1;/*even though we don't want to re-compute this solution's energy, we still want to allow this best solution to mate*/
                        	new_pop[*new_pop_count].fitness = pop[tmp_pop[count].index].fitness;
                        	(*new_pop_count)++;
                        	pop[tmp_pop[count].index].copied = 1;/*don't allow for migration but only crossover*/
                        	IDX = count;
				count++;
                        	i++;
			}
			else {
				count++;
			}
		}
		else {
			memmove(new_pop[*new_pop_count].genes,pop[tmp_pop[count].index].genes,sizeof(int)*ga.N);
			memmove(new_pop[*new_pop_count].pw_vdw,pop[tmp_pop[count].index].pw_vdw,sizeof(VDWCONF)*ga.N);
			new_pop[*new_pop_count].flag = 1;/*even though we don't want to re-compute this solution's energy, we still want to allow this best solution to mate*/
        		new_pop[*new_pop_count].fitness = pop[tmp_pop[count].index].fitness;
        		(*new_pop_count)++;
        		pop[tmp_pop[count].index].copied = 1;/*don't allow for migration but only crossover*/
			IDX = count;
			count++;
			i++;
		}
	}
	free(tmp_pop);
};

/*main one*/
inline int condition(CHROMOSOME *pop, int parent1, int parent2) {
	return ( ((pop[parent1].mated==1) && (pop[parent2].mated==1)) || (parent1 == parent2) );
}

/*example of another condition ga operator*/
inline int condition3(CHROMOSOME *pop, int parent1, int parent2) {
        return ( (parent1 == parent2) );
}

int more_to_insert(int *cuts, int c) {
	int i;
	int sum=0;
	for(i=0; i<(c-1); i++) {
		sum += abs(cuts[i]-cuts[i+1]-1);
	}
	if(sum < ga.N/2) return 1;
	return 0;
};

void Random_N_Point_Crossover(CHROMOSOME *pop, CHROMOSOME *new_pop, int *new_pop_count, int parent1_idx, int parent2_idx, int max_cut) {
	/*establish number of N-points for the crossover*/
	int nb_cuts = random_number_1_M(max_cut);
	int *cuts = malloc(nb_cuts*sizeof(int));/*each cut index of the chromosome*/
	int ins_cuts=0;/*number of currently inserted cut points*/
	int i, j;
	int tmp_cut;
	int good, p_sel;
	int parents[2];
	parents[0] = parent1_idx; parents[1] = parent2_idx;
	/*each cut needs to have at least 1 gene spaced out with the previous closest and next closest cut*/
	while(ins_cuts<nb_cuts && more_to_insert(cuts, ins_cuts)) {
		tmp_cut = random_number_0_M(ga.N);
		good = 1;
		/*check against all other cuts and find whether the closest larger and the closest smaller have at least distance of 2 between each other*/
		for(i=0; i<ins_cuts; i++) {
			if( !(abs(cuts[i] - tmp_cut) >= 2) ) {
				good = 0; break;
			}
		}
		if(good) {
			cuts = realloc(cuts,sizeof(int)*(ins_cuts+1));
			cuts[ins_cuts] = tmp_cut; ins_cuts++;
			qsort(cuts, ins_cuts, sizeof(int), compareINT);
		}
	}
	
	/*proceed on with the crossover for offspring1*/
	p_sel = 0;
	for(i=0; i<ins_cuts; i++) {
		if(i==0) {
			memmove(&new_pop[*new_pop_count].genes[0],&pop[parents[p_sel]].genes[0],sizeof(int)*(cuts[i]));
			memmove(&new_pop[*new_pop_count].pw_vdw[0],&pop[parents[p_sel]].pw_vdw[0],sizeof(VDWCONF)*(cuts[i]));
		}
		else /*if(i<(ins_cuts-1))*/ {
			memmove(&new_pop[*new_pop_count].genes[cuts[i-1]],&pop[parents[p_sel]].genes[cuts[i-1]],sizeof(int)*(cuts[i]-cuts[i-1]));
			memmove(&new_pop[*new_pop_count].pw_vdw[cuts[i-1]],&pop[parents[p_sel]].pw_vdw[cuts[i-1]],sizeof(VDWCONF)*(cuts[i]-cuts[i-1]));
		}
		p_sel ^= 1;/*swap from parent 1 to parent 2*/
	}
	memmove(&new_pop[*new_pop_count].genes[cuts[i-1]],&pop[parents[p_sel]].genes[cuts[i-1]],sizeof(int)*(ga.N-cuts[i-1]));
	memmove(&new_pop[*new_pop_count].pw_vdw[cuts[i-1]],&pop[parents[p_sel]].pw_vdw[cuts[i-1]],sizeof(VDWCONF)*(ga.N-cuts[i-1]));
	new_pop[*new_pop_count].flag = 0;
	(*new_pop_count)++;
	pop[parent1_idx].mated = 1;
	pop[parent2_idx].mated = 1;

	if((*new_pop_count)<ga.pop_size) {
        	/*proceed on with the crossover for offspring2*/
        	p_sel = 1;
        	for(i=0; i<ins_cuts; i++) {
                	if(i==0) {
                        	memmove(&new_pop[*new_pop_count].genes[0],&pop[parents[p_sel]].genes[0],sizeof(int)*(cuts[i]));
                		memmove(&new_pop[*new_pop_count].pw_vdw[0],&pop[parents[p_sel]].pw_vdw[0],sizeof(VDWCONF)*(cuts[i]));
			}
                	else {
                        	memmove(&new_pop[*new_pop_count].genes[cuts[i-1]],&pop[parents[p_sel]].genes[cuts[i-1]],sizeof(int)*(cuts[i]-cuts[i-1]));
                		memmove(&new_pop[*new_pop_count].pw_vdw[cuts[i-1]],&pop[parents[p_sel]].pw_vdw[cuts[i-1]],sizeof(VDWCONF)*(cuts[i]-cuts[i-1]));
			}
                	p_sel ^= 1;/*swap from parent 1 to parent 2*/
        	}
		memmove(&new_pop[*new_pop_count].genes[cuts[i-1]],&pop[parents[p_sel]].genes[cuts[i-1]],sizeof(int)*(ga.N-cuts[i-1]));
		memmove(&new_pop[*new_pop_count].pw_vdw[cuts[i-1]],&pop[parents[p_sel]].pw_vdw[cuts[i-1]],sizeof(VDWCONF)*(ga.N-cuts[i-1]));
        	new_pop[*new_pop_count].flag = 0;
        	(*new_pop_count)++;
	}
	free(cuts);
};

__inline void copy_population_over(CHROMOSOME *new_population, CHROMOSOME *population) {
	int i;
	#pragma omp parallel for
	for(i=0; i<ga.pop_size; i++) {
        	memmove(population[i].genes, new_population[i].genes,sizeof(int)*ga.N);
		memmove(population[i].pw_vdw, new_population[i].pw_vdw,sizeof(VDWCONF)*ga.N);
                population[i].flag = 0;
                population[i].mated = 0;
		population[i].copied = 0;
		population[i].fitness = new_population[i].fitness;
	}
};

int compute_nb_rotamers(PROT prot, CHROMOSOME *pop) {
	int i,j,p,count, **res;
	res = calloc(prot.n_res,sizeof(int*));
	count=0;
	for(i=0; i<prot.n_res; i++) {
		res[i] = calloc(prot.res[i].n_conf,sizeof(int));
		for(j=0; j<prot.res[i].n_conf; j++) {
			res[i][j] = 0;
		}
	}

	for(p=0; p<ga.pop_size; p++) {
		for(i=0; i<ga.N; i++) {
			if(res[ga.list[i]][pop[p].genes[i]] !=1) {
				res[ga.list[i]][pop[p].genes[i]] = 1;
				count++;
			}
		}
	}
	for(i=0; i<prot.n_res; i++) {
		free(res[i]);
	}
	free(res);
	return count;
};

int final_rotamers(PROT prot, CHROMOSOME *final_pop, int fcount) {
        int p;
        int i,j;
        int **res;
        res = calloc(prot.n_res,sizeof(int*));
        int count=0;
	for(i=0; i<prot.n_res; i++) {
                res[i] = calloc(prot.res[i].n_conf,sizeof(int));
                for(j=0; j<prot.res[i].n_conf; j++) {
                        res[i][j] = 0;
                }
        }

        for(p=0; p<fcount; p++) {
                for(i=0; i<ga.N; i++) {
			res[ga.list[i]][final_pop[p].genes[i]]++;
                }
        }
	for(i=0; i<prot.n_res; i++) {
		for(j=0; j<prot.res[i].n_conf; j++) {
			if(res[i][j] > 0) {count++;}
		}
	}
        for(i=0; i<prot.n_res; i++) {
                free(res[i]);
        }
        free(res);
        return count;
};

__inline float compute_population_average(CHROMOSOME *pop) {
	int p;
	double average = 0.0f;
	#pragma omp parallel for private(p) reduction(+:average)
	for(p=0; p<ga.pop_size; p++) {
		average += pop[p].fitness;
	}
	return (average/ga.pop_size);
}

void write_backbone_atoms(FILE *crd, PROT prot) {
        int i, j, k;
        for (i=0; i<prot.n_res; i++) {
                for (j=0; j<1; j++) {
                        for (k=0; k<prot.res[i].conf[j].n_atom; k++) {
                                if (!(prot.res[i].conf[j].atom[k].on)) continue;
				fprintf(crd,"%4.8f %4.8f %4.8f\n",prot.res[i].conf[j].atom[k].xyz.x, prot.res[i].conf[j].atom[k].xyz.y, prot.res[i].conf[j].atom[k].xyz.z);
                        }
                }
        }
};

__inline int genetic_difference(int *genes1, int *genes2) {
	int i;
	for(i=0; i<ga.N; i++) {
		if(genes1[i]!=genes2[i]) {return 1;}
	}
	return 0;
}

__inline int genetic_difference_stat(int *genes1, int *genes2) {
        int diff=0;
        int i;
        for(i=0; i<ga.N; i++) {
                if(genes1[i]!=genes2[i]) {diff += 1;}
        }
        return diff;
}
/*
void write_best_solution_as_pdb(PROT prot, CHROMOSOME *new_population, CHROMOSOME *pop, CHROMOSOME *original, CHROMOSOME *final_population, int final_pop_count, float KT) {
	int i, j, k, iConf, c;
        FILE* stream;

        stream = fopen("GA_BEST_out.pdb","wt");
        c = 0;
        for (i=0; i<prot.n_res; i++) {
                iConf = 0;
                for (j=0; j<1; j++) {
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
                }
        }
        
	for (i=0; i<ga.N; i++) {
                iConf = 1;
                for (k=0; k<prot.res[ga.list[i]].conf[new_population[0].genes[i]].n_atom; k++) {
                        if (!(prot.res[ga.list[i]].conf[new_population[0].genes[i]].atom[k].on)) continue;
                        if (c<99999) c++;
                        fprintf(stream, "ATOM  %5d %4s%c%3s %c%04d%c%03d%8.3f%8.3f%8.3f %7.3f      %6.3f      %-11s\n",
                        c, prot.res[ga.list[i]].conf[new_population[0].genes[i]].atom[k].name,
                        prot.res[ga.list[i]].conf[new_population[0].genes[i]].altLoc,
                        prot.res[ga.list[i]].resName,
                        prot.res[ga.list[i]].chainID,
                        prot.res[ga.list[i]].resSeq,
                        prot.res[ga.list[i]].iCode,
                        iConf,
                        prot.res[ga.list[i]].conf[new_population[0].genes[i]].atom[k].xyz.x,
                        prot.res[ga.list[i]].conf[new_population[0].genes[i]].atom[k].xyz.y,
                        prot.res[ga.list[i]].conf[new_population[0].genes[i]].atom[k].xyz.z,
                        prot.res[ga.list[i]].conf[new_population[0].genes[i]].atom[k].rad,
                        prot.res[ga.list[i]].conf[new_population[0].genes[i]].atom[k].crg,
                        prot.res[ga.list[i]].conf[new_population[0].genes[i]].history);
                }
                iConf++;
        }
        fclose(stream);
	
	//write Chromosomes as .CDR coordinate file
	FILE* crd;
	crd = fopen("GA_CRD.crd","wt");
	int nb = ga.pop_size;
	qsort(pop, ga.pop_size, sizeof(CHROMOSOME), compareFLT);
	ATOM *atm;
	//skip the first one as it will be the best one which was already outputted as .pdb file
	fprintf(crd,"rdparm modified trajectory\n");
	for(i=1; i<nb; i++) {
		write_backbone_atoms(crd,prot);
		for(j=0; j<ga.N; j++) {
			for(k=0; k<prot.res[ga.list[j]].conf[pop[i].genes[j]].on_atoms; k++) {
				atm = &prot.res[ga.list[j]].conf[pop[i].genes[j]].atom[prot.res[ga.list[j]].conf[pop[i].genes[j]].on_atm_idx[k]];
				fprintf(crd,"%4.8f %4.8f %4.8f\n", atm->xyz.x, atm->xyz.y, atm->xyz.z);
			}
		}
	}
	fclose(crd);
	printf("CRD Coordinate file written.\n");	

	crd = fopen("Fitness_GA_FULL.txt","wt");
	for(i=0; i<nb; i++) {
		fprintf(crd,"[%04i] Orig. Fit: %15.7f Best Fit: %15.7f Current Fit: %15.7f Genetic Difference: %5.2f\n", i, original->fitness, new_population[0].fitness, fabs(pop[i].fitness-new_population[0].fitness), genetic_difference_stat(new_population[0].genes, pop[i].genes)/(float)ga.N*100.0f);
	}
	fclose(crd);
	printf("Fitness File written.\n\n");

        float min, max;

	printf("Generating GA's Final Population Histogram JPEG and DATA...\n");
	//init ouput to jpeg file
	plsdev("jpeg");
	plsfnam("GA_Histogram.jpeg");
	
	plscol0(0, 255,255,255);
	plinit();
	//generate histogram

	PLINT n = ga.pop_size;

	PLFLT *data = calloc(n, sizeof(PLFLT));
	FILE* out;
	out = fopen("GA_Histogram.txt","wt");
	fprintf(out,"%i\n",n);
	for(i=0; i<n; i++) {
		data[i] = fabs(pop[i].fitness-new_population[0].fitness);
		//if(data[i] > 100.0f) {n=i; break;}
		fprintf(out,"%f\n",data[i]);
	}
	fclose(out);

	PLFLT datmin,datmax;
	PLINT nbin;
	PLINT opt;
        min = FLT_MAX; max = FLT_MIN;
        for(i=0; i<n; i++) {
                if(data[i]>max) {max=data[i];}
                if(data[i]<min) {min=data[i];}
        }
        datmin = min; datmax = max;

	//datmax = 5.0f;
	
        nbin = 70;
	plcol0(8);
	plhist (n, data , datmin , datmax , nbin , 0 );
	plcol0(1);
	pllab("#frEnergy Difference from Min", "#frFrequency", "GA Population Histogram");

	//generate histogram
	plend();//exit call for plplot
	printf("done.\n");
	free(data);
	
	printf("Generating Bucket Population's Histogram within %f Kcal/mol from %f JPEG...\n", ga.dist_center, new_population[0].fitness);
	qsort(final_population, final_pop_count, sizeof(CHROMOSOME), compareFLT);
        plsdev("jpeg");
        plsfnam("BUCKET_Histogram.jpeg");
        
	plscol0(0, 255,255,255);
	plinit();
	
	out = fopen("BUCKET_Histogram.txt","wt");
	n = final_pop_count;
	data = calloc(n,sizeof(PLFLT));
	for(i=0; i<n; i++) {
		data[i] = fabs(final_population[i].fitness-new_population[0].fitness);
		if(data[i] > 100.0f) {n=i; break;}
		fprintf(out,"%f\n",data[i]);
	}
	fclose(out);
	
        min = FLT_MAX; max = FLT_MIN;
        for(i=0; i<n; i++) {
                if(data[i]>max) {max=data[i];}
                if(data[i]<min) {min=data[i];}
        }
        datmin = min; datmax = max;
	nbin=70;
	plcol0(8);
	plhist (n, data , datmin , datmax , nbin , 0 );
	plcol0(1);
	pllab("#frEnergy Difference from Min", "#frFrequency", "Bucket Population Histogram within KT of the MIN");
	plend();
	printf("done.\n");
	free(data);
	
	mkdir("./Histograms",S_IRWXU);
	chdir("./Histograms");
	char name[64];
	float *tmp;
	nbin = 200;
	for(i=0; i<ga.N; i++) {
		min = FLT_MAX;
		max = FLT_MIN;
		plsdev("jpeg");
		sprintf(name,"Res_%s%04i.jpeg",prot.res[ga.list[i]].resName, ga.list[i]);
		plsfnam(name);
		plscol0(0, 255,255,255);
		plinit();
		tmp = calloc(final_pop_count,sizeof(float));
		data = calloc(final_pop_count,sizeof(PLFLT));
		for(j=0; j<final_pop_count; j++) {
			tmp[j] = fabs(final_population[j].pw_vdw[i].vdw - final_population[j].pw_vdw[i].min);
		}
		
		for(j=0; j<final_pop_count; j++) {
			data[j] = (PLFLT)tmp[j];
			if(tmp[j]>max) {max=tmp[j];}
			if(tmp[j]<min) {min=tmp[j];}
		}
	        
		datmin = min;
	        datmax = max;
	        if( fabs(datmin-datmax) > 0.00001 ) {
			plcol0(8);
        		plhist (final_pop_count, data , datmin , datmax , nbin , 0 );
        		plcol0(1);
			sprintf(name,"%s%04i Min.E:%0.3f Fit:%0.3f",prot.res[ga.list[i]].resName, ga.list[i], final_population[0].pw_vdw[i].min, final_population[0].pw_vdw[i].occ_fit);
			pllab("#frEnergy Difference from Min", "#frFrequency", name);
		}
		plend();
        	free(data);
		free(tmp);
	}
	chdir("../");
};
*/

__inline float deltaE(CHROMOSOME *pop, float E_min) {
	int i;
	double d_e=0.0f;
	int max = (int)(ga.pop_size*ga.dE);/*IMPORTANT PARAMETER!! ONLY 70% of pop considered for delta_E (empirically determined by optimal)*/
	SOBJ *tmp_pop = malloc(sizeof(SOBJ)*ga.pop_size);

	#pragma omp parallel for	
	for(i=0; i<ga.pop_size; i++) {
		tmp_pop[i].index = i; tmp_pop[i].fitness = pop[i].fitness;
	}
	qsort(tmp_pop, ga.pop_size, sizeof(SOBJ), compareSOBJ);

	#pragma omp parallel for reduction(+:d_e)	
	for(i=0; i<max; i++) {
		d_e += (E_min - tmp_pop[i].fitness);
	}

	free(tmp_pop);
	return fabs(d_e/max);
};

__inline float find_Emin(CHROMOSOME *pop) {
	int i;
	float min = FLT_MAX;
	for(i=0; i<ga.pop_size; i++) {
		if (pop[i].fitness < min) {
			min = pop[i].fitness;
		}
	}
	return min;
};

void similarity_check_slow(CHROMOSOME *new_pop) {
	int i, j, k, new_value;
	for(i=0; i<ga.pop_size; i++) {
		#pragma omp parallel for private(i,j,k,new_value)
		for(j=i+1; j<ga.pop_size; j++) {
			if( genetic_difference(new_pop[i].genes, new_pop[j].genes) == 0)  {
				for(k=0; k<ga.N; k++) {
                                       	if(ga.values[k]==2) continue;/*if can't mutate this residue because it only have 1 state (its initial conformation from xray structure), then continue*/
                                       	do {
                                               	new_value = random_number_1_M(ga.values[k]-1);
                                       	} while(new_value==new_pop[j].genes[k]);/*check to ensure mutation generates a different value than the current one*/
                                       	new_pop[j].genes[k] = new_value;
                        	}
                        	new_pop[j].flag = 0; new_pop[j].fitness = 0.0f;
			}
		}
	}
};

void similarity_check_fast(CHROMOSOME *new_pop) {
        int i, j, k;
        int new_value;
	/*this algorithm is faster than the other one*/
	/*first flag in parallel those chromosomes that are similar using critical section*/
	/*then in parallel, mutate 100% those that need to be changed*/
	#pragma omp parallel for private(i,j)
        for(i=0; i<ga.pop_size; i++) {
                for(j=i+1; j<ga.pop_size; j++) {
                        if( genetic_difference(new_pop[i].genes, new_pop[j].genes) == 0)  {
                                #pragma omp critical
				{
				new_pop[j].flag = 5; new_pop[j].fitness = 0.0f;
				}
                        }
                }
        }
	
	#pragma omp parallel for private(i,j,k,new_value)
	for(i=0; i<ga.pop_size; i++) {
		if(new_pop[i].flag==5) {
			for(k=0; k<ga.N; k++) {
                        	if(ga.values[k]==2) continue;/*if can't mutate this residue because it only have 1 state (its initial conformation from xray structure), then continue*/
                                do {
                                	new_value = random_number_1_M(ga.values[k]-1);
                                } while(new_value==new_pop[i].genes[k]);/*check to ensure mutation generates a different value than the current one*/
                                new_pop[i].genes[k] = new_value;
                        }
                        new_pop[i].flag = 0;//reset flag to proper value
		}
	}
};

__inline void clear_sumvdw(CHROMOSOME *pop) {
	int i,j;
	#pragma omp parallel for private(i,j)
	for(i=0; i<ga.pop_size; i++) {
		if(pop[i].flag==1) {continue;}
		for(j=0; j<ga.N; j++) {
			pop[i].pw_vdw[j].vdw = 0.0;
		}
	}
};

__inline void set_vdwmin(CHROMOSOME *pop) {
	int i, j;
	float min;
	for(j=0; j<ga.N; j++) {
		min = FLT_MAX;
		for(i=0; i<ga.pop_size; i++) {
			if(pop[i].pw_vdw[j].vdw < min) {
				min = pop[i].pw_vdw[j].vdw;
			}
		}
		for(i=0; i<ga.pop_size; i++) {
			if(min < pop[i].pw_vdw[j].min) {
				pop[i].pw_vdw[j].min = min;
				pop[i].pw_vdw[j].occ_fit = pop[i].fitness;
			}
		}
	}
};


void RUN_GA(PROT *prot, int patch_nb, int eval_p) {
	
	int i,j,k,parent1_idx,parent2_idx;
	int new_pop_count;
	struct timeb start, end, start2, end2;
	double elapsed, elap;
	int nb_rotamers;
	float pop_average;       
	float d_e, E_min;
	int converged, diverged;
	int (*Boltzmann_Selection)(CHROMOSOME*, float, float)=NULL;
	void (*eval_population)(PROT prot, CHROMOSOME *pop)=NULL;

        int check_count=0;
        int break_count=0;
        int phase = 10;
        int shift = 5;
        int total_count = 0;
        int max_check_counter = 50;
        float old_best = 0;

        CHROMOSOME *population = malloc(ga.pop_size*sizeof(CHROMOSOME));
        CHROMOSOME *new_population = malloc(ga.pop_size*sizeof(CHROMOSOME));
        CHROMOSOME *final_population = malloc(ga.pop_bucket*sizeof(CHROMOSOME));

	int final_pop_count;	
	char filename_ga[64]; sprintf(filename_ga,"GA_OUTPUT%05i.txt",patch_nb);
	FILE *GA_OUTPUT;
	
	if(env.output==0) {GA_OUTPUT = stdout;fprintf(stdout,"GA: Outputting GA information to display screen\n");}
	else {GA_OUTPUT = fopen(filename_ga,"wt"); fprintf(stdout,"GA: Outputting GA information to file 'GA_OUTPUT%05i.txt'\n",patch_nb);}
	
	fprintf(GA_OUTPUT,"GA: Chromosome length is %i (residues)\n",ga.N);
	
	if(ga.rand_cut_points>= ga.N/2) {
		fprintf(GA_OUTPUT,"GA: Number of random cutting points for cross-over operator exceeds chromosome {length / 2} of %i.\nAssigning default value of %i\n\n",ga.N/2,(int)(ga.N*0.4));
		ga.rand_cut_points = (int)(ga.N*0.4);
	}

	if(ga.seed==-1) {srand(time(0));}
	else {srand(ga.seed);}

	eval_population = evaluate_population;
	if(eval_p==1) {eval_population = evaluate_population2;}

	fprintf(GA_OUTPUT,"GA: Initializing data, please wait...\n\n");	
        /*sampling bucket*/
        for(i=0; i<ga.pop_bucket; i++) {
                final_population[i].genes = calloc(ga.N, sizeof(int));
                final_population[i].pw_vdw = calloc(ga.N, sizeof(VDWCONF));
        }
	fprintf(GA_OUTPUT,"GA: Setting up initial GA population...\n");
	initialize_populations(population, new_population);
	fprintf(GA_OUTPUT,"GA: Evaluating initial GA population...\n");
	eval_population((*prot), population);/*evaluate initial population*/
	fprintf(GA_OUTPUT,"GA: Initialization finished.\n"); 

     	CHROMOSOME orig;
        orig.genes = calloc(ga.N,sizeof(int));
	orig.pw_vdw = calloc(ga.N,sizeof(VDWCONF));

        /*setting up the solution to the initial conformation for each residue*/
	for(i=0; i<ga.N; i++) {
		orig.genes[i] = 0;
		orig.pw_vdw[i].vdw = 0.0f;
		orig.pw_vdw[i].min = FLT_MAX;
		orig.genes[i] = 1;
        }
	evaluate_ind((*prot),&orig);
	fprintf(GA_OUTPUT,"GA: Original Configuration's Energy [%08.7f] NB Residues[%06i]\n",orig.fitness,ga.N);

	memmove(population[ga.pop_size-1].genes, orig.genes, sizeof(int)*ga.N);
	memmove(population[ga.pop_size-1].pw_vdw, orig.pw_vdw, sizeof(VDWCONF)*ga.N);
	population[ga.pop_size-1].fitness = orig.fitness;
	population[ga.pop_size-1].flag = 0;
	population[ga.pop_size-1].mated = 0;
	population[ga.pop_size-1].copied = 0;
	set_vdwmin(population);

	Boltzmann_Selection = Boltzmann_Pressure_Converge;
	converged = diverged = 0;
	final_pop_count = 0;
	
	fprintf(GA_OUTPUT,"GA: Running Genetic Algorithm: INDIVIDUALS[%i] GENERATIONS[AUTOMATIC] CROSSOVER[%3.1f%%] MUTATION[%3.1f%%] MIGRATION[%3.1f%%]\n\n", ga.pop_size, ga.xover_rate*100, ga.mutation_rate*100, ga.migration_rate*100);


	for(i=0; i>=0; i++) {
		ftime(&start);
		new_pop_count = 0;
		/* elitism; usually only 1% but here, only 1 out of ga.pop_size*/
		elitism(&new_pop_count,population,new_population);
		pop_average = compute_population_average(population);
		E_min = new_population[0].fitness;//elite's fitness

		d_e = deltaE(population, E_min);
	        if(!converged) {
                        if(break_count==phase) {
                                if(check_count<shift) {
                                        if(fabs(E_min-old_best) < 0.00001 ) {/*same best*/
                                                total_count++;
                                                check_count++;
                                        }
                                        else {
                                                total_count = 0;
                                                check_count = 0;
                                                break_count = 0;
                                        }
                                }
                                else {
                                        break_count = 0;/*reset counter*/
                                        check_count = 0;/*reset counter*/
                                }
                                old_best = E_min;
                        }
                        else {
                                break_count++;/*increment counter*/
                        }
                }
	
		fprintf(GA_OUTPUT,"GA: Converged[%6.2f%%] Gen[%07i] BestFit[%10.3f] AveFit[%10.3f] Dist.Center[%6.1f] dE[%6.2f] ", total_count/(float)max_check_counter*100,i, E_min, pop_average, ga.dist_center, d_e);
		fflush(GA_OUTPUT);

		if(!converged) {
			fprintf(GA_OUTPUT, "Converging\t"); fflush(GA_OUTPUT);
		}
		else {
			if(!diverged) {
				if( d_e < ga.dist_center ) {
					fprintf(GA_OUTPUT,"Evol.Sampling: Diverging\t"); fflush(GA_OUTPUT);
					Boltzmann_Selection = Boltzmann_Pressure_Diverge2;
				}
				else {
					fprintf(GA_OUTPUT, "Evol.Sampling: DIVERGED   \t"); fflush(GA_OUTPUT);
					diverged = 1;
				}
			}
			else {
				fprintf(GA_OUTPUT, "Evol.Sampling: Converging\t"); fflush(GA_OUTPUT);
				Boltzmann_Selection = Boltzmann_Pressure_Converge2;
				if( d_e < ga.dist_center ) diverged = 0;
			}
		}
	
		while(new_pop_count<ga.pop_size) {
			/* Cross-over operator*/
			if(random_number_0_1() < ga.xover_rate) {
				do {
					parent1_idx = Boltzmann_Selection(population, d_e, E_min); parent2_idx = Boltzmann_Selection(population, d_e, E_min);
                        	} while( condition(population,parent1_idx,parent2_idx) );
				Random_N_Point_Crossover(population, new_population, &new_pop_count, parent1_idx, parent2_idx, ga.rand_cut_points);
			}
			/*Crossover operator*/

			/*Migration operator*/
			if(new_pop_count<ga.pop_size) {
				if(random_number_0_1() < ga.migration_rate) {
					do {
						parent1_idx = Boltzmann_Selection(population, d_e, new_population[0].fitness);
					} while ( population[parent1_idx].copied==1 );
					memmove(new_population[new_pop_count].genes,population[parent1_idx].genes,sizeof(int)*ga.N);
					memmove(new_population[new_pop_count].pw_vdw,population[parent1_idx].pw_vdw,sizeof(VDWCONF)*ga.N);
					new_population[new_pop_count].flag = 1;
					new_population[new_pop_count].fitness = population[parent1_idx].fitness;
					new_pop_count++;
					population[parent1_idx].copied = 1;
				}
			}
			/*Migration operator*/
		}

		/*mutate the whole new population, except the elitist ones*/
		mutate_population(new_population);

		/*check for similar solutions and discourage niching by applying forced 100% gene mutation*/
		/*this also guarantees the introduction of new genetic baggage*/
		/*provided there are exact and too convergent solutions in the population*/
		/*stresses variability of solutions and prevents stagnation of solutions*/
		/*in other words, we prevent 'd_e' from ever reaching a value of 0*/
		similarity_check_fast(new_population);

		/*re-evaluate individuals; only those with .flag=0 will be computed.*/
		/*At this point, the whole population is flagged .flag=0 anyway,*/
  		/*except for the newly introduced mutants*/
		clear_sumvdw(new_population);
		eval_population((*prot), new_population);
		set_vdwmin(new_population);
		copy_population_over(new_population,population);

                if ( (!converged) && (total_count>=max_check_counter) ) {
                        converged = 1;
                }

		if(converged && (final_pop_count < ga.pop_bucket) ) {
			/*filter population and keep those within FIXED_BAR*EPS of the min*/
			#pragma omp parallel for
			for(k=0; k<ga.pop_size; k++) {
				if( fabs(population[k].fitness - E_min) < (ga.dist_center_eps*ga.dist_center) ) {
					#pragma omp critical
					{
						if(final_pop_count < ga.pop_bucket) {
							memmove(final_population[final_pop_count].genes, population[k].genes, sizeof(int)*ga.N);
							memmove(final_population[final_pop_count].pw_vdw, population[k].pw_vdw, sizeof(VDWCONF)*ga.N);
							final_population[final_pop_count].fitness = population[k].fitness;
							final_pop_count++;
						}
						else {
							i = -1;
						}
					}
				}
			}
			if(i==-1) {goto out;}
		}

		ftime(&end);
                elap = (1000.0*(end.time-start.time)+(end.millitm-start.millitm))/1000;
                fprintf(GA_OUTPUT," FinalPop[%i] Time[%8.3f]\n",final_pop_count, elap); fflush(GA_OUTPUT);
                elapsed += elap;
	}
	out:
	/*called elitism to get the last best solution*/
	new_pop_count = 0;
	elitism(&new_pop_count,population,new_population);
	nb_rotamers = final_rotamers((*prot), final_population, final_pop_count);//compute_nb_rotamers((*prot), population);
	pop_average = compute_population_average(population);
	E_min = new_population[0].fitness;

	d_e = deltaE(population, E_min);
	fprintf(GA_OUTPUT,"GA: BestFit[%10.3f] AveFit[%10.3f] UniqueRotamers[%i] DeltaE[%6.2f]\n", E_min, pop_average, nb_rotamers, d_e);
	fprintf(GA_OUTPUT,"GA: Total Time[%8.3f seconds]\n\n",elapsed); fflush(GA_OUTPUT);

	/*if(env.gen_histogram==1) {
		write_best_solution_as_pdb((*prot),new_population, population, &orig, final_population, final_pop_count, ga.dist_center);
	}*/
	
	/*Now starting post-GA filtering*/
	/*Keep track of the rotamers to keep/remove in 2d structure 'res[i][j]'*/
        /*use 'tmp_flag' of CONF structure to keep track of the conformers*/
        /*we need to keep for MCCE step3 and step4*/

	int **res;
        res = calloc(prot->n_res,sizeof(int*));
        int count=0;
	int conf_count=0; int conf_del=0; int tot_conf=0;
	float tmp_min, occ_fit;
       
	//initialize 'res' 
	#pragma omp parallel for private(i,j)
	for(i=0; i<prot->n_res; i++) {
                res[i] = calloc(prot->res[i].n_conf,sizeof(int));
                for(j=0; j<prot->res[i].n_conf; j++) {res[i][j] = 0;}
        }
	
	//flag all the rotamers found in the sampling bucket
	for(i=0; i<final_pop_count; i++) {
		for(j=0; j<ga.N; j++) {
				res[ga.list[j]][final_population[i].genes[j]]++;
		}
	}

        for(i=0; i<final_pop_count; i++) {
		for(j=0; j<ga.N; j++) {
			final_population[i].pw_vdw[j].min = FLT_MAX;
		}
        }

	/*update global minium pairwise VDW to BUCKET pairwise VDW minium*/
	/*(more presentative than an absolute minimum from entire GA run)*/
	for(i=0; i<ga.N; i++) {
		tmp_min = FLT_MAX;
		//#pragma omp parallel for
		for(j=0; j<final_pop_count; j++) {
      			if(final_population[j].pw_vdw[i].vdw < tmp_min) {
				tmp_min = final_population[j].pw_vdw[i].vdw;
				occ_fit = final_population[j].fitness;
			}
			
        	}
		for(j=0; j<final_pop_count; j++) {
			final_population[j].pw_vdw[i].min = tmp_min;
			final_population[j].pw_vdw[i].occ_fit = occ_fit;
		}
	}
	
	/*removing rotamers that fall outside the range*/
	if(ga.occupancy>0) {
		fprintf(GA_OUTPUT,"GA: Occupancy cutoff %f with local residue cutoff %f\n",ga.occupancy,ga.residue_check); fflush(GA_OUTPUT);
		/*initialize occupancy array*/
		typedef struct {
			float val;
			int count;
		}OCC;
		OCC **oc = calloc(prot->n_res,sizeof(OCC*));
		for(i=0; i<prot->n_res; i++) {
			oc[i] = calloc(prot->res[i].n_conf,sizeof(OCC));
			for(j=0; j<prot->res[i].n_conf; j++) {oc[i][j].val=0.0f; oc[i][j].count=0;}
		}
		/*count the number of times a particular rotamer appears in the solutions*/
		/*keep track of the number of times it falls within the acceptance level*/
		for(i=0; i<final_pop_count; i++) {
			for(j=0; j<ga.N; j++) {
				if( fabs(final_population[i].pw_vdw[j].vdw - final_population[i].pw_vdw[j].min) <= ga.residue_check ) {
					oc[ga.list[j]][final_population[i].genes[j]].val += 1.0f;
				}
				oc[ga.list[j]][final_population[i].genes[j]].count += 1;
			}
		}
		/*express occupancies as percentage*/
		for(i=0; i<prot->n_res; i++) {
			for(j=0; j<prot->res[i].n_conf; j++) {
				if(oc[i][j].count>0) {
					oc[i][j].val /= oc[i][j].count;
					if(oc[i][j].val<ga.occupancy) {
						res[i][j] = 0;
					}
				}
			}
		}
		for(i=0; i<prot->n_res; i++) {free(oc[i]);}
		free(oc);
	}
	else {
		fprintf(GA_OUTPUT,"GA: No Occupancy cutoff with local residue cutoff %f\n",ga.residue_check); fflush(GA_OUTPUT);
		for(i=0; i<final_pop_count; i++) {
        	        for(j=0; j<ga.N; j++) {
                	        if( fabs(final_population[i].pw_vdw[j].vdw - final_population[i].pw_vdw[j].min) > ga.residue_check ) {
					res[ga.list[j]][final_population[i].genes[j]] = 0;
                        	}
	                }
        	}
	}
        
	for(i=0; i<prot->n_res; i++) {
		prot->res[i].conf[0].tmp_flag = 1;/*always keep the backbone*/
		if(prot->res[i].n_conf>1) {
			prot->res[i].conf[1].tmp_flag =1;/*and always keep initial conformation*/
			conf_count++;
		}
		tot_conf += prot->res[i].n_conf-1;
		for(j=2; j<prot->res[i].n_conf; j++) {
			if(res[i][j] > 0) {prot->res[i].conf[j].tmp_flag = 1; conf_count++;}
			else {prot->res[i].conf[j].tmp_flag = 0; conf_del++;}
		}
	}
	fprintf(GA_OUTPUT,"GA: Total Conformers: %i Conformers kept: %i Conformers deleted: %i\n",tot_conf, conf_count, conf_del); fflush(GA_OUTPUT);


	/*write out which rotamers to keep for each residue*/	
	char dir[64];
	sprintf(dir,"./pdb_patches/ga_output_patch%04i", patch_nb);
        FILE *out = fopen(dir,"wt");
        for(j=0; j<prot->n_res; j++) {
 	       for(k=0; k<prot->res[j].n_conf; k++) {
        	       fprintf(out,"%i\n", res[j][k]);
               }
        }
        fclose(out);

	/*remove unused rotamers from protein structure*/	
	for(i=0; i<prot->n_res; i++) {
		del:
		for(j=2; j<prot->res[i].n_conf; j++) {//don't delete backbone or initial conformation
			if(prot->res[i].conf[j].tmp_flag==0) {
				del_conf(&prot->res[i], j);
				goto del;
			}
		}
	}
	
	/*reset the confID of each conformer*/
	for(i=0; i<prot->n_res; i++) {
		for(j=2; j<prot->res[i].n_conf; j++) {
			sprintf(prot->res[i].conf[j].confID,"%08i\0",j);
			for(k=0; k<prot->res[i].conf[j].n_atom; k++) {
				sprintf(prot->res[i].conf[j].atom[k].confID,"%08i\0",j);
			}
		}
	}

	fprintf(GA_OUTPUT,"GA: GA Finished successfully\n"); fflush(GA_OUTPUT);
	if(env.output) {fclose(GA_OUTPUT);}

	/*free dynamically allocated memory*/
	for(i=0; i<ga.pop_size; i++) {
		free(population[i].genes); free(population[i].pw_vdw);
		population[i].genes = NULL; population[i].pw_vdw = NULL;
		free(new_population[i].genes);free(new_population[i].pw_vdw);
		new_population[i].genes = NULL; new_population[i].pw_vdw = NULL;
	}
	free(population); population = NULL;
	free(new_population); new_population = NULL;

        for(i=0; i<final_pop_count; i++) {
                free(final_population[i].genes); final_population[i].genes = NULL;
                free(final_population[i].pw_vdw); final_population[i].pw_vdw = NULL;
        }
        free(final_population); final_population = NULL;

	free(orig.genes); orig.genes=NULL;
	free(orig.pw_vdw); orig.pw_vdw=NULL;

	if(eval_p==1) {//if we pre-calculated pairwise VDW then we need to free memory
		for(i=0; i<ga.tot_conf; i++) {
			if(ga.conf_nrg[i]!=NULL) {
				free(ga.conf_nrg[i]);
				ga.conf_nrg[i] = NULL;
			}
		}
		free(ga.conf_nrg); ga.conf_nrg = NULL;
	}

	for(i=0; i<prot->n_res; i++) {
		free(res[i]);
		res[i] = NULL;
	}
	free(res); res=NULL;
	return;//return res;/*elements of 2-d array corresponds to res[res_index][n_conf]; if value is >0 then keep as conformer, otherwise deleted*/
};

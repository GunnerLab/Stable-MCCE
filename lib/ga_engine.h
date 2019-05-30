/* Pascal Comte, Brock University, St. Catharines, Ontario, Canada - 2010
 * Master's Thesis Project: GA sidechain packing optimization and sampling technique
 * 
 * Read thesis chapter 4 before attempting to modify code
 * One needs to clear understanding of genetic algorithm & simulated annealing theories,
 * Boltzmann & exponential laws.
 * 
 * IMPORTANT:
 * This Algorithm for sidechain packing is MEANT to execute on a system with multiple CORES
 * and operates using openMP.  We suggest a system with at least a QUAD-CORE CPU.
 * 
 * Work that uses this algorithm needs to cite the paper reference: []
 *  
 */

typedef struct {
        float x;
        float y;
        float z;
        float vdw_eps;
        float vdw_rad;
        float sas;
} B2_ATM;

typedef struct {
        float min;
        double vdw;/*double for summation*/
	float occ_fit;
} VDWCONF;

typedef struct {
        int *genes;/*index of the rotamer to use for this chromosome*/
       	VDWCONF *pw_vdw; 
	char flag;
        char mated;
        char copied;
        double fitness;/*vdw energy value; double for summation*/
} CHROMOSOME;

typedef struct {
        float fitness;
        int index;
} SOBJ;

typedef struct {
        int *list;/*index of residues to include in the ga*/
        int *values;/*number of rotamers for each residue (1-1 mapping)*/
        int N;/*chromosome length;*/
        int pop_size;
        float mutation_rate;
        float xover_rate;
        float migration_rate;
        int seed;
        int k_select;
        int rand_cut_points;
	int generations;
        //int b2_nb_atoms;
        //B2_ATM *b2_atoms;
        //float **eps_atm;
        float elitism;

	float dE;
	int phase;
	int shift;
	int sphere_focus_resid;
	float sphere_focus_probe_radius;

	float occupancy;	
	float dist_center;/*within KT from the min*/
	float dist_center_eps;/*used for filling the bucket*/
	float epsilon;/*convergence criteria*/
	int pop_bucket;/*how many chromosomes (protein conformations) in final filling bucket (duplicates allowed to create a distribution)*/
	float residue_check;/*used to check how far from the minimum occupied energy to keep solutions from*/
	float **conf_nrg;
	int tot_conf;
	
} GA_STRUCTURE;

int compareSOBJ(const void *A, const void *B);
int compareFLT(const void *A, const void *B);
int compareINT(const void *A, const void *B);

__inline int random_number_1_M(int M);
__inline int random_number_0_M(int M);
__inline float random_number_0_1(void);

__inline void initialize_populations(CHROMOSOME *pop, CHROMOSOME *new_pop);
void initialize_eps_matrix(PROT prot);
void evaluate_ind(PROT prot, CHROMOSOME *ind);
void evaluate_population(PROT prot, CHROMOSOME *pop);

__inline void mutate_population(CHROMOSOME *new_pop);
__inline int Boltzmann_Pressure_Converge(CHROMOSOME *pop, float dE, float E_min);
__inline int Boltzmann_Pressure_Diverge(CHROMOSOME *pop, float dE, float E_min);

__inline int Boltzmann_Pressure_Converge2(CHROMOSOME *pop, float dE, float E_min);
__inline int Boltzmann_Pressure_Diverge2(CHROMOSOME *pop, float dE, float E_min);


__inline void elitism(int *new_pop_count, CHROMOSOME *pop, CHROMOSOME *new_pop);
__inline int condition(CHROMOSOME *pop, int parent1, int parent2);

void Random_N_Point_Crossover(CHROMOSOME *pop, CHROMOSOME *new_pop, int *new_pop_count, int parent1_idx, int parent2_idx, int max_cut);
void copy_population_over(CHROMOSOME *new_population, CHROMOSOME *population);
int compute_nb_rotamers(PROT prot, CHROMOSOME *pop);
__inline float compute_population_average(CHROMOSOME *pop);
void write_backbone_atoms(FILE *crd, PROT prot);

int genetic_difference(int *genes1, int *genes2);
void write_best_solution_as_pdb(PROT prot, CHROMOSOME *new_population, CHROMOSOME *pop, CHROMOSOME *original, CHROMOSOME *final_population, int final_pop_count, float KT);

__inline float deltaE(CHROMOSOME *pop, float E_min);

void GA_SETUP(PROT *prot, GA_STRUCTURE gs, int patch_nb);
void RUN_GA(PROT *prot, int patch_nb, int eval_p);

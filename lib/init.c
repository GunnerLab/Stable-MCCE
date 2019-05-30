#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

ENV env;

int init()
{   FILE *fp;
    float kscale;
    time_t now;

    printf("   Load run control file \"%s\"...\n", FN_RUNPRM); fflush(stdout);
    if (get_env()) {printf("   FATAL: init(): \"failed initializing.\"\n"); return USERERR;}
    else {printf("   Done\n\n"); fflush(stdout);}

    printf("   Tentatively load local param file \"%s\"...", env.new_tpl); fflush(stdout);
    if (env.do_premcce) remove(env.new_tpl);
    if ((fp=fopen(env.new_tpl, "r"))) {
        fclose(fp);
        if (load_param(env.new_tpl)) {
            printf("\n   FATAL: init(): Failed loading file \"%s\".\n", env.new_tpl);
            return USERERR;
        }
        printf("   File loaded.\n");
    }
    else printf("   No such file, ignore.\n");
    printf("   Done\n\n");
    fflush(stdout);
    
    printf("   Load parameters from directory \"%s\" ... \n", env.param); fflush(stdout);
    if (load_all_param(env.param)) {printf("   FATAL: init(): \"failed.\"\n"); return USERERR;}
    else {printf("   Done\n\n"); fflush(stdout);}

    printf("   Load linear free energy correction parameters from \"%s\"...", env.extra);fflush(stdout);
    if ((fp=fopen(env.extra, "r"))) {
        printf("%s\n", env.extra);
        fclose(fp);
        if (load_param(env.extra)) {
            printf("\n   FATAL: init(): Failed loading file \"%s\".\n", env.extra);
            return USERERR;
        }
        printf("   File loaded.\n");
    }
    else printf("   No such file, ignore.\n");
    printf("   Done\n\n");
    fflush(stdout);

    /* Order of scale values (highest has priority):
     * 1. values from (SCALE_VDW), (SCALE_VDW0), (SCALE_VDW1) options in run.prm
     * 2. tpl file (usually extra.tpl)
     * 3. default value from env.epsilon_prot
     */
    kscale = 1.0; /*/env.epsilon_prot;*/ /* scaling factor based on dielectric constant */
    if (env.scale_vdw0 < 0) {
    	if (param_get("SCALING", "VDW0", "", &env.scale_vdw0))  env.scale_vdw0   = 1.0*kscale;
    }
    if (env.scale_vdw1 < 0) {
    	if (param_get("SCALING", "VDW1",  "", &env.scale_vdw1)) env.scale_vdw1   = 1.0*kscale;
    }
    if (env.scale_vdw < 0) {
    	if (param_get("SCALING", "VDW",   "", &env.scale_vdw)) env.scale_vdw     = 1.0*kscale;
    }
    if (param_get("SCALING", "TORS",  "", &env.scale_tor)) env.scale_tor     = 1.0*kscale;
    if (param_get("SCALING", "ELE",   "", &env.scale_ele)) env.scale_ele     = 1.0;
    if (param_get("SCALING", "DSOLV", "", &env.scale_dsolv)) env.scale_dsolv = 1.0;

    remove(env.debug_log);
    remove(env.progress_log);
	
	now = time(NULL);
    if (env.test_seed < 0) srand(now); //allows random numbers to be fixed for testing
    else srand(env.test_seed);

    return 0;
}


int get_env()
{   FILE *fp;
    char sbuff[256];
    char *str1;
    /* vars related to run.trace */
    time_t now;
    struct tm *time_ptr;
    int dotrace = 0;
    FILE *tr;
    char trbuff[256];
    FILE *ent;

    memset(&env, 0, sizeof(ENV));

    /* Default values */
    env.test_seed = -1;
    env.minimize_size = 0;
    env.PI                = 4.*atan(1.);
    env.d2r               = env.PI/180.;

    strcpy(env.debug_log,    "debug.log");
    strcpy(env.new_tpl,      "new.tpl");
    strcpy(env.progress_log, "progress.log");
    strcpy(env.extra, "extra.tpl");
    env.reassign     = 0;
    env.pbe_start = 0;
    env.pbe_end   = 999999;
    strcpy(env.pbe_solver, "delphi");
    strcpy(env.rxn_method, "self");		/*surface or self energies*/
    env.rot_specif   = 0;
    env.prune_thr = 0.01;
    env.ngh_vdw_thr  = 0.1;
    env.repack_e_thr_exposed  = 0.5;
    env.repack_e_thr_buried  = 4.;
    env.repack_fav_vdw_off   = 0;
    env.nconf_limit       =    0;
    env.n_hv_conf_limit   =   20;
    env.relax_wat         =    1;
    env.trans_dist        =  0.5;

    env.hdirected         =    0;
    env.hdirdiff          =  1.0;
    env.hdirlimt          =   36;

    env.water_relax_thr   =  2.4;

    env.n_initial_relax     =  0;
    //env.initial_relax_rebuild = 0;

    env.hv_relax_ncycle     =  0;
    env.hv_relax_niter      = 50;
    env.hv_relax_vdw_thr    =  5;
    env.hv_relax_hv_vdw_thr =  5;
    env.hv_relax_elec_thr   =  -2.0;
    env.hv_relax_elec_crg_thr =  0.1;
    env.hv_relax_elec_dist_thr = 2.4;
    env.hv_relax_dt         =  1;
    env.hv_tors_scale       =  1;
    env.hv_relax_constraint =  1.;
    env.hv_relax_constraint_frc = 10.;
    env.hv_relax_n_shake    =  3000;
    env.hv_relax_shake_tol =  1e-4;  /* Ratio to constraint distance */
    env.hv_relax_include_ngh    =  0;
    env.hv_relax_ngh_thr    =  4.;
    env.prune_rmsd        = 2.0;
    env.prune_ele         = 2.0;
    env.prune_vdw         = 2.0;


    env.relax_n_hyd       =    6;
    env.relax_clash_thr   =  10.;

    env.recalc_tors     = 0;

    env.default_radius = 1.7;
    env.factor_14lj = 0.5;
    env.epsilon_coulomb = 6.;
    
    env.sas2vdw = -0.06;
    
    env.warn_pairwise     = 20.0;
    env.big_pairwise      = 5.0;

    env.monte_adv_opt     =    0;
    env.anneal_temp_start = ROOMT;
    env.anneal_nstep      =    1;
    env.monte_tsx         =    0;
    env.anneal_niter_step =   30;
    env.monte_niter_max   =   -1;
    env.adding_conf       =    0;
    env.monte_old_input   =    0;
    env.monte_niter_chk   =  100;
    env.monte_converge    = 1e-4;
    env.monte_do_energy   =    0;
    env.monte_print_nonzero =  1;
    strcpy(env.pbe_folder, "/tmp");
    env.delphi_clean      =  1;
    env.ionrad = 0.0;
    env.salt =0.00;

    /* default value for IPECE */
    memset(&env.ipece,0,sizeof(IPECE));
    
    env.ipece.grid_space = 1.0;
    
    env.ipece.mem_position_defined = 0;
    env.ipece.probe_radius = 1.40;
    env.ipece.surface_exp_dist = 5.;
    
    env.ipece.boundary_extention = 5.;
    env.ipece.half_mem_thickness = 15.;
    env.ipece.mem_separation = 3.;
    env.ipece.membrane_size = 10.;
    strcpy(env.ipece.mem_resName, "MEM");
    env.ipece.mem_chainID = 'X';
    env.ipece.mem_atom_radius = 1.7;

    env.ipece.beta = 0.2;
    env.ipece.n_iteration = 500;
    env.ipece.translation_max = 3.;
    env.ipece.rotation_max = 10.*env.d2r;

    env.mfe_cutoff = 0.0;
    env.mfe_point = 0.0;
    env.mfe_flag = 0;
    env.scale_vdw0=-1.0;
    env.scale_vdw1=-1.0;
    env.scale_vdw=-1.0;
    /* apbs */

    env.fg_scale = 1.5;
    env.grids_apbs = 129;
    strcpy(env.apbs_method, "mg-auto");
    strcpy(env.srfm, "mol");     /*model used to construct the dielectric and ion-accessibility coefficients*/
    strcpy(env.chgm, "spl0");    /*method by which the biomolecular point charges are mapped to the grid*/
    strcpy(env.bcfl, "sdh");     /*specify the type of boundry condition used -"Single Debye-H√ºckel" boundary condition*/
    strcpy(env.apbs_exe, "apbs");
    
    env.ignore_input_h = 1;
    env.do_corrections = 1;
     
    /* open "run.prm" to read in mcce environment variables */
    if ((fp=fopen(FN_RUNPRM, "r")) == NULL) {
        printf("   FATAL: get_env(): \"No run control file %s.\"\n", FN_RUNPRM);
        return USERERR;
    }

    /* user values */
    while (fgets(sbuff, sizeof(sbuff), fp)) {
    	if (sbuff[0] == '#') {
             continue;
        }
    	if (strstr(sbuff, "(DO_TRACE)")) {
        	str1 = strtok(sbuff, " ");
        	if (str1[0] == 't') dotrace = 1;
        }
        else if (strstr(sbuff, "(TEST_SEED)")) {
            env.test_seed = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(IGNORE_INPUT_H)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.ignore_input_h = 1;
            else env.ignore_input_h = 0;
        }
        else if (strstr(sbuff, "(INPDB)")) {
            strcpy(env.inpdb, strtok(sbuff, " "));
        }
        
        else if (strstr(sbuff, "(MCCE_HOME)")) {
            strcpy(env.mcce_home, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(DEBUG_LOG)")) {
            strcpy(env.debug_log, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PROGRESS_LOG)")) {
            strcpy(env.progress_log, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(NEWTPL)")) {
            strcpy(env.new_tpl, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(EXTRA)")) {
            strcpy(env.extra, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MINIMIZE_SIZE)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.minimize_size = 1;
            else env.minimize_size = 0;
        }

        
        else if (strstr(sbuff, "(DO_PREMCCE)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.do_premcce = 1;
            }
            else env.do_premcce = 0;
        }
        else if (strstr(sbuff, "(TERMINALS)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.terminals = 1;
            }
            else env.terminals = 0;
        }
        else if (strstr(sbuff, "(DO_ROTAMERS)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.do_rotamers = 1;
            else env.do_rotamers = 0;
        }
        else if (strstr(sbuff, "(DO_ENERGY)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.do_energies = 1;
            else env.do_energies = 0;
        }
        else if (strstr(sbuff, "(DO_MONTE)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.do_monte = 1;
            else env.do_monte = 0;
        }

        else if (strstr(sbuff, "(REBUILD_SC)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.rebuild_sc = 1;
            else env.rebuild_sc = 0;
        }
        else if (strstr(sbuff, "(ROT_SWAP)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.rot_swap = 1;
            else env.rot_swap = 0;
        }
        else if (strstr(sbuff, "(DEFAULT_RADIUS)")) {
            env.default_radius = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(FACTOR_14LJ)")) {
            env.factor_14lj = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(EPSILON_COULOMB)")) {
            env.epsilon_coulomb = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(SAS2VDW)")) {
            env.sas2vdw = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(CLASH_DISTANCE)")) {
            env.clash_distance = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(H2O_SASCUTOFF)")) {
            env.h2o_sascutoff = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RENAME_RULES)")) {
            strcpy(env.rename_rules, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(ROT_SPECIF)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.rot_specif = 1;
            else env.rot_specif = 0;
        }
        else if (strstr(sbuff, "(SWING)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.swing = 1;
            else env.swing = 0;
        }
        else if (strstr(sbuff, "(HDIRECTED)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.hdirected = 1;
            else env.hdirected = 0;
        }
        else if (strstr(sbuff, "(HDIRDIFF)")) {
            env.hdirdiff = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HDIRLIMT)")) {
            env.hdirlimt = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PACK)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.pack = 1;
            else env.pack = 0;
        }
        else if (strstr(sbuff, "(N_TRANS)")) {
            env.n_trans = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(TRANS_DIST)")) {
            env.trans_dist = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PRUNE_THR)")) {
            env.prune_thr = atof(strtok(sbuff, " "));
        }

        else if (strstr(sbuff, "(SAS_CUTOFF)")) {
            env.sas_cutoff = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(VDW_CUTOFF)")) {
            env.vdw_cutoff = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(REPACK_CUTOFF)")) {
            env.repack_cutoff = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(REPACK_FAV_VDW_OFF)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.repack_fav_vdw_off = 1;
            else env.repack_fav_vdw_off = 0;
        }
        else if (strstr(sbuff, "(REPACK_E_THR_EXPOSED)")) {
            env.repack_e_thr_exposed = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(REPACK_E_THR_BURIED)")) {
            env.repack_e_thr_buried = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(NGH_VDW_THR)")) {
            env.ngh_vdw_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PHI_SWING)")) {
            env.phi_swing = atof(strtok(sbuff, " ")) * env.d2r;
        }
        else if (strstr(sbuff, "(ROTATIONS)")) {
            env.rotations = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(ROTAMER_LIMIT)")) {
            env.rotamer_limit = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(REPACKS)")) {
            env.repacks = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(NCONF_LIMIT)")) {
            env.nconf_limit = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(N_HV_CONF_LIMIT)")) {
            env.n_hv_conf_limit = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RELAX_WAT)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.relax_wat = 1;
            else env.relax_wat = 0;
        }
        else if (strstr(sbuff, "(WATER_RELAX_THR)")) {
            env.water_relax_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(N_INITIAL_RELAX)")) {
            env.n_initial_relax = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_NCYCLE)")) {
            env.hv_relax_ncycle = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_NITER)")) {
            env.hv_relax_niter = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_VDW_THR)")) {
            env.hv_relax_vdw_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_HV_VDW_THR)")) {
            env.hv_relax_hv_vdw_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_ELEC_THR)")) {
            env.hv_relax_elec_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_ELEC_CRG_THR)")) {
            env.hv_relax_elec_crg_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_ELEC_DIST_THR)")) {
            env.hv_relax_elec_dist_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_DT)")) {
            env.hv_relax_dt = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_TORS_SCALE)")) {
            env.hv_tors_scale = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_CONSTRAINT)")) {
            env.hv_relax_constraint = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_CONSTRAINT_FRC)")) {
            env.hv_relax_constraint_frc = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_N_SHAKE)")) {
            env.hv_relax_n_shake = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_SHAKE_TOL)")) {
            env.hv_relax_shake_tol = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(HV_RELAX_INCLUDE_NGH)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.hv_relax_include_ngh = 1;
            else env.hv_relax_include_ngh = 0;
        }
        else if (strstr(sbuff, "(HV_RELAX_NGH_THR)")) {
            env.hv_relax_ngh_thr = atof(strtok(sbuff, " "));
        }

        else if (strstr(sbuff, "(RELAX_H)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.relax_h = 1;
            else env.relax_h = 0;
        }
        else if (strstr(sbuff, "(RELAX_E_THR)")) {
            env.relax_e_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RELAX_NSTATES)")) {
            env.relax_nstates = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RELAX_N_HYD)")) {
            env.relax_n_hyd = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RELAX_CLASH_THR)")) {
            env.relax_clash_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RELAX_PHI)")) {
            env.relax_phi = 3.1415926/180.0 * atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RELAX_NITER)")) {
            env.relax_niter = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RELAX_TORQ_THR)")) {
            env.relax_torq_thr = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PRUNE_RMSD)")) {
            env.prune_rmsd = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PRUNE_ELE)")) {
            env.prune_ele = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PRUNE_VDW)")) {
            env.prune_vdw = atof(strtok(sbuff, " "));
        }


        else if (strstr(sbuff, "(REASSIGN)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.reassign = 1;
            }
            else env.average_pairwise = 0;
        }
        else if (strstr(sbuff, "(EPSILON_PROT)")) {
            env.epsilon_prot = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(EPSILON_SOLV)")) {
            env.epsilon_solv = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GRIDS_DELPHI)")) {
            env.grids_delphi = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GRIDS_PER_ANG)")) {
            env.grids_per_ang = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(RADIUS_PROBE)")) {
            env.radius_probe = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(IONRAD)")) {
            env.ionrad = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(SALT)")) {
            env.salt = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(DELPHI_EXE)")) {
            strcpy(env.delphi_exe, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(DELPHI_FAILS)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 'd' || str1[0] == 'D') {
                env.delphi_fails = 'd';
            }
            else  env.delphi_fails = 's';
        }
        else if (strstr(sbuff, "(RECALC_TORS)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.recalc_tors = 1;
            }
            else env.recalc_tors = 0;
        }

        else if (strstr(sbuff, "(AVERAGE_PAIRWISE)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.average_pairwise = 1;
            }
            else env.average_pairwise = 0;
        }
        else if (strstr(sbuff, "(WARN_PAIRWISE)")) {
            env.warn_pairwise = atof(strtok(sbuff, " "));
        }
        
        else if (strstr(sbuff, "(MONTE_ADV_OPT)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.monte_adv_opt = 1;
            else env.monte_adv_opt = 0;
        }
        else if (strstr(sbuff, "(MONTE_SEED)")) {
            env.monte_seed = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_T)")) {
            env.monte_temp = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_RUNS)")) {
            env.monte_runs = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_NSTART)")) {
            env.monte_nstart = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_NEQ)")) {
            env.monte_neq = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_NITER)")) {
            env.monte_niter = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_FLIPS)")) {
            env.monte_flips = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_TRACE)")) {
            env.monte_trace = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_TSX)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.monte_tsx = 1;
            }
            else env.monte_tsx = 0;
        }
        else if (strstr(sbuff, "(MONTE_N_REDUCE)")) {
            env.monte_n_red = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_REDUCE)")) {
            env.monte_reduce = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(NSTATE_MAX)")) {
            env.nstate_max = atoi(strtok(sbuff, " "));
        }
        
        else if (strstr(sbuff, "(ADDING_CONF)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.adding_conf = 1;
            }
            else env.adding_conf = 0;
        }
        else if (strstr(sbuff, "(MONTE_OLD_INPUT)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.monte_old_input = 1;
            }
            else env.monte_old_input = 0;
        }
        else if (strstr(sbuff, "(MONTE_NITER_CYCLE)")) {
            env.monte_niter_cycle = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_NITER_MIN)")) {
            env.monte_niter_min = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_NITER_MAX)")) {
            env.monte_niter_max = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_NITER_CHK)")) {
            env.monte_niter_chk = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_CONVERGE)")) {
            env.monte_converge = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MONTE_DO_ENERGY)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.monte_do_energy = 1;
            }
            else env.monte_do_energy = 0;
        }
        else if (strstr(sbuff, "(MONTE_PRINT_NONZERO)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') {
                env.monte_print_nonzero = 1;
            }
            else env.monte_print_nonzero = 0;
        }
        
        else if (strstr(sbuff, "(ANNEAL_TEMP_START)")) {
            env.anneal_temp_start = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(ANNEAL_NSTEP)")) {
            env.anneal_nstep = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(ANNEAL_NITER_STEP)")) {
            env.anneal_niter_step = atoi(strtok(sbuff, " "));
        }
        
        else if (strstr(sbuff, "(TITR_TYPE)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 'p' || str1[0] == 'P') env.titr_type = 'p';
            else env.titr_type = 'e';
        }
        else if (strstr(sbuff, "(TITR_PH0)")) {
            env.titr_ph0 = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(TITR_EH0)")) {
            env.titr_eh0 = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(TITR_PHD)")) {
            env.titr_phd = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(TITR_EHD)")) {
            env.titr_ehd = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(TITR_STEPS)")) {
            env.titr_steps = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(BIG_PAIRWISE)")) {
            env.big_pairwise = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(VDWF1)")) {
            env.vdwf1 = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(VDWF2)")) {
            env.vdwf2 = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(ROTATE_RES)")) {
            strcpy(env.rotate_res, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PBE_START)")) {
            env.pbe_start = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(PBE_END)")) {
            env.pbe_end = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(SCALE_VDW)")) {
            env.scale_vdw = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(SCALE_VDW0)")) {
            env.scale_vdw0 = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(SCALE_VDW1)")) {
            env.scale_vdw1 = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(SCALE_ELE)")) {
            env.scale_ele = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(SKIP_ELE)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.skip_ele = 1;
            else env.skip_ele = 0;
        }
        else if (strstr(sbuff, "(PBE_FOLDER)")) {
            strcpy(env.pbe_folder, strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(DELPHI_CLEAN)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.delphi_clean = 1;
            else env.delphi_clean = 0;
        }
        else if (strstr(sbuff, "(PBE_SOLVER)")) {
             strcpy(sbuff, strtok(sbuff, " "));
             if (strstr(sbuff, "apbs") || strstr(sbuff, "APBS")) {
                strcpy(env.pbe_solver, "apbs");
             }
             else if (strstr(sbuff, "delphi") || strstr(sbuff, "DELPHI")) {
                strcpy(env.pbe_solver, "delphi");
 	    	 }
 	    	 else if (strstr(sbuff, "zap") || strstr(sbuff, "ZAP")) {
                strcpy(env.pbe_solver, "zap");
 	    	 }
             else {
                 printf("\n   Not known PBE solver: \"%s\". Using delphi in step 3...\n", sbuff);
                 strcpy(env.pbe_solver, "delphi");
             }
         }
 		else if (strstr(sbuff, "(RXN_METHOD)")) {
             strcpy(sbuff, strtok(sbuff, " "));
             if (strstr(sbuff, "self") || strstr(sbuff, "SELF")) {
                strcpy(env.rxn_method, "self");
             }
             else if (strstr(sbuff, "ntsurface") || strstr(sbuff, "NTSURFACE")) {
                strcpy(env.rxn_method, "ntsurface");
 	    	 }
             else if (strstr(sbuff, "surface") || strstr(sbuff, "SURFACE")) {
                strcpy(env.rxn_method, "surface");
 	    	 }
             else {
                 printf("\n   Not known RXN method: \"%s\". Using self energies in step 3...\n", sbuff);
                 strcpy(env.rxn_method, "self");
             }
         }
 	/* apbs*/
 	else if (strstr(sbuff, "(GRIDS_APBS)")) {
             env.grids_apbs = atoi(strtok(sbuff, " "));
         }
 	else if (strstr(sbuff, "(SURFACE_APBS)")) {
             strcpy(env.srfm, strtok(sbuff, " "));
         }
 	else if (strstr(sbuff, "(CHARGES_APBS)")) {
             strcpy(env.chgm, strtok(sbuff, " "));
         }
 	else if (strstr(sbuff, "(BOUND_COND_APBS)")) {
             strcpy(env.bcfl, strtok(sbuff, " "));
         }
 	else if (strstr(sbuff, "(FINE_SCALE)")) {
             env.fg_scale = atof(strtok(sbuff, " "));
         }
 	else if (strstr(sbuff, "(APBS_EXE)")) {
             strcpy(env.apbs_exe, strtok(sbuff, " "));
         }


        /* IPECE */
        else if (strstr(sbuff, "(IPECE_GRID_SPACE)")) {
            env.ipece.grid_space = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(IPECE_ADD_MEM)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 't' || str1[0] == 'T') env.ipece.add_mem = 1;
            else env.ipece.add_mem = 0;
        }
        else if (strstr(sbuff, "(IPECE_MEM_THICKNESS)")) {
            env.ipece.half_mem_thickness = 0.5 * atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(IPECE_MEM_CHAINID)")) {
            env.ipece.mem_chainID = *strtok(sbuff, " ");    /* first character of the line */
        }
        else if (strstr(sbuff, "(IPECE_MEM_ATOM_RADIUS)")) {
            env.ipece.grid_space = atof(strtok(sbuff, " "));
        }
        /*Pascal's MSC Ga Parameters*/
        else if (strstr(sbuff, "(GA_POP_SIZE")) {
                env.pop_size = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_MUTATION_RATE)")) {
                env.mutation_rate = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_MIGRATION_RATE)")) {
                env.migration_rate = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_CROSSOVER_RATE)")) {
                env.xover_rate = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_RANDOM_CUT_POINTS)")) {
                env.rand_cut_points = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_SEED)")) {
                env.ga_seed = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_GENERATIONS)")) {
                env.generations = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_PHASE)")) {
                env.ga_phase = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_SHIFT)")) {
                env.ga_shift = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_DIST_CENTER)")) {
                env.ga_dist_center = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_DIST_CENTER_EPS)")) {
                env.ga_dist_center_eps = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_MAX_BUCKET_POP)")) {
                env.pop_bucket = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_RESIDUE_MIN_ENERGY_CUTOFF)")) {
                env.residue_check = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_OCCUPANCY_CUTOFF)")) {
                env.occupancy = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_DELTA_E)")) {
                env.ga_deltaE = atof(strtok(sbuff, " "));
        }
        else if(strstr(sbuff,"(SIDECHAIN_OPT)")) {
                env.sidechain_opt = atoi(strtok(sbuff, " "));
                if (env.sidechain_opt == 2) {
                	env.sidechain_opt = 1;
                	env.rot_mhd_prune = 1;
                }
        }
        else if(strstr(sbuff,"(GA_OUTPUT)")) {
                env.output = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_SPHERE_FOCUS_RESID)")) {
                env.ga_focus_resid = atoi(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(GA_SPHERE_PROBE_RADIUS)")) {
                env.ga_focus_probe_radius = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MFE_CUTOFF)")) {
                env.mfe_cutoff = atof(strtok(sbuff, " "));
        }
        else if (strstr(sbuff, "(MFE_POINT)")) {
                if (!(strchr(strtok(sbuff, " "), 'f'))) {
                    env.mfe_flag = 1; 
                    env.mfe_point = atof(strtok(sbuff, " "));
                }
        }
        else if (strstr(sbuff, "(DO_CORRECTIONS)")) {
            str1 = strtok(sbuff, " ");
            if (str1[0] == 'f' || str1[0] == 'F') {
                env.do_corrections = 0;
            }
            else env.do_corrections = 1;
        }
    }

    fclose(fp);
    if (env.monte_niter_cycle <= 0) env.monte_niter_cycle = 1000;
    if (env.repack_e_thr_exposed <= 0) env.repack_e_thr_exposed = 1e-4;
    if (env.repack_e_thr_buried <= 0) env.repack_e_thr_buried = 1e-4;
    if (env.ipece.boundary_extention <= env.ipece.probe_radius*4.)
        env.ipece.boundary_extention = env.ipece.probe_radius*4.;

    /* round the dielectric constant to the nearest integer number */
    sprintf(env.param, "%s/param%02d", env.mcce_home, (int) (env.epsilon_prot));
    
    /* sets env.ga_seed and env.monte_seed equal to env.test_seed if env.test_seed has been set to a 
    non random value, eliminating all randomness from MCCE */
    if (env.test_seed >= 1) {
    	env.monte_seed = env.test_seed;
    	env.ga_seed = env.test_seed;
    }
    	
	
	/*adds in the default behavior of the run.trace output*/
	if (dotrace == 1) {
		
		if ((fp=fopen(FN_RUNPRM, "r")) == NULL) {
        	printf("   FATAL: get_env(): \"No run control file %s.\"\n", FN_RUNPRM);
        	return USERERR;
        }
        if ((tr=fopen("run.trace", "a")) == NULL) {
        	printf("   FATAL: get_env(): \"Cannot create run.trace file.\"\n");
        	return USERERR;
        }
        
        /* this block grabs the subversion revision number from .svn/entries */
        if ((ent=fopen("/home/mcce/mcce2.5.1/.svn/entries", "r")) != NULL) {
        	int i;
        	for (i = 0; i < 4; ++i) {
        		fgets(trbuff, sizeof(trbuff), ent);
        	}
        	fprintf(tr, "The MCCE directory revision number is: %s", trbuff);
        	fclose(ent);
        }
        else {
        	fprintf(tr, "%s\n", "can't open entries file");
        }
        
        /* this block prints a timestamp */
        now = time(NULL);
        time_ptr = localtime(&now);
        fprintf(tr, "%s\n", asctime(time_ptr));
        
        /* this block prints headers for all of the columns */
        fprintf(tr, "%31s%31s\t%31s%31s\n", "run.prm var name", "run.prm var value", "init.c var name", "init.c var value");
        
        while (fgets(trbuff, sizeof(trbuff), fp)) {
			if (strstr(trbuff, "(ADDING_CONF)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(ADDING_CONF)", strtok(trbuff, " "), "env.adding_conf", env.adding_conf);
			}
			else if (strstr(trbuff, "(ANNEAL_NITER_STEP)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(ANNEAL_NITER_STEP)", strtok(trbuff, " "), "env.anneal_niter_step", env.anneal_niter_step);
			}
			else if (strstr(trbuff, "(ANNEAL_NSTEP)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(ANNEAL_NSTEP)", strtok(trbuff, " "), "env.anneal_nstep", env.anneal_nstep);
			}
			else if (strstr(trbuff, "(ANNEAL_TEMP_START)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(ANNEAL_TEMP_START)", strtok(trbuff, " "), "env.anneal_temp_start", env.anneal_temp_start);
			}
			else if (strstr(trbuff, "(APBS_EXE)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(APBS_EXE)", strtok(trbuff, " "), "env.apbs_exe", env.apbs_exe);
			}
			else if (strstr(trbuff, "(AVERAGE_PAIRWISE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(AVERAGE_PAIRWISE)", strtok(trbuff, " "), "env.average_pairwise", env.average_pairwise);
			}
			else if (strstr(trbuff, "(BIG_PAIRWISE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(BIG_PAIRWISE)", strtok(trbuff, " "), "env.big_pairwise", env.big_pairwise);
			}
			else if (strstr(trbuff, "(BOUND_COND_APBS)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(BOUND_COND_APBS)", strtok(trbuff, " "), "env.bcfl", env.bcfl);
			}
			else if (strstr(trbuff, "(CHARGES_APBS)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(CHARGES_APBS)", strtok(trbuff, " "), "env.chgm", env.chgm);
			}
			else if (strstr(trbuff, "(CLASH_DISTANCE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(CLASH_DISTANCE)", strtok(trbuff, " "), "env.clash_distance", env.clash_distance);
			}
			else if (strstr(trbuff, "(DEBUG_LOG)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(DEBUG_LOG)", strtok(trbuff, " "), "env.debug_log", env.debug_log);
			}
			else if (strstr(trbuff, "(DEFAULT_RADIUS)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(DEFAULT_RADIUS)", strtok(trbuff, " "), "env.default_radius", env.default_radius);
			}
			else if (strstr(trbuff, "(DELPHI_CLEAN)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(DELPHI_CLEAN)", strtok(trbuff, " "), "env.delphi_clean", env.delphi_clean);
			}
			else if (strstr(trbuff, "(DELPHI_EXE)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(DELPHI_EXE)", strtok(trbuff, " "), "env.delphi_exe", env.delphi_exe);
			}
			else if (strstr(trbuff, "(DELPHI_FAILS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(DELPHI_FAILS)", strtok(trbuff, " "), "env.delphi_fails", env.delphi_fails);
			}
			else if (strstr(trbuff, "(DO_ENERGY)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(DO_ENERGY)", strtok(trbuff, " "), "env.do_energies", env.do_energies);
			}
			else if (strstr(trbuff, "(DO_MONTE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(DO_MONTE)", strtok(trbuff, " "), "env.do_monte", env.do_monte);
			}
			else if (strstr(trbuff, "(DO_PREMCCE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(DO_PREMCCE)", strtok(trbuff, " "), "env.do_premcce", env.do_premcce);
			}
			else if (strstr(trbuff, "(DO_ROTAMERS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(DO_ROTAMERS)", strtok(trbuff, " "), "env.do_rotamers", env.do_rotamers);
			}
			else if (strstr(trbuff, "(EPSILON_SOLV)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(EPSILON_SOLV)", strtok(trbuff, " "), "env.epsilon_solv", env.epsilon_solv);
			}
			else if (strstr(trbuff, "(EPSILON_COULOMB)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(EPSILON_COULOMB)", strtok(trbuff, " "), "env.epsilon_coulomb", env.epsilon_coulomb);
			}
			else if (strstr(trbuff, "(EPSILON_PROT)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(EPSILON_PROT)", strtok(trbuff, " "), "env.epsilon_prot", env.epsilon_prot);
			}
			else if (strstr(trbuff, "(EXTRA)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(EXTRA)", strtok(trbuff, " "), "env.extra", env.extra);
			}
			else if (strstr(trbuff, "(FACTOR_14LJ)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(FACTOR_14LJ)", strtok(trbuff, " "), "env.factor_14lj", env.factor_14lj);
			}
			else if (strstr(trbuff, "(FINE_SCALE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(FINE_SCALE)", strtok(trbuff, " "), "env.fg_scale", env.fg_scale);
			}
			else if (strstr(trbuff, "(GA_CROSSOVER_RATE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_CROSSOVER_RATE)", strtok(trbuff, " "), "env.xover_rate", env.xover_rate);
			}
			else if (strstr(trbuff, "(GA_DELTA_E)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_DELTA_E)", strtok(trbuff, " "), "env.ga_deltaE", env.ga_deltaE);
			}
			else if (strstr(trbuff, "(GA_DIST_CENTER)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_DIST_CENTER)", strtok(trbuff, " "), "env.ga_dist_center", env.ga_dist_center);
			}
			else if (strstr(trbuff, "(GA_DIST_CENTER_EPS)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_DIST_CENTER_EPS)", strtok(trbuff, " "), "env.ga_dist_center_eps", env.ga_dist_center_eps);
			}
			else if (strstr(trbuff, "(GA_GENERATIONS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GA_GENERATIONS)", strtok(trbuff, " "), "env.generations", env.generations);
			}
			else if (strstr(trbuff, "(GA_MAX_BUCKET_POP)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GA_MAX_BUCKET_POP)", strtok(trbuff, " "), "env.pop_bucket", env.pop_bucket);
			}
			else if (strstr(trbuff, "(GA_MIGRATION_RATE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_MIGRATION_RATE)", strtok(trbuff, " "), "env.migration_rate", env.migration_rate);
			}
			else if (strstr(trbuff, "(GA_MUTATION_RATE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_MUTATION_RATE)", strtok(trbuff, " "), "env.mutation_rate", env.mutation_rate);
			}
			else if (strstr(trbuff, "(GA_OCCUPANCY_CUTOFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_OCCUPANCY_CUTOFF)", strtok(trbuff, " "), "env.occupancy", env.occupancy);
			}
			else if (strstr(trbuff, "(GA_OUTPUT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GA_OUTPUT)", strtok(trbuff, " "), "env.output", env.output);
			}
			else if (strstr(trbuff, "(GA_PHASE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GA_PHASE)", strtok(trbuff, " "), "env.ga_phase", env.ga_phase);
			}
			else if (strstr(trbuff, "(GA_RANDOM_CUT_POINTS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GA_RANDOM_CUT_POINTS)", strtok(trbuff, " "), "env.rand_cut_points", env.rand_cut_points);
			}
			else if (strstr(trbuff, "(GA_RESIDUE_MIN_ENERGY_CUTOFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_RESIDUE_MIN_ENERGY_CUTOFF)", strtok(trbuff, " "), "env.residue_check", env.residue_check);
			}
			else if (strstr(trbuff, "(GA_SEED)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GA_SEED)", strtok(trbuff, " "), "env.ga_seed", env.ga_seed);
			}
			else if (strstr(trbuff, "(GA_SHIFT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GA_SHIFT)", strtok(trbuff, " "), "env.ga_shift", env.ga_shift);
			}
			else if (strstr(trbuff, "(GA_SPHERE_FOCUS_RESID)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GA_SPHERE_FOCUS_RESID)", strtok(trbuff, " "), "env.ga_focus_resid", env.ga_focus_resid);
			}
			else if (strstr(trbuff, "(GA_SPHERE_PROBE_RADIUS)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GA_SPHERE_PROBE_RADIUS)", strtok(trbuff, " "), "env.ga_focus_probe_radius", env.ga_focus_probe_radius);
			}
			else if (strstr(trbuff, "(GRIDS_APBS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GRIDS_APBS)", strtok(trbuff, " "), "env.grids_apbs", env.grids_apbs);
			}
			else if (strstr(trbuff, "(GRIDS_DELPHI)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(GRIDS_DELPHI)", strtok(trbuff, " "), "env.grids_delphi", env.grids_delphi);
			}
			else if (strstr(trbuff, "(GRIDS_PER_ANG)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(GRIDS_PER_ANG)", strtok(trbuff, " "), "env.grids_per_ang", env.grids_per_ang);
			}
			else if (strstr(trbuff, "(H2O_SASCUTOFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(H2O_SASCUTOFF)", strtok(trbuff, " "), "env.h2o_sascutoff", env.h2o_sascutoff);
			}
			else if (strstr(trbuff, "(HDIRDIFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HDIRDIFF)", strtok(trbuff, " "), "env.hdirdiff", env.hdirdiff);
			}
			else if (strstr(trbuff, "(HDIRECTED)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(HDIRECTED)", strtok(trbuff, " "), "env.hdirected", env.hdirected);
			}
			else if (strstr(trbuff, "(HDIRLIMT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(HDIRLIMT)", strtok(trbuff, " "), "env.hdirlimt", env.hdirlimt);
			}
			else if (strstr(trbuff, "(HV_RELAX_ELEC_DIST_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_ELEC_DIST_THR)", strtok(trbuff, " "), "env.hv_relax_elec_dist_thr", env.hv_relax_elec_dist_thr);
			}
			else if (strstr(trbuff, "(HV_RELAX_HV_VDW_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_HV_VDW_THR)", strtok(trbuff, " "), "env.hv_relax_hv_vdw_thr", env.hv_relax_hv_vdw_thr);
			}
			else if (strstr(trbuff, "(HV_RELAX_NGH_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_NGH_THR)", strtok(trbuff, " "), "env.hv_relax_ngh_thr", env.hv_relax_ngh_thr);
			}
			else if (strstr(trbuff, "(HV_RELAX_INCLUDE_NGH)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(HV_RELAX_INCLUDE_NGH)", strtok(trbuff, " "), "env.hv_relax_include_ngh", env.hv_relax_include_ngh);
			}
			else if (strstr(trbuff, "(HV_RELAX_ELEC_CRG_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_ELEC_CRG_THR)", strtok(trbuff, " "), "env.hv_relax_elec_crg_thr", env.hv_relax_elec_crg_thr);
			}
			else if (strstr(trbuff, "(HV_RELAX_ELEC_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_ELEC_THR)", strtok(trbuff, " "), "env.hv_relax_elec_thr", env.hv_relax_elec_thr);
			}
			else if (strstr(trbuff, "(HV_RELAX_NCYCLE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(HV_RELAX_NCYCLE)", strtok(trbuff, " "), "env.hv_relax_ncycle", env.hv_relax_ncycle);
			}
			else if (strstr(trbuff, "(HV_RELAX_N_SHAKE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(HV_RELAX_N_SHAKE)", strtok(trbuff, " "), "env.hv_relax_n_shake", env.hv_relax_n_shake);
			}
			else if (strstr(trbuff, "(HV_RELAX_CONSTRAINT)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_CONSTRAINT)", strtok(trbuff, " "), "env.hv_relax_constraint", env.hv_relax_constraint);
			}
			else if (strstr(trbuff, "(HV_RELAX_SHAKE_TOL)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_SHAKE_TOL)", strtok(trbuff, " "), "env.hv_relax_shake_tol", env.hv_relax_shake_tol);
			}
			else if (strstr(trbuff, "(HV_RELAX_CONSTRAINT_FRC)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_CONSTRAINT_FRC)", strtok(trbuff, " "), "env.hv_relax_constraint_frc", env.hv_relax_constraint_frc);
			}
			else if (strstr(trbuff, "(HV_RELAX_NITER)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(HV_RELAX_NITER)", strtok(trbuff, " "), "env.hv_relax_niter", env.hv_relax_niter);
			}
			else if (strstr(trbuff, "(HV_RELAX_DT)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_DT)", strtok(trbuff, " "), "env.hv_relax_dt", env.hv_relax_dt);
			}
			else if (strstr(trbuff, "(HV_RELAX_VDW_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_RELAX_VDW_THR)", strtok(trbuff, " "), "env.hv_relax_vdw_thr", env.hv_relax_vdw_thr);
			}
			else if (strstr(trbuff, "(HV_TORS_SCALE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(HV_TORS_SCALE)", strtok(trbuff, " "), "env.hv_tors_scale", env.hv_tors_scale);
			}
			else if (strstr(trbuff, "(INPDB)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(INPDB)", strtok(trbuff, " "), "env.inpdb", env.inpdb);
			}
			else if (strstr(trbuff, "(IONRAD)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(IONRAD)", strtok(trbuff, " "), "env.ionrad", env.ionrad);
			}
			else if (strstr(trbuff, "(IPECE_ADD_MEM)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(IPECE_ADD_MEM)", strtok(trbuff, " "), "env.ipece", env.ipece);
			}
			else if (strstr(trbuff, "(IPECE_GRID_SPACE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(IPECE_GRID_SPACE)", strtok(trbuff, " "), "env.ipece", env.ipece);
			}
			else if (strstr(trbuff, "(IPECE_MEM_CHAINID)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(IPECE_MEM_CHAINID)", strtok(trbuff, " "), "env.ipece", env.ipece);
			}
			else if (strstr(trbuff, "(IPECE_MEM_THICKNESS)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(IPECE_MEM_THICKNESS)", strtok(trbuff, " "), "env.ipece", env.ipece);
			}
			else if (strstr(trbuff, "(IPECE_MEM_ATOM_RADIUS)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(IPECE_MEM_ATOM_RADIUS)", strtok(trbuff, " "), "env.ipece", env.ipece);
			}
			else if (strstr(trbuff, "(MCCE_HOME)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(MCCE_HOME)", strtok(trbuff, " "), "env.mcce_home", env.mcce_home);
			}
			else if (strstr(trbuff, "(MFE_CUTOFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(MFE_CUTOFF)", strtok(trbuff, " "), "env.mfe_cutoff", env.mfe_cutoff);
			}
			else if (strstr(trbuff, "(MFE_POINT)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(MFE_POINT)", strtok(trbuff, " "), "env.mfe_flag", env.mfe_flag);
			}
			else if (strstr(trbuff, "(MINIMIZE_SIZE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MINIMIZE_SIZE)", strtok(trbuff, " "), "env.minimize_size", env.minimize_size);
			}
			else if (strstr(trbuff, "(MONTE_ADV_OPT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_ADV_OPT)", strtok(trbuff, " "), "env.monte_adv_opt", env.monte_adv_opt);
			}
			else if (strstr(trbuff, "(MONTE_CONVERGE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(MONTE_CONVERGE)", strtok(trbuff, " "), "env.monte_converge", env.monte_converge);
			}
			else if (strstr(trbuff, "(MONTE_DO_ENERGY)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_DO_ENERGY)", strtok(trbuff, " "), "env.monte_do_energy", env.monte_do_energy);
			}
			else if (strstr(trbuff, "(MONTE_FLIPS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_FLIPS)", strtok(trbuff, " "), "env.monte_flips", env.monte_flips);
			}
			else if (strstr(trbuff, "(MONTE_NEQ)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_NEQ)", strtok(trbuff, " "), "env.monte_neq", env.monte_neq);
			}
			else if (strstr(trbuff, "(MONTE_NITER_MAX)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_NITER_MAX)", strtok(trbuff, " "), "env.monte_niter_max", env.monte_niter_max);
			}
			else if (strstr(trbuff, "(MONTE_NITER)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_NITER)", strtok(trbuff, " "), "env.monte_niter", env.monte_niter);
			}
			else if (strstr(trbuff, "(MONTE_NITER_MIN)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_NITER_MIN)", strtok(trbuff, " "), "env.monte_niter_min", env.monte_niter_min);
			}
			else if (strstr(trbuff, "(MONTE_NITER_CHK)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_NITER_CHK)", strtok(trbuff, " "), "env.monte_niter_chk", env.monte_niter_chk);
			}
			else if (strstr(trbuff, "(MONTE_NITER_CYCLE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_NITER_CYCLE)", strtok(trbuff, " "), "env.monte_niter_cycle", env.monte_niter_cycle);
			}
			else if (strstr(trbuff, "(MONTE_NSTART)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_NSTART)", strtok(trbuff, " "), "env.monte_nstart", env.monte_nstart);
			}
			else if (strstr(trbuff, "(MONTE_N_REDUCE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_N_REDUCE)", strtok(trbuff, " "), "env.monte_n_red", env.monte_n_red);
			}
			else if (strstr(trbuff, "(MONTE_OLD_INPUT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_OLD_INPUT)", strtok(trbuff, " "), "env.monte_old_input", env.monte_old_input);
			}
			else if (strstr(trbuff, "(MONTE_PRINT_NONZERO)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_PRINT_NONZERO)", strtok(trbuff, " "), "env.monte_print_nonzero", env.monte_print_nonzero);
			}
			else if (strstr(trbuff, "(MONTE_REDUCE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(MONTE_REDUCE)", strtok(trbuff, " "), "env.monte_reduce", env.monte_reduce);
			}
			else if (strstr(trbuff, "(MONTE_RUNS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_RUNS)", strtok(trbuff, " "), "env.monte_runs", env.monte_runs);
			}
			else if (strstr(trbuff, "(MONTE_SEED)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_SEED)", strtok(trbuff, " "), "env.monte_seed", env.monte_seed);
			}
			else if (strstr(trbuff, "(MONTE_T)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(MONTE_T)", strtok(trbuff, " "), "env.monte_temp", env.monte_temp);
			}
			else if (strstr(trbuff, "(MONTE_TRACE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_TRACE)", strtok(trbuff, " "), "env.monte_trace", env.monte_trace);
			}
			else if (strstr(trbuff, "(MONTE_TSX)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(MONTE_TSX)", strtok(trbuff, " "), "env.monte_tsx", env.monte_tsx);
			}
			else if (strstr(trbuff, "(NCONF_LIMIT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(NCONF_LIMIT)", strtok(trbuff, " "), "env.nconf_limit", env.nconf_limit);
			}
			else if (strstr(trbuff, "(NEWTPL)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(NEWTPL)", strtok(trbuff, " "), "env.new_tpl", env.new_tpl);
			}
			else if (strstr(trbuff, "(NGH_VDW_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(NGH_VDW_THR)", strtok(trbuff, " "), "env.ngh_vdw_thr", env.ngh_vdw_thr);
			}
			else if (strstr(trbuff, "(NSTATE_MAX)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(NSTATE_MAX)", strtok(trbuff, " "), "env.nstate_max", env.nstate_max);
			}
			else if (strstr(trbuff, "(N_HV_CONF_LIMIT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(N_HV_CONF_LIMIT)", strtok(trbuff, " "), "env.n_hv_conf_limit", env.n_hv_conf_limit);
			}
			else if (strstr(trbuff, "(N_INITIAL_RELAX)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(N_INITIAL_RELAX)", strtok(trbuff, " "), "env.n_initial_relax", env.n_initial_relax);
			}
			else if (strstr(trbuff, "(N_TRANS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(N_TRANS)", strtok(trbuff, " "), "env.n_trans", env.n_trans);
			}
			else if (strstr(trbuff, "(PACK)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(PACK)", strtok(trbuff, " "), "env.pack", env.pack);
			}
			else if (strstr(trbuff, "(PBE_END)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(PBE_END)", strtok(trbuff, " "), "env.pbe_end", env.pbe_end);
			}
			else if (strstr(trbuff, "(PBE_FOLDER)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(PBE_FOLDER)", strtok(trbuff, " "), "env.pbe_folder", env.pbe_folder);
			}
			else if (strstr(trbuff, "(PBE_SOLVER)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(PBE_SOLVER)", strtok(trbuff, " "), "env.pbe_solver", env.pbe_solver);
			}
			else if (strstr(trbuff, "(PBE_START)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(PBE_START)", strtok(trbuff, " "), "env.pbe_start", env.pbe_start);
			}
			else if (strstr(trbuff, "(PHI_SWING)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(PHI_SWING)", strtok(trbuff, " "), "env.phi_swing", env.phi_swing);
			}
			else if (strstr(trbuff, "(PROGRESS_LOG)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(PROGRESS_LOG)", strtok(trbuff, " "), "env.progress_log", env.progress_log);
			}
			else if (strstr(trbuff, "(PRUNE_ELE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(PRUNE_ELE)", strtok(trbuff, " "), "env.prune_ele", env.prune_ele);
			}
			else if (strstr(trbuff, "(PRUNE_RMSD)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(PRUNE_RMSD)", strtok(trbuff, " "), "env.prune_rmsd", env.prune_rmsd);
			}
			else if (strstr(trbuff, "(PRUNE_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(PRUNE_THR)", strtok(trbuff, " "), "env.prune_thr", env.prune_thr);
			}
			else if (strstr(trbuff, "(PRUNE_VDW)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(PRUNE_VDW)", strtok(trbuff, " "), "env.prune_vdw", env.prune_vdw);
			}
			else if (strstr(trbuff, "(RADIUS_PROBE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(RADIUS_PROBE)", strtok(trbuff, " "), "env.radius_probe", env.radius_probe);
			}
			else if (strstr(trbuff, "(REASSIGN)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(REASSIGN)", strtok(trbuff, " "), "env.reassign", env.reassign);
			}
			else if (strstr(trbuff, "(REBUILD_SC)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(REBUILD_SC)", strtok(trbuff, " "), "env.rebuild_sc", env.rebuild_sc);
			}
			else if (strstr(trbuff, "(RECALC_TORS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(RECALC_TORS)", strtok(trbuff, " "), "env.recalc_tors", env.recalc_tors);
			}
			else if (strstr(trbuff, "(RELAX_CLASH_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(RELAX_CLASH_THR)", strtok(trbuff, " "), "env.relax_clash_thr", env.relax_clash_thr);
			}
			else if (strstr(trbuff, "(RELAX_E_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(RELAX_E_THR)", strtok(trbuff, " "), "env.relax_e_thr", env.relax_e_thr);
			}
			else if (strstr(trbuff, "(RELAX_H)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(RELAX_H)", strtok(trbuff, " "), "env.relax_h", env.relax_h);
			}
			else if (strstr(trbuff, "(RELAX_NITER)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(RELAX_NITER)", strtok(trbuff, " "), "env.relax_niter", env.relax_niter);
			}
			else if (strstr(trbuff, "(RELAX_NSTATES)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(RELAX_NSTATES)", strtok(trbuff, " "), "env.relax_nstates", env.relax_nstates);
			}
			else if (strstr(trbuff, "(RELAX_N_HYD)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(RELAX_N_HYD)", strtok(trbuff, " "), "env.relax_n_hyd", env.relax_n_hyd);
			}
			else if (strstr(trbuff, "(RELAX_PHI)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(RELAX_PHI)", strtok(trbuff, " "), "env.relax_phi", env.relax_phi);
			}
			else if (strstr(trbuff, "(RELAX_TORQ_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(RELAX_TORQ_THR)", strtok(trbuff, " "), "env.relax_torq_thr", env.relax_torq_thr);
			}
			else if (strstr(trbuff, "(RELAX_WAT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(RELAX_WAT)", strtok(trbuff, " "), "env.relax_wat", env.relax_wat);
			}
			else if (strstr(trbuff, "(RENAME_RULES)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(RENAME_RULES)", strtok(trbuff, " "), "env.rename_rules", env.rename_rules);
			}
			else if (strstr(trbuff, "(REPACKS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(REPACKS)", strtok(trbuff, " "), "env.repacks", env.repacks);
			}
			else if (strstr(trbuff, "(REPACK_CUTOFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(REPACK_CUTOFF)", strtok(trbuff, " "), "env.repack_cutoff", env.repack_cutoff);
			}
			else if (strstr(trbuff, "(REPACK_E_THR_BURIED)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(REPACK_E_THR_BURIED)", strtok(trbuff, " "), "env.repack_e_thr_buried", env.repack_e_thr_buried);
			}
			else if (strstr(trbuff, "(REPACK_E_THR_EXPOSED)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(REPACK_E_THR_EXPOSED)", strtok(trbuff, " "), "env.repack_e_thr_exposed", env.repack_e_thr_exposed);
			}
			else if (strstr(trbuff, "(REPACK_FAV_VDW_OFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(REPACK_FAV_VDW_OFF)", strtok(trbuff, " "), "env.repack_fav_vdw_off", env.repack_fav_vdw_off);
			}
			else if (strstr(trbuff, "(ROTAMER_LIMIT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(ROTAMER_LIMIT)", strtok(trbuff, " "), "env.rotamer_limit", env.rotamer_limit);
			}
			else if (strstr(trbuff, "(ROTATE_RES)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(ROTATE_RES)", strtok(trbuff, " "), "env.rotate_res", env.rotate_res);
			}
			else if (strstr(trbuff, "(ROTATIONS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(ROTATIONS)", strtok(trbuff, " "), "env.rotations", env.rotations);
			}
			else if (strstr(trbuff, "(ROT_SPECIF)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(ROT_SPECIF)", strtok(trbuff, " "), "env.rot_specif", env.rot_specif);
			}
			else if (strstr(trbuff, "(ROT_SWAP)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(ROT_SWAP)", strtok(trbuff, " "), "env.rot_swap", env.rot_swap);
			}
			else if (strstr(trbuff, "(RXN_METHOD)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(RXN_METHOD)", strtok(trbuff, " "), "env.rxn_method", env.rxn_method);
			}
			else if (strstr(trbuff, "(SALT)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(SALT)", strtok(trbuff, " "), "env.salt", env.salt);
			}
			else if (strstr(trbuff, "(SAS2VDW)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(SAS2VDW)", strtok(trbuff, " "), "env.sas2vdw", env.sas2vdw);
			}
			else if (strstr(trbuff, "(SAS_CUTOFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(SAS_CUTOFF)", strtok(trbuff, " "), "env.sas_cutoff", env.sas_cutoff);
			}
			else if (strstr(trbuff, "(SCALE_ELE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(SCALE_ELE)", strtok(trbuff, " "), "env.scale_ele", env.scale_ele);
			}
			else if (strstr(trbuff, "(SCALE_VDW)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(SCALE_VDW)", strtok(trbuff, " "), "env.scale_vdw", env.scale_vdw);
			}
			else if (strstr(trbuff, "(SIDECHAIN_OPT)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(SIDECHAIN_OPT)", strtok(trbuff, " "), "env.sidechain_opt", env.sidechain_opt);
			}
			else if (strstr(trbuff, "(SKIP_ELE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(SKIP_ELE)", strtok(trbuff, " "), "env.skip_ele", env.skip_ele);
			}
			else if (strstr(trbuff, "(SURFACE_APBS)")) {
				fprintf(tr, "%31s%31s\t%31s%31s\n", "(SURFACE_APBS)", strtok(trbuff, " "), "env.srfm", env.srfm);
			}
			else if (strstr(trbuff, "(SWING)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(SWING)", strtok(trbuff, " "), "env.swing", env.swing);
			}
			else if (strstr(trbuff, "(TERMINALS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(TERMINALS)", strtok(trbuff, " "), "env.terminals", env.terminals);
			}
			else if (strstr(trbuff, "(TEST_SEED)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(TEST_SEED)", strtok(trbuff, " "), "env.test_seed", env.test_seed);
			}
			else if (strstr(trbuff, "(TITR_EH0)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(TITR_EH0)", strtok(trbuff, " "), "env.titr_eh0", env.titr_eh0);
			}
			else if (strstr(trbuff, "(TITR_EHD)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(TITR_EHD)", strtok(trbuff, " "), "env.titr_ehd", env.titr_ehd);
			}
			else if (strstr(trbuff, "(TITR_PH0)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(TITR_PH0)", strtok(trbuff, " "), "env.titr_ph0", env.titr_ph0);
			}
			else if (strstr(trbuff, "(TITR_PHD)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(TITR_PHD)", strtok(trbuff, " "), "env.titr_phd", env.titr_phd);
			}
			else if (strstr(trbuff, "(TITR_STEPS)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(TITR_STEPS)", strtok(trbuff, " "), "env.titr_steps", env.titr_steps);
			}
			else if (strstr(trbuff, "(TITR_TYPE)")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(TITR_TYPE)", strtok(trbuff, " "), "env.titr_type", env.titr_type);
			}
			else if (strstr(trbuff, "(TRANS_DIST)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(TRANS_DIST)", strtok(trbuff, " "), "env.trans_dist", env.trans_dist);
			}
			else if (strstr(trbuff, "(VDWF1)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(VDWF1)", strtok(trbuff, " "), "env.vdwf1", env.vdwf1);
			}
			else if (strstr(trbuff, "(VDWF2)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(VDWF2)", strtok(trbuff, " "), "env.vdwf2", env.vdwf2);
			}
			else if (strstr(trbuff, "(VDW_CUTOFF)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(VDW_CUTOFF)", strtok(trbuff, " "), "env.vdw_cutoff", env.vdw_cutoff);
			}
			else if (strstr(trbuff, "(WARN_PAIRWISE)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(WARN_PAIRWISE)", strtok(trbuff, " "), "env.warn_pairwise", env.warn_pairwise);
			}
			else if (strstr(trbuff, "(WATER_RELAX_THR)")) {
				fprintf(tr, "%31s%31s\t%31s%31f\n", "(WATER_RELAX_THR)", strtok(trbuff, " "), "env.water_relax_thr", env.water_relax_thr);
			}
			else if (strstr(trbuff, "(DO_CORRECTIONS")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(DO_CORRECTIONS)", strtok(trbuff, " "), "env.do_corrections", env.do_corrections);
			}
			else if (strstr(trbuff, "(IGNORE_INPUT_H")) {
				fprintf(tr, "%31s%31s\t%31s%31d\n", "(IGNORE_INPUT_H)", strtok(trbuff, " "), "env.ignore_input_h", env.ignore_input_h);
			}
			
		}
    	fclose(fp);
    	fclose(tr);
    }
    return 0;
}


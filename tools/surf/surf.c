#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "mcce.h"
extern ENV env;

PROT load_pdb_no_param(FILE *fp);
int atom2pdbline(char *line, ATOM atom);
float estimate_rxn(CONF *conf_p, PROT *boundary_prot);

float shell_size;

int main(int argc, char *argv[]) {
    FILE *pdb_fp;
    PROT prot;
    int i_res, i_conf;
    time_t   timer_start, timer_end;
    
    if (argc < 2) {
        printf("   Usage: load pdbfile\n");
        return -1;
    }
    
    db_open();
    if (get_env()) {
        printf("   No run.prm, use default settings\n");
    }
    env.ipece.grid_space = 0.5;
    env.ipece.boundary_extention = env.ipece.probe_radius*5.;
    shell_size = env.ipece.probe_radius*2.5;
    
    timer_start = time(NULL);

    /* load pdb file */
    if ((pdb_fp=fopen(argv[1], "r"))) {
        prot = load_pdb_no_param(pdb_fp);
        fclose(pdb_fp);
        
        if (prot.n_res == 0) {
            printf("   Fail to load pdb file %s\n",argv[1]);
            return -1;
        }
    }
    else {
        printf("   Fail to open pdb file %s\n",argv[1]);
        return -1;
    }
    
    id_conf(prot);
    setup_vdw_fast(prot);


    for (i_res=0; i_res<prot.n_res; i_res++) {
        PROT rxn_calc_prot;
        int k_res,j_res;

        /* calc. rxn in protein */
        rxn_calc_prot = new_prot();
        
        /* add the backbone of the current residue */
        k_res = ins_res(&rxn_calc_prot, rxn_calc_prot.n_res);
        ins_conf(&rxn_calc_prot.res[k_res], 0, prot.res[i_res].conf[0].n_atom);
        cpy_conf(&rxn_calc_prot.res[k_res].conf[0], &prot.res[i_res].conf[0]);
        
        /* add backbone and first conformer of nearby residues */
        for (j_res=0; j_res<prot.n_res; j_res++) {
            float d2 = (2.+shell_size)*(2.+shell_size);
            if (i_res == j_res) continue;
            if (out_of_range(prot.res[i_res].r_min, prot.res[i_res].r_max, prot.res[j_res].r_min, prot.res[j_res].r_max, d2)) continue;
            
            k_res = ins_res(&rxn_calc_prot, rxn_calc_prot.n_res);
            ins_conf(&rxn_calc_prot.res[k_res], 0, prot.res[j_res].conf[0].n_atom);
            cpy_conf(&rxn_calc_prot.res[k_res].conf[0], &prot.res[j_res].conf[0]);
            
            if (prot.res[j_res].n_conf >1) {
                ins_conf(&rxn_calc_prot.res[k_res], 1, prot.res[j_res].conf[1].n_atom);
                cpy_conf(&rxn_calc_prot.res[k_res].conf[1], &prot.res[j_res].conf[1]);
            }
        }
        
        for (i_conf = 1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            if (rxn_calc_prot.res[0].n_conf > 1) {
                del_conf(&rxn_calc_prot.res[0], 1);
            }
            
            ins_conf(&rxn_calc_prot.res[0], 1, prot.res[i_res].conf[i_conf].n_atom);
            cpy_conf(&rxn_calc_prot.res[0].conf[1], &prot.res[i_res].conf[i_conf]);
            
            prot.res[i_res].conf[i_conf].E_rxn = estimate_rxn(&rxn_calc_prot.res[0].conf[1], &rxn_calc_prot);
            //printf("%s %8.3f\n",prot.res[i_res].conf[i_conf].uniqID, estimate_rxn(&rxn_calc_prot.res[0].conf[1], &rxn_calc_prot));
        }
        
        del_prot(&rxn_calc_prot);
        
        
        /* calculate solution (reference) rxn */
        rxn_calc_prot = new_prot();

        /* add the backbone of the current residue */
        k_res = ins_res(&rxn_calc_prot, rxn_calc_prot.n_res);
        ins_conf(&rxn_calc_prot.res[k_res], 0, prot.res[i_res].conf[0].n_atom);
        cpy_conf(&rxn_calc_prot.res[k_res].conf[0], &prot.res[i_res].conf[0]);

        for (i_conf = 1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            float ref_rxn;
            while (rxn_calc_prot.res[0].n_conf > 1) {
                del_conf(&rxn_calc_prot.res[0], 1);
            }
            
            ins_conf(&rxn_calc_prot.res[0], 1, prot.res[i_res].conf[i_conf].n_atom);
            cpy_conf(&rxn_calc_prot.res[0].conf[1], &prot.res[i_res].conf[i_conf]);
            
            //prot.res[i_res].conf[i_conf].E_rxn = estimate_rxn(&rxn_calc_prot.res[0].conf[1], &rxn_calc_prot);
            ref_rxn = estimate_rxn(&rxn_calc_prot.res[0].conf[1], &rxn_calc_prot);
            printf("%s rxn=%8.3f ref=%8.3f desolv=%8.3f\n", prot.res[i_res].conf[i_conf].uniqID, prot.res[i_res].conf[i_conf].E_rxn, ref_rxn, prot.res[i_res].conf[i_conf].E_rxn-ref_rxn);
        }
        
        del_prot(&rxn_calc_prot);
        
    }
    timer_end = time(NULL);
    printf("Time used: %ld seconds.\n", timer_end-timer_start); fflush(stdout);

    db_close();
    
    return 0;
}

float estimate_rxn(CONF *conf_p, PROT *boundary_prot) {
    INT_VECT k_grid;
    int i_atom;

    probe(*boundary_prot, &env.ipece);
    
    /* initialize the image crg array */
    env.ipece.image_crg = malloc(env.ipece.n_grid.x * sizeof(float **));
    for (k_grid.x = 0; k_grid.x<env.ipece.n_grid.x; k_grid.x++) {
        env.ipece.image_crg[k_grid.x] = malloc(env.ipece.n_grid.y * sizeof(float *));
        for (k_grid.y = 0; k_grid.y<env.ipece.n_grid.y; k_grid.y++) {
            env.ipece.image_crg[k_grid.x][k_grid.y] = malloc(env.ipece.n_grid.z * sizeof(float));
            for (k_grid.z = 0; k_grid.z<env.ipece.n_grid.z; k_grid.z++) {
                env.ipece.image_crg[k_grid.x][k_grid.y][k_grid.z] = 0.;
            }
        }
    }
    
    /* assign image charges */
    for (i_atom=0; i_atom<conf_p->n_atom; i_atom++) {
        INT_VECT corner1, corner2;
        VECTOR sum_rad;
        
        if (!conf_p->atom[i_atom].on) continue;
        
        sum_rad.x = sum_rad.y = sum_rad.z = conf_p->atom[i_atom].rad+env.ipece.probe_radius+shell_size;
        corner1 = coor2grid(vector_vminusv(conf_p->atom[i_atom].xyz, sum_rad), &env.ipece);
        corner2 = coor2grid(vector_vplusv(conf_p->atom[i_atom].xyz, sum_rad), &env.ipece);
        
        float n_grid_surf = 4./3.*env.PI*(pow(sum_rad.x,3)-pow(conf_p->atom[i_atom].rad+env.ipece.probe_radius,3)) / pow(env.ipece.grid_space,3);
        //float n_grid_surf = 4./3.*env.PI*(pow(sum_rad.x,3)-pow(conf_p->atom[i_atom].rad+env.ipece.probe_radius,3)) * (conf_p->atom[i_atom].rad+env.ipece.probe_radius) / pow(env.ipece.grid_space,3);
        
        for (k_grid.x = corner1.x; k_grid.x < corner2.x; k_grid.x++) {
            for (k_grid.y = corner1.y; k_grid.y < corner2.y; k_grid.y++) {
                for (k_grid.z = corner1.z; k_grid.z < corner2.z; k_grid.z++) {
                    if (*reach_label(k_grid, &env.ipece) == 'o' || *reach_label(k_grid, &env.ipece) == 'c') {
                        float d;
                        d = dvv(conf_p->atom[i_atom].xyz, grid2coor(k_grid, &env.ipece));
                        if (d>sum_rad.x) continue;
                        
                        env.ipece.image_crg[k_grid.x-env.ipece.grid_boundary_lower.x][k_grid.y-env.ipece.grid_boundary_lower.y][k_grid.z-env.ipece.grid_boundary_lower.z]
                        -= conf_p->atom[i_atom].crg/(n_grid_surf*d);
                    }
                }
            }
        }
        
    }
    
    /* calc. charge-image interactions */
    float sum_e = 0.;
    for (i_atom=0; i_atom<conf_p->n_atom; i_atom++) {
        float e = 0.;
        if (fabs(conf_p->atom[i_atom].crg) < 1e-4) continue;
        for (k_grid.x = env.ipece.grid_boundary_lower.x; k_grid.x<env.ipece.grid_boundary_higher.x; k_grid.x++) {
            for (k_grid.y = env.ipece.grid_boundary_lower.y; k_grid.y<env.ipece.grid_boundary_higher.y; k_grid.y++) {
                for (k_grid.z = env.ipece.grid_boundary_lower.z; k_grid.z<env.ipece.grid_boundary_higher.z; k_grid.z++) {
                    float d;
                    if (*reach_label(k_grid, &env.ipece) != 'o' && *reach_label(k_grid, &env.ipece) != 'c') continue;
                    
                    d = dvv(conf_p->atom[i_atom].xyz, grid2coor(k_grid, &env.ipece));
                    
                    e += 331.5*conf_p->atom[i_atom].crg*
                    env.ipece.image_crg[k_grid.x-env.ipece.grid_boundary_lower.x][k_grid.y-env.ipece.grid_boundary_lower.y][k_grid.z-env.ipece.grid_boundary_lower.z]/(0.5259 * d);
                    
                    
                }
            }
        }
        
        sum_e += e;
    }
    
    /* free up memory */
    free_probe(&env.ipece);
    for (k_grid.x = 0; k_grid.x<env.ipece.n_grid.x; k_grid.x++) {
        for (k_grid.y = 0; k_grid.y<env.ipece.n_grid.y; k_grid.y++) {
            free(env.ipece.image_crg[k_grid.x][k_grid.y]);
        }
        free(env.ipece.image_crg[k_grid.x]);
    }
    free(env.ipece.image_crg);
    
    return sum_e;

}

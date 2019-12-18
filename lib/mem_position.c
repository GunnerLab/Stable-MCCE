#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mcce.h"

#define PROBE_RADIUS 1.4

int get_scored_atoms(PROT prot, IPECE *ipece);
int free_scored_atoms(IPECE *ipece);
extern long idum;

/* Move the protein to the position so the membrane is between
z=+/- half_membrane_thichness */
int mem_position(PROT *prot_p, IPECE *ipece)
{
    VECTOR *orig_r;
    ATOM **all_atoms; /* a list of pointers pointing to all atoms in prot */
    int i_res, i_conf, i_atom, na;
    ATOM *atom_p;

    VECTOR vec_orig, vec_i, vec_j, vec_k;
    vec_orig.x = 0.;    vec_orig.y = 0.;    vec_orig.z = 0.;
    vec_i.x=1.;         vec_i.y=0.;         vec_i.z=0.;
    vec_j.x=0.;         vec_j.y=1.;         vec_j.z=0.;
    vec_k.x=0.;         vec_k.y=0.;         vec_k.z=1.;
    
    if (env.test_seed < 0) idum = time(NULL); //allows random numbers to be fixed for testing
    else idum = env.test_seed;
    
    /* set up. add a residue to contain anchor vectors,
    initialize all_atoms array, backup all atom coordinates to orig_r */
    ins_res(prot_p, prot_p->n_res);
    ins_conf(&prot_p->res[prot_p->n_res-1],0,4);
    for (i_atom=0;i_atom<4;i_atom++) {
        prot_p->res[prot_p->n_res-1].conf[0].atom[i_atom].on = 1;
    }
    prot_p->res[prot_p->n_res-1].conf[0].atom[0].xyz = vec_orig;
    prot_p->res[prot_p->n_res-1].conf[0].atom[1].xyz = vector_rescale(vec_i,100.);
    prot_p->res[prot_p->n_res-1].conf[0].atom[2].xyz = vector_rescale(vec_j,100.);
    prot_p->res[prot_p->n_res-1].conf[0].atom[3].xyz = vector_rescale(vec_k,100.);
    
    na = 0;
    all_atoms = NULL;
    orig_r = NULL;
    for (i_res=0; i_res<prot_p->n_res; i_res++) {
        for (i_conf=0; i_conf<prot_p->res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot_p->res[i_res].conf[i_conf].n_atom; i_atom++) {
                atom_p = &prot_p->res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                na++;
                orig_r = realloc(orig_r, na*sizeof(VECTOR));
                orig_r[na-1] = atom_p->xyz;
                all_atoms = realloc(all_atoms, na*sizeof(ATOM *));
                all_atoms[na-1] = atom_p;
                atom_p->i_atom_prot = na-1;
            }
        }
    }
    
    /* move the center of the protein to the origin */
    VECTOR center;
    center = vec_orig;
    int ia;
    for (ia=0; ia<na; ia++) {
        center = vector_vplusv(center, all_atoms[ia]->xyz);
    }
    center = vector_rescale(center, 1./na);
    for (ia=0; ia<na; ia++) {
        all_atoms[ia]->xyz = vector_vminusv(all_atoms[ia]->xyz, center);
    }
    
    float score, min_score;
    int best_pos;
    /* get the initial score for the protein */
    probe(*prot_p, ipece);
    /* collect atoms used for calculating scores */
    get_scored_atoms(*prot_p, ipece);
    min_score = calc_score(ipece);
    best_pos = 1;
    
    /* rotate around x axis by 90¡ to calculate score again */
    rotate_atoms(vec_orig, vec_i, env.PI/2., ipece->n_scored_atom, ipece->scored_atoms);
    score = calc_score(ipece);
    if (score < min_score) {
        min_score = score;
        best_pos = 2;
    }
    rotate_atoms(vec_orig, vec_i, -env.PI/2., ipece->n_scored_atom, ipece->scored_atoms);
    
    /* rotate around y axis by 90¡ to calculate score again */
    rotate_atoms(vec_orig, vec_j, env.PI/2., ipece->n_scored_atom, ipece->scored_atoms);
    score = calc_score(ipece);
    if (score < min_score) {
        min_score = score;
        best_pos = 3;
    }
    rotate_atoms(vec_orig, vec_j, -env.PI/2., ipece->n_scored_atom, ipece->scored_atoms);
    
    /* choose the lowest scored position to start with */
    if (best_pos == 2) {
        rotate_atoms(vec_orig, vec_i, env.PI/2., na, all_atoms);
        rotate_atoms(vec_orig, vec_i, env.PI/2., ipece->n_scored_atom, ipece->scored_atoms);
    }
    if (best_pos == 3) {
        rotate_atoms(vec_orig, vec_j, env.PI/2., na, all_atoms);
        rotate_atoms(vec_orig, vec_j, env.PI/2., ipece->n_scored_atom, ipece->scored_atoms);
    }
    
    /* randomly translate the protein in z direction or rotation around an axis
    on x-y plane to find the lowest scored position.
    sampling the positions in a way similar as monte carlo */
    VECTOR vec_trans = vec_orig;
    VECTOR vec_axis = vec_orig;
    
    score = calc_score(ipece);
    min_score = score;
    for (ia=0;ia<4;ia++) {
        ipece->membrane_position[ia] = prot_p->res[prot_p->n_res-1].conf[0].atom[ia].xyz;
    }
    
    int iter;
    for (iter=0;iter<ipece->n_iteration;iter++) {
        float new_score, delta_score;
        /* translation */
        if (ran2(&idum) < 0.5) {
            vec_trans.z =  2.*(ran2(&idum) - 0.5)*ipece->translation_max;
            
            translate_atoms(vec_trans, ipece->n_scored_atom, ipece->scored_atoms);
            new_score = calc_score(ipece);
            delta_score = new_score - score;
            if (ran2(&idum) < exp(-ipece->beta * delta_score)) {
                translate_atoms(vec_trans, na, all_atoms);
                score = new_score;
            }
            else {
                vec_trans.z = -vec_trans.z;
                translate_atoms(vec_trans, ipece->n_scored_atom, ipece->scored_atoms);
            }
        }
        /* rotation */
        else {
            float phi, theta;
            phi = ran2(&idum)*2.*env.PI;
            theta = ran2(&idum)*ipece->rotation_max;    /* rotation_max has been
                                                        converted to radian in
                                                        parameter reading process */
            vec_axis.x = cos(phi);
            vec_axis.y = sin(phi);
            rotate_atoms(vec_orig, vec_axis, theta, ipece->n_scored_atom, ipece->scored_atoms);
            
            new_score = calc_score(ipece);
            delta_score = new_score - score;
            if (ran2(&idum) < exp(-ipece->beta * delta_score)) {
                rotate_atoms(vec_orig, vec_axis, theta, na, all_atoms);
                score = new_score;
            }
            else {
                vec_trans.z = -vec_trans.z;
                rotate_atoms(vec_orig, vec_axis, -theta, ipece->n_scored_atom, ipece->scored_atoms);
            }

        }
        
        if (score < min_score) {
            min_score = score;
            for (ia=0;ia<4;ia++) {
                ipece->membrane_position[ia] = prot_p->res[prot_p->n_res-1].conf[0].atom[ia].xyz;
            }
        }
        
        /*
        printf("\n");
        for (i=0;i<4;i++) {
            printf("       %10.6f %10.6f %10.6f\n",
            prot_p->res[prot_p->n_res-1].conf[0].atom[i].xyz.x,
            prot_p->res[prot_p->n_res-1].conf[0].atom[i].xyz.y,
            prot_p->res[prot_p->n_res-1].conf[0].atom[i].xyz.z);
        }
        */
    }
    
    ipece->mem_position_defined = 1;
    
    /* move the protein back to orginal position */
    for (ia=0; ia<na; ia++) {
        all_atoms[ia]->xyz = orig_r[ia];
    }

    /* free memory */
    free_scored_atoms(ipece);
    free_probe(ipece);
    free(all_atoms);
    free(orig_r);
    na = 0;
    
    /* delete the extra anchor residue */
    del_res(prot_p, prot_p->n_res-1);
    return 0;
}

int get_scored_atoms(PROT prot, IPECE *ipece)
{
    int i_res, i_conf, i_atom;
    int exp_dist_ngrid;  /* number of grids to be searched for surface, see the description in the later code */
    
    /* initialize scored atom list */
    ipece->n_scored_atom = 0;
    
    /* convert surface_exp_dist to exp_dist_ngrid */
    exp_dist_ngrid = (int)(ipece->surface_exp_dist/ipece->grid_space + 0.5);
    
    surfw(prot, PROBE_RADIUS);
    
    for (i_res=0; i_res<prot.n_res; i_res++) {
        /* only score one sidechain conformer */
        int n_conf = 2;
        if (prot.res[i_res].n_conf<n_conf)
            n_conf=prot.res[i_res].n_conf;
        for (i_conf=0; i_conf<n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                float score;
                ATOM *atom_p;
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                
                if ( !param_get("IPECE_SC",prot.res[i_res].resName, atom_p->name, &score) ) {
                    /* only put the atom into the score list if it is within
                    surface_exp_dist to the solvent. this condition is added for
                    the outer membrane proteins, such as porin, where there are
                    many exposed ionizable residues in the transmembrane region
                    facing the pore side. These residues will not contribute to
                    the score if the protein is aligned in the right orientation */
                    
                    INT_VECT atom_grid;
                    atom_grid = coor2grid(atom_p->xyz, ipece);
                    int added = 0;  /* flag to indicate whether atom_p has been added to the scored atom list */
                    INT_VECT probe_grid;
                    for (probe_grid.x = atom_grid.x-exp_dist_ngrid;
                    probe_grid.x <= atom_grid.x+exp_dist_ngrid;
                    probe_grid.x++) {
                        for (probe_grid.y = atom_grid.y-exp_dist_ngrid;
                        probe_grid.y <= atom_grid.y+exp_dist_ngrid;
                        probe_grid.y++) {
                            for (probe_grid.z = atom_grid.z-exp_dist_ngrid;
                            probe_grid.z <= atom_grid.z+exp_dist_ngrid;
                            probe_grid.z++) {
                                if ( *reach_label(probe_grid, ipece) == 's' ||
                                *reach_label(probe_grid, ipece) == 'o' ) {
                                    ipece->n_scored_atom++;
                                    
                                    ipece->scored_atoms = realloc(ipece->scored_atoms, ipece->n_scored_atom*sizeof(ATOM *));
                                    ipece->scored_atoms[ipece->n_scored_atom-1] = malloc(sizeof(ATOM));
                                    memcpy(ipece->scored_atoms[ipece->n_scored_atom-1], atom_p, sizeof(ATOM));
                                    
                                    ipece->atom_scores = realloc(ipece->atom_scores, ipece->n_scored_atom*sizeof(float));
                                    ipece->atom_scores[ipece->n_scored_atom-1] = score*atom_p->sas;
                                    
                                    added = 1;
                                    break;
                                }
                            }
                            if (added) break;
                        }
                        if (added) break;
                    }
                }
                else
                    continue;
            }
        }
    }
    
    /*
    int ia;
    printf("    List of atoms being scored:\n");
    for (ia=0; ia<ipece->n_scored_atom; ia++) {
        printf("%s acc=%8.3f\n",
        ipece->scored_atoms[ia]->name,ipece->scored_atoms[ia]->sas);
    }
    */
    return 0;
}

int free_scored_atoms(IPECE *ipece)
{
    int ia;
    for (ia=0; ia<ipece->n_scored_atom; ia++) {
        free(ipece->scored_atoms[ia]);
    }
    free(ipece->scored_atoms);
    ipece->n_scored_atom = 0;
    return 0;
}

float calc_score(IPECE *ipece)
{
    float score = 0.;
    int ia;
    for (ia=0; ia<ipece->n_scored_atom; ia++) {
        if (fabs(ipece->scored_atoms[ia]->xyz.z)<ipece->half_mem_thickness)
            score += ipece->atom_scores[ia];
    }
    
    return score;
}

int translate_atoms(VECTOR vec, int na, ATOM **atoms_p)
{
    int ia;
    for (ia=0; ia<na; ia++) {
        atoms_p[ia]->xyz = vector_vplusv(atoms_p[ia]->xyz, vec);
    }
    return 0;
}

int rotate_atoms(VECTOR v0, VECTOR v1, float angle, int na, ATOM **atoms_p)
{
    LINE axis = line_2v(v0, v1);
    GEOM op;
    int ia;
    geom_reset(&op);
    geom_roll(&op, angle, axis);
    for (ia=0; ia<na; ia++) {
        geom_apply(op, &(atoms_p[ia]->xyz));
    }
    return 0;
}

    

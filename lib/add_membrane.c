#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mcce.h"

int anchor2defined(PROT *prot_p, VECTOR *pos);

int add_membrane(PROT *prot_p, IPECE *ipece)
{
    int i_res, i_conf, i_atom;

    /* vectors for origin, and in x, y, z direction. */
    VECTOR vec_anchor[4];
    vec_anchor[0].x = 0.;    vec_anchor[0].y = 0.;    vec_anchor[0].z = 0.;
    vec_anchor[1].x = 100.;  vec_anchor[1].y = 0.;    vec_anchor[1].z = 0.;
    vec_anchor[2].x = 0.;    vec_anchor[2].y = 100.;  vec_anchor[2].z = 0.;
    vec_anchor[3].x = 0.;    vec_anchor[3].y = 0.;    vec_anchor[3].z = 100.;

    /* set up. add a residue to contain anchor vectors,
    initialize all_atoms array, backup all atom coordinates to orig_r */
    ins_res(prot_p, prot_p->n_res);
    ins_conf(&prot_p->res[prot_p->n_res-1],0,4);
    for (i_atom=0;i_atom<4;i_atom++) {
        prot_p->res[prot_p->n_res-1].conf[0].atom[i_atom].on = 1;
        prot_p->res[prot_p->n_res-1].conf[0].atom[i_atom].xyz = vec_anchor[i_atom];
    }
    
    /* if membrane postion is pre-defined in the mem_pos file, move the protein
    so that the membrane can be placed on the defined position */
    FILE *fp;
    fp = fopen(MEM_POS,"r");
    if (!fp) {
        ipece->mem_position_defined = 0;
    }
    else {
        char sbuff[MAXCHAR_LINE];              /* string buffer */
        ipece->mem_position_defined = 1;
        int i = 0;
        while (fgets(sbuff, MAXCHAR_LINE*sizeof(char), fp)) {
            strip(sbuff, sbuff);
            if (strlen(sbuff) == 0) continue;
            
            sscanf(sbuff, "%lf %lf %lf",
            &ipece->membrane_position[i].x,
            &ipece->membrane_position[i].y,
            &ipece->membrane_position[i].z);
            
            i++;
            if (i >=4) break;
        }
        if (i <= 3) ipece->mem_position_defined = 0;
        fclose(fp);
    }
    
    if (!ipece->mem_position_defined) {
        printf("   Error in mem_pos file, no membrane added.\n");
        return -1;
    }
    
    /* move the protein */
    anchor2defined(prot_p, ipece->membrane_position);
    
    /* print out the current position of the anchors */
    ATOM *prot_anchor;
    
    prot_anchor = prot_p->res[prot_p->n_res-1].conf[0].atom;
    printf("    The protein is moved to the current position:\n");
    int i;
    for (i = 0; i<4; i++) {
        printf("       %10.6f %10.6f %10.6f\n",
        prot_anchor[i].xyz.x,
        prot_anchor[i].xyz.y,
        prot_anchor[i].xyz.z);
    }
    
    float save_boundary_extention = ipece->boundary_extention;
    ipece->boundary_extention = ipece->membrane_size;
    probe(*prot_p, ipece);
    ipece->boundary_extention = save_boundary_extention;
    
    /* define the boundary for membrane atoms: in z direction, membrane is located
    from z=-half_membrane_thichness to z=half_membrane_thichness; in x and y direction
    a box is made extending from the boundary of the protein by the size of the membrane
    similar as the subroutine create_grid_box, however the boundary of transmembrane
    region is used instead of the whole protein */

    VECTOR  mem_coor_boundary_lower;    /* lower  boundary of the membrane box defined in coordinates */
    VECTOR  mem_coor_boundary_higher;   /* higher boundary of the membrane box defined in coordinates */
    INT_VECT  mem_grid_boundary_lower;  /* lower  boundary of the membrane box defined in grid index */
    INT_VECT  mem_grid_boundary_higher; /* higher boundary of the membrane box defined in grid index */
    /* Initialize the variables for lower and higher boundary of the box */
    int initialized = 0;
    for (i_res=0; i_res<prot_p->n_res-1; i_res++) { /* the last residue is anchor residue, is not counted */
        for (i_conf=0; i_conf<prot_p->res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot_p->res[i_res].conf[i_conf].n_atom; i_atom++) {
                ATOM *atom_p = &prot_p->res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                if (fabs(atom_p->xyz.z) > ipece->half_mem_thickness) continue;
                mem_coor_boundary_lower = atom_p->xyz;
                mem_coor_boundary_higher = atom_p->xyz;

                initialized = 1;
                break;
            }
            if (initialized) break;
        }
        if (initialized) break;
    }
    
    /* Find the size of the box to contain the protein */
    /* Define lower and higher boundary of the box */
    mem_coor_boundary_lower.z = -ipece->half_mem_thickness;
    mem_coor_boundary_higher.z = ipece->half_mem_thickness;
    for (i_res=0; i_res<prot_p->n_res-1; i_res++) { /* the last residue is anchor residue, is not counted */
        for (i_conf=0; i_conf<prot_p->res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot_p->res[i_res].conf[i_conf].n_atom; i_atom++) {
                ATOM *atom_p = &prot_p->res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                if (fabs(atom_p->xyz.z) > ipece->half_mem_thickness) continue;
                
                if (atom_p->xyz.x < mem_coor_boundary_lower.x) mem_coor_boundary_lower.x = atom_p->xyz.x;
                if (atom_p->xyz.y < mem_coor_boundary_lower.y) mem_coor_boundary_lower.y = atom_p->xyz.y;
                
                if (atom_p->xyz.x > mem_coor_boundary_higher.x) mem_coor_boundary_higher.x = atom_p->xyz.x;
                if (atom_p->xyz.y > mem_coor_boundary_higher.y) mem_coor_boundary_higher.y = atom_p->xyz.y;
            }
        }
    }
    
    /* Expand the boundary of the box a little */
    mem_coor_boundary_lower.x  -= ipece->membrane_size;
    mem_coor_boundary_lower.y  -= ipece->membrane_size;
    
    mem_coor_boundary_higher.x += ipece->membrane_size;
    mem_coor_boundary_higher.y += ipece->membrane_size;
    
    /* Define lower and higher boundary of the grid box in grid index*/
    mem_grid_boundary_lower = coor2grid(mem_coor_boundary_lower,ipece);
    mem_grid_boundary_higher = coor2grid(mem_coor_boundary_higher,ipece);;
    
    /* initialize residue and atom index for membrane atoms: */
    int k_res = prot_p->n_res-1;    /* n_res-1 is the anchor residue, insert membrane before it */
    int n_mem_res = 1;
    int k_atom = 0;
    
    /* loop over all grid points in transmembrane region */
    int incre_grid = (int)(ipece->mem_separation/ipece->grid_space);
    INT_VECT mem_grid;
    for (mem_grid.x = mem_grid_boundary_lower.x;
    mem_grid.x <= mem_grid_boundary_higher.x;
    mem_grid.x+=incre_grid) {
        for (mem_grid.y = mem_grid_boundary_lower.y;
        mem_grid.y <= mem_grid_boundary_higher.y;
        mem_grid.y+=incre_grid) {
            for (mem_grid.z = mem_grid_boundary_lower.z;
            mem_grid.z <= mem_grid_boundary_higher.z;
            mem_grid.z+=incre_grid) {
                /* only add membrane to grids marked 'o' */
                if (*reach_label(mem_grid, ipece) != 'o') {
                    continue;
                }

                /* add a residue with 1000 atoms in the backbone conformer
                if this is a new residue */
                if (k_res>=prot_p->n_res-1) {
                    k_res = prot_p->n_res-1;    /* insert membrane residue before the anchor residue */
                    ins_res(prot_p, k_res);
                    prot_p->res[k_res].iCode = '_';
                    strcpy(prot_p->res[k_res].resName, ipece->mem_resName);
                    prot_p->res[k_res].chainID = ipece->mem_chainID;
                    prot_p->res[k_res].resSeq = n_mem_res;
                    
                    ins_conf(&prot_p->res[k_res],0,1000);
                    strcpy(prot_p->res[k_res].conf[0].confName, ipece->mem_resName);
                    strcpy(prot_p->res[k_res].conf[0].confName+3, "BK");
                    prot_p->res[k_res].conf[0].altLoc = ' ';
                    strcpy(prot_p->res[k_res].conf[0].history, "BK________");
                }
                
                /* add membrane atom */
                prot_p->res[k_res].conf[0].atom[k_atom].on = 1;
                /* name the atom as 0C00, 0C01, ... 1C00, 1C01, ... 9C99 */
                sprintf(prot_p->res[k_res].conf[0].atom[k_atom].name+1,"%03i",k_atom);
                prot_p->res[k_res].conf[0].atom[k_atom].name[0]
                = prot_p->res[k_res].conf[0].atom[k_atom].name[1];
                prot_p->res[k_res].conf[0].atom[k_atom].name[1] = 'C';
                prot_p->res[k_res].conf[0].atom[k_atom].name[4] = '\0';

                prot_p->res[k_res].conf[0].atom[k_atom].xyz = grid2coor(mem_grid, ipece);
                prot_p->res[k_res].conf[0].atom[k_atom].rad = ipece->mem_atom_radius;

                /* update residue and atom index numbers */
                k_atom++;
                if (k_atom >=1000) {
                    k_res++;
                    n_mem_res++;
                    k_atom = 0;
                }
            }
        }
    }
    
    /* move the protein back */
    anchor2defined(prot_p, vec_anchor);
    
    /* free memory for grids */
    free_probe(ipece);
    
    /* delete the extra anchor residue */
    del_res(prot_p, prot_p->n_res-1);
    return 0;
}

/* move protein so that the anchor is aligned to a defined orientation */
int anchor2defined(PROT *prot_p, VECTOR *pos) {
    float cos_theta, theta;

    int i_res, i_conf, i_atom, na, ia;
    ATOM **all_atoms; /* a list of pointers pointing to all atoms in prot */
    ATOM *prot_anchor;

    prot_anchor = prot_p->res[prot_p->n_res-1].conf[0].atom;
    na = 0;
    all_atoms = NULL;
    for (i_res=0; i_res<prot_p->n_res; i_res++) {
        for (i_conf=0; i_conf<prot_p->res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot_p->res[i_res].conf[i_conf].n_atom; i_atom++) {
                ATOM *atom_p = &prot_p->res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                na++;
                all_atoms = realloc(all_atoms, na*sizeof(ATOM *));
                all_atoms[na-1] = atom_p;
            }
        }
    }
    
    /* move anchor[0] to pos[0] */
    VECTOR shift = vector_vminusv(pos[0], prot_anchor[0].xyz);
    for (ia=0; ia<na; ia++) {
        all_atoms[ia]->xyz = vector_vplusv(all_atoms[ia]->xyz, shift);
    }
    
    /* rotate anchor[1] to be aligned with pos[1] */
    /* get i direction given by the defined orientaion (pos) and the current
    protein position (prot_anchor) */
    VECTOR vec_axis;
    VECTOR pos_i, pos_j, pos_k;
    VECTOR anchor_i, anchor_j;
    
    pos_i = vector_normalize(vector_vminusv(pos[1],pos[0]));
    anchor_i = vector_normalize(vector_vminusv(prot_anchor[1].xyz, prot_anchor[0].xyz));

    cos_theta = vdotv(anchor_i, pos_i);
    theta = acos(cos_theta);
    
    vec_axis = vector_vplusv(pos[0],vector_vxv(anchor_i, pos_i));
    
    rotate_atoms(pos[0], vec_axis, theta, na, all_atoms);
    
    /* rotate anchor[2] to be aligned with pos[2] */
    pos_j = vector_normalize(vector_vminusv(pos[2],pos[0]));
    anchor_j = vector_normalize(vector_vminusv(prot_anchor[2].xyz, prot_anchor[0].xyz));

    cos_theta = vdotv(anchor_j, pos_j);
    theta = acos(cos_theta);
    
    pos_k = vector_normalize(vector_vminusv(pos[3],pos[0]));
    if (vdotv(vector_vxv(anchor_j,pos_j),pos_i) <0.)
        theta = 2.*env.PI - theta;

    rotate_atoms(pos[0], pos[1], theta, na, all_atoms);

    free(all_atoms);
    na = 0;
    return 0;
}


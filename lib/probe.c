#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

int probe(PROT prot, IPECE *ipece)
{
    int i_res, i_conf, i_atom;
    int i_probe, label_updated;
    INT_VECT probe_grid;
    ATOM *atom_p;
    
    /* Make a grid box using the protein */
    create_grid_box(prot, ipece);
    
    /* The grid points within radius of an atom are labeled 'p', those further
    than radius of an atom but within (atom radius + probe radius) are labeled
    's' for surface*/
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                int box_size;
                int i,j,k;
                VECTOR corner_coord, probe_coord;
                INT_VECT corner_grid, probe_grid;
                
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                if (atom_p->rad < 0.001) continue;
                
                /* a small grid box around each atom is probed */
                /* the number of grids covered by this box is box_size */
                box_size = (int)(2.*(atom_p->rad + ipece->probe_radius)/ipece->grid_space) + 1;
                
                /* define a corner of the box */
                corner_coord.x = atom_p->xyz.x - atom_p->rad - ipece->probe_radius - 0.5;
                corner_coord.y = atom_p->xyz.y - atom_p->rad - ipece->probe_radius - 0.5;
                corner_coord.z = atom_p->xyz.z - atom_p->rad - ipece->probe_radius - 0.5;
                
                /* convert to grid space */
                corner_grid = coor2grid(corner_coord, ipece);
                
                /* two r^2 are calculated for comparison later */
                float radsq = atom_p->rad*atom_p->rad;
                float sumrad_sq = (atom_p->rad + ipece->probe_radius)*(atom_p->rad + ipece->probe_radius);
                
                /* probe the small box in three dimensions */
                for (i=0;i<box_size;i++) {
                    probe_grid.x = corner_grid.x+i;
                    for (j=0;j<box_size;j++) {
                        probe_grid.y = corner_grid.y+j;
                        for (k=0;k<box_size;k++) {
                            probe_grid.z = corner_grid.z+k;
                            
                            probe_coord = grid2coor(probe_grid, ipece);
                            
                            /* label grid point 'p' if its distance to the center
                            of an atom is within the radius */
                            if (ddvv(atom_p->xyz, probe_coord) < radsq) {
                                *reach_label(probe_grid, ipece) = 'p';
                            }
                            /* label grid point 's' if its distance to the center
                            of an atom is longer than the radius but smaller than
                            the sum of probe and atom radii
                            's' is only labeled if the grid point is not labeled yet */
                            else if (ddvv(atom_p->xyz, probe_coord) < sumrad_sq) {
                                if (*reach_label(probe_grid, ipece) == '\0') {
                                    *reach_label(probe_grid, ipece) = 's';
                                }
                            }
                        }
                    }
                }
                
            }
        }
    }
    
    /* Starting from [0,0,0] corner, label the grids that are easy to access from outside.
    Loop over x, y, and z direction. No break in x and y loops so that all grids
    in the x-y plane are covered. In z direction, the loop is stopped when the first
    non-outside grid is reached. */
    for (probe_grid.x = ipece->grid_boundary_lower.x;
    probe_grid.x <= ipece->grid_boundary_higher.x;
    probe_grid.x++) {
        for (probe_grid.y = ipece->grid_boundary_lower.y;
        probe_grid.y <= ipece->grid_boundary_higher.y;
        probe_grid.y++) {
            for (probe_grid.z = ipece->grid_boundary_lower.z;
            probe_grid.z <= ipece->grid_boundary_higher.z;
            probe_grid.z++) {
                if (*reach_label(probe_grid, ipece) == 'o') continue;
                /* if the grid is not labeled yet, label it as outside */
                if (*reach_label(probe_grid, ipece) == '\0') {
                    *reach_label(probe_grid, ipece) = 'o';
                }
                else {
                    break;
                }
            }
            /* no break in x and y direction in order to cover the entire x-y area */
        }
    }
    /* do the same to in x direction for grids on y-z plane */
    for (probe_grid.y = ipece->grid_boundary_lower.y;
    probe_grid.y <= ipece->grid_boundary_higher.y;
    probe_grid.y++) {
        for (probe_grid.z = ipece->grid_boundary_lower.z;
        probe_grid.z <= ipece->grid_boundary_higher.z;
        probe_grid.z++) {
            for (probe_grid.x = ipece->grid_boundary_lower.x;
            probe_grid.x <= ipece->grid_boundary_higher.x;
            probe_grid.x++) {
                if (*reach_label(probe_grid, ipece) == 'o') continue;
                if (*reach_label(probe_grid, ipece) == '\0') {
                    *reach_label(probe_grid, ipece) = 'o';
                }
                else {
                    break;
                }
            }
        }
    }
    /* do the same to in y direction for grids on z-x plane */
    for (probe_grid.z = ipece->grid_boundary_lower.z;
    probe_grid.z <= ipece->grid_boundary_higher.z;
    probe_grid.z++) {
        for (probe_grid.x = ipece->grid_boundary_lower.x;
        probe_grid.x <= ipece->grid_boundary_higher.x;
        probe_grid.x++) {
            for (probe_grid.y = ipece->grid_boundary_lower.y;
            probe_grid.y <= ipece->grid_boundary_higher.y;
            probe_grid.y++) {
                if (*reach_label(probe_grid, ipece) == 'o') continue;
                if (*reach_label(probe_grid, ipece) == '\0') {
                    *reach_label(probe_grid, ipece) = 'o';
                }
                else {
                    break;
                }
            }
        }
    }
    /* Now starting from [1,1,1] corner and do the same as earlier */
    for (probe_grid.x = ipece->grid_boundary_higher.x;
    probe_grid.x >= ipece->grid_boundary_lower.x;
    probe_grid.x--) {
        for (probe_grid.y = ipece->grid_boundary_higher.y;
        probe_grid.y >= ipece->grid_boundary_lower.y;
        probe_grid.y--) {
            for (probe_grid.z = ipece->grid_boundary_higher.z;
            probe_grid.z >= ipece->grid_boundary_lower.z;
            probe_grid.z--) {
                if (*reach_label(probe_grid, ipece) == 'o') continue;
                if (*reach_label(probe_grid, ipece) == '\0') {
                    *reach_label(probe_grid, ipece) = 'o';
                }
                else {
                    break;
                }
            }
        }
    }
    for (probe_grid.z = ipece->grid_boundary_higher.z;
    probe_grid.z >= ipece->grid_boundary_lower.z;
    probe_grid.z--) {
        for (probe_grid.x = ipece->grid_boundary_higher.x;
        probe_grid.x >= ipece->grid_boundary_lower.x;
        probe_grid.x--) {
            for (probe_grid.y = ipece->grid_boundary_higher.y;
            probe_grid.y >= ipece->grid_boundary_lower.y;
            probe_grid.y--) {
                if (*reach_label(probe_grid, ipece) == 'o') continue;
                if (*reach_label(probe_grid, ipece) == '\0') {
                    *reach_label(probe_grid, ipece) = 'o';
                }
                else {
                    break;
                }
            }
        }
    }
    for (probe_grid.y = ipece->grid_boundary_higher.y;
    probe_grid.y >= ipece->grid_boundary_lower.y;
    probe_grid.y--) {
        for (probe_grid.z = ipece->grid_boundary_higher.z;
        probe_grid.z >= ipece->grid_boundary_lower.z;
        probe_grid.z--) {
            for (probe_grid.x = ipece->grid_boundary_higher.x;
            probe_grid.x >= ipece->grid_boundary_lower.x;
            probe_grid.x--) {
                if (*reach_label(probe_grid, ipece) == 'o') continue;
                if (*reach_label(probe_grid, ipece) == '\0') {
                    *reach_label(probe_grid, ipece) = 'o';
                }
                else {
                    break;
                }
            }
        }
    }
    
    /* Collect all the grids that are not labeled yet */
    ipece->n_probe = 0;
    ipece->probes = NULL;
    for (probe_grid.x = ipece->grid_boundary_lower.x;
    probe_grid.x <= ipece->grid_boundary_higher.x;
    probe_grid.x++) {
        for (probe_grid.y = ipece->grid_boundary_lower.y;
        probe_grid.y <= ipece->grid_boundary_higher.y;
        probe_grid.y++) {
            for (probe_grid.z = ipece->grid_boundary_lower.z;
            probe_grid.z <= ipece->grid_boundary_higher.z;
            probe_grid.z++) {
                if (*reach_label(probe_grid, ipece) == '\0') {
                    ipece->n_probe++;
                    ipece->probes = realloc(ipece->probes, ipece->n_probe*sizeof(INT_VECT));
                    ipece->probes[ipece->n_probe-1] = probe_grid;
                }
            }
        }
    }
    
    /* Look for the grids in the list adjacent to an exposed grid, and label it
    as 'o' */
    label_updated = 1;  /* a flag indicates that new exposed grid are found in one cycle */
    while (label_updated) {
        label_updated = 0;
        for (i_probe = ipece->n_probe-1; i_probe>=0; i_probe--) {
            INT_VECT adj_grid;
            
            /* going through nearby grids in three directions */
            for (adj_grid.x = ipece->probes[i_probe].x-1;
            adj_grid.x <= ipece->probes[i_probe].x+1;
            adj_grid.x++) {
                for (adj_grid.y = ipece->probes[i_probe].y-1;
                adj_grid.y <= ipece->probes[i_probe].y+1;
                adj_grid.y++) {
                    for (adj_grid.z = ipece->probes[i_probe].z-1;
                    adj_grid.z <= ipece->probes[i_probe].z+1;
                    adj_grid.z++) {
                        /* if the adjacent grid is 'o', label the grid in the list as 'o' */
                        if (*reach_label(adj_grid, ipece) == 'o') {
                            *reach_label(ipece->probes[i_probe], ipece) = 'o';
                            break;
                        }
                    }
                    if (*reach_label(ipece->probes[i_probe], ipece) == 'o') break;
                }
                if (*reach_label(ipece->probes[i_probe], ipece) == 'o') break;
            }
            
            /* if the grid is labeled as 'o', remove it from the list */
            if (*reach_label(ipece->probes[i_probe], ipece) == 'o') {
                
                memmove(&ipece->probes[i_probe], &ipece->probes[i_probe+1],
                (ipece->n_probe-i_probe-1)*sizeof(INT_VECT));
                
                ipece->n_probe--;
                label_updated = 1;
            }
        }
        ipece->probes = realloc(ipece->probes, (ipece->n_probe)*sizeof(INT_VECT));
    }
    
    /* All the remaining grids in the list now are not connected to outside of
    the protein, so labeled as 'c' for cavity */
    for (i_probe = 0; i_probe < ipece->n_probe; i_probe++) {
        *reach_label(ipece->probes[i_probe], ipece) = 'c';
    }
    
    ipece->n_probe = 0;
    free(ipece->probes);

    return 0;
}

/* Initialize a grid box using the current position of the protein */
int create_grid_box(PROT prot, IPECE *ipece)
{
    int i_grid, j_grid, k_grid;
    int i_res, i_conf, i_atom;
    VECTOR vec1;
    ATOM *atom_p;
    int initialized = 0;
    
    /* Initialize the variables for lower and higher boundary of the box */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                ipece->coor_boundary_lower = atom_p->xyz;
                ipece->coor_boundary_higher = atom_p->xyz;

                initialized = 1;
                break;
            }
            if (initialized) break;
        }
        if (initialized) break;
    }
    
    /* Find the size of the box to contain the protein */
    /* Define lower and higher boundary of the box */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!atom_p->on) continue;
                if (atom_p->xyz.x < ipece->coor_boundary_lower.x) ipece->coor_boundary_lower.x = atom_p->xyz.x;
                if (atom_p->xyz.y < ipece->coor_boundary_lower.y) ipece->coor_boundary_lower.y = atom_p->xyz.y;
                if (atom_p->xyz.z < ipece->coor_boundary_lower.z) ipece->coor_boundary_lower.z = atom_p->xyz.z;
                
                if (atom_p->xyz.x > ipece->coor_boundary_higher.x) ipece->coor_boundary_higher.x = atom_p->xyz.x;
                if (atom_p->xyz.y > ipece->coor_boundary_higher.y) ipece->coor_boundary_higher.y = atom_p->xyz.y;
                if (atom_p->xyz.z > ipece->coor_boundary_higher.z) ipece->coor_boundary_higher.z = atom_p->xyz.z;
            }
        }
    }
    
    /* Expand the boundary of the box a little */
    ipece->coor_boundary_lower.x  -= ipece->boundary_extention;
    ipece->coor_boundary_lower.y  -= ipece->boundary_extention;
    ipece->coor_boundary_lower.z  -= ipece->boundary_extention;
    
    ipece->coor_boundary_higher.x += ipece->boundary_extention;
    ipece->coor_boundary_higher.y += ipece->boundary_extention;
    ipece->coor_boundary_higher.z += ipece->boundary_extention;
    
    /* Define lower and higher boundary of the grid box in grid index*/
    vec1 = vector_rescale(ipece->coor_boundary_lower, 1./ipece->grid_space);
    ipece->grid_boundary_lower.x = (int)vec1.x - 1;
    ipece->grid_boundary_lower.y = (int)vec1.y - 1;
    ipece->grid_boundary_lower.z = (int)vec1.z - 1;
    
    vec1 = vector_rescale(ipece->coor_boundary_higher, 1./ipece->grid_space);
    ipece->grid_boundary_higher.x = (int)vec1.x + 1;
    ipece->grid_boundary_higher.y = (int)vec1.y + 1;
    ipece->grid_boundary_higher.z = (int)vec1.z + 1;
    
    /* number of grids in each direction */
    ipece->n_grid.x = ipece->grid_boundary_higher.x - ipece->grid_boundary_lower.x + 1;
    ipece->n_grid.y = ipece->grid_boundary_higher.y - ipece->grid_boundary_lower.y + 1;
    ipece->n_grid.z = ipece->grid_boundary_higher.z - ipece->grid_boundary_lower.z + 1;
    
    /* initialize the label array */
    ipece->label = malloc(ipece->n_grid.x * sizeof(char **));
    for (i_grid = 0; i_grid<ipece->n_grid.x; i_grid++) {
        ipece->label[i_grid] = malloc(ipece->n_grid.y * sizeof(char *));
        for (j_grid = 0; j_grid<ipece->n_grid.y; j_grid++) {
            ipece->label[i_grid][j_grid] = malloc(ipece->n_grid.z * sizeof(char));
            for (k_grid = 0; k_grid<ipece->n_grid.z; k_grid++) {
                ipece->label[i_grid][j_grid][k_grid] = '\0';
            }
        }
    }
    
    return 0;
}

INT_VECT coor2grid(VECTOR r, IPECE *ipece)
{
    INT_VECT grid;
    grid.x = (int)(r.x / ipece->grid_space + 0.5);
    grid.y = (int)(r.y / ipece->grid_space + 0.5);
    grid.z = (int)(r.z / ipece->grid_space + 0.5);
    return grid;
}

VECTOR grid2coor(INT_VECT grid, IPECE *ipece)
{
    VECTOR r;
    r.x = (float)grid.x * ipece->grid_space + 0.5;
    r.y = (float)grid.y * ipece->grid_space + 0.5;
    r.z = (float)grid.z * ipece->grid_space + 0.5;
    return r;
}

char *reach_label(INT_VECT grid, IPECE *ipece)
{
    INT_VECT index = grid;
    index.x -= ipece->grid_boundary_lower.x;
    index.y -= ipece->grid_boundary_lower.y;
    index.z -= ipece->grid_boundary_lower.z;
    return &(ipece->label[index.x][index.y][index.z]);
}

int free_probe(IPECE *ipece)
{
    int i_grid, j_grid;
    for (i_grid = 0; i_grid<ipece->n_grid.x; i_grid++) {
        for (j_grid = 0; j_grid<ipece->n_grid.y; j_grid++) {
            free(ipece->label[i_grid][j_grid]);
        }
        free(ipece->label[i_grid]);
    }
    free(ipece->label);

    return 0;
}

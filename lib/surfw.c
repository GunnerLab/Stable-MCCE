#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"
#define  ATOM_RAD 2.0
#define  TRUE     1
#define  FALSE    0

typedef struct {
    INT_VECT gsize;
    //INT_VECT gindex;
    VECTOR xyz_min, xyz_max;
    float  PROBE_RAD;
    float  grid_interval;
    CONF   ***grid;
} GRID;  /* using a structure to pass grid space data instead of using global variables - Yifan */

int  insAtomToConf(CONF *conf, int ins);
void calc_gsize();
void alloc_3d_array (int x_size, int y_size, int z_size, GRID* grid_p);
void free_3d_array(int x_size, int y_size, int z_size, GRID* grid_p);
int  fill_grid(PROT prot, GRID* grid_p);
int  mkacc(ATOM *atom, GRID* grid_p);
float get_sas_res(ATOM atom, RES res);
void get_pdb_size(PROT prot, GRID* grid_p);
void extend_grid(GRID* grid_p);
void set_vdw_rad(PROT prot, float probe_rad); 
void reset_atom_rad(PROT prot, float probe_rad);
void trim_conf(PROT prot);
void copy_sas(PROT src, PROT target);
int surfw_res(PROT prot, int ir, float probe_rad);

/* preset level 1: i = j = 2, 122 uniformly distributed points on a sphere */
#define   num_pts 122.
float point_preset[][3] = {
    { 0.580799102783,      0.760826885700,      0.289507985115},
    { 0.525731086731,      0.850650787354,      0.000000000000},
    { 0.000000000000,      0.525731086731,     -0.850650787354},
    {-0.525731086731,      0.850650787354,      0.000000000000},
    { 0.000000000000,      0.525731086731,      0.850650787354},
    { 0.850650787354,      0.000000000000,      0.525731086731},
    {-0.850650787354,      0.000000000000,      0.525731086731},
    { 0.000000000000,     -0.525731086731,      0.850650787354},
    {-0.850650787354,      0.000000000000,     -0.525731086731},
    {-0.525731086731,     -0.850650787354,      0.000000000000},
    { 0.850650787354,      0.000000000000,     -0.525731086731},
    { 0.525731086731,     -0.850650787354,      0.000000000000},
    { 0.000000000000,     -0.525731086731,     -0.850650787354},
    { 0.577350258827,      0.577350258827,      0.577350258827},
    { 0.356822103262,      0.000000000000,     -0.934172332287},
    {-0.577350258827,      0.577350258827,     -0.577350258827},
    {-0.356822103262,      0.000000000000,      0.934172332287},
    { 0.934172332287,      0.356822103262,      0.000000000000},
    { 0.000000000000,      0.934172332287,     -0.356822103262},
    {-0.577350258827,      0.577350258827,      0.577350258827},
    { 0.577350258827,     -0.577350258827,      0.577350258827},
    {-0.934172332287,     -0.356822103262,      0.000000000000},
    { 0.577350258827,      0.577350258827,     -0.577350258827},
    { 0.000000000000,      0.934172332287,      0.356822103262},
    {-0.577350258827,     -0.577350258827,      0.577350258827},
    { 0.934172332287,     -0.356822103262,      0.000000000000},
    { 0.356822103262,      0.000000000000,      0.934172332287},
    {-0.356822103262,      0.000000000000,     -0.934172332287},
    { 0.000000000000,     -0.934172332287,     -0.356822103262},
    { 0.577350258827,     -0.577350258827,     -0.577350258827},
    { 0.000000000000,     -0.934172332287,      0.356822103262},
    {-0.934172332287,      0.356822103262,      0.000000000000},
    {-0.577350258827,     -0.577350258827,     -0.577350258827},
    { 1.000000000000,      0.000000000000,      0.000000000000},
    {-0.500000000000,     -0.309017002583,     -0.809017002583},
    {-1.000000000000,      0.000000000000,      0.000000000000},
    { 0.500000000000,     -0.309017002583,      0.809017002583},
    { 0.500000000000,      0.309017002583,     -0.809017002583},
    {-0.500000000000,      0.309017002583,     -0.809017002583},
    {-0.500000000000,     -0.309017002583,      0.809017002583},
    { 0.809017002583,     -0.500000000000,     -0.309017002583},
    {-0.309017002583,     -0.809017002583,      0.500000000000},
    {-0.309017002583,      0.809017002583,     -0.500000000000},
    {-0.809017002583,      0.500000000000,      0.309017002583},
    { 0.000000000000,     -1.000000000000,      0.000000000000},
    { 0.809017002583,      0.500000000000,     -0.309017002583},
    { 0.000000000000,      1.000000000000,      0.000000000000},
    {-0.809017002583,     -0.500000000000,      0.309017002583},
    { 0.809017002583,      0.500000000000,      0.309017002583},
    { 0.309017002583,      0.809017002583,     -0.500000000000},
    {-0.309017002583,     -0.809017002583,     -0.500000000000},
    { 0.000000000000,      0.000000000000,      1.000000000000},
    { 0.309017002583,      0.809017002583,      0.500000000000},
    { 0.000000000000,      0.000000000000,     -1.000000000000},
    { 0.309017002583,     -0.809017002583,      0.500000000000},
    {-0.309017002583,      0.809017002583,      0.500000000000},
    { 0.309017002583,     -0.809017002583,     -0.500000000000},
    {-0.809017002583,      0.500000000000,     -0.309017002583},
    {-0.809017002583,     -0.500000000000,     -0.309017002583},
    { 0.809017002583,     -0.500000000000,      0.309017002583},
    { 0.500000000000,      0.309017002583,      0.809017002583},
    { 0.500000000000,     -0.309017002583,     -0.809017002583},
    {-0.500000000000,      0.309017002583,      0.809017002583},
    { 0.178925782442,      0.291291087866,     -0.939752638340},
    {-0.580799102783,      0.760826885700,     -0.289507985115},
    {-0.178925782442,      0.291291087866,      0.939752638340},
    { 0.759724855423,      0.650244653225,     -0.000000000141},
    {-0.291291087866,      0.939752638340,     -0.178925782442},
    {-0.289507985115,      0.580799102783,      0.760826885700},
    { 0.760826885700,     -0.289507985115,      0.580799102783},
    {-0.939752638340,     -0.178925782442,      0.291291087866},
    { 0.580799102783,      0.760826885700,     -0.289507985115},
    {-0.000000000141,      0.759724855423,      0.650244653225},
    {-0.289507985115,     -0.580799102783,      0.760826885700},
    { 0.939752638340,     -0.178925782442,      0.291291087866},
    { 0.289507985115,      0.580799102783,      0.760826885700},
    {-0.178925782442,     -0.291291087866,      0.939752638340},
    { 0.178925782442,      0.291291087866,      0.939752638340},
    {-0.291291087866,      0.939752638340,      0.178925782442},
    {-0.650244653225,     -0.000000000141,     -0.759724855423},
    {-0.760826885700,     -0.289507985115,      0.580799102783},
    {-0.580799102783,      0.760826885700,      0.289507985115},
    {-0.760826885700,      0.289507985115,     -0.580799102783},
    {-0.291291087866,     -0.939752638340,     -0.178925782442},
    { 0.650244653225,     -0.000000000141,      0.759724855423},
    { 0.760826885700,     -0.289507985115,     -0.580799102783},
    { 0.291291087866,      0.939752638340,     -0.178925782442},
    { 0.760826885700,      0.289507985115,      0.580799102783},
    { 0.939752638340,     -0.178925782442,     -0.291291087866},
    {-0.178925782442,      0.291291087866,     -0.939752638340},
    { 0.291291087866,     -0.939752638340,      0.178925782442},
    {-0.939752638340,      0.178925782442,      0.291291087866},
    { 0.000000000141,     -0.759724855423,      0.650244653225},
    {-0.760826885700,      0.289507985115,      0.580799102783},
    {-0.580799102783,     -0.760826885700,     -0.289507985115},
    { 0.289507985115,     -0.580799102783,      0.760826885700},
    {-0.650244653225,      0.000000000141,      0.759724855423},
    { 0.291291087866,     -0.939752638340,     -0.178925782442},
    { 0.939752638340,      0.178925782442,      0.291291087866},
    { 0.178925782442,     -0.291291087866,      0.939752638340},
    { 0.650244653225,      0.000000000141,     -0.759724855423},
    {-0.759724855423,     -0.650244653225,     -0.000000000141},
    { 0.580799102783,     -0.760826885700,      0.289507985115},
    {-0.289507985115,     -0.580799102783,     -0.760826885700},
    { 0.289507985115,      0.580799102783,     -0.760826885700},
    {-0.760826885700,     -0.289507985115,     -0.580799102783},
    {-0.759724855423,      0.650244653225,      0.000000000141},
    { 0.000000000141,      0.759724855423,     -0.650244653225},
    { 0.289507985115,     -0.580799102783,     -0.760826885700},
    {-0.939752638340,     -0.178925782442,     -0.291291087866},
    {-0.289507985115,      0.580799102783,     -0.760826885700},
    { 0.939752638340,      0.178925782442,     -0.291291087866},
    {-0.000000000141,     -0.759724855423,     -0.650244653225},
    { 0.760826885700,      0.289507985115,     -0.580799102783},
    {-0.939752638340,      0.178925782442,     -0.291291087866},
    { 0.178925782442,     -0.291291087866,     -0.939752638340},
    { 0.580799102783,     -0.760826885700,     -0.289507985115},
    {-0.580799102783,     -0.760826885700,      0.289507985115},
    {-0.291291087866,     -0.939752638340,      0.178925782442},
    { 0.291291087866,      0.939752638340,      0.178925782442},
    {-0.178925782442,     -0.291291087866,     -0.939752638340},
    { 0.759724855423,     -0.650244653225,      0.000000000141}
};
float area_coeff = 4. * 3.1415926 / num_pts;

int surfw(PROT protein, float probe_rad)
{
  int   i, j, k;
  int i_res, i_conf;
  float res_sas;                /* sas per residue (disregard all other residues) */
  float atom_sas;               /* sas per atom */
  PROT  prot;
  GRID  grid;
  
  prot = new_prot();
  cpy_prot(&prot, &protein);
  trim_conf(prot);
  for (i_res=0; i_res<prot.n_res; i_res++) {
      for (i_conf=0; i_conf<prot.res[i_res].n_conf; i_conf++) {
          prot.res[i_res].conf[i_conf].on  = 1;
      }
  }
  
  grid.PROBE_RAD  = probe_rad;       /* get parameters */
  grid.grid_interval = grid.PROBE_RAD + ATOM_RAD;

  get_pdb_size(prot, &grid);

  extend_grid(&grid);

  calc_gsize(&grid);

  //grid.gindex.x = grid.gsize.x - 1;
  //grid.gindex.y = grid.gsize.y - 1;
  //grid.gindex.z = grid.gsize.z - 1;

  set_vdw_rad(prot, grid.PROBE_RAD);

  alloc_3d_array(grid.gsize.x, grid.gsize.y, grid.gsize.z, &grid);
  /* end drawing the grid */

  fill_grid(prot, &grid);

  /* make the surface for each atom */
  for (i = 0 ; i < prot.n_res; i++) {
      prot.res[i].sas = 0;
      for (j = 0; j < prot.res[i].n_conf; j++) {
          for (k = 0; k < prot.res[i].conf[j].n_atom; k++) {
              if (!(prot.res[i].conf[j].atom[k].on)) {
                  /* printf("%s, %d, %s\n", prot.res[i].resName, prot.res[i].resSeq, prot.res[i].conf[j].atom[k].name); */
                  continue;
              }
              if (prot.res[i].conf[j].atom[k].vdw_rad < 1e-4) continue; /* ignore 0 radius atoms -Yifan*/
              mkacc(&(prot.res[i].conf[j].atom[k]), &grid);
              prot.res[i].sas += prot.res[i].conf[j].atom[k].sas;
          }
      }
  }

  /* compute the % of solvent-accessible surface of each residue */
  for (i = 0 ; i < prot.n_res; i++) {
      res_sas = 0;
      for (j = 0; j < prot.res[i].n_conf; j++) {
          for (k = 0; k < prot.res[i].conf[j].n_atom; k++) {
              if (!(prot.res[i].conf[j].atom[k].on)) continue;
              if (prot.res[i].conf[j].atom[k].vdw_rad < 1e-4) continue; /* ignore 0 radius atoms -Yifan*/
              atom_sas = get_sas_res(prot.res[i].conf[j].atom[k], prot.res[i]);
              res_sas += atom_sas;
          }
      }
      prot.res[i].sas /= res_sas;
  }

  reset_atom_rad(prot, grid.PROBE_RAD);
  
  /* print the surface area per atom */
  /*
  for (i = 0; i < prot.n_res; i++) {
      for (j = 0; j < prot.res[i].n_conf; j++) {
          for (k = 0; k < prot.res[i].conf[j].n_atom; k++) {
              if (prot.res[i].conf[j].atom[k].on) {
                  printf("%4s  %4s %c %3d%8.3f\n",
                      prot.res[i].conf[j].atom[k].name, prot.res[i].conf[j].atom[k].resName,
                      prot.res[i].conf[j].atom[k].chainID, prot.res[i].conf[j].atom[k].resSeq,
                      prot.res[i].conf[j].atom[k].sas);
              }
          }
      }
  }
  */

  /* print the surface area per residue */
  /*
  for (i = 0; i < prot.n_res; i++) {
    printf("%4s %c %4d %8.3f\n", prot.res[i].resName, prot.res[i].chainID, prot.res[i].resSeq,
            prot.res[i].sas);
  }
  */

  free_3d_array(grid.gsize.x, grid.gsize.y, grid.gsize.z, &grid);
  copy_sas(prot, protein);
  del_prot(&prot);
  return 0;

}

/* remove all conformers other than conform 0 and 1 */
void trim_conf(PROT prot) {
  int i, j;

  for (i = 0; i < prot.n_res; i++) {
    if (prot.res[i].n_conf <= 2) continue;
    for (j = prot.res[i].n_conf; j > 2; j--) {
      del_conf( &(prot.res[i]), j-1 );
    }
  }

  return;
}

/* calculate ASA of the given atom using the given res as the boundary condition - Yifan*/
/* used to estimate ASA of a residue/conformer when it is fully exposed */
/* for multi-conformation, conf[].on flag is used to control which conformer is used for the boundary condition */
float get_sas_res(ATOM atom, RES res)
{
    int j, k, m;
    int count;
    int status;
    VECTOR point;
    float distance2;
    float surface_area;
    
    count = num_pts;
    
    for (m = 0; m < num_pts; m++) {
        status = TRUE;
        point.x = point_preset[m][0] * atom.vdw_rad + atom.xyz.x;
        point.y = point_preset[m][1] * atom.vdw_rad + atom.xyz.y;
        point.z = point_preset[m][2] * atom.vdw_rad + atom.xyz.z;
        
        /* go through all atoms in this residue, and excluding point from total count if overlaps with any atom -cmt. by Yifan */
        for (j = 0; j < res.n_conf; j++) {
            if (status == FALSE) break;
            if (j >= 1) { /* only using backbone (j=0) and the conformer that is turned on -added Yifan */
                if (!res.conf[j].on) continue;
            }
            for (k = 0; k < res.conf[j].n_atom; k++) {
                if (!res.conf[j].atom[k].on) continue;
                if (res.conf[j].atom[k].vdw_rad < 1e-4) continue; /* ignore 0 radius atoms -Yifan*/
                if (status != TRUE) continue;
                if (!memcmp(&atom, &(res.conf[j].atom[k]), sizeof(ATOM))) continue;
                
                distance2 = ddvv(point, res.conf[j].atom[k].xyz);
                if (distance2 < res.conf[j].atom[k].vdw_rad * res.conf[j].atom[k].vdw_rad) {
                    count--;
                    status = FALSE;
                }
            }
        }
    }
    
    surface_area = area_coeff * atom.vdw_rad * atom.vdw_rad * (float)count;
    return surface_area;
}

/* copy the SAS value back to the original pdb */
void copy_sas(PROT src, PROT target)
{
  int i, j, k;

  for (i=0; i<src.n_res; i++) {
    target.res[i].sas = src.res[i].sas;
    for (j=0; j<src.res[i].n_conf; j++) {
      for (k=0; k<src.res[i].conf[j].n_atom; k++) {
        if (!(src.res[i].conf[j].atom[k].on)) continue;
        target.res[i].conf[j].atom[k].sas = src.res[i].conf[j].atom[k].sas;
      }
    }
  }

  return;
}

/* calculate ASA of all conformers, assuming the rest of the protein occupied in the first conformer -Yifan */
int sas_native(PROT prot)
{
    int ir, ic;
    
    for (ir=0; ir<prot.n_res; ir++) {
        prot.res[ir].conf[0].on = 1;
        if (prot.res[ir].n_conf > 1) prot.res[ir].conf[1].on = 1;
        for (ic=2; ic<prot.res[ir].n_conf; ic++) prot.res[ir].conf[ic].on = 0;
    }
    
    /* calculate ASA for all residues */
    for (ir=0; ir<prot.n_res; ir++) {
        surfw_res(prot, ir, env.radius_probe);
    }
    
    return 0;
}

/* calculate conf[].sas and atom[].sas (absolute value) for all conformers in k_res -Yifan */
int surfw_res(PROT prot, int i_res, float probe_rad)
{
    int   i_conf, j_conf, i_atom;
    GRID  grid;
    
    grid.PROBE_RAD  = probe_rad;       /* get parameters */
    grid.grid_interval = grid.PROBE_RAD + ATOM_RAD;
    
    set_vdw_rad(prot, grid.PROBE_RAD);
    
    /* make the surface for each atom */
    for (i_conf = 1; i_conf < prot.res[i_res].n_conf; i_conf++) {
        float solution_sas;
        /* turn on i_conf, and turn off the rest */
        for (j_conf = 1; j_conf < prot.res[i_res].n_conf; j_conf++) {
            if (i_conf==j_conf) prot.res[i_res].conf[j_conf].on = 1;
            else prot.res[i_res].conf[j_conf].on = 0;
        }
        
        get_pdb_size(prot, &grid);
        extend_grid(&grid);
        calc_gsize(&grid);
        //grid.gindex.x = grid.gsize.x - 1;
        //grid.gindex.y = grid.gsize.y - 1;
        //grid.gindex.z = grid.gsize.z - 1;
        alloc_3d_array(grid.gsize.x, grid.gsize.y, grid.gsize.z, &grid);
        fill_grid(prot, &grid);
        
        prot.res[i_res].conf[i_conf].sas = 0.;
        for (i_atom = 0; i_atom < prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            if (!(prot.res[i_res].conf[i_conf].atom[i_atom].on)) continue;
            if (prot.res[i_res].conf[i_conf].atom[i_atom].vdw_rad < 1e-4) continue; /* ignore 0 radius atoms -Yifan*/
            mkacc(&(prot.res[i_res].conf[i_conf].atom[i_atom]), &grid);
            prot.res[i_res].conf[i_conf].sas += prot.res[i_res].conf[i_conf].atom[i_atom].sas;
        }
        free_3d_array(grid.gsize.x, grid.gsize.y, grid.gsize.z, &grid);
        
        /* calculate the faction of exposure - added Yifan */
        solution_sas = 0.;
        for (i_atom = 0; i_atom < prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
            if (!(prot.res[i_res].conf[i_conf].atom[i_atom].on)) continue;
            if (prot.res[i_res].conf[i_conf].atom[i_atom].vdw_rad < 1e-4) continue; /* ignore 0 radius atoms -Yifan*/
            solution_sas += get_sas_res(prot.res[i_res].conf[i_conf].atom[i_atom], prot.res[i_res]);
        }
        prot.res[i_res].conf[i_conf].sas_fraction = prot.res[i_res].conf[i_conf].sas/solution_sas;
    }
    
    /* reset control flags in this residue */
    if (prot.res[i_res].n_conf > 1) prot.res[i_res].conf[1].on = 1;
    for (i_conf=2; i_conf<prot.res[i_res].n_conf; i_conf++) prot.res[i_res].conf[i_conf].on = 0;
    
    reset_atom_rad(prot, grid.PROBE_RAD);
    
    return 0;
    
}

/* calculate number of grids in the x, y, z axis */
void calc_gsize(GRID* grid_p)
{
  grid_p->gsize.x = (int)((grid_p->xyz_max.x - grid_p->xyz_min.x) / grid_p->grid_interval + 1);
  grid_p->gsize.y = (int)((grid_p->xyz_max.y - grid_p->xyz_min.y) / grid_p->grid_interval + 1);
  grid_p->gsize.z = (int)((grid_p->xyz_max.z - grid_p->xyz_min.z) / grid_p->grid_interval + 1);

  return;
}

int mkacc(ATOM *atom, GRID* grid_p)
{
    VECTOR point;
    int i, j, k, l, m;
    float d_square;
    int count;
    int ix, iy, iz;
    int ax, ay, az;         /* the grid index of this atom */
    int status;
    
    count = num_pts;
    
    ax = (int) ((atom -> xyz.x - grid_p->xyz_min.x) / grid_p->grid_interval);
    ay = (int) ((atom -> xyz.y - grid_p->xyz_min.y) / grid_p->grid_interval);
    az = (int) ((atom -> xyz.z - grid_p->xyz_min.z) / grid_p->grid_interval);
    
    /* loop over each point */
    
    for (m = 0; m < num_pts; m++) {
        status = TRUE;
        
        point.x = point_preset[m][0] * atom ->vdw_rad + atom -> xyz.x;
        point.y = point_preset[m][1] * atom ->vdw_rad + atom -> xyz.y;
        point.z = point_preset[m][2] * atom ->vdw_rad + atom -> xyz.z;
        
        /* find out which grid this point belongs to */
        ix = (int) ((point.x - grid_p->xyz_min.x) / grid_p->grid_interval);
        iy = (int) ((point.y - grid_p->xyz_min.y) / grid_p->grid_interval);
        iz = (int) ((point.z - grid_p->xyz_min.z) / grid_p->grid_interval);
        
        if (grid_p->grid[ix][iy][iz].n_atom >= 1) {
            for (l = 0; l < grid_p->grid[ix][iy][iz].n_atom; l++) {
                if (grid_p->grid[ix][iy][iz].atom[l].on &&
                    status == TRUE &&
                    (ix != ax || iy != ay || iz != az ||
                        memcmp(atom, &(grid_p->grid[ix][iy][iz].atom[l]), sizeof(ATOM)))) {
                d_square = ddvv(point, grid_p->grid[ix][iy][iz].atom[l].xyz);
                if (d_square < grid_p->grid[ix][iy][iz].atom[l].vdw_rad * grid_p->grid[ix][iy][iz].atom[l].vdw_rad) {
                    status = FALSE;
                    count --;
                }
                        }
            }
        }
        
        for (i = ix - 1; i <= ix + 1; i++) {
            if (i<0) continue;
            if (i>=grid_p->gsize.x) continue;
            if (status == FALSE) break;
            for (j = iy - 1; j <= iy + 1; j++) {
                if (j<0) continue;
                if (j>=grid_p->gsize.y) continue;
                if (status == FALSE) break;
                for (k = iz - 1; k <= iz + 1; k++) {
                    if (k<0) continue;
                    if (k>=grid_p->gsize.z) continue;
                    if (status == FALSE) break;
                    for (l = 0; l < grid_p->grid[i][j][k].n_atom; l++) {
                        if (grid_p->grid[i][j][k].atom[l].on && status == TRUE) {
                            if (i != ax || j != ay || k != az || memcmp(atom, &(grid_p->grid[i][j][k].atom[l]), sizeof(ATOM))) {
                                d_square = ddvv(point, grid_p->grid[i][j][k].atom[l].xyz);
                                if (d_square < grid_p->grid[i][j][k].atom[l].vdw_rad * grid_p->grid[i][j][k].atom[l].vdw_rad) {
                                    status = FALSE;
                                    count --;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    atom->sas = area_coeff * atom->vdw_rad * atom->vdw_rad * (float) count;
    //atom->sas = (float) count/num_pts;
    //printf("%s %s%4d %4d %8.3f %8.3f\n",atom->name,atom->resName,atom->resSeq,count,atom->vdw_rad,atom->sas);
    return 0;
}

/* find minimum and maximum x,y,z of the whole prot */
void get_pdb_size(PROT prot, GRID* grid_p)
{
    int i, j, k;
    
    grid_p->xyz_min.x =  1.0E20;
    grid_p->xyz_max.x = -1.0E20;
    grid_p->xyz_min.y =  1.0E20;
    grid_p->xyz_max.y = -1.0E20;
    grid_p->xyz_min.z =  1.0E20;
    grid_p->xyz_max.z = -1.0E20;
    
    for (i = 0; i < prot.n_res; i++) {
        for (j = 0; j < prot.res[i].n_conf; j++) {
            if (!prot.res[i].conf[j].on) continue;
            for (k = 0; k < prot.res[i].conf[j].n_atom; k++) {
                if (!(prot.res[i].conf[j].atom[k].on)) continue;
                
                if      (grid_p->xyz_min.x > prot.res[i].conf[j].atom[k].xyz.x)
                    grid_p->xyz_min.x = prot.res[i].conf[j].atom[k].xyz.x;
                if       (grid_p->xyz_max.x < prot.res[i].conf[j].atom[k].xyz.x)
                    grid_p->xyz_max.x = prot.res[i].conf[j].atom[k].xyz.x;
                
                if      (grid_p->xyz_min.y > prot.res[i].conf[j].atom[k].xyz.y)
                    grid_p->xyz_min.y = prot.res[i].conf[j].atom[k].xyz.y;
                if      (grid_p->xyz_max.y < prot.res[i].conf[j].atom[k].xyz.y)
                    grid_p->xyz_max.y = prot.res[i].conf[j].atom[k].xyz.y;
                
                if      (grid_p->xyz_min.z > prot.res[i].conf[j].atom[k].xyz.z)
                    grid_p->xyz_min.z = prot.res[i].conf[j].atom[k].xyz.z;
                if      (grid_p->xyz_max.z < prot.res[i].conf[j].atom[k].xyz.z)
                    grid_p->xyz_max.z = prot.res[i].conf[j].atom[k].xyz.z;
            }
        }
    }
    
    /*
    printf("min: %8.3f %8.3f %8.3f\n", xyz_min.x, xyz_min.y, xyz_min.z);
    printf("max: %8.3f %8.3f %8.3f\n", xyz_max.x, xyz_max.y, xyz_max.z);
    */
    
    return;
}

/* malloc a 3-D array of size x_size * y_size * z_size, and set array[i][j][k] to 0 */
void alloc_3d_array (int x_size, int y_size, int z_size, GRID* grid_p)
{
    int i, j, k;
    
    grid_p->grid = (CONF ***)calloc(x_size, sizeof(CONF **)); /* cast "void *" to "CONF ***" */
    
    for (i = 0; i < x_size; i++) {
        grid_p->grid[i] = (CONF **) calloc( y_size, sizeof(CONF *) );
    }
    for (i = 0; i < x_size; i++) {
        for (j = 0; j < y_size; j++) {
            grid_p->grid[i][j] = (CONF *) calloc(z_size, sizeof(CONF));
        }
    }
    
    for (i = 0; i < x_size; i++) {
        for (j = 0; j < y_size; j++) {
            for (k = 0; k < z_size; k++) {
                /*
                grid[i][j][k].n_atom = 0;
                grid[i][j][k].atom = NULL;
                */
                memset(&grid_p->grid[i][j][k], 0, sizeof(CONF));    
            }
        }
    }
    return;
}

void free_3d_array(int x_size, int y_size, int z_size, GRID* grid_p)
{
    int i, j, k;
    
    for (i = 0; i < x_size; i++) {
        for (j = 0; j < y_size; j++) {
            for (k = 0; k < z_size; k++) {
                if (grid_p->grid[i][j][k].atom)
                    free(grid_p->grid[i][j][k].atom);
            }
            free(grid_p->grid[i][j]);
        }
        free(grid_p->grid[i]);
    }
    free(grid_p->grid);
    
    return;
}

int fill_grid(PROT prot, GRID* grid_p)
{
    int ix, iy, iz;                  /* x, y, z grid index for an atom */
    int i, j, k;
    //int iConf;
    int ins_pos;
    
    /* put each atom into approriate grid box */
    for (i = 0; i < prot.n_res; i++) {
        //iConf = 0;
        for (j = 0; j < prot.res[i].n_conf; j++) {
            if (!prot.res[i].conf[j].on) continue;
            for (k = 0; k < prot.res[i].conf[j].n_atom; k++) {
                if (!prot.res[i].conf[j].atom[k].on) continue;
                if (prot.res[i].conf[j].atom[k].vdw_rad < 1e-4) continue; /* ignore 0 radius atoms -Yifan*/
                ix = (int) ((prot.res[i].conf[j].atom[k].xyz.x - grid_p->xyz_min.x) / grid_p->grid_interval);
                iy = (int) ((prot.res[i].conf[j].atom[k].xyz.y - grid_p->xyz_min.y) / grid_p->grid_interval);
                iz = (int) ((prot.res[i].conf[j].atom[k].xyz.z - grid_p->xyz_min.z) / grid_p->grid_interval);
                
                ins_pos = insAtomToConf(&(grid_p->grid[ix][iy][iz]), grid_p->grid[ix][iy][iz].n_atom);
                //grid[ix][iy][iz].atom[ins_pos] = prot.res[i].conf[j].atom[k];
                memcpy( &(grid_p->grid[ix][iy][iz].atom[ins_pos]), &(prot.res[i].conf[j].atom[k]), sizeof(ATOM) );
            }
            //iConf++;
        }
    }
    
    return 0;
}

void extend_grid(GRID* grid_p)
{
    grid_p->xyz_min.x -= 2 * grid_p->grid_interval+1.0;
    grid_p->xyz_min.y -= 2 * grid_p->grid_interval+1.0;
    grid_p->xyz_min.z -= 2 * grid_p->grid_interval+1.0;
    grid_p->xyz_max.x += 2 * grid_p->grid_interval+1.0;
    grid_p->xyz_max.y += 2 * grid_p->grid_interval+1.0;
    grid_p->xyz_max.z += 2 * grid_p->grid_interval+1.0;
    
    return;
}

void set_vdw_rad(PROT prot, float probe_rad)
{
    
    int i, j, k;
    
    for (i = 0; i < prot.n_res; i++) {
        for (j = 0; j < prot.res[i].n_conf; j++) {
            for (k = 0; k < prot.res[i].conf[j].n_atom; k++) {
                //printf("atom=%s rad= %8.3f\n",prot.res[i].conf[j].atom[k].name,prot.res[i].conf[j].atom[k].vdw_rad);
                if (prot.res[i].conf[j].atom[k].vdw_rad < 1e-4) continue; /* ignore 0 radius atoms -Yifan*/
                prot.res[i].conf[j].atom[k].vdw_rad += probe_rad;
            }
        }
    }
    return;
}

void reset_atom_rad(PROT prot, float probe_rad)
{
    int i, j, k;
    
    for (i = 0; i < prot.n_res; i++) {
        for (j = 0; j < prot.res[i].n_conf; j++) {
            for (k = 0; k < prot.res[i].conf[j].n_atom; k++) {
                if (prot.res[i].conf[j].atom[k].vdw_rad < 1e-4) continue; /* ignore 0 radius atoms -Yifan*/
                prot.res[i].conf[j].atom[k].vdw_rad -= probe_rad;
            }
        }
    }
    
    return;
}

int insAtomToConf(CONF *conf, int ins)
{ if (ins > conf->n_atom) {
    printf("insAtomToConf(): off range insertion.\n");
    return USERERR;
    }
    
    /* resize the memory */
    if (!(conf->atom = (ATOM *) realloc(conf->atom, (conf->n_atom + 1) * sizeof(ATOM)))) {
        printf("insAtomToConf(): Fails resizing memory.\n");
        return USERERR;
    }
    conf->n_atom ++;
    
    /* move contents after position "ins" forward by 1 */
    memmove(conf->atom+ins+1, conf->atom+ins, (conf->n_atom - ins - 1) * sizeof(ATOM));
    
    /* reset the new slot */
    memset(conf->atom+ins, 0, sizeof(ATOM));
    
    return ins;
}

/* calculation ASA of terminal N and O atoms of ionizable residues - added Yifan*/
int sas_ionizable(PROT prot, float probe_rad)
{
    int i_res, i_conf, i_atom;
    float score;
    GRID grid;
    grid.PROBE_RAD  = probe_rad;       /* get parameters */
    grid.grid_interval = grid.PROBE_RAD + ATOM_RAD;
    
    delete_h(prot);
    set_vdw_rad(prot, grid.PROBE_RAD);
    
    /* initialize conf[].on */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        prot.res[i_res].conf[0].on = 1;
        if (prot.res[i_res].n_conf > 1) prot.res[i_res].conf[1].on = 1;
        for (i_conf=2; i_conf<prot.res[i_res].n_conf; i_conf++) prot.res[i_res].conf[i_conf].on = 0;
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            prot.res[i_res].conf[i_conf].sas = 0.;
            prot.res[i_res].conf[i_conf].sas_fraction = 0.;
        }
    }
    
    /* calculate ASA for all residues */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            int j_conf;
            float solution_sas;
            int use_ipece; /* a flag to decide whether ipece information is used to determine which atom counts as polar */
            if (i_conf == 1) { /* if the first conformer, check to see if this is an ionizable residue */
                for (i_atom = 0; i_atom < prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                    ATOM* atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                    if (!(atom_p->on)) continue;
                    if ( !param_get("IPECE_SC",prot.res[i_res].resName, atom_p->name, &score) ) break;
                }
                if (i_atom >= prot.res[i_res].conf[i_conf].n_atom) {
                    /* cannot find an atom matches parameter, move to the next residue */
                    use_ipece = 0;
                }
                else {
                    use_ipece = 1;
                }
            }
            /* ionizable residue */
            
            /* turn on i_conf, and turn off the rest */
            for (j_conf = 1; j_conf < prot.res[i_res].n_conf; j_conf++) {
                if (i_conf==j_conf) prot.res[i_res].conf[j_conf].on = 1;
                else prot.res[i_res].conf[j_conf].on = 0;
            }
            
            /* make the surface for each atom */
            get_pdb_size(prot, &grid);
            extend_grid(&grid);
            calc_gsize(&grid);
            alloc_3d_array(grid.gsize.x, grid.gsize.y, grid.gsize.z, &grid);
            fill_grid(prot, &grid);
            
            prot.res[i_res].conf[i_conf].sas = 0.;
            solution_sas = 0.;
            for (i_atom = 0; i_atom < prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                ATOM* atom_p = &prot.res[i_res].conf[i_conf].atom[i_atom];
                if (!(atom_p->on)) continue;
                if (use_ipece) {
                    if ( param_get("IPECE_SC",prot.res[i_res].resName, atom_p->name, &score) != 0) continue; /* skip atoms not in parameter file */
                }
                else {
                    if ( atom_p->name[1] != 'O' && atom_p->name[1] != 'N' ) continue; /* skip atoms not O or N */
                }
                mkacc(atom_p, &grid);
                prot.res[i_res].conf[i_conf].sas += prot.res[i_res].conf[i_conf].atom[i_atom].sas;
                solution_sas += get_sas_res(prot.res[i_res].conf[i_conf].atom[i_atom], prot.res[i_res]);
            }
            if (solution_sas > 1e-3) {
                prot.res[i_res].conf[i_conf].sas_fraction = prot.res[i_res].conf[i_conf].sas/solution_sas;
            }
            free_3d_array(grid.gsize.x, grid.gsize.y, grid.gsize.z, &grid);
        }
        
        /* reset control flags in this residue */
        if (prot.res[i_res].n_conf > 1) prot.res[i_res].conf[1].on = 1;
        for (i_conf=2; i_conf<prot.res[i_res].n_conf; i_conf++) prot.res[i_res].conf[i_conf].on = 0;
        
    }
    
    reset_atom_rad(prot, grid.PROBE_RAD);
    
    return 0;
}


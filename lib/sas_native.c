#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"
#define  ATOM_RAD 2.0
#define  TRUE     1
#define  FALSE    0

typedef struct {
  int x, y, z;
} INT_VECTOR;

static int  insAtomToConf(CONF *conf, int ins);
static void calc_gsize();
static void alloc_3d_array(int x_size, int y_size, int z_size);
static void free_3d_array(int x_size, int y_size, int z_size);
static int  fill_grid(PROT prot);
static int  mkacc(ATOM *atom);
static void get_pdb_size(PROT prot);
static void extend_grid();
static void set_vdw_rad(PROT prot, float probe_rad);
static void reset_atom_rad(PROT prot, float probe_rad);

static INT_VECTOR gsize;
static INT_VECTOR gindex;
static VECTOR xyz_min, xyz_max;
static float  PROBE_RAD;
static float  grid_interval;
static CONF   ***grid;
static int    num_pts;





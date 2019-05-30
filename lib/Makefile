CC      = gcc
#CC	= icc
CFLAGS  = -g -Wall `pkg-config --cflags --libs glib-2.0`
#CFLAGS	= -xT -openmp -I/opt/local/include -I/usr/local/include -I/opt/GDBM/include
#CFLAGS  = -O3 -g -Wall -I/opt/GDBM/include -fopenmp

LIB     = mcce.a
AR      = ar
ARFLAGS = rvs


SRC     =  all.c       ddvv.c          geom_3v_onto_3v.c  ins_res.c         strip.c\
app.c       del_conf.c      geom_apply.c       line_2v.c         vdotv.c\
avv.c       del_prot.c      geom_inverse.c     load_all_param.c  vector_normalize.c\
avvvv.c     del_res.c       geom_move.c        load_param.c      vector_vminusv.c\
cpy_conf.c  det3.c          geom_reset.c       mxm4.c            vector_vplusv.c\
cpy_prot.c  det4.c          geom_roll.c        new_prot.c        vector_vxv.c\
cpy_res.c   dll.c           get_files.c        param_get.c	param_exist.c\
db_close.c  dvv.c           iatom.c            param_sav.c\
db_open.c   free_strings.c  ins_conf.c         plane_3v.c pdbline2atom.c\
premcce.c   init.c          load_pdb.c write_pdb.c   rotamers.c assign_rad.c get_connect12.c\
surfw.c   vdw.c vdw_conf.c  shuffle_n.c  cmp_conf.c  sort_conf.c	sort_res.c	 id_conf.c\
energies.c  assign_crg.c    coulomb.c   coulomb_conf.c\
get_vdw0.c  get_vdw1.c      relax_water.c      relax_h.c monte.c monte2.c  ran2.c\
relaxation.c collect_connect.c  torsion.c	vdw_fast.c hbond_extra.c swap.c quick_e.c\
check_tpl.c zip.c del_dir.c make_matrices.c\
mem_position.c probe.c add_membrane.c    load_pdb_no_param.c ga_engine.c rotamers_ga.c compute_patches.c

OBJ     = $(SRC:.c=.o)

HEADER  = mcce.h

$(LIB): $(OBJ)
	$(AR) $(ARFLAGS) $(LIB) $(OBJ)

$(OBJ): $(HEADER)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -f *.o mcce.a

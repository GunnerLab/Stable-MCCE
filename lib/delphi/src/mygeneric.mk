SRC1=qdiff4v.f objvertices.f distobj.f elb.f up.f \
 phintp4.f scaler4.f setbc4.f conplt.f itit4j.f \
 wrteps4b.f relfac4b.f react2.f clb.f setcrg.f setfcrg.f nitit.f \
 cputime.f setrc4d.f  rdiv.f chkcrg.f expand4.f \
 rdlog4.f wrtprm.f  off4.f extrm.f wrtatd.f \
 distrTOpoint.f wrtvisual.f extrmobjects.f\
 nlener.f timef.f \
 grdatm.f crgarr.f dbsfd.f  mkdbsf.f phicon.f wrtphi.f wrtsit4.f \
 encalc4.f wrtgcrg.f  namlen3.f watput.f qinttot.f cfind4.f rfind4.f \
 radass4.f crgass4.f ichash4.f irhash4.f \
 form2.f  watpte.f omalt.f rforce.f  ts.f setout.f epsmak.f
#
SRC2=anagrd4.f anasurf.f cent4.f rdhcrg.f rdhrad.f rent4.f wrt4.f \
 getatm2.f
#
SRC3=vwtms2.f scale.f indver.f sas.f  cube.f  msrf.f \
     mkvtl.f ex.f fxhl.f wrtsurf.f
#
SRC4= callerF.f
#
SRC5=memalloc.c
#
SRC6=creapdb.c
#----------------------------------------------------
OBJ1=$(SRC1:.f=.o)
OBJ2=$(SRC2:.f=.o)
OBJ3=$(SRC3:.f=.o)
OBJ4=$(SRC4:.f=.o)
OBJ5=$(SRC5:.c=.o)
OBJ6=$(SRC6:.c=.o)
#----------------------------------------------------
.f.o:
	$(FC) $(FLAGS) -c  $(VPATH)/$*.f
.c.o:
	$(CC) $(CFLAGS) -c  $(VPATH)/$*.c
#----------------------------------------------------
default: delphi
#----------------------------------------------------
libdelphi:$(OBJ1) $(OBJ3) $(OBJ4) $(OBJ5) $(SRC1) $(SRC3) $(SRC4) $(SRC5)
	$(AR) $(ARFLAGS) delphi.a $(OBJ1) $(OBJ3) $(OBJ4) $(OBJ5)
	$(RANLIB) delphi.a
#----------------------------------------------------
clean:
	rm -f *.o delphi.a
#----------------------------------------------------
delphi:$(OBJ1) $(OBJ2) $(OBJ3) $(OBJ5) $(OBJ6) $(SRC1) $(SRC2) $(SRC3) $(SRC5) $(SRC6)
	$(FC) $(LFLAGS) -o $@ $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ5) $(OBJ6)
#----------------------------------------------------
delphiroutine:$(OBJ1) $(OBJ3) $(OBJ4) $(OBJ5) $(SRC1) $(SRC3) $(SRC4) $(SRC5)
	$(FC) $(LFLAGS) -o $@ $(OBJ1) $(OBJ3) $(OBJ4) $(OBJ5)
#----------------------------------------------------
$(OBJ1) : $(VPATH)/qdiffpar5.h $(VPATH)/qlog.h $(VPATH)/pointer.h
$(OBJ2) : $(VPATH)/qdiffpar5.h $(VPATH)/qlog.h $(VPATH)/pointer.h
$(OBJ3) : $(VPATH)/acc2.h $(VPATH)/pointer.h
$(OBJ4) : $(VPATH)/qdiffpar5.h $(VPATH)/qlog.h $(VPATH)/pointer.h
#----------------------------------------------------

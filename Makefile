CC      = gcc -g 
LDIR = lib
DEPS    = $(LDIR)/mcce.h
LIB     = $(LDIR)/mcce.a
AR      = ar
ARFLAGS = rvs

SRC = $(wildcard $(LDIR)/*.c)
DELPHI = bin/delphi

bin/mcce: mcce.c $(LIB) $(DEPS) $(DELPHI)
	$(CC) -o bin/mcce mcce.c $(LIB) -lm 

$(DELPHI): $(LDIR)/delphi/delphi
	$(MAKE) -C $(LDIR)/delphi
	cp $(LDIR)/delphi/delphi $(DELPHI)

OBJ     = $(SRC:.c=.o)

$(LIB): $(OBJ)
	cd $(LDIR)
	$(AR) $(ARFLAGS) $(LIB) $(OBJ)

$(LDIR)/%.o: $(LDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -I/usr/include/glib-2.0


clean:
	-rm -f bin/mcce bin/delphi $(LIB) $(LDIR)/*.o
	$(MAKE) clean -C $(LDIR)/delphi

cleanbin/mcce:
	-rm -f bin/mcce bin/delphi $(LIB) $(LDIR)/*.o
	$(MAKE) clean -C $(LDIR)/delphi


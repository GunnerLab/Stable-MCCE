VPATH=src
#LINUX INTEL
  FC=gfortran -O3 -static -fcray-pointer -ffixed-line-length-none
  CC=gcc -O3

# Double precision performance
#RECFLAGS = -O3 -I$(VPATH)/ -DLINUX -DIFC -DDP -fdefault-real-8
#FLAGS = -O3 -I$(VPATH)/   -DLINUX -DIFC  -DDP -fdefault-real-8
RECFLAGS = -O3 -I$(VPATH)/ -DLINUX -DIFC -DDP 
FLAGS = -O3 -I$(VPATH)/   -DLINUX -DIFC  -DDP
CFLAGS = -O3  -DLINUX 
LFLAGS = -O3  -L/usr/lib/ 
#RECFLAGS = -g -I$(VPATH)/ -DLINUX -DIFC -DDP 
#FLAGS = -g -I$(VPATH)/   -DLINUX -DIFC  -DDP
#CFLAGS = -g  -DLINUX 
#LFLAGS = -g  -L/usr/lib/ 


include $(VPATH)/genericLinux.mk

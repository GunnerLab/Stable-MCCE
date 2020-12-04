VPATH=src
#LINUX INTEL
#  FC=gfortran -O3 -fcray-pointer -ffixed-line-length-none
#  CC=gcc -O3
#  FC=gfortran -g -static -fcray-pointer -ffixed-line-length-none
  FC=x86_64-apple-darwin13.4.0-gfortran -fcray-pointer -ffixed-line-length-none
  CC=x86_64-apple-darwin13.4.0-clang

# Double precision performance
#RECFLAGS = -O3 -I$(VPATH)/ -DLINUX -DIFC -DDP -fdefault-real-8
#FLAGS = -O3 -I$(VPATH)/   -DLINUX -DIFC  -DDP -fdefault-real-8

#RECFLAGS = -O3 -I$(VPATH)/ -DLINUX -DIFC -DDP 
#FLAGS = -O3 -I$(VPATH)/   -DLINUX -DIFC  -DDP
#CFLAGS = -O3  -DLINUX 
#LFLAGS = -O3  -L/usr/lib/ 

RECFLAGS = -g -I$(VPATH)/ -DLINUX -DIFC -DDP 
FLAGS = -g -I$(VPATH)/   -DLINUX -DIFC  -DDP
CFLAGS = -g  -DLINUX 
LFLAGS = -g 


include $(VPATH)/genericLinux.mk

# makefile for R
SHELL=/bin/sh

SRC=checkbox.c ds11.f fastcov.c gamrob.f glmrob.f lmrobbi.f \
    lmrobml.f lmrobmm.f mmprnt.c rb.c regar.f rocovm.f rqbr.f rqfn.f \
    tdmean.f tmlfor.f wblrob.f xlapack.f Rcompatibility.c \
		fastmcd.f sfastmcd.c

OBJ=checkbox.o ds11.o fastcov.o gamrob.o glmrob.o lmrobbi.o \
    lmrobml.o lmrobmm.o mmprnt.o rb.o regar.o rocovm.o rqbr.o rqfn.o \
    tdmean.o tmlfor.o wblrob.o xlapack.o Rcompatibility.o \
		fastmcd.o sfastmcd.o

NOTBUILT=lmrobfs.f

all : $(OBJ)
	$(R_HOME)/bin/R CMD SHLIB -o robust$(SHLIB_EXT) $(OBJ)

clean :
	rm -f robust$(SHLIB_EXT) $(OBJ)



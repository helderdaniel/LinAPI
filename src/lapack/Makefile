#Makefile to build static and shared LinAPI libraries
#Lapack version

#==== ATLAS and LAPACK Libs =====

#include ATLAS BLAS Makefile base definitions
#check if it is really needed
#include /shared/opt/atlas-3.8.2/Make.inc

#LIB dirs
ATLASdir = /shared/opt/atlas-3.8.2
LAPACKdir = /shared/opt/lapack-3.1.1

#ATLAS original LAPACK lib (do not include ATLAS BLAS optimized functions)
#To compile with this accessing F77 intreface LAPACKLLIB+BLASLIB must be linked
#LAPACKLIB = -llapack
#or
#LAPACKLIB = $(LAPACKdir)/lapack_LINUX.a

#ATLAS optimized BLAS + LAPACK interfaces in C and F77
#CLAPACKATLLIB is full ATLAS POWERED LAPACK F77 LIB + C ATLAS interface to LAPA
#So for C77 or C (limited) interfaces CLAPACKATLLIB+BLASLIB+CBLASLIB
CLAPACKATLLIB = $(ATLASdir)/lib/liblapack.a
BLASLIB = -L$(ATLASdir)/lib/ -lf77blas -latlas
CBLASLIB = -L$(ATLASdir)/lib/ -lcblas -latlas
#To use threadded lybraries for SMP: -lptf77blas or -lptcblas

LAPLIBS=$(CBLASLIB) $(CLAPACKATLLIB) $(BLASLIB) $(LAPACKLIB)

#==== end ATLAS and LAPACK Libs =====

CC = g++
CEXT = cpp

SRCDIR  = .
CMNSRCDIR = ..
# := assignement copy value of vars when issued, so if $SRCDIR changed in 
# future its value on $SRC will be the one when assigned.
# = assignemet sets reference to vars, so if $SRCDIR change, $SRC will
# change accordingly
#generate list with all <SRCDIR>/*.cpp (files on cur dir)
SPCSRC := $(wildcard $(SRCDIR)/*.$(CEXT))

#generate list with all <CMNSRCDIR>/*.cpp (files on parent dir)
CMNSRC := $(wildcard $(CMNSRCDIR)/*.$(CEXT))

SRC := $(SPCSRC) $(CMNSRC)

#Eliminate from $(SRC) *.cpp of templates (which are included in *.h
#until gcc support for export keyword avaiable) so that it will not 
#be compiled again
#To remove when gcc support for export keyword is available
SRC := $(patsubst ./operator2d%.cpp,,$(SRC))
SRC := $(patsubst ./linmatrix%.cpp,,$(SRC))
SRC := $(patsubst ./linvector%.cpp,,$(SRC))

#generate list equal to $(SCR) but with extension *.o and ../ subst as ./
OBJ := $(patsubst %.$(CEXT),%.o,$(SRC))
OBJ := $(patsubst ../%,./%,$(OBJ))

#Library location and type
LIBTYPE = lapack
LIBIDMACRO = _LAPACK_

#static version data
LIBNAME   = linapi
LIBFILEst = lib$(LIBNAME).a

#Shared version data
Major     = 1
Minor     = 0
SubMinor  = 1
COMPLIB   = lib$(LIBNAME).so
RUNLIB    = $(COMPLIB).$(Major)
LIBFILEsh = $(COMPLIB).$(Major).$(Minor).$(SubMinor)

CCFLAGS  += -Wall -D$(LIBIDMACRO) -I./ -I../

all: static shared

static:
	@#-M tells compiler to look on source file to find all #includes and generating
	@#a rule like: main.o: main.cpp defs.h stdlib.h ....
	@$(CC) -M $(SRC) $(CCFLAGS)
	$(CC) -static -c $(SRC) $(CCFLAGS) $(ICCFLAGS) -DATL_USEPTHREADS $(LAPLIBS)
	ar rcs $(LIBFILEst) $(OBJ)
	#To include lapack libs on LinApi archive
	#ar x $(CBLASLIB)
	#ar x $(CLAPACKATLLIB)
	#ar x $(BLASLIB)
	#ar x $(LAPACKLIB)
	#ar rcs $(LIBFILEst) *.o
	
shared:
	@#-M tells compiler to look on source file to find #includes generatin
	@#a rule like: main.o: main.cpp defs.h stdlib.h ....
	@$(CC) -M $(SRC) $(CCFLAGS)
	$(CC) -c -fPIC $(SRC) $(CCFLAGS) $(ICCFLAGS) -DATL_USEPTHREADS $(LAPLIBS)
	$(CC) -shared -Wl,-soname,$(RUNLIB) -o $(LIBFILEsh) $(OBJ)
	#To include lapack libs on LinApi archive
	#ar x $(CBLASLIB)
	#ar x $(CLAPACKATLLIB)
	#ar x $(BLASLIB)
	#ar x $(LAPACKLIB)
	#$(CC) -shared -Wl,-soname,$(RUNLIB) -o $(LIBFILEsh) *.o
	@-ln -s $(LIBFILEsh) $(COMPLIB)
	@-ln -s $(LIBFILEsh) $(RUNLIB)
                                                        
clean:
        #-rm ignores rm cmd if file do not exist
	@-rm *~ *.o 

cleanall: 
	$(MAKE) clean
	#clean static lib
	@-rm $(LIBFILEst)
	#clean shared (dynamic) lib
	@-rm $(LIBFILEsh) $(COMPLIB) $(RUNLIB)
	@echo $(CCFLAGS)

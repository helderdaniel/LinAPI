#Makefile to build static and shared LinAPI libraries
#Scalapack version

#include BLACS+MPI+LAPACK+SCALAPACK Makefile definitions
include /shared/opt/scalapack-1.8.0/SLmake.inc

#Enable C++ compiler wrapped for mpich2 should not be required
#since $(RLIBS) or $(PBLIBS) already include MPICH2 LIB. g++ should be enough
#however includes for MPI (mpi.h for instance) are not set
#CC = g++
CC = /shared/opt/mpich2/bin/mpicxx
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
LIBTYPE = scalapack
LIBIDMACRO = _SCALAPACK_

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
	$(CC) -static -c $(SRC) $(CCFLAGS) 
	ar rcs $(LIBFILEst) $(OBJ) 
	
shared:
	@#-M tells compiler to look on source file to find #includes generatin
	@#a rule like: main.o: main.cpp defs.h stdlib.h ....
	@$(CC) -M $(SRC) $(CCFLAGS)
	$(CC) -c -fPIC $(SRC) $(CCFLAGS)
	$(CC) -shared -Wl,-soname,$(RUNLIB) -o $(LIBFILEsh) $(OBJ)
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

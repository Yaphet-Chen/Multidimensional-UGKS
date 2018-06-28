SHELL = /bin/bash -O extglob

#--------------------------------------------------
#variables
#--------------------------------------------------
#compilation flag
OMP = yes
# FC = ifort
FC = gfortran
FCFLAGS = -module $(BIN) -O3
OMPFLAG = -openmp -parallel -fpp

ifeq ($(FC),gfortran)
    # FCFLAGS = -J$(BIN) -O3 -ffree-line-length-none
	FCFLAGS = -J$(BIN) -O3 -ffree-line-length-none -march=native
    OMPFLAG = -fopenmp
	# gfortran -ffree-line-length-none -march=native -funroll-loops -flto -pipe -O3
endif

ifeq ($(OMP),yes)
    FCFLAGS += $(OMPFLAG)
endif

#directory for executables
BIN = bin

#--------------------------------------------------
#compiling
#--------------------------------------------------
all: checkdir StationaryShockStructure

#mkdir
checkdir:
	mkdir -p $(BIN)

#build executables
StationaryShockStructure: checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/StationaryShockStructure src/StationaryShockStructure.f90

#clean
clean:
	rm -f bin/*

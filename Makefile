SHELL = /bin/bash -O extglob

#--------------------------------------------------
#variables
#--------------------------------------------------
#compilation flag
OMP = yes
# FC = ifort
FC = gfortran
FCFLAGS = -module $(BIN) -O3 -fpconstant
OMPFLAG = -openmp -parallel -fpp

ifeq ($(FC),gfortran)
	FCFLAGS = -J$(BIN) -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -march=native 
    OMPFLAG = -fopenmp
	# gfortran -ffree-line-length-none -march=native -funroll-loops -flto -pipe -O3
	# -Wall -Wextra -fcheck=all -fbacktrace
endif

ifeq ($(OMP),yes)
    FCFLAGS += $(OMPFLAG)
endif

#directory for executables
BIN = bin

#--------------------------------------------------
#compiling
#--------------------------------------------------
all: checkdir StationaryShockStructure SodShockTube Cavity BoundaryLayer Cavity_LocalTimeStepping Oscillatory_Cavity

#mkdir
checkdir:
	mkdir -p $(BIN)

#build executables
StationaryShockStructure: checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/StationaryShockStructure src/StationaryShockStructure.f90

SodShockTube: checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/SodShockTube src/SodShockTube.f90

Cavity: checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/Cavity src/Cavity.f90

BoundaryLayer: checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/BoundaryLayer src/BoundaryLayer.f90

Cavity_LocalTimeStepping: checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/Cavity_LocalTimeStepping src/Cavity_LocalTimeStepping.f90

Oscillatory_Cavity: checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/Oscillatory_Cavity src/Oscillatory_Cavity.f90

#clean
clean:
	rm -f bin/*

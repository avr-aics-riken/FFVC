###################################################################
#
# FFV : Frontflow / violet Cartesian
#
# Copyright (c) 2012-2013  All right reserved.
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

########################
# INTEL (gcc, build openmpi with gcc)
########################

AR          = ar cru
RANLIB      = ranlib
RM          = \rm -f
MPI_DIR	    = /opt/openmpi
CPM_DIR     = /usr/local/CPMlib
TP_DIR      = /usr/local/TextParser
POLYLIB_DIR = ../Polylib-2.4
CUTLIB_DIR  = ../Cutlib-2.0.5
PMLIB_DIR   = ../PMlib-1.7
CIOLIB_DIR  = ../CIOlib-1.0
OMP_FLAGS   = -fopenmp 
UDEF_OPT    =
CC          = mpicc
CFLAGS      = -O3
CXX         = mpicxx
CXXFLAGS    = -O3 $(OMP_FLAGS)
FC          = mpif90
FCFLAGS     = -O3 -ffree-form
F90         = mpif90
F90FLAGS    = -O3 -fpp $(OMP_FLAGS) -ffree-form
LDFLAGS     = 
LIBS        =
#LIBS        = -lgfortran
UDEF_LIB_PATH_SPEC =
UDEF_LIBS_SPEC     =

## iff double
#CFLAGS     += -D_REAL_IS_DOUBLE_
#CXXFLAGS   += -D_REAL_IS_DOUBLE_
#FCFLAGS    += 
#F90FLAGS   += 

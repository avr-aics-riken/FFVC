###################################################################
#
# CIOlib - Cartesian Input / Output library
#
# Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
# All right reserved.
#
###################################################################

#  If CIOlib is placed inside the directory of FFVC solver, this make_setting is ignored.


########################
# FX10, K
########################

AR          = ar cr
RANLIB      = ranlib
RM          = \rm -f
MPI_DIR     =
TP_DIR      = /usr/local/TextParser
OMP_FLAGS   =
CC          = mpifccpx
CFLAGS      = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,openmp -Nsrc
CXX         = mpiFCCpx
CXXFLAGS    = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,openmp,optmsg=2 -V -Nsrc
FC          = mpifrtpx
FCFLAGS     = -Cpp -Kfast,ocl,preex,simd=2,uxsimd,array_private,auto,parallel,openmp -Qt
F90         = mpifrtpx
F90FLAGS    = -Cpp -Kfast,ocl,preex,simd=2,uxsimd,array_private,auto,parallel,openmp -Qt
LDFLAGS     = --linkfortran
UDEF_INC_PATH = -I$(TP_DIR)/include -I$(MPI_DIR)/include

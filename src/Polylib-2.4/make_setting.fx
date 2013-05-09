###################################################################
#
# Polylib - Polygon Management library
#
# Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
# All right reserved.
#
###################################################################

#  If Polylib is placed inside the directory of FFVC solver, this make_setting is ignored.


########################
# FX10, K
########################

AR          = ar cr
RANLIB      = ranlib
RM          = \rm -f
MPI_DIR     =
TP_DIR      = /usr/local/TextParser
OMP_FLAGS   =
CXX         = mpiFCCpx
CXXFLAGS    = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,openmp,optmsg=2 -V -Nsrc
UDEF_INC_PATH = -I$(TP_DIR)/include -I$(MPI_DIR)/include

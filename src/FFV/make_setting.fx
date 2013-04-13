###################################################################
#
# FFV : Frontflow / violet Cartesian
#
# Copyright (c) 2012-2013  All right reserved.
#
# Institute of Industrial Science, University of Tokyo, Japan. 
#
###################################################################

########################
# FX10, K
########################

AR          = ar cr
RANLIB      = ranlib
RM          = \rm -f
MPI_DIR     =
CPM_DIR     = /hoge
TP_DIR      = /funya
POLYLIB_DIR = ../Polylib-2.3
CUTLIB_DIR  = ../Cutlib-2.0.5
PMLIB_DIR   = ../PMlib-1.6
CIOLIB_DIR  = ../CIOlib-1.0
OMP_FLAGS   =
UDEF_OPT    = -D__K_FPCOLL -D__ARCH_FX
CC          =  mpifccpx
CFLAGS      = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,openmp -Nsrc
CXX         = mpiFCCpx
CXXFLAGS    = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,openmp,optmsg=2 -V -Nsrc
FC          = mpifrtpx
FCFLAGS     = -Cpp -Kfast,ocl,preex,simd=2,uxsimd,array_private,auto,parallel,openmp -Qt
F90         = mpifrtpx
F90FLAGS    = -Cpp -Kfast,ocl,preex,simd=2,uxsimd,array_private,auto,parallel,openmp -Qt
LDFLAGS     = --linkfortran
LIBS        =
UDEF_LIB_PATH_SPEC =
UDEF_LIBS_SPEC     =

## iff double
#CFLAGS     += -D_REAL_IS_DOUBLE_
#CXXFLAGS   += -D_REAL_IS_DOUBLE_
#FCFLAGS    += -CcdRR8
#F90FLAGS   += -CcdRR8

## iff different restart with staging
#CFLAGS     += -D_STAGING_
#CXXFLAGS   += -D_STAGING_
#FCFLAGS    += 
#F90FLAGS   += 

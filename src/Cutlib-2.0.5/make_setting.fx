###################################################################
#
# Cutlib - Cut Information Library
#
# Copyright (c) 2010-2013 Advanced Institute for Computational Science, RIKEN.
# All right reserved.
#
###################################################################

#  If Cutlib is placed inside the directory of FFVC solver, this make_setting is ignored.


########################
# FX10, K
########################

AR          = ar cr
RANLIB      = ranlib
RM          = \rm -f
TP_DIR      = /usr/local/TextParser
OMP_FLAGS   =
CXX         = mpiFCCpx
CXXFLAGS    = -Kfast,ocl,preex,simd=2, array_private,parallel,openmp,optmsg=2 -V -Nsrc
POLYLIB_DIR = ../Polylib-2.4


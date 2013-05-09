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
# GNU
########################

AR          = ar cru
RANLIB      = ranlib
RM          = \rm -f
TP_DIR      = /usr/local/TextParser
OMP_FLAGS   =
CXX         = mpicxx
CXXFLAGS    = -O3 $(OMP_FLAGS) 
POLYLIB_DIR = ../Polylib-2.4

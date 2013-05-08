###################################################################
#
# PMlib - Performance Monitor library
# 
# Copyright (c) 2010-2013 Advanced Institute for Computational Science, RIKEN.
# All right reserved.
#
###################################################################

#  If PMlib is placed inside the directory of FFVC solver, this make_setting is ignored.


##########
# GNU
##########

AR          = ar cru
RANLIB      = ranlib
RM          = \rm -f
MPI_DIR	    = /opt/openmpi
CXX         = mpicxx
CXXFLAGS    = -O3 -Wall
UDEF_INC_PATH = -I$(MPI_DIR)/include

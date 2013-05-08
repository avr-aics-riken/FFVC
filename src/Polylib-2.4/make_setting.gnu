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
# GNU
########################

AR          = ar cru
RANLIB      = ranlib
RM          = \rm -f
MPI_DIR	    = /opt/openmpi
TP_DIR      = /usr/local/TextParser
OMP_FLAGS   =
CXX         = mpicxx
CXXFLAGS    = -O3 $(OMP_FLAGS) 
UDEF_INC_PATH = -I$(TP_DIR)/include -I$(MPI_DIR)/include -I../include

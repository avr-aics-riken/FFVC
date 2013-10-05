###################################################################################
#
# FFV-C : Frontflow / violet Cartesian
#
# Copyright (c) 2007-2011 VCAD System Research Program, RIKEN. 
# All rights reserved.
#
# Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
# All rights reserved. 
#
###################################################################################

noinst_lib_LIBRARIES = libPLOT3D.a

noinst_libdir = ../lib

libPLOT3D_a_CXXFLAGS = \
  -I. \
  -I../FB \
  -I../F_CORE \
  @CPM_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@ \
  @TP_CFLAGS@ \
  @FREALOPT@ 

libPLOT3D_a_SOURCES = \
  dfi.C \
  dfi.h \
  PLOT3D.C \
  PLOT3D.h \
  PLOT3D_read.C \
  PLOT3D_read.h \
  PLOT3D_write.C \
  PLOT3D_write.h \
  ffv_PLOT3D.C \
  ffv_PLOT3D.h \
  PLT3D.f90 \
  PLT3D_read.f90 \
  PLT3D_write.f90


EXTRA_DIST = Makefile_hand depend.inc ReleaseNote_PLT3D.txt

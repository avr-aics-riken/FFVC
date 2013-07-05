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

bin_PROGRAMS = combsph

combsph_CXXFLAGS = \
  -I. \
  -I../FB \
  -I../PLOT3D \
  @CPM_CFLAGS@ \
  @CIO_CFLAGS@ \
  @TP_CFLAGS@ \
  @PM_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@

combsph_SOURCES = \
  COMB_Define.h \
  FileIO_read_sph.C \
  FileIO_sph.C \
  FileIO_sph.h \
  comb.C \
  comb.h \
  comb_plot3d.C \
  comb_sph.C \
  endianUtil.h \
  main.C \
  util_Path.h


EXTRA_DIST = Makefile_hand depend.inc README

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
  @TP_CFLAGS@ \
  @PM_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@ \
  @FREALOPT@

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

combsph_LDADD = \
  @CPM_LDFLAGS@ \
  @TP_LDFLAGS@ \
  @PM_LDFLAGS@ \
  @MPICH_LDFLAGS@ \
  @MPICH_LIBS@ \
  @MPICH_FCLIBS@ \
  @FFVC_LIBS@ \
  @FFVC_FC_LD@ \
  -L../FB -lFB \
  -L../PLOT3D -lPLOT3D


EXTRA_DIST = Makefile_hand depend.inc README

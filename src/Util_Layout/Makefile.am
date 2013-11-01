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

bin_PROGRAMS = layout



layout_CXXFLAGS = \
  -I. \
  -I../FB \
  @CPM_CFLAGS@ \
  @TP_CFLAGS@ \
  @PM_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@ \
  @FREALOPT@

layout_SOURCES = \
  LAYOUT_Define.h \
  layout.C \
  layout.h \
  main.C

layout_LDADD = \
  -L../FB -lFB \
  -L../PLOT3D -lPLOT3D \
  @CPM_LDFLAGS@ \
  @TP_LDFLAGS@ \
  @PM_LDFLAGS@ \
  @MPICH_LDFLAGS@ \
  @MPICH_LIBS@ \
  @MPICH_FCLIBS@ \
  @FFVC_LIBS@ \
  @FFVC_FC_LD@


EXTRA_DIST = Makefile_hand depend.inc README

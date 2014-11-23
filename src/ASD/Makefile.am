###################################################################################
#
# FFV-C : Frontflow / violet Cartesian
#
# Copyright (c) 2007-2011 VCAD System Research Program, RIKEN. 
# All rights reserved.
#
# Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
# All rights reserved. 
#
###################################################################################

noinst_lib_LIBRARIES = libASD.a

noinst_libdir = ../lib

libASD_a_CXXFLAGS = \
  -I. \
  -I../FB \
  -I../Geometry \
  @CPM_CFLAGS@ \
  @CUT_CFLAGS@ \
  @PL_CFLAGS@ \
  @TP_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@

libASD_a_SOURCES = \
  SubDomain.C \
  SubDomain.h \
  ASDmodule.C \
  ASDmodule.h

EXTRA_DIST = Makefile_hand depend.inc


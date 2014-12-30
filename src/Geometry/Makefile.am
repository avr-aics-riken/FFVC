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

noinst_lib_LIBRARIES = libGEOM.a

noinst_libdir = ../lib

libGEOM_a_CXXFLAGS = \
  -I. \
  -I../FB \
  @CPM_CFLAGS@ \
  @PL_CFLAGS@ \
  @TP_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@

libGEOM_a_SOURCES = \
  Geometry.C \
  Geometry.h \
  CompoFraction.C \
  CompoFraction.h \
  fill_bid_naive.h \
  Glyph.C \
  Glyph.h

EXTRA_DIST = Makefile_hand depend.inc fill_bid.h


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

noinst_lib_LIBRARIES = libFLS.a

noinst_libdir = ../lib

libFLS_a_FCFLAGS = \
  -I. \
  -I../FB \
  @CPM_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @FREALOPT@ 

libFLS_a_SOURCES = \
  ffv_LSfunc.h \
  ffv_blas.f90 \
  ffv_bicg.f90 \
  ffv_poisson.f90 \
  ffv_poisson2.f90

EXTRA_DIST = Makefile_hand depend.inc ffv_poisson_cds.f90

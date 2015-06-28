###################################################################################
#
# FFV-C : Frontflow / violet Cartesian
#
# Copyright (c) 2007-2011 VCAD System Research Program, RIKEN. 
# All rights reserved.
#
# Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved. 
#
###################################################################################

noinst_lib_LIBRARIES = libFLS.a

noinst_libdir = ../lib

libFLS_a_FCFLAGS = \
  -I$(top_srcdir)/src/F_LS \
  -I$(top_srcdir)/src/FB \
  @CPM_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @FREALOPT@ 

libFLS_a_SOURCES = \
  ffv_LSfunc.h \
  ffv_blas.f90 \
  ffv_SOR.f90 \
  core_psor.h \
  core_sor2sma.h


EXTRA_DIST = Makefile_hand depend.inc ffv_poisson_cds.f90 ffv_poisson2.f90 ffv_rc.f90

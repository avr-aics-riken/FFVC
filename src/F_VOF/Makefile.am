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

lib_LIBRARIES = libFCORE.a

libFCORE_a_FCFLAGS = \
  -I. \
  -I../FB \
  @MPICH_CFLAGS@ \
  @FREALOPT@ 

libFCORE_a_SOURCES = \
  d_o_o_p.h \
  ffv_Ffunc.h \
  ffv_blas.f90 \
  ffv_forcing.f90 \
  ffv_pbc.f90 \
  ffv_poisson.f90 \
  ffv_poisson2.f90 \
  ffv_poisson_cds.f90 \
  ffv_pscalar.f90 \
  ffv_utility.f90 \
  ffv_vbc_inner.f90 \
  ffv_vbc_outer.f90 \
  ffv_velocity_binary.f90 \
  ffv_velocity_cds.f90 \
  force.h \
  load_var_stencil5.h \
  muscl.h

libFCORE_LDADD = \
    @CPM_LDFLAGS@ \
    @MPICH_LDFLAGS@ \
    @MPICH_FCLIBS@

EXTRA_DIST = Makefile_hand depend.inc unused.f90

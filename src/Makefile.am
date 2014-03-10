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

bin_PROGRAMS = ffvc

ffvc_CXXFLAGS = \
  -I. \
  -I./FB \
  -I./F_CORE \
  -I./IP \
  -I./ASD \
  -I./FFV \
  @CPM_CFLAGS@ \
  @CIO_CFLAGS@ \
  @TP_CFLAGS@ \
  @PM_CFLAGS@ \
  @PL_CFLAGS@ \
  @CUT_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@

ffvc_SOURCES = main.C

ffvc_LDADD = \
  -L./FB -lFB \
  -L./IP -lIP \
  -L./F_CORE -lFCORE \
  -L./ASD -lASD \
  -L./FFV -lFFV \
  @CPM_LDFLAGS@ \
  @CIO_LDFLAGS@ \
  @TP_LDFLAGS@ \
  @PM_LDFLAGS@ \
  @PL_LDFLAGS@ \
  @CUT_LDFLAGS@ \
  @MPICH_LDFLAGS@ \
  @MPICH_LIBS@ \
  @MPICH_FCLIBS@ \
  @FFVC_LIBS@

EXTRA_DIST = Makefile_hand Makefile_ffvc depend.inc make_setting

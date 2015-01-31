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

bin_PROGRAMS = ffvc

ffvc_CXXFLAGS = \
  -I. \
  -I./FFV \
  -I./F_CORE \
  -I./F_LS \
  -I./FILE_IO \
  -I./ASD \
  -I./IP \
  -I./FB \
  -I./Geometry \
  @CDM_CFLAGS@ \
  @CPM_CFLAGS@ \
  @PL_CFLAGS@ \
  @PM_CFLAGS@ \
  @TP_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@

ffvc_SOURCES = main.C

ffvc_LDADD = \
  -L./FFV -lFFV \
  -L./FILE_IO -lFIO \
  -L./F_CORE -lFCORE \
  -L./F_LS -lFLS \
  -L./ASD -lASD \
  -L./Geometry -lGEOM \
  -L./IP -lIP \
  -L./FB -lFB \
  @CDM_LDFLAGS@ \
  @CPM_LDFLAGS@ \
  @PL_LDFLAGS@ \
  @PM_LDFLAGS@ \
  @TP_LDFLAGS@ \
  @MPICH_LDFLAGS@ \
  @MPICH_LIBS@ \
  @MPICH_FCLIBS@ \
  @FFVC_LIBS@

EXTRA_DIST = Makefile_hand Makefile_ffvc depend.inc make_setting

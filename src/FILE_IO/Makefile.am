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

noinst_lib_LIBRARIES = libFIO.a

noinst_libdir = ../lib

libFIO_a_CXXFLAGS = \
  -I$(top_srcdir)/src/FILE_IO \
  -I$(top_srcdir)/src/FB \
  -I$(top_srcdir)/src/F_CORE \
  -I$(top_srcdir)/src/F_LS \
  -I$(top_srcdir)/src/FFV \
  @CDM_CFLAGS@ \
  @CPM_CFLAGS@ \
  @PL_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @PM_CFLAGS@ \
  @TP_CFLAGS@ \
  @REALOPT@


libFIO_a_SOURCES = \
  ffv_io_base.h \
  ffv_io_base.C \
  ffv_sph.h \
  ffv_sph.C \
  ffv_plot3d.h \
  ffv_plot3d.C \
  BlockSaver.C \
  BlockSaver.h \
  BitVoxel.C \
  BitVoxel.h \
  RLE.h \
  type.h \
  FileSystemUtil.C \
  FileSystemUtil.h \
  FileCommon.h


EXTRA_DIST = Makefile_hand depend.inc

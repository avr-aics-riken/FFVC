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

noinst_lib_LIBRARIES = libFFV.a

noinst_libdir = ../lib


libFFV_a_CXXFLAGS = \
  -I$(top_srcdir)/src/FFV \
  -I$(top_srcdir)/src/F_CORE \
  -I$(top_srcdir)/src/F_LS \
  -I$(top_srcdir)/src/IP \
  -I$(top_srcdir)/src/FB \
  -I$(top_srcdir)/src/Geometry \
  -I$(top_srcdir)/src/FILE_IO \
  @CDM_CFLAGS@ \
  @CPM_CFLAGS@ \
  @PL_CFLAGS@ \
  @PM_CFLAGS@ \
  @TP_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@

libFFV_a_SOURCES = \
  NS_FS_E_Binary.C \
  NS_FS_E_CDS.C \
  PS_Binary.C \
  ffv.C \
  ffv.h \
  ffv_Alloc.C \
  ffv_Alloc.h \
  ffv_Define.h \
  ffv_Filter.C \
  ffv_Heat.C \
  ffv_Initialize.C \
  ffv_LS.C \
  ffv_LS.h \
  ffv_Loop.C \
  ffv_Post.C \
  ffv_SetBC.C \
  ffv_SetBC.h \
  ffv_TerminateCtrl.C \
  ffv_TerminateCtrl.h \
  ffv_Version.h

EXTRA_DIST = Makefile_hand depend.inc

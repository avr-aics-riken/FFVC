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

noinst_lib_LIBRARIES = libFFV.a

noinst_libdir = ../lib


libFFV_a_CXXFLAGS = \
  -I. \
  -I../F_CORE \
  -I../F_LS \
  -I../IP \
  -I../FB \
  @CDM_CFLAGS@ \
  @CPM_CFLAGS@ \
  @CUT_CFLAGS@ \
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
  ffv_Geometry.C \
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

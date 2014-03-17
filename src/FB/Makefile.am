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

noinst_lib_LIBRARIES = libFB.a

noinst_libdir = ../lib

libFB_a_CXXFLAGS = \
  -I. \
  @CIO_CFLAGS@ \
  @CPM_CFLAGS@ \
  @TP_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@ \
  @FREALOPT@

libFB_a_SOURCES = \
  Alloc.C \
  Alloc.h \
  BndOuter.C \
  BndOuter.h \
  CompoFraction.C \
  CompoFraction.h \
  Component.C \
  Component.h \
  Control.C \
  Control.h \
  DataHolder.C \
  DataHolder.h \
  DomainInfo.h \
  FBUtility.C \
  FBUtility.h \
  FB_Define.h \
  FB_Ffunc.h \
  FB_util.f90 \
  FileIO.C \
  FileIO.h \
  FindexS3D.h \
  History.C \
  History.h \
  IntervalManager.h \
  Intrinsic.C \
  Intrinsic.h \
  IterationControl.C \
  IterationControl.h \
  Medium.h \
  MonCompo.C \
  MonCompo.h \
  Monitor.C \
  Monitor.h \
  ParseBC.C \
  ParseBC.h \
  ParseMat.C \
  ParseMat.h \
  Sampling.C \
  Sampling.h \
  SetBC.C \
  SetBC.h \
  VoxInfo.C \
  VoxInfo.h \
  basic_func.h \
  dbwrite.f90 \
  ffv_f_params.h \
  fill_bid.h \
  mydebug.h \
  util_Path.h \
  Vec3.h \
  Glyph.C \
  Glyph.h

EXTRA_DIST = Makefile_hand depend.inc ReleaseNote_FB.txt


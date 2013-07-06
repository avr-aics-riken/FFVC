###################################################################################
#
# FFV-C : Frontflow / violet Cartesian
#
# Copyright (c) 2007-2011 VCAD System Research Program, RIKEN. 
# All rights reserved.
#
# Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
# All rights reserved. 
#
###################################################################################

noinst_lib_LIBRARIES = libIP.a

noinst_libdir = ../lib

libIP_a_CXXFLAGS = \
  -I. \
  -I../FB \
  @CPM_CFLAGS@ \
  @TP_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@

libIP_a_SOURCES = \
  IP_Cylinder.C \
  IP_Cylinder.h \
  IP_Define.h \
  IP_Duct.C \
  IP_Duct.h \
  IP_Jet.C \
  IP_Jet.h \
  IP_PMT.C \
  IP_PMT.h \
  IP_PPLT2D.C \
  IP_PPLT2D.h \
  IP_Polygon.C \
  IP_Polygon.h \
  IP_Rect.C \
  IP_Rect.h \
  IP_SHC1D.C \
  IP_SHC1D.h \
  IP_Sphere.C \
  IP_Sphere.h \
  IP_Step.C \
  IP_Step.h

EXTRA_DIST = Makefile_hand depend.inc ReleaseNote_IP.txt IP_RSP.C IP_RSP.h

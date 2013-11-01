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

bin_PROGRAMS = ffvc


ffvc_CXXFLAGS = \
  -I. \
  -I../FB \
  -I../F_CORE \
  -I../IP \
  @CPM_CFLAGS@ \
  @CIO_CFLAGS@ \
  @TP_CFLAGS@ \
  @PM_CFLAGS@ \
  @PL_CFLAGS@ \
  @CUT_CFLAGS@ \
  @MPICH_CFLAGS@ \
  @REALOPT@

ffvc_SOURCES = \
  NS_FS_E_Binary.C \
  NS_FS_E_CDS.C \
  PS_Binary.C \
  ffv.C \
  ffv.h \
  ffv_Alloc.C \
  ffv_Define.h \
  ffv_Heat.C \
  ffv_Initialize.C \
  ffv_LS.C \
  ffv_Loop.C \
  ffv_Post.C \
  ffv_Restart.C \
  ffv_SetBC.C \
  ffv_SetBC.h \
  ffv_TerminateCtrl.C \
  ffv_TerminateCtrl.h \
  main.C \
  ffv_Version.h

ffvc_LDADD = \
  -L../FB -lFB \
  -L../IP -lIP \
  -L../F_CORE -lFCORE \
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


EXTRA_DIST = Makefile_hand depend.inc

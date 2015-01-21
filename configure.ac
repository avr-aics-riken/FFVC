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
#
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.57])
AC_INIT([FFVC], [2.0.5], [keno@riken.jp], [FFVC])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/main.C])
AC_CONFIG_HEADERS([config.h])

#
# revision No.
#
FFVC_REVISION="20150117_1800"
AC_SUBST(FFVC_REVISION)


#
# FORTRAN
#

FC=$F90
FCFLAGS=$F90FLAGS


#
# compiler
#
AC_ARG_WITH(comp, [AC_HELP_STRING([--with-comp=(INTEL|FJ|GNU)],[Specify Compiler type])], , with_comp=none)


#
# REAL type
#
AC_ARG_WITH(precision, [AC_HELP_STRING([--with-precision=(single|double)],[Specify REAL type [single]])], , with_precision=single)
AC_SUBST(REALOPT)
AC_SUBST(FREALOPT)

if test x"$with_precision" = x"double" ; then
  REALOPT=-D_REAL_IS_DOUBLE_
  case "$with_comp" in
    INTEL) FREALOPT="-r8" ;;
    FJ)    FREALOPT="-CcdRR8" ;;
    GNU)   FREALOPT="-fdefault-real-8" ;;
  esac
else
  REALOPT=
  FREALOPT=
fi



#
# MPI setting. wrappers like mpicxx, mpicc, mpiCC, mpic++, mpiFCCpx are taking into consideration, if not add in case.
#

enable_mpi="none"


FFVC_CXX="$CXX"
AC_SUBST(FFVC_CXX)

case "$FFVC_CXX" in
  mpicxx)   enable_mpi="wrappermpi" ;;
  mpicc)    enable_mpi="wrappermpi" ;;
  mpiCC)    enable_mpi="wrappermpi" ;;
  mpic++)   enable_mpi="wrappermpi" ;;
  mpiFCCpx) enable_mpi="wrappermpi" ;;
  CC)       enable_mpi="wrappermpi" ;;
esac

# mpiFCCpx ; Fujitsu FX
# CC       ; NERSC Hopper

#
# MPICH environment
#
AC_SUBST(MPICH_DIR)
AC_SUBST(MPICH_CFLAGS)
AC_SUBST(MPICH_LDFLAGS)
AC_SUBST(MPICH_LIBS)
AC_SUBST(MPICH_FCLIBS)



#
# if wrapper compiler is specified, skip this part
#

if test "$enable_mpi" = "none" ; then

# OpenMPI // ompが指定されていれば変数をセット
  AC_ARG_WITH(ompi, [AC_HELP_STRING([--with-ompi=dir],[Specify OpenMPI install directory])], , with_ompi=none)

  if test "$with_ompi" != "none" ; then
    MPICH_DIR=$with_ompi;
    MPICH_CFLAGS=-I$MPICH_DIR/include
    MPICH_LDFLAGS=-L$MPICH_DIR/lib
    enable_mpi="ompi"
  fi

fi



#
# MPI Library for Fortran
#

if test "$with_comp" = "FJ" ; then
  MPICH_LIBS=""
  case "$host" in
    sparc*) MPICH_FCLIBS="$MPICH_LIBS"" -lmpi_f77 -lmpi_f90" ;;
  esac
elif test "$enable_mpi" = "ompi" ; then
  MPICH_LIBS="-lmpi -lmpi_cxx"
  MPICH_FCLIBS="$MPICH_LIBS"" -lmpi_f77 -lmpi_f90"
fi



#
# warning
#

if test "$enable_mpi" != "wrappermpi" ; then
  if test "$with_ompi" = "none" ; then
    echo "  Warning: OpenMPI library is not specified."
    echo "           If compilation error occurs about MPI, specify --with-ompi option."
  fi
fi



FFVC_FC="$FC"
AC_SUBST(FFVC_FC)

FFVC_LIBS=""
AC_SUBST(FFVC_LIBS)

FFVC_FC_LD=""
AC_SUBST(FFVC_FC_LD)

#if test "$with_comp" = "FJ"; then
#  FFVC_FC_LD="$CXX"" --linkfortran"
#  FFVC_LIBS=" --linkfortran"
#elif test "$with_comp" = "INTEL" ; then
#  FFVC_LIBS="$FFVC_LIBS"" -lstdc++ -lifcore"
#elif test "$with_comp" = "GNU" ; then
#  FFVC_LIBS="$FFVC_LIBS"" -lstdc++ -lgfortran"
#else
#  FFVC_LIBS="$FFVC_LIBS"" -lstdc++"
#fi

case "$with_comp" in
  INTEL) FFVC_LIBS="$FFVC_LIBS"" -lstdc++ -lifcore" ;;
  FJ)    FFVC_FC_LD="$CXX"" --linkfortran" \
         FFVC_LIBS=" --linkfortran" ;;
  GNU)   FFVC_LIBS="$FFVC_LIBS"" -lstdc++ -lgfortran" ;;
esac

FFVC_LDFLAGS=""
AC_SUBST(FFVC_LDFLAGS)


#
# default prefix setting
#
AC_PREFIX_DEFAULT([/usr/local/FFVC])
if test x"$prefix" = x"NONE" ; then
  FFVC_INST_DIR=/usr/local/FFVC
else
  FFVC_INST_DIR="$prefix"
fi



#
# CPMlib
#
AC_ARG_WITH(cpm, [AC_HELP_STRING([--with-cpm=dir],[Specify CPMlib installed directory])], , with_cpm=none)
AC_SUBST(CPM_DIR)
AC_SUBST(CPM_CFLAGS)
AC_SUBST(CPM_LDFLAGS)
if test x"$with_cpm" != x"none" ; then
  CPM_DIR=$with_cpm;
  CPM_CFLAGS=`$CPM_DIR/bin/cpm-config --cflags`
  CPM_LDFLAGS=`$CPM_DIR/bin/cpm-config --libs`
else
  echo "  Error: CPM library is not specified."
  echo "         Specify --with-cpm option."
  exit 1
fi



#
# CDMlib
#
AC_ARG_WITH(cdm, [AC_HELP_STRING([--with-cdm=dir],[Specify CDMlib installed directory])], , with_cdm=none)
AC_SUBST(CDM_DIR)
AC_SUBST(CDM_CFLAGS)
AC_SUBST(CDM_LDFLAGS)
if test x"$with_cdm" != x"none" ; then
  CDM_DIR=$with_cdm;
  CDM_CFLAGS=`$CDM_DIR/bin/cdm-config --cflags`
  CDM_LDFLAGS=`$CDM_DIR/bin/cdm-config --libs`
else
  echo "  Error: CDM library is not specified."
  echo "         Specify --with-cdm option."
  exit 1
fi




#
# TextParser
#
AC_ARG_WITH(parser, [AC_HELP_STRING([--with-parser=dir],[Specify TextParser install directory])], , with_parser=none)
AC_SUBST(TP_DIR)
AC_SUBST(TP_CFLAGS)
AC_SUBST(TP_LDFLAGS)
if test x"$with_parser" != x"none" ; then
  TP_DIR=$with_parser;
  TP_CFLAGS=`$TP_DIR/bin/tp-config --cflags`
  TP_LDFLAGS=`$TP_DIR/bin/tp-config --libs`
else
  echo "  Error: TextParser library is not specified."
  echo "         Specify --with-parser option."
  exit 1
fi



#
# PMlib
#
AC_SUBST(PM_DIR)
AC_SUBST(PM_CFLAGS)
AC_ARG_WITH(pm, [AC_HELP_STRING([--with-pm=dir],[Specify PMlib install directory])], , with_pm=none)
AC_SUBST(PM_LDFLAGS)

if test x"$with_pm" != x"none" ; then
  PM_DIR=$with_pm;
  PM_CFLAGS=`$PM_DIR/bin/pm-config --cflags`
  PM_LDFLAGS=`$PM_DIR/bin/pm-config --libs`
else
  echo "  Error: PM library is not specified."
  echo "         Specify --with-pm option."
  exit 1
fi



#
# Polylib
#
AC_ARG_WITH(polylib, [AC_HELP_STRING([--with-polylib=dir],[Specify Polylib install directory])], , with_polylib=none)
AC_SUBST(PL_DIR)
AC_SUBST(PL_CFLAGS)
AC_SUBST(PL_LDFLAGS)

if test x"$with_polylib" != x"none" ; then
  PL_DIR=$with_polylib;
  PL_CFLAGS=`$PL_DIR/bin/polylib-config --cflags`
  PL_LDFLAGS=`$PL_DIR/bin/polylib-config --libs`
else
  echo "  Error: Polylib is not specified."
  echo "         Specify --with-polylib option."
  exit 1
fi




#
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src/IP src/FB src/F_CORE src/F_LS src/FILE_IO src/ASD src/Geometry src/FFV src src/Util_Combsph doc"



#
# Checks for programs.
#
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
AC_PROG_FC
AC_PROG_RANLIB


#
# Checks for libraries.
#



#
# Checks for header files.
#
AC_CHECK_HEADERS([fcntl.h float.h stdlib.h string.h strings.h sys/param.h unistd.h stdarg.h])



#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_SIZE_T
#AC_TYPE_INT32_T


#
# Checks for library functions.
#
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset mkdir pow realpath sqrt strcasecmp strchr strdup strerror va_start])



AC_CONFIG_FILES([Makefile \
                 doc/Makefile \
                 src/F_CORE/Makefile \
                 src/F_LS/Makefile \
                 src/FILE_IO/Makefile \
                 src/FB/Makefile \
                 src/IP/Makefile \
                 src/ASD/Makefile \
                 src/Geometry/Makefile \
                 src/FFV/Makefile \
                 src/Makefile \
                 src/Util_Combsph/Makefile \
                 src/FFV/ffv_Version.h \
                ])

AC_OUTPUT

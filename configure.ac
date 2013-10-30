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
#
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FFVC], [1.4.8], [keno@riken.jp], [FFVC])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/FB/Alloc.C])
AC_CONFIG_HEADERS([config.h])

#
# revision No.
#
FFVC_REVISION="20131030_1830"
AC_SUBST(FFVC_REVISION)


#
# default prefix setting
#
AC_PREFIX_DEFAULT([/usr/local/FFVC])
if test "$prefix" = "NONE" ; then
  FFVC_INST_DIR=/usr/local/FFVC
else
  FFVC_INST_DIR="$prefix"
fi


#
# platform
#
TGT_DEVICE="$($SHELL ./chk-uname)"
AC_SUBST(TGT_DEVICE)


#
# FORTRAN f90 compiler check
#

if test "_$FC" = "_"; then
  if test "_$F90" != "_"; then
    FC=$F90
  fi
fi
if test "_$FCFLAGS" = "_"; then
  if test "_$F90FLAGS" != "_"; then
    FCFLAGS=$F90FLAGS
  fi
fi



#
# compiler
#
AC_ARG_WITH(comp, [AC_HELP_STRING([--with-comp=(INTEL|FJ)],[Specify Compiler type])], , with_comp=none)



#
# REAL type
#
AC_ARG_WITH(precision, [AC_HELP_STRING([--with-precision=(single|double)],[Specify REAL type [single]])], , with_precision=single)
AC_SUBST(REALOPT)
AC_SUBST(FREALOPT)
if test "$with_precision" = "double" ; then
  REALOPT=-D_REAL_IS_DOUBLE_
  if test "$with_comp" = "INTEL" ; then
    FREALOPT="-r8"
  elif test "$with_comp" = "FJ" ; then
    FREALOPT="-CcdRR8"
  fi
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
esac


#
# MPI environment
#
AC_SUBST(MPI_DIR)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_FCLIBS)



#
# if wrapper compiler is specified, skip this part
#

if test "$enable_mpi" = "none" ; then

# MPICH
  AC_ARG_WITH(mpich, [AC_HELP_STRING([--with-mpich=dir],[Specify MPICH install directory])], , with_mpich=none)

  if test "$with_mpich" != "none" ; then
    MPI_DIR=$with_mpich;
    MPI_CFLAGS=-I$MPI_DIR/include
    MPI_LDFLAGS=-L$MPI_DIR/lib
    enable_mpi="mpich"
  fi

# OpenMPI // mpichが指定されておらず，かつ，ompが指定されていれば変数をセット
  if test "$with_mpich" = "none" ; then
    AC_ARG_WITH(ompi, [AC_HELP_STRING([--with-ompi=dir],[Specify OpenMPI install directory])], , with_ompi=none)

    if test "$with_ompi" != "none" ; then
      MPI_DIR=$with_ompi;
      MPI_CFLAGS=-I$MPI_DIR/include
      MPI_LDFLAGS=-L$MPI_DIR/lib
      enable_mpi="ompi"
    fi
  fi

fi



#
# MPI Library for Fortran
#

if test "$with_comp" = "FJ" ; then
  MPI_LIBS=""
  case "$host" in
    sparc*) MPI_FCLIBS="$MPI_LIBS"" -lmpi_f77 -lmpi_f90" ;;
  esac
elif test "$enable_mpi" = "mpich" ; then
  MPI_LIBS="-lmpich"
  MPI_FCLIBS="$MPI_LIBS"" -lfmpich"
  case "$TGT_DEVICE" in
    "Lion"|"Snow_Leopard"|"Leopard"|"Mountain_Lion"|"Marvericks") ;;
  esac
elif test "$enable_mpi" = "ompi" ; then
  MPI_LIBS="-lmpi -lmpi_cxx"
  MPI_FCLIBS="$MPI_LIBS"" -lmpi_f77 -lmpi_f90"
  case "$TGT_DEVICE" in
    "Lion"|"Snow_Leopard"|"Leopard"|"Mountain_Lion"|"Marvericks") ;;
  esac
fi



# warning. if a wrapper compiler is used, "enable_mpi" is empty, but parallel.
#

if test "$enable_mpi" = "none" ; then
  echo "  Warning: Neither MPICH nor OpenMPI library is not specified."
  echo "           If compilation error occurs about MPI, specify --with-mpich or --with-ompi option."
fi


#
# FFVC special flags
#


FFVC_FC="$FC"
AC_SUBST(FFVC_FC)

FFVC_LIBS=""
AC_SUBST(FFVC_LIBS)

FFVC_LDFLAGS=""
AC_SUBST(FFVC_LDFLAGS)

FFVC_FC_LD="$FC"
AC_SUBST(FFVC_FC_LD)

if test "$with_comp" = "FJ"; then
  FFVC_FC_LD="$CXX"" --linkfortran"
  FFVC_LIBS=" --linkfortran"
elif test "$with_comp" = "INTEL" ; then
  FFVC_LIBS="$FFVC_LIBS"" -lstdc++ -lifcore"
else
  FFVC_LIBS="$FFVC_LIBS"" -lstdc++"
fi





#
# CPMlib
#
AC_ARG_WITH(cpm, [AC_HELP_STRING([--with-cpm=dir],[Specify CPMlib installed directory])], , with_cpm=none)
AC_SUBST(CPM_DIR)
AC_SUBST(CPM_CFLAGS)
AC_SUBST(CPM_LDFLAGS)
if test "$with_cpm" != "none" ; then
  CPM_DIR=$with_cpm;
  CPM_CFLAGS=`$CPM_DIR/bin/cpm-config --cflags`
  CPM_LDFLAGS=`$CPM_DIR/bin/cpm-config --libs`
else
  echo "  Error: CPM library is not specified."
  echo "         Specify --with-cpm option."
  exit 1
fi



#
# CIOlib
#
AC_ARG_WITH(cio, [AC_HELP_STRING([--with-cio=dir],[Specify CIOlib installed directory])], , with_cio=none)
AC_SUBST(CIO_DIR)
AC_SUBST(CIO_CFLAGS)
AC_SUBST(CIO_LDFLAGS)
if test "$with_cio" != "none" ; then
  CIO_DIR=$with_cio;
  CIO_CFLAGS=`$CIO_DIR/bin/cio-config --cflags`
  CIO_LDFLAGS=`$CIO_DIR/bin/cio-config --libs`
else
  echo "  Error: CIO library is not specified."
  echo "         Specify --with-cio option."
  exit 1
fi




#
# TextParser
#
AC_ARG_WITH(parser, [AC_HELP_STRING([--with-parser=dir],[Specify TextParser install directory])], , with_parser=none)
AC_SUBST(TP_DIR)
AC_SUBST(TP_CFLAGS)
AC_SUBST(TP_LDFLAGS)
if test "$with_parser" != "none" ; then
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

if test "$with_pm" != "none" ; then
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

if test "$with_polylib" != "none" ; then
  PL_DIR=$with_polylib;
  PL_CFLAGS=`$PL_DIR/bin/polylib-config --cflags`
  PL_LDFLAGS=`$PL_DIR/bin/polylib-config --libs`
else
  echo "  Error: Polylib is not specified."
  echo "         Specify --with-polylib option."
  exit 1
fi



#
# Cutlib
#
AC_ARG_WITH(cut, [AC_HELP_STRING([--with-cut=dir],[Specify Cutlib install directory])], , with_cut=none)
AC_SUBST(CUT_DIR)
AC_SUBST(CUT_CFLAGS)
AC_SUBST(CUT_LDFLAGS)

if test "$with_cut" != "none" ; then
  CUT_DIR=$with_cut;
  CUT_CFLAGS=`$CUT_DIR/bin/cut-config --cflags`
  CUT_LDFLAGS=`$CUT_DIR/bin/cut-config --libs`
else
  echo "  Error: Cutlib is not specified."
  echo "         Specify --with-cut option."
  exit 1
fi




#
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src/IP src/FB src/F_CORE src/FFV src/Util_Combsph doc"



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
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_MODE_T
AC_TYPE_SIZE_T


#
# Checks for library functions.
#
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset mkdir pow realpath sqrt strcasecmp strchr strdup strerror va_start])



AC_CONFIG_FILES([Makefile \
                 doc/Makefile \
                 src/F_CORE/Makefile \
                 src/FB/Makefile \
                 src/IP/Makefile \
                 src/FFV/Makefile \
                 src/Util_Combsph/Makefile \
                 ffvc-config \
                 src/FFV/ffv_Version.h \
                ])

AC_OUTPUT

chmod +x ./ffvc-config

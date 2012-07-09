###
###########################################
# CPM - Cartesian Partition Manager       #
#                                         #
# Copyright (C) 2012 University of Tokyo. #
###########################################
###

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(CPMlib,1.0.4,,CPMlib)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/cpm_ParaManager.cpp])
AM_CONFIG_HEADER([config.h])

#
# default prefix setting
#
AC_PREFIX_DEFAULT([/usr/local/CPMlib])
if test "$prefix" = "NONE" ; then
  CPM_INST_DIR=/usr/local/CPMlib
else
  CPM_INST_DIR="$prefix"
fi

#
# platform
#
CPM_DEVICE="$($SHELL ./cpm-uname)"
AC_SUBST(CPM_DEVICE)

#
# compiler libs
#
AC_ARG_WITH(comp, [  --with-comp=(INTEL|FJ)      Specify Compiler.], , with_comp=none)
AC_SUBST(COMPILER_LIBS)
case "$with_comp" in
"INTEL")
  COMPILER_LIBS="-lifport -lifcore"
  ;;
"FJ")
  COMPILER_LIBS="--linkfortran"
  ;;
esac

enable_mpi="none"

#
# MPICH environment
#
AC_SUBST(MPICH_DIR)
AC_SUBST(MPICH_CFLAGS)
AC_SUBST(MPICH_LDFLAGS)
AC_SUBST(MPICH_LIBS)
AC_SUBST(MPICH_FCLIBS)
AC_ARG_WITH(mpich, [  --with-mpich=dir            Specify MPICH install directory.], , with_mpich=none)
if test "$with_mpich" != "none" ; then
  MPICH_DIR=$with_mpich;
  MPICH_CFLAGS=-I$MPICH_DIR/include
  MPICH_LDFLAGS=-L$MPICH_DIR/lib
  enable_mpi="mpich"
fi

#
# OpenMPI environment
#
if test "$with_mpich" = "none" ; then
  AC_ARG_WITH( ompi, [  --with-ompi=dir             Specify OpenMPI install directory.], , with_ompi=none)
  if test "$with_ompi" != "none" ; then
    MPICH_DIR=$with_ompi;
    MPICH_CFLAGS=-I$MPICH_DIR/include
    MPICH_LDFLAGS=-L$MPICH_DIR/lib
    enable_mpi="ompi"
  fi
fi

#
# TextParser environment
#
AC_ARG_WITH(parser, [  --with-parser=dir           Specify TextParser install directory.], , with_parser=none)
AC_SUBST(TP_DIR)
AC_SUBST(TP_CFLAGS)
AC_SUBST(TP_LDFLAGS)
AC_SUBST(TP_LIBS)
if test "$with_parser" != "none" ; then
  TP_DIR=$with_parser;
  TP_CFLAGS=-I$TP_DIR/include
  TP_LDFLAGS=-L$TP_DIR/lib
#  TP_LIBS=-lTextParser
  TP_LIBS=$TP_DIR/lib/libTextParser.a
else
  echo "  Error: TextParser library not specified."
  echo "         set --with-parser option"
  exit 1
fi

#
# REAL environment
#
AC_ARG_WITH(real, [  --with-real=(float|double)  Specify REAL_TYPE size. [float]], , with_real=float)
AC_SUBST(REALOPT)
AC_SUBST(FREALOPT)
if test "$with_real" = "double" ; then
  REALOPT=-D_REAL_IS_DOUBLE_
  if test "$with_comp" = "INTEL" ; then
    FREALOPT="-r8"
  elif test "$with_comp" = "FJ" ; then
    FREALOPT="-CcdRR8"
  fi
else
  REALOPT=
fi

#
# MPI Library
#
if test "$with_comp" = "FJ" ; then
  MPICH_LIBS=""
  case "$host" in
  sparc*)
    MPICH_FCLIBS="$MPICH_LIBS"" -lmpi_f77 -lmpi_f90"
    ;;
  esac
elif test "$enable_mpi" = "mpich" ; then
  MPICH_LIBS="-lmpich"
  MPICH_FCLIBS="$MPICH_LIBS"" -lfmpich"
  case "$CPM_DEVICE" in
  "Lion"|"Snow_Leopard"|"Leopard"|"Tiger_Intel"|"Tiger_PPC"|"Panther")
#    MPICH_LIBS="$MPICH_LIBS"" -lpmpich"
    ;;
  esac
elif test "$enable_mpi" = "ompi" ; then
  MPICH_LIBS="-lmpi"
  MPICH_FCLIBS="$MPICH_LIBS"" -lmpi_f77 -lmpi_f90"
  case "$CPM_DEVICE" in
  "Lion"|"Snow_Leopard"|"Leopard"|"Tiger_Intel"|"Tiger_PPC"|"Panther")
#    MPICH_LIBS="$MPICH_LIBS"" -lmpi_cxx"
    ;;
  esac
fi


#
# f90 example
#
AC_ARG_WITH(f90example, [  --with-f90example=(yes|no)  make f90 example.], , with_f90example=yes)


#
# f90 compiler check
#
if test "_$FC" = "_"; then
  if test "_$F90" = "_"; then
    with_f90example="no"
  else
    FC=$F90
  fi
fi
if test "_$FCFLAGS" = "_"; then
  if test "_$F90FLAGS" != "_"; then
    FCFLAGS=$F90FLAGS
  fi
fi

#
# CPM special flags
#
CPM_CFLAGS="-I$CPM_INST_DIR/include"
AC_SUBST(CPM_CFLAGS)
CPM_LDFLAGS="-L$CPM_INST_DIR/lib"
AC_SUBST(CPM_LDFLAGS)
CPM_LIBS="-lcpmlib"
AC_SUBST(CPM_LIBS)
CPM_CXX="$CXX"
AC_SUBST(CPM_CXX)
CPM_FC="$FC"
AC_SUBST(CPM_FC)
CPM_FC_LD="$FC"
if test "$with_comp" = "FJ"; then
  CPM_FC_LD="$CXX"" --linkfortran"
else
  CPM_LIBS="$CPM_LIBS"" -lstdc++"
fi
AC_SUBST(CPM_FC_LD)


#
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src doc Examples/cxx"
if test "$with_f90example" = "yes" ; then
  MAKE_SUB_DIRS="$MAKE_SUB_DIRS"" Examples/f90"
fi

#
# Checks for programs.
#
AC_PROG_CXX
AC_PROG_FC
AC_PROG_RANLIB

#
# Checks for header files.
#
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h \
                  sys/time.h\
                 ])

#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME

#
# Checks for library functions.
#
AC_CHECK_FUNCS([gettimeofday memset])

AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 doc/Makefile \
                 cpm-config \
                 include/cpm_Version.h \
                 Examples/cxx/Makefile \
                 Examples/f90/Makefile \
                ])

AC_OUTPUT

chmod +x ./cpm-config


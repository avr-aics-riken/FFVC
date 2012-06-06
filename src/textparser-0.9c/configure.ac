#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([TextParser],[0.9c_devel2],[tba@tba.tba])

AC_CONFIG_SRCDIR([include/TextParser.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
#AC_CONFIG_MACRO_DIR([m4 libltdl/m4])
#AC_CONFIG_AUX_DIR([libltdl/config])

AM_INIT_AUTOMAKE([foreign])

dnl mpi option
AC_MSG_CHECKING(whether to enable mpi feature)
AC_ARG_ENABLE(mpi,AS_HELP_STRING([--enable-mpi],
	[turn on mpi feature [default no]]),
	[enable_mpi="$enableval"],[enable_mpi="no"])
if test $enable_mpi != "no"; then
   AC_MSG_RESULT(yes)
   AC_DEFINE(ENABLE_MPI, ,
   Enable mpi. Please set CXX to mpi enabled c++ compiler.)

#   AC_SUBST(enablempi,[yes])
else
   AC_MSG_RESULT(no)
#   AC_SUBST(enablempi,[no])
fi
   AM_CONDITIONAL(MPI,test $enable_mpi != "no")	

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
AC_PROG_FC
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL
AM_PROG_LIBTOOL

#AC_LIBLTDL_CONVENIENCE
dnl AC_LIBLTDL_INSTALLABLE
#AC_LIBTOOL_DLOPEN

#AC_CONFIG_SUBDIRS(libltdl)

#LT_CONFIG_LTDL_DIR([libltdl])
#AC_LIB_LTDL
#AC_SUBST(INCLTDL)

#AC_SUBST(LIBLTDL)




# Checks for libraries.
# FIXME: Replace `main' with a function in `-ltextparser':
AC_CHECK_LIB([TextParser], [TextParser::get_instance()])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.

AC_CONFIG_FILES([Examples/Makefile
                 Makefile
                 include/Makefile
                 src/Makefile])
AC_OUTPUT

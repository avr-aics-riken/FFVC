###
###########################################
# CPM - Cartesian Partition Manager       #
#                                         #
# Copyright (C) 2012 University of Tokyo. #
###########################################
###
noinst_PROGRAMS=mconvp mconvp_CPM

DISTCLEANFILES=*~
CLEANFILES=*.mod data/*.sph

mconvp_SOURCES = \
    mconvp.f90

mconvp_CPM_SOURCES = \
    mconvp_CPM.f90 \
    README \
    ./data/in.ftcs \
    ./data/in.jacobi \
    ./data/run_ftcs.sh \
    ./data/run_jacobi.sh


mconvp_CPM_FCFLAGS = \
    -I$(top_builddir)/include \
    @MPICH_CFLAGS@ \
    @FREALOPT@

mconvp_CPM_LINK = @CPM_FC_LD@ $(FCFLAGS) \
  $(AM_LDFLAGS) $(LDFLAGS) -o $@

mconvp_CPM_LDADD = \
    -L$(top_builddir)/src \
    @CPM_LIBS@ \
    @MPICH_LDFLAGS@ \
    @MPICH_FCLIBS@ \
    @TP_LDFLAGS@ \
    @TP_LIBS@

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)


###
###########################################
# CPM - Cartesian Partition Manager       #
#                                         #
# Copyright (C) 2012 University of Tokyo. #
###########################################
###
noinst_PROGRAMS=exampleCXX

DISTCLEANFILES=*~
CLEANFILES=data/*.log

exampleCXX_SOURCES = \
    main.C \
    ./include/_debug.h \
    ./include/_mpi_test.h \
    README \
    ./data/input1.txt \
    ./data/input2.txt \
    ./data/run.sh

exampleCXX_CXXFLAGS = \
    -I$(top_builddir)/include \
    @MPICH_CFLAGS@ \
    @TP_CFLAGS@ \
    @REALOPT@

exampleCXX_LDADD = \
    -L$(top_builddir)/src \
    @CPM_LIBS@ \
    @MPICH_LDFLAGS@ \
    @MPICH_LIBS@ \
    @TP_LDFLAGS@ \
    @TP_LIBS@

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)


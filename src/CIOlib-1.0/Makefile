###################################################################
#
# CIOlib - Cartesian Input / Output library
#
# Copyright (c) AICS, RIKEN. All right reserved. 2013
#
# Advanced Institute for Computational Science (AICS), RIKEN
#
###################################################################
#
#  At first, edit a make_setting file, then make
#  If CIOlib is placed inside the directory of FFVC solver, this make file is ignored.


include make_setting

all:
	( \
	cd src; \
	make \
		CXX='$(CXX)' \
    CXXFLAGS='$(CXXFLAGS)' \
    F90='$(F90)' \
    F90FLAGS='$(F90FLAGS)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
	)

clean:
	(cd src; make clean)

depend:
	(cd src; make depend)

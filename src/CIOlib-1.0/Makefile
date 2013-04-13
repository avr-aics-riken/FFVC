###################################################################
#
# CIOlib - Cartesian Input / Output library
#
# Copyright (c) AICS, RIKEN. All right reserved. 2013
#
###################################################################
#
#  At first, editing a make_setting file, then make
#  If CIOlib is placed inside the directory of FFVC solver, this make file is ignored.
#

include make_setting

all:
	( \
	cd src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		MPI_DIR='$(MPI_DIR)' \
	)

clean:
	(cd src; make clean)

allclean: clean
	(cd lib; $(RM) $(TARGET) )

depend:
	(cd src; make depend)

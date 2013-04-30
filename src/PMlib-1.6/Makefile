###################################################################
#
# PMlib - Performance Monitor library
# 
# Copyright (c) 2010-2013 Advanced Institute for Computational Science, RIKEN.
# All right reserved.
#
###################################################################
#
#  At first, editing a make_setting file, then make
#  If PMlib is placed inside the directory of FFVC solver, this make file is ignored.


include make_setting

all:
	( \
	cd src; \
	make \
		CXX='$(CXX)' \
    CXXFLAGS='$(CXXFLAGS)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
	)

clean:
	(cd src; make clean)

depend:
	(cd src; make depend)

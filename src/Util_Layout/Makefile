# #################################################################
#
# output layoutfiles
#
# Copyright (c) All right reserved. 2012-2013
#
# Institute of Industrial Science, University of Tokyo, Japan.
#
# #################################################################

include ../FFV/make_setting

all:

	( \
	cd ../FB; \
	make \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		CC='$(CC)' \
		CFLAGS='$(CFLAGS)' \
		FC='$(FC)' \
		FCFLAGS='$(FCFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
	)

	( \
	make -f MakeLAYOUT \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		LDFLAGS='$(LDFLAGS)' \
		LIBS='$(LIBS)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
		UDEF_LIB_PATH='$(UDEF_LIB_PATH)' \
		UDEF_LIBS='$(UDEF_LIBS_UTIL)' \
		RM='$(RM)' \
	)

clean: 
	(make -f MakeLAYOUT clean)
  
allclean: clean
	(cd ../FB; make clean)

depend:
	(cd ../FB; make depend)
	(make -f MakeLAYOUT depend)



# DO NOT DELETE

main.o: layout.h ../FB/FB_Define.h ../FB/mydebug.h ../FB/TPControl.h
main.o: ../FB/dfiinfo.h ../FB/DomainInfo.h LAYOUT_Define.h
layout.o: layout.h ../FB/FB_Define.h ../FB/mydebug.h ../FB/TPControl.h
layout.o: ../FB/dfiinfo.h ../FB/DomainInfo.h LAYOUT_Define.h

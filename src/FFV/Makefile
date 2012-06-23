###################################################################
#
# CAERU Library
#
# Copyright (c) All right reserved. 2012
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

BASEDIR  = $(PWD)
TARGET   = ffv

include make_setting

all:
	( \
	cd ../IP; \
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
	make -f MakeFFV \
		TARGET='$(TARGET)' \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		LDFLAGS='$(LDFLAGS)' \
		LIBS='$(LIBS)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
		UDEF_LIB_PATH='$(UDEF_LIB_PATH)' \
		UDEF_LIBS='$(UDEF_LIBS)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
	)

clean:
	(cd ../IP; make clean)
	(cd ../FB; make clean)
	(make -f MakeFFV clean)

allclean: clean
	(cd ../../bin; $(RM) $(TARGET) )

depend:
	(cd ../IP; make depend)
	(cd ../FB; make depend)
	(make -f MakeFFV depend)

# DO NOT DELETE

ffv.o: ffv.h ../FB/FB_Define.h ../FB/mydebug.h ffv_Define.h ../FB/FBUtility.h
ffv.o: ../FB/FB_Define.h ../FB/Control.h ../FB/Medium.h ../FB/Component.h
ffv.o: ../FB/BndOuter.h ../FB/Interval_Mngr.h ../FB/TPControl.h
ffv_Initialize.o: ffv.h ../FB/FB_Define.h ../FB/mydebug.h ffv_Define.h
ffv_Initialize.o: ../FB/FBUtility.h ../FB/FB_Define.h ../FB/Control.h
ffv_Initialize.o: ../FB/Medium.h ../FB/Component.h ../FB/BndOuter.h
ffv_Initialize.o: ../FB/Interval_Mngr.h ../FB/TPControl.h
main.o: ffv.h ../FB/FB_Define.h ../FB/mydebug.h ffv_Define.h
main.o: ../FB/FBUtility.h ../FB/FB_Define.h ../FB/Control.h ../FB/Medium.h
main.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
main.o: ../FB/TPControl.h

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
ffv.o: ../FB/Control.h ../FB/Medium.h ../FB/Component.h ../FB/BndOuter.h
ffv.o: ../FB/Interval_Mngr.h ../FB/TPControl.h ../FB/Alloc.h ../FB/FileIO.h
ffv.o: ../FB/FB_Ffunc.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h ../IP/IP_Duct.h
ffv.o: ../IP/IP_Define.h ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
ffv.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h ../IP/IP_Polygon.h
ffv.o: ../IP/IP_Sphere.h
ffv_Initialize.o: ffv.h ../FB/FB_Define.h ../FB/mydebug.h ffv_Define.h
ffv_Initialize.o: ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Initialize.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Initialize.o: ../FB/TPControl.h ../FB/Alloc.h ../FB/FileIO.h
ffv_Initialize.o: ../FB/FB_Ffunc.h ../FB/ParseBC.h ../FB/vec3.h
ffv_Initialize.o: ../FB/vec3_func.h ../FB/vec3f_func.h ../FB/Intrinsic.h
ffv_Initialize.o: ../FB/ParseMat.h ../IP/IP_Duct.h ../IP/IP_Define.h
ffv_Initialize.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
ffv_Initialize.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
ffv_Initialize.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_SetBC.o: ffv_SetBC.h ../FB/SetBC.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_SetBC.o: ../FB/FBUtility.h ../FB/BndOuter.h ../FB/Control.h
ffv_SetBC.o: ../FB/Medium.h ../FB/Component.h ../FB/Interval_Mngr.h
ffv_SetBC.o: ../FB/TPControl.h ../FB/Intrinsic.h ../FB/vec3.h
ffv_SetBC.o: ../FB/vec3_func.h ../FB/vec3f_func.h ffv_FortranFunc.h
ffv_SetBC.o: ffv_Define.h
main.o: ffv.h ../FB/FB_Define.h ../FB/mydebug.h ffv_Define.h
main.o: ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h ../FB/Component.h
main.o: ../FB/BndOuter.h ../FB/Interval_Mngr.h ../FB/TPControl.h
main.o: ../FB/Alloc.h ../FB/FileIO.h ../FB/FB_Ffunc.h ../FB/ParseBC.h
main.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h ../FB/Intrinsic.h
main.o: ../FB/ParseMat.h ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h
main.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
main.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h

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
	cd ../F_CORE; \
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
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
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
	(cd ../F_CORE; make clean)
	(make -f MakeFFV clean)

allclean: clean
	(cd ../../bin; $(RM) $(TARGET) )

depend:
	(cd ../IP; make depend)
	(cd ../FB; make depend)
	(cd ../F_CORE; make depend)
	(make -f MakeFFV depend)

# DO NOT DELETE

ffv.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv.o: ../FB/TPControl.h ../FB/Alloc.h ../FB/FileIO.h ../FB/FB_Ffunc.h
ffv.o: ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h
ffv.o: ../FB/Intrinsic.h ../FB/ParseMat.h ../FB/VoxInfo.h ../FB/SetBC.h
ffv.o: ffv_SetBC.h ffv_Ffunc.h ../IP/IP_Duct.h ../IP/IP_Define.h
ffv.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h
ffv.o: ../IP/IP_Step.h ../IP/IP_Cylinder.h ../IP/IP_Polygon.h
ffv.o: ../IP/IP_Sphere.h
ffv_Alloc.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Alloc.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Alloc.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Alloc.o: ../FB/TPControl.h ../FB/Alloc.h ../FB/FileIO.h ../FB/FB_Ffunc.h
ffv_Alloc.o: ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Alloc.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Alloc.o: ../FB/VoxInfo.h ../FB/SetBC.h ffv_SetBC.h ffv_Ffunc.h
ffv_Alloc.o: ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h
ffv_Alloc.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Alloc.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Initialize.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Initialize.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h
ffv_Initialize.o: ../FB/Medium.h ../FB/Component.h ../FB/BndOuter.h
ffv_Initialize.o: ../FB/Interval_Mngr.h ../FB/TPControl.h ../FB/Alloc.h
ffv_Initialize.o: ../FB/FileIO.h ../FB/FB_Ffunc.h ../FB/ParseBC.h
ffv_Initialize.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h
ffv_Initialize.o: ../FB/Intrinsic.h ../FB/ParseMat.h ../FB/VoxInfo.h
ffv_Initialize.o: ../FB/SetBC.h ffv_SetBC.h ffv_Ffunc.h ../IP/IP_Duct.h
ffv_Initialize.o: ../IP/IP_Define.h ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h
ffv_Initialize.o: ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Initialize.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_SetBC.o: ffv_SetBC.h ../FB/SetBC.h ../FB/DomainInfo.h ../FB/FB_Define.h
ffv_SetBC.o: ../FB/mydebug.h ../FB/FBUtility.h ../FB/BndOuter.h
ffv_SetBC.o: ../FB/Control.h ../FB/Medium.h ../FB/Component.h
ffv_SetBC.o: ../FB/Interval_Mngr.h ../FB/TPControl.h ../FB/Intrinsic.h
ffv_SetBC.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h ffv_Define.h
ffv_SetBC.o: ffv_Ffunc.h
main.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
main.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
main.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
main.o: ../FB/TPControl.h ../FB/Alloc.h ../FB/FileIO.h ../FB/FB_Ffunc.h
main.o: ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h
main.o: ../FB/Intrinsic.h ../FB/ParseMat.h ../FB/VoxInfo.h ../FB/SetBC.h
main.o: ffv_SetBC.h ffv_Ffunc.h ../IP/IP_Duct.h ../IP/IP_Define.h
main.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h
main.o: ../IP/IP_Step.h ../IP/IP_Cylinder.h ../IP/IP_Polygon.h
main.o: ../IP/IP_Sphere.h

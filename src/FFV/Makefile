###################################################################
#
# FFV : Frontflow / violet Cartesian
#
# Copyright (c) 2012 All right reserved.
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

BASEDIR  = $(PWD)
TARGET   = ffvc

include make_setting

all:
	( \
	cd ../Polylib-2.1.1/src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		MPI_DIR='$(MPI_DIR)' \
		TP_DIR='$(TP_DIR)' \
	)

	( \
	cd ../Cutlib-2.0.4/src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		TP_DIR='$(TP_DIR)' \
	)

	( \
	cd ../PMlib-1.4/src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		MPI_DIR='$(MPI_DIR)' \
	)

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
		CPM_DIR='$(CPM_DIR)' \
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
		UDEF_LIB_PATH_SPEC='$(UDEF_LIB_PATH_SPEC)' \
	)

clean: 
	(cd ../IP; make clean)
	(cd ../FB; make clean)
	(cd ../F_CORE; make clean)
	(make -f MakeFFV clean)
  
mostlyclean: clean
	(cd ../Polylib-2.1.1/src; make clean)
	(cd ../Cutlib-2.0.4/src; make clean)
	(cd ../PMlib-1.4/src; make clean)

allclean: mostlyclean
	(cd ../../bin; $(RM) $(TARGET) )

depend:
	(cd ../Polylib-2.1.1/src; make depend)
	(cd ../Cutlib-2.0.4/src; make depend)
	(cd ../PMlib-1.4/src; make depend)
	(cd ../IP; make depend)
	(cd ../FB; make depend)
	(cd ../F_CORE; make depend)
	(make -f MakeFFV depend)

# DO NOT DELETE

main.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
main.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
main.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
main.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h ../FB/FB_Ffunc.h
main.o: ../FB/PLOT3D_write.h ../FB/Alloc.h ../FB/FileIO.h ../FB/ParseBC.h
main.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h ../FB/Intrinsic.h
main.o: ../FB/ParseMat.h ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h
main.o: ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h
main.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
main.o: ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h
main.o: ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
main.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h ../FB/FB_Ffunc.h
ffv.o: ../FB/PLOT3D_write.h ../FB/Alloc.h ../FB/FileIO.h ../FB/ParseBC.h
ffv.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h ../FB/Intrinsic.h
ffv.o: ../FB/ParseMat.h ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h
ffv.o: ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h
ffv.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv.o: ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h
ffv.o: ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
ffv.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Initialize.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Initialize.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h
ffv_Initialize.o: ../FB/Medium.h ../FB/Component.h ../FB/BndOuter.h
ffv_Initialize.o: ../FB/Interval_Mngr.h ../FB/TPControl.h ../FB/PLOT3D_read.h
ffv_Initialize.o: ../FB/PLOT3D.h ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h
ffv_Initialize.o: ../FB/Alloc.h ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h
ffv_Initialize.o: ../FB/vec3_func.h ../FB/vec3f_func.h ../FB/Intrinsic.h
ffv_Initialize.o: ../FB/ParseMat.h ../FB/VoxInfo.h ../FB/SetBC.h
ffv_Initialize.o: ../FB/CompoFraction.h ffv_SetBC.h ../F_CORE/ffv_Ffunc.h
ffv_Initialize.o: ../FB/dfi.h ../FB/History.h ../FB/Monitor.h
ffv_Initialize.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv_Initialize.o: ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h
ffv_Initialize.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h
ffv_Initialize.o: ../IP/IP_Step.h ../IP/IP_Cylinder.h ../IP/IP_Polygon.h
ffv_Initialize.o: ../IP/IP_Sphere.h
ffv_Alloc.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Alloc.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Alloc.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Alloc.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
ffv_Alloc.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h
ffv_Alloc.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Alloc.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Alloc.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Alloc.o: ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h
ffv_Alloc.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
ffv_Alloc.o: ../FB/Sampling.h ../IP/IP_Duct.h ../IP/IP_Define.h
ffv_Alloc.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
ffv_Alloc.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
ffv_Alloc.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_SetBC.o: ffv_SetBC.h ../FB/SetBC.h ../FB/DomainInfo.h ../FB/FB_Define.h
ffv_SetBC.o: ../FB/mydebug.h ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h
ffv_SetBC.o: ../FB/Component.h ../FB/FBUtility.h ../FB/Interval_Mngr.h
ffv_SetBC.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
ffv_SetBC.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Intrinsic.h
ffv_SetBC.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h ffv_Define.h
ffv_SetBC.o: ../F_CORE/ffv_Ffunc.h
ffv_Restart.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Restart.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Restart.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Restart.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
ffv_Restart.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h
ffv_Restart.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Restart.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Restart.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h
ffv_Restart.o: ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h
ffv_Restart.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
ffv_Restart.o: ../FB/Sampling.h ../IP/IP_Duct.h ../IP/IP_Define.h
ffv_Restart.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
ffv_Restart.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
ffv_Restart.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Loop.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Loop.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Loop.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Loop.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
ffv_Loop.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h
ffv_Loop.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Loop.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Loop.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Loop.o: ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h ../FB/Monitor.h
ffv_Loop.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv_Loop.o: ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h
ffv_Loop.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Loop.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_LS.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_LS.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_LS.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_LS.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
ffv_LS.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h ../FB/FileIO.h
ffv_LS.o: ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h
ffv_LS.o: ../FB/Intrinsic.h ../FB/ParseMat.h ../FB/VoxInfo.h ../FB/SetBC.h
ffv_LS.o: ../FB/CompoFraction.h ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../FB/dfi.h
ffv_LS.o: ../FB/History.h ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
ffv_LS.o: ../FB/Sampling.h ../IP/IP_Duct.h ../IP/IP_Define.h
ffv_LS.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h
ffv_LS.o: ../IP/IP_Step.h ../IP/IP_Cylinder.h ../IP/IP_Polygon.h
ffv_LS.o: ../IP/IP_Sphere.h
ffv_Post.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Post.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Post.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Post.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
ffv_Post.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h
ffv_Post.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Post.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Post.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Post.o: ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h ../FB/Monitor.h
ffv_Post.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv_Post.o: ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h
ffv_Post.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Post.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Plot3d.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Plot3d.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Plot3d.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Plot3d.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
ffv_Plot3d.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h
ffv_Plot3d.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Plot3d.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Plot3d.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Plot3d.o: ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h
ffv_Plot3d.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
ffv_Plot3d.o: ../FB/Sampling.h ../IP/IP_Duct.h ../IP/IP_Define.h
ffv_Plot3d.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
ffv_Plot3d.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
ffv_Plot3d.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Heat.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Heat.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Heat.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Heat.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
ffv_Heat.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h
ffv_Heat.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Heat.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Heat.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Heat.o: ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h ../FB/Monitor.h
ffv_Heat.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv_Heat.o: ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h
ffv_Heat.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Heat.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Heat.o: ../FB/FindexS3D.h
NS_FS_E_Binary.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
NS_FS_E_Binary.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h
NS_FS_E_Binary.o: ../FB/Medium.h ../FB/Component.h ../FB/BndOuter.h
NS_FS_E_Binary.o: ../FB/Interval_Mngr.h ../FB/TPControl.h ../FB/PLOT3D_read.h
NS_FS_E_Binary.o: ../FB/PLOT3D.h ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h
NS_FS_E_Binary.o: ../FB/Alloc.h ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h
NS_FS_E_Binary.o: ../FB/vec3_func.h ../FB/vec3f_func.h ../FB/Intrinsic.h
NS_FS_E_Binary.o: ../FB/ParseMat.h ../FB/VoxInfo.h ../FB/SetBC.h
NS_FS_E_Binary.o: ../FB/CompoFraction.h ffv_SetBC.h ../F_CORE/ffv_Ffunc.h
NS_FS_E_Binary.o: ../FB/dfi.h ../FB/History.h ../FB/Monitor.h
NS_FS_E_Binary.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
NS_FS_E_Binary.o: ../IP/IP_Duct.h ../IP/IP_Define.h ../IP/IP_PPLT2D.h
NS_FS_E_Binary.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h
NS_FS_E_Binary.o: ../IP/IP_Step.h ../IP/IP_Cylinder.h ../IP/IP_Polygon.h
NS_FS_E_Binary.o: ../IP/IP_Sphere.h
NS_FS_E_CDS.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
NS_FS_E_CDS.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
NS_FS_E_CDS.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
NS_FS_E_CDS.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
NS_FS_E_CDS.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h
NS_FS_E_CDS.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
NS_FS_E_CDS.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
NS_FS_E_CDS.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h
NS_FS_E_CDS.o: ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h
NS_FS_E_CDS.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
NS_FS_E_CDS.o: ../FB/Sampling.h ../IP/IP_Duct.h ../IP/IP_Define.h
NS_FS_E_CDS.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
NS_FS_E_CDS.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
NS_FS_E_CDS.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h
PS_Binary.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
PS_Binary.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
PS_Binary.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
PS_Binary.o: ../FB/TPControl.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
PS_Binary.o: ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h ../FB/Alloc.h
PS_Binary.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
PS_Binary.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
PS_Binary.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
PS_Binary.o: ../F_CORE/ffv_Ffunc.h ../FB/dfi.h ../FB/History.h
PS_Binary.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
PS_Binary.o: ../FB/Sampling.h ../IP/IP_Duct.h ../IP/IP_Define.h
PS_Binary.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
PS_Binary.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
PS_Binary.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h

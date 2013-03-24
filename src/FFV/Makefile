###################################################################
#
# FFV : Frontflow / violet Cartesian
#
# Copyright (c) 2012-2013  All right reserved.
#
# Institute of Industrial Science, University of Tokyo, Japan. 
#
###################################################################

# Default environment
FFVC=intel

### compiler detection

ifeq ($(FFVC),gnu)
include make_setting.gnu
endif

ifeq ($(FFVC),fx)
include make_setting.fx
endif

ifeq ($(FFVC),ibm)
include make_setting.ibm
endif

ifeq ($(FFVC),intel)
include make_setting.intel
endif



#####



UDEF_INC_PATH=-I$(CUTLIB_DIR)/include \
              -I$(POLYLIB_DIR)/include \
              -I$(PMLIB_DIR)/include \
              -I$(TP_DIR)/include \
              -I$(CPM_DIR)/include \
              -I../FB \
              -I../IP \
              -I../F_CORE \
              -I../PLOT3D \
              -I$(MPI_DIR)/include
              
UDEF_LIB_PATH=-L$(POLYLIB_DIR)/lib \
              -L$(CUTLIB_DIR)/lib \
              -L$(PMLIB_DIR)/lib \
              -L$(TP_DIR)/lib \
              -L$(CPM_DIR)/lib \
              -L../FB \
              -L../IP \
              -L../F_CORE \
              -L../PLOT3D \
              -L$(MPI_DIR)/lib
              
UDEF_LIBS = -lmpi -lPLOT3D -lFB -lIP -lFFVcore -lMPIPolylib -lcutlib -lpmlib -lTextParser_mpi -lcpmlib \
            $(UDEF_LIBS_SPEC)

UDEF_LIBS_UTIL = -lmpi -lPLOT3D -lFB -lpmlib -lTextParser_mpi -lcpmlib $(UDEF_LIBS_SPEC)



#####



all:
	( \
	cd $(POLYLIB_DIR)/src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		MPI_DIR='$(MPI_DIR)' \
		TP_DIR='$(TP_DIR)' \
	)

	( \
	cd $(CUTLIB_DIR)/src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		TP_DIR='$(TP_DIR)' \
		POLYLIB_DIR='$(POLYLIB_DIR)' \
	)

	( \
	cd $(PMLIB_DIR)/src; \
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
	cd ../PLOT3D; \
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
  
	( \
	cd ../Util_Combsph; \
	make -f MakeCOMB \
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
		UDEF_LIB_PATH_SPEC='$(UDEF_LIB_PATH_SPEC)' \
	)
  
	( \
	cd ../Util_Layout; \
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
		UDEF_LIB_PATH_SPEC='$(UDEF_LIB_PATH_SPEC)' \
	)
  

clean: 
	(cd ../IP; make clean)
	(cd ../FB; make clean)
	(cd ../F_CORE; make clean)
	(cd ../PLOT3D; make clean)
	(make -f MakeFFV clean)
  
mostlyclean: clean
	(cd $(POLYLIB_DIR)/src; make clean)
	(cd $(CUTLIB_DIR)/src; make clean)
	(cd $(PMLIB_DIR)/src; make clean)

allclean: mostlyclean
	(cd ../Util_Combsph; make -f MakeCOMB clean)
	(cd ../Util_Layout; make -f MakeLAYOUT clean)

depend:
	(cd $(POLYLIB_DIR)/src; make depend)
	(cd $(CUTLIB_DIR)/src; make depend)
	(cd $(PMLIB_DIR)/src; make depend)
	(cd ../IP; make depend)
	(cd ../FB; make depend)
	(cd ../F_CORE; make depend)
	(cd ../PLOT3D; make depend)
	(make -f MakeFFV depend)
	(cd ../Util_Combsph; make -f MakeCOMB depend)
	(cd ../Util_Layout; make -f MakeLAYOUT depend)

# DO NOT DELETE

main.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
main.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
main.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
main.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
main.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
main.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
main.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h ../FB/VoxInfo.h
main.o: ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h ../F_CORE/ffv_Ffunc.h
main.o: ../IP/IP_Jet.h ../IP/IP_Define.h ../FB/dfi.h ../FB/dfiinfo.h
main.o: ../FB/History.h ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
main.o: ../FB/Sampling.h ffv_TerminateCtrl.h ../IP/IP_Duct.h
main.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h
main.o: ../IP/IP_Step.h ../IP/IP_Cylinder.h ../IP/IP_Polygon.h
main.o: ../IP/IP_Sphere.h ../IP/IP_Jet.h ../PLOT3D/ffv_PLOT3D.h
ffv.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
ffv.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h ../FB/FileIO.h
ffv.o: ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h
ffv.o: ../FB/Intrinsic.h ../FB/ParseMat.h ../FB/VoxInfo.h ../FB/SetBC.h
ffv.o: ../FB/CompoFraction.h ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h
ffv.o: ../IP/IP_Define.h ../FB/dfi.h ../FB/dfiinfo.h ../FB/History.h
ffv.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv.o: ffv_TerminateCtrl.h ../IP/IP_Duct.h ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h
ffv.o: ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
ffv.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h ../IP/IP_Jet.h
ffv.o: ../PLOT3D/ffv_PLOT3D.h
ffv_Initialize.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Initialize.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h
ffv_Initialize.o: ../FB/Medium.h ../FB/Component.h ../FB/BndOuter.h
ffv_Initialize.o: ../FB/Interval_Mngr.h ../FB/TPControl.h
ffv_Initialize.o: ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h ../FB/FB_Ffunc.h
ffv_Initialize.o: ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h ../FB/FileIO.h
ffv_Initialize.o: ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Initialize.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Initialize.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h
ffv_Initialize.o: ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h
ffv_Initialize.o: ../IP/IP_Define.h ../FB/dfi.h ../FB/dfiinfo.h
ffv_Initialize.o: ../FB/History.h ../FB/Monitor.h ../FB/MonCompo.h
ffv_Initialize.o: ../FB/basic_func.h ../FB/Sampling.h ffv_TerminateCtrl.h
ffv_Initialize.o: ../IP/IP_Duct.h ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h
ffv_Initialize.o: ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Initialize.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Initialize.o: ../IP/IP_Jet.h ../PLOT3D/ffv_PLOT3D.h
ffv_Alloc.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Alloc.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Alloc.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Alloc.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
ffv_Alloc.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
ffv_Alloc.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Alloc.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Alloc.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Alloc.o: ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h ../IP/IP_Define.h
ffv_Alloc.o: ../FB/dfi.h ../FB/dfiinfo.h ../FB/History.h ../FB/Monitor.h
ffv_Alloc.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv_Alloc.o: ffv_TerminateCtrl.h ../IP/IP_Duct.h ../IP/IP_PPLT2D.h
ffv_Alloc.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Alloc.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Alloc.o: ../IP/IP_Jet.h ../PLOT3D/ffv_PLOT3D.h
ffv_SetBC.o: ffv_SetBC.h ../FB/SetBC.h ../FB/DomainInfo.h ../FB/FB_Define.h
ffv_SetBC.o: ../FB/mydebug.h ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h
ffv_SetBC.o: ../FB/Component.h ../FB/FBUtility.h ../FB/Interval_Mngr.h
ffv_SetBC.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
ffv_SetBC.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Intrinsic.h
ffv_SetBC.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h ffv_Define.h
ffv_SetBC.o: ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h ../IP/IP_Define.h
ffv_Restart.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Restart.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Restart.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Restart.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
ffv_Restart.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
ffv_Restart.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Restart.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Restart.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h
ffv_Restart.o: ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h
ffv_Restart.o: ../IP/IP_Define.h ../FB/dfi.h ../FB/dfiinfo.h ../FB/History.h
ffv_Restart.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
ffv_Restart.o: ../FB/Sampling.h ffv_TerminateCtrl.h ../IP/IP_Duct.h
ffv_Restart.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
ffv_Restart.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
ffv_Restart.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h ../IP/IP_Jet.h
ffv_Restart.o: ../PLOT3D/ffv_PLOT3D.h
ffv_Loop.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Loop.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Loop.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Loop.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
ffv_Loop.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
ffv_Loop.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Loop.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Loop.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Loop.o: ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h ../IP/IP_Define.h
ffv_Loop.o: ../FB/dfi.h ../FB/dfiinfo.h ../FB/History.h ../FB/Monitor.h
ffv_Loop.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv_Loop.o: ffv_TerminateCtrl.h ../IP/IP_Duct.h ../IP/IP_PPLT2D.h
ffv_Loop.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Loop.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Loop.o: ../IP/IP_Jet.h ../PLOT3D/ffv_PLOT3D.h
ffv_LS.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_LS.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_LS.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_LS.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
ffv_LS.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
ffv_LS.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_LS.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_LS.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_LS.o: ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h ../IP/IP_Define.h ../FB/dfi.h
ffv_LS.o: ../FB/dfiinfo.h ../FB/History.h ../FB/Monitor.h ../FB/MonCompo.h
ffv_LS.o: ../FB/basic_func.h ../FB/Sampling.h ffv_TerminateCtrl.h
ffv_LS.o: ../IP/IP_Duct.h ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
ffv_LS.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
ffv_LS.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h ../IP/IP_Jet.h
ffv_LS.o: ../PLOT3D/ffv_PLOT3D.h
ffv_Post.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Post.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Post.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Post.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
ffv_Post.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
ffv_Post.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Post.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Post.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Post.o: ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h ../IP/IP_Define.h
ffv_Post.o: ../FB/dfi.h ../FB/dfiinfo.h ../FB/History.h ../FB/Monitor.h
ffv_Post.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv_Post.o: ffv_TerminateCtrl.h ../IP/IP_Duct.h ../IP/IP_PPLT2D.h
ffv_Post.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Post.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Post.o: ../IP/IP_Jet.h ../PLOT3D/ffv_PLOT3D.h
ffv_Heat.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
ffv_Heat.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
ffv_Heat.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
ffv_Heat.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
ffv_Heat.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
ffv_Heat.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
ffv_Heat.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
ffv_Heat.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
ffv_Heat.o: ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h ../IP/IP_Define.h
ffv_Heat.o: ../FB/dfi.h ../FB/dfiinfo.h ../FB/History.h ../FB/Monitor.h
ffv_Heat.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
ffv_Heat.o: ffv_TerminateCtrl.h ../IP/IP_Duct.h ../IP/IP_PPLT2D.h
ffv_Heat.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
ffv_Heat.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
ffv_Heat.o: ../IP/IP_Jet.h ../PLOT3D/ffv_PLOT3D.h ../FB/FindexS3D.h
ffv_TerminateCtrl.o: ffv_TerminateCtrl.h
NS_FS_E_Binary.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
NS_FS_E_Binary.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h
NS_FS_E_Binary.o: ../FB/Medium.h ../FB/Component.h ../FB/BndOuter.h
NS_FS_E_Binary.o: ../FB/Interval_Mngr.h ../FB/TPControl.h
NS_FS_E_Binary.o: ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h ../FB/FB_Ffunc.h
NS_FS_E_Binary.o: ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h ../FB/FileIO.h
NS_FS_E_Binary.o: ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
NS_FS_E_Binary.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
NS_FS_E_Binary.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h
NS_FS_E_Binary.o: ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h
NS_FS_E_Binary.o: ../IP/IP_Define.h ../FB/dfi.h ../FB/dfiinfo.h
NS_FS_E_Binary.o: ../FB/History.h ../FB/Monitor.h ../FB/MonCompo.h
NS_FS_E_Binary.o: ../FB/basic_func.h ../FB/Sampling.h ffv_TerminateCtrl.h
NS_FS_E_Binary.o: ../IP/IP_Duct.h ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h
NS_FS_E_Binary.o: ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
NS_FS_E_Binary.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
NS_FS_E_Binary.o: ../IP/IP_Jet.h ../PLOT3D/ffv_PLOT3D.h
NS_FS_E_CDS.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
NS_FS_E_CDS.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
NS_FS_E_CDS.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
NS_FS_E_CDS.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
NS_FS_E_CDS.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
NS_FS_E_CDS.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
NS_FS_E_CDS.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
NS_FS_E_CDS.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h
NS_FS_E_CDS.o: ffv_SetBC.h ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h
NS_FS_E_CDS.o: ../IP/IP_Define.h ../FB/dfi.h ../FB/dfiinfo.h ../FB/History.h
NS_FS_E_CDS.o: ../FB/Monitor.h ../FB/MonCompo.h ../FB/basic_func.h
NS_FS_E_CDS.o: ../FB/Sampling.h ffv_TerminateCtrl.h ../IP/IP_Duct.h
NS_FS_E_CDS.o: ../IP/IP_PPLT2D.h ../IP/IP_SHC1D.h ../IP/IP_PMT.h
NS_FS_E_CDS.o: ../IP/IP_Rect.h ../IP/IP_Step.h ../IP/IP_Cylinder.h
NS_FS_E_CDS.o: ../IP/IP_Polygon.h ../IP/IP_Sphere.h ../IP/IP_Jet.h
NS_FS_E_CDS.o: ../PLOT3D/ffv_PLOT3D.h
PS_Binary.o: ffv.h ../FB/DomainInfo.h ../FB/FB_Define.h ../FB/mydebug.h
PS_Binary.o: ffv_Define.h ../FB/FBUtility.h ../FB/Control.h ../FB/Medium.h
PS_Binary.o: ../FB/Component.h ../FB/BndOuter.h ../FB/Interval_Mngr.h
PS_Binary.o: ../FB/TPControl.h ../PLOT3D/PLOT3D_read.h ../PLOT3D/PLOT3D.h
PS_Binary.o: ../FB/FB_Ffunc.h ../PLOT3D/PLOT3D_write.h ../FB/Alloc.h
PS_Binary.o: ../FB/FileIO.h ../FB/ParseBC.h ../FB/vec3.h ../FB/vec3_func.h
PS_Binary.o: ../FB/vec3f_func.h ../FB/Intrinsic.h ../FB/ParseMat.h
PS_Binary.o: ../FB/VoxInfo.h ../FB/SetBC.h ../FB/CompoFraction.h ffv_SetBC.h
PS_Binary.o: ../F_CORE/ffv_Ffunc.h ../IP/IP_Jet.h ../IP/IP_Define.h
PS_Binary.o: ../FB/dfi.h ../FB/dfiinfo.h ../FB/History.h ../FB/Monitor.h
PS_Binary.o: ../FB/MonCompo.h ../FB/basic_func.h ../FB/Sampling.h
PS_Binary.o: ffv_TerminateCtrl.h ../IP/IP_Duct.h ../IP/IP_PPLT2D.h
PS_Binary.o: ../IP/IP_SHC1D.h ../IP/IP_PMT.h ../IP/IP_Rect.h ../IP/IP_Step.h
PS_Binary.o: ../IP/IP_Cylinder.h ../IP/IP_Polygon.h ../IP/IP_Sphere.h
PS_Binary.o: ../IP/IP_Jet.h ../PLOT3D/ffv_PLOT3D.h

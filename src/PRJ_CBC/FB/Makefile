#
# SPHERE - Skeleton for PHysicalSvCube and Engineering REsearch
#
# Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
#

TARGET = libsphnsolvFB.a

CSRCS = \


CXXSRCS = \
       ../../FB/Alloc.C \
       ../../FB/VoxInfo.C \
       ../../FB/BndOuter.C \
       ../../FB/Component.C \
       ../../FB/CompoFraction.C \
       ../../FB/Control.C \
       ../../FB/FBUtility.C \
       ../../FB/FileIO.C \
       ../../FB/History.C \
       ../../FB/Intrinsic.C \
       ../../FB/ParseBC.C \
       ../../FB/ParseMat.C \
       ../../FB/SetBC.C \
       ../../FB/Monitor.C \
       ../../FB/MonCompo.C \
       ../../FB/Interval_Mngr.C \
       ../../FB/Sampling.C \
       ../../FB/DataHolder.C \
       ../../FB/dfi.C \
       ../../FB/TPControl.C \


FCSRCS = \


F90SRCS = \
       ../../FB/FB_util.f90 \


SRCS  = $(CSRCS) $(CXXSRCS) $(FCSRCS) $(F90SRCS)

.SUFFIXES: .o .c .C .f .f90
COBJS = $(CSRCS:.c=.o)
CXXOBJS = $(CXXSRCS:.C=.o)
FCOBJS = $(FCSRCS:.f=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS  = $(COBJS) $(CXXOBJS) $(FCOBJS) $(F90OBJS)

$(TARGET):$(OBJS)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.c.o:
	$(CC) $(CFLAGS) $(XML2FLAGS) $(SPHERE_CFLAGS) $(REALOPT) -o $@ -c $<

.C.o:
	$(CXX) $(CXXFLAGS) $(XML2FLAGS) $(SPHERE_CFLAGS) $(REALOPT) -o $@ -c $<

.f.o:
	$(FC) $(FCFLAGS) -I$(SPHEREDIR)/include -o $@ -c $<

.f90.o:
	$(F90) $(F90FLAGS) -I$(SPHEREDIR)/include -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -I. -I../app -I$(SPHEREDIR)/include -Y $(CXXSRCS) $(CSRCS)

# DO NOT DELETE

../../FB/Alloc.o: ../../FB/Alloc.h ../../FB/FB_Define.h ../../FB/mydebug.h
../../FB/Alloc.o: /usr/local/sphere/include/Skl.h
../../FB/Alloc.o: /usr/local/sphere/include/endianUtil.h
../../FB/Alloc.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/Alloc.o: /usr/local/sphere/include/SklBase.h
../../FB/Alloc.o: /usr/local/sphere/include/SklDefine.h
../../FB/Alloc.o: /usr/local/sphere/include/SklReservWord.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/Alloc.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/Alloc.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/Alloc.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/Alloc.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/Alloc.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/Alloc.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/Alloc.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/Alloc.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/Alloc.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/Alloc.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/Alloc.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/Alloc.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/Alloc.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/Alloc.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/Alloc.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/Alloc.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/Alloc.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/Alloc.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/Alloc.o: ../../FB/FB_Ffunc.h
../../FB/VoxInfo.o: ../../FB/VoxInfo.h /usr/local/sphere/include/Skl.h
../../FB/VoxInfo.o: /usr/local/sphere/include/endianUtil.h
../../FB/VoxInfo.o: ../../FB/FBUtility.h ../../FB/FB_Define.h
../../FB/VoxInfo.o: ../../FB/mydebug.h ../../FB/Component.h ../../FB/Medium.h
../../FB/VoxInfo.o: ../../FB/Parallel_node.h ../../FB/IDtable.h
../../FB/VoxInfo.o: ../../FB/SetBC.h
../../FB/VoxInfo.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/VoxInfo.o: /usr/local/sphere/include/SklBase.h
../../FB/VoxInfo.o: /usr/local/sphere/include/SklDefine.h
../../FB/VoxInfo.o: /usr/local/sphere/include/SklReservWord.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/VoxInfo.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/VoxInfo.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/VoxInfo.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/VoxInfo.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/VoxInfo.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/VoxInfo.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/VoxInfo.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/VoxInfo.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/VoxInfo.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/VoxInfo.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/VoxInfo.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/VoxInfo.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/VoxInfo.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/VoxInfo.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/VoxInfo.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/VoxInfo.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/VoxInfo.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/VoxInfo.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/VoxInfo.o: ../../FB/BndOuter.h ../../FB/Control.h
../../FB/VoxInfo.o: /usr/local/sphere/include/SklUtil.h ../../FB/Monitor.h
../../FB/VoxInfo.o: ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/VoxInfo.o: ../../FB/vec3f_func.h ../../FB/MonCompo.h
../../FB/VoxInfo.o: ../../FB/basic_func.h ../../FB/Sampling.h
../../FB/VoxInfo.o: ../../FB/Interval_Mngr.h ../../FB/TPControl.h
../../FB/VoxInfo.o: ../../FB/Intrinsic.h
../../FB/BndOuter.o: ../../FB/BndOuter.h ../../FB/FB_Define.h
../../FB/BndOuter.o: ../../FB/mydebug.h
../../FB/Component.o: ../../FB/Component.h ../../FB/FB_Define.h
../../FB/Component.o: ../../FB/mydebug.h
../../FB/CompoFraction.o: ../../FB/CompoFraction.h ../../FB/FB_Define.h
../../FB/CompoFraction.o: ../../FB/mydebug.h ../../FB/FBUtility.h
../../FB/CompoFraction.o: ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/CompoFraction.o: ../../FB/vec3f_func.h
../../FB/Control.o: ../../FB/Control.h ../../FB/FB_Define.h
../../FB/Control.o: ../../FB/mydebug.h /usr/local/sphere/include/Skl.h
../../FB/Control.o: /usr/local/sphere/include/endianUtil.h
../../FB/Control.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/Control.o: /usr/local/sphere/include/SklBase.h
../../FB/Control.o: /usr/local/sphere/include/SklDefine.h
../../FB/Control.o: /usr/local/sphere/include/SklReservWord.h
../../FB/Control.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/Control.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/Control.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/Control.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/Control.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/Control.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/Control.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/Control.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/Control.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/Control.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/Control.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/Control.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/Control.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/Control.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/Control.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/Control.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/Control.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/Control.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/Control.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/Control.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/Control.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/Control.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/Control.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/Control.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../FB/Control.o: ../../FB/Component.h ../../FB/FBUtility.h
../../FB/Control.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/Control.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../FB/Control.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../FB/Control.o: ../../FB/Sampling.h ../../FB/BndOuter.h
../../FB/Control.o: ../../FB/Interval_Mngr.h ../../FB/TPControl.h
../../FB/FBUtility.o: ../../FB/FBUtility.h ../../FB/FB_Define.h
../../FB/FBUtility.o: ../../FB/mydebug.h
../../FB/FileIO.o: ../../FB/FileIO.h /usr/local/sphere/include/Skl.h
../../FB/FileIO.o: /usr/local/sphere/include/endianUtil.h
../../FB/FileIO.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/FileIO.o: /usr/local/sphere/include/SklBase.h
../../FB/FileIO.o: /usr/local/sphere/include/SklDefine.h
../../FB/FileIO.o: /usr/local/sphere/include/SklReservWord.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/FileIO.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/FileIO.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/FileIO.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/FileIO.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/FileIO.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/FileIO.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/FileIO.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/FileIO.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/FileIO.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/FileIO.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/FileIO.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/FileIO.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/FileIO.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/FileIO.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/FileIO.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/FileIO.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/FileIO.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/FileIO.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/FileIO.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../FB/FileIO.o: /usr/local/sphere/include/SklUtil.h
../../FB/FileIO.o: ../../FB/Parallel_node.h
../../FB/FileIO.o: /usr/local/sphere/include/fileio/SklSbxDataSet.h
../../FB/FileIO.o: ../../FB/FBUtility.h ../../FB/FB_Ffunc.h
../../FB/History.o: ../../FB/History.h ../../FB/Control.h
../../FB/History.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../FB/History.o: /usr/local/sphere/include/Skl.h
../../FB/History.o: /usr/local/sphere/include/endianUtil.h
../../FB/History.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/History.o: /usr/local/sphere/include/SklBase.h
../../FB/History.o: /usr/local/sphere/include/SklDefine.h
../../FB/History.o: /usr/local/sphere/include/SklReservWord.h
../../FB/History.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/History.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/History.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/History.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/History.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/History.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/History.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/History.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/History.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/History.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/History.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/History.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/History.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/History.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/History.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/History.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/History.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/History.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/History.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/History.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/History.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/History.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/History.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/History.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../FB/History.o: ../../FB/Component.h ../../FB/FBUtility.h
../../FB/History.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/History.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../FB/History.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../FB/History.o: ../../FB/Sampling.h ../../FB/BndOuter.h
../../FB/History.o: ../../FB/Interval_Mngr.h ../../FB/TPControl.h
../../FB/Intrinsic.o: ../../FB/Intrinsic.h ../../FB/FB_Define.h
../../FB/Intrinsic.o: ../../FB/mydebug.h ../../FB/BndOuter.h
../../FB/Intrinsic.o: ../../FB/Control.h /usr/local/sphere/include/Skl.h
../../FB/Intrinsic.o: /usr/local/sphere/include/endianUtil.h
../../FB/Intrinsic.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/Intrinsic.o: /usr/local/sphere/include/SklBase.h
../../FB/Intrinsic.o: /usr/local/sphere/include/SklDefine.h
../../FB/Intrinsic.o: /usr/local/sphere/include/SklReservWord.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/Intrinsic.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/Intrinsic.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/Intrinsic.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/Intrinsic.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/Intrinsic.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/Intrinsic.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/Intrinsic.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/Intrinsic.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/Intrinsic.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/Intrinsic.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/Intrinsic.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/Intrinsic.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/Intrinsic.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/Intrinsic.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/Intrinsic.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/Intrinsic.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/Intrinsic.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/Intrinsic.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/Intrinsic.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../FB/Intrinsic.o: ../../FB/Component.h ../../FB/FBUtility.h
../../FB/Intrinsic.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/Intrinsic.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../FB/Intrinsic.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../FB/Intrinsic.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../FB/Intrinsic.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../FB/ParseBC.o: ../../FB/ParseBC.h /usr/local/sphere/include/Skl.h
../../FB/ParseBC.o: /usr/local/sphere/include/endianUtil.h
../../FB/ParseBC.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/ParseBC.o: /usr/local/sphere/include/SklBase.h
../../FB/ParseBC.o: /usr/local/sphere/include/SklDefine.h
../../FB/ParseBC.o: /usr/local/sphere/include/SklReservWord.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/ParseBC.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/ParseBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/ParseBC.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/ParseBC.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/ParseBC.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/ParseBC.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/ParseBC.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/ParseBC.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/ParseBC.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/ParseBC.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/ParseBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/ParseBC.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/ParseBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/ParseBC.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/ParseBC.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/ParseBC.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/ParseBC.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/ParseBC.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/ParseBC.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../FB/ParseBC.o: ../../FB/FBUtility.h ../../FB/BndOuter.h
../../FB/ParseBC.o: ../../FB/Control.h /usr/local/sphere/include/SklUtil.h
../../FB/ParseBC.o: ../../FB/Medium.h ../../FB/Component.h ../../FB/Monitor.h
../../FB/ParseBC.o: ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/ParseBC.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../FB/ParseBC.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../FB/ParseBC.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../FB/ParseBC.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../FB/ParseBC.o: ../../FB/Intrinsic.h
../../FB/ParseMat.o: ../../FB/ParseMat.h ../../FB/FB_Define.h
../../FB/ParseMat.o: ../../FB/mydebug.h ../../FB/Control.h
../../FB/ParseMat.o: /usr/local/sphere/include/Skl.h
../../FB/ParseMat.o: /usr/local/sphere/include/endianUtil.h
../../FB/ParseMat.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/ParseMat.o: /usr/local/sphere/include/SklBase.h
../../FB/ParseMat.o: /usr/local/sphere/include/SklDefine.h
../../FB/ParseMat.o: /usr/local/sphere/include/SklReservWord.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/ParseMat.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/ParseMat.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/ParseMat.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/ParseMat.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/ParseMat.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/ParseMat.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/ParseMat.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/ParseMat.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/ParseMat.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/ParseMat.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/ParseMat.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/ParseMat.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/ParseMat.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/ParseMat.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/ParseMat.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/ParseMat.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/ParseMat.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/ParseMat.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/ParseMat.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../FB/ParseMat.o: ../../FB/Component.h ../../FB/FBUtility.h
../../FB/ParseMat.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/ParseMat.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../FB/ParseMat.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../FB/ParseMat.o: ../../FB/Sampling.h ../../FB/BndOuter.h
../../FB/ParseMat.o: ../../FB/Interval_Mngr.h ../../FB/TPControl.h
../../FB/SetBC.o: ../../FB/SetBC.h /usr/local/sphere/include/Skl.h
../../FB/SetBC.o: /usr/local/sphere/include/endianUtil.h
../../FB/SetBC.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/SetBC.o: /usr/local/sphere/include/SklBase.h
../../FB/SetBC.o: /usr/local/sphere/include/SklDefine.h
../../FB/SetBC.o: /usr/local/sphere/include/SklReservWord.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/SetBC.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/SetBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/SetBC.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/SetBC.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/SetBC.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/SetBC.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/SetBC.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/SetBC.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/SetBC.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/SetBC.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/SetBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/SetBC.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/SetBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/SetBC.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/SetBC.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/SetBC.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/SetBC.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/SetBC.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/SetBC.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../FB/SetBC.o: ../../FB/FBUtility.h ../../FB/BndOuter.h ../../FB/Control.h
../../FB/SetBC.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../FB/SetBC.o: ../../FB/Component.h ../../FB/Monitor.h ../../FB/vec3.h
../../FB/SetBC.o: ../../FB/vec3_func.h ../../FB/vec3f_func.h
../../FB/SetBC.o: ../../FB/Parallel_node.h ../../FB/MonCompo.h
../../FB/SetBC.o: ../../FB/basic_func.h ../../FB/Sampling.h
../../FB/SetBC.o: ../../FB/Interval_Mngr.h ../../FB/TPControl.h
../../FB/SetBC.o: ../../FB/Intrinsic.h ../../FB/IDtable.h
../../FB/Monitor.o: ../../FB/Monitor.h ../../FB/FB_Define.h
../../FB/Monitor.o: ../../FB/mydebug.h ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/Monitor.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../FB/Monitor.o: ../../FB/MonCompo.h /usr/local/sphere/include/SklUtil.h
../../FB/Monitor.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/Monitor.o: /usr/local/sphere/include/Skl.h
../../FB/Monitor.o: /usr/local/sphere/include/endianUtil.h
../../FB/Monitor.o: /usr/local/sphere/include/SklBase.h
../../FB/Monitor.o: /usr/local/sphere/include/SklDefine.h
../../FB/Monitor.o: /usr/local/sphere/include/SklReservWord.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/Monitor.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/Monitor.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/Monitor.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/Monitor.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/Monitor.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/Monitor.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/Monitor.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/Monitor.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/Monitor.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/Monitor.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/Monitor.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/Monitor.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/Monitor.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/Monitor.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/Monitor.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/Monitor.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/Monitor.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/Monitor.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/Monitor.o: ../../FB/basic_func.h ../../FB/FBUtility.h
../../FB/Monitor.o: ../../FB/Component.h ../../FB/Sampling.h
../../FB/MonCompo.o: ../../FB/MonCompo.h ../../FB/FB_Define.h
../../FB/MonCompo.o: ../../FB/mydebug.h /usr/local/sphere/include/SklUtil.h
../../FB/MonCompo.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/MonCompo.o: /usr/local/sphere/include/Skl.h
../../FB/MonCompo.o: /usr/local/sphere/include/endianUtil.h
../../FB/MonCompo.o: /usr/local/sphere/include/SklBase.h
../../FB/MonCompo.o: /usr/local/sphere/include/SklDefine.h
../../FB/MonCompo.o: /usr/local/sphere/include/SklReservWord.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/MonCompo.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/MonCompo.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/MonCompo.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/MonCompo.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/MonCompo.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/MonCompo.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/MonCompo.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/MonCompo.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/MonCompo.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/MonCompo.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/MonCompo.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/MonCompo.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/MonCompo.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/MonCompo.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/MonCompo.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/MonCompo.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/MonCompo.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/MonCompo.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/MonCompo.o: ../../FB/vec3.h ../../FB/vec3_func.h
../../FB/MonCompo.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../FB/MonCompo.o: ../../FB/basic_func.h ../../FB/FBUtility.h
../../FB/MonCompo.o: ../../FB/Component.h ../../FB/Sampling.h
../../FB/Interval_Mngr.o: ../../FB/Interval_Mngr.h ../../FB/FB_Define.h
../../FB/Interval_Mngr.o: ../../FB/mydebug.h
../../FB/Sampling.o: ../../FB/Sampling.h ../../FB/FB_Define.h
../../FB/Sampling.o: ../../FB/mydebug.h ../../FB/FBUtility.h ../../FB/vec3.h
../../FB/Sampling.o: ../../FB/vec3_func.h ../../FB/vec3f_func.h
../../FB/Sampling.o: ../../FB/basic_func.h
../../FB/DataHolder.o: ../../FB/DataHolder.h ../../FB/FB_Define.h
../../FB/DataHolder.o: ../../FB/mydebug.h /usr/local/sphere/include/Skl.h
../../FB/DataHolder.o: /usr/local/sphere/include/endianUtil.h
../../FB/DataHolder.o: /usr/local/sphere/include/SklSolverBase.h
../../FB/DataHolder.o: /usr/local/sphere/include/SklBase.h
../../FB/DataHolder.o: /usr/local/sphere/include/SklDefine.h
../../FB/DataHolder.o: /usr/local/sphere/include/SklReservWord.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgBase.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../FB/DataHolder.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../FB/DataHolder.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../FB/DataHolder.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../FB/DataHolder.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../FB/DataHolder.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../FB/DataHolder.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklAray.h
../../FB/DataHolder.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../FB/DataHolder.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../FB/DataHolder.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../FB/DataHolder.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../FB/DataHolder.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../FB/DataHolder.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../FB/DataHolder.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../FB/DataHolder.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../FB/DataHolder.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../FB/DataHolder.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../FB/DataHolder.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../FB/DataHolder.o: /usr/local/sphere/include/coord/SklCoord.h
../../FB/DataHolder.o: ../../FB/Parallel_node.h
../../FB/dfi.o: ../../FB/dfi.h ../../FB/FB_Define.h ../../FB/mydebug.h
../../FB/dfi.o: ../../FB/util_Path.h
../../FB/TPControl.o: ../../FB/TPControl.h ../../FB/FB_Define.h
../../FB/TPControl.o: ../../FB/mydebug.h

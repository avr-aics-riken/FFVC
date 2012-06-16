#
# SPHERE - Skeleton for PHysicalSvCube and Engineering REsearch
#
# Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
#

TARGET = libsphnsolvIP.a

CSRCS = \


CXXSRCS = \
	../../IP/IP_PMT.C \
	../../IP/IP_Rect.C \
	../../IP/IP_Polygon.C \
	../../IP/IP_Cylinder.C \
	../../IP/IP_Step.C \
	../../IP/IP_PPLT2D.C \
	../../IP/IP_Duct.C \
	../../IP/IP_SHC1D.C \
	../../IP/IP_Sphere.C \


FCSRCS = \


F90SRCS = \


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

../../IP/IP_PMT.o: ../../IP/IP_PMT.h ../../FB/Intrinsic.h
../../IP/IP_PMT.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_PMT.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_PMT.o: /usr/local/sphere/include/Skl.h
../../IP/IP_PMT.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_PMT.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_PMT.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_PMT.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_PMT.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_PMT.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_PMT.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_PMT.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_PMT.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_PMT.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_PMT.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_PMT.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_PMT.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_PMT.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_PMT.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_PMT.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_PMT.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_PMT.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_PMT.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_PMT.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_PMT.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_PMT.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_PMT.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_PMT.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_PMT.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_PMT.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../IP/IP_PMT.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../IP/IP_PMT.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../IP/IP_PMT.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../IP/IP_PMT.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../IP/IP_PMT.o: ../../IP/IP_Define.h
../../IP/IP_Rect.o: ../../IP/IP_Rect.h ../../FB/Intrinsic.h
../../IP/IP_Rect.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_Rect.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_Rect.o: /usr/local/sphere/include/Skl.h
../../IP/IP_Rect.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_Rect.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_Rect.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_Rect.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_Rect.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Rect.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Rect.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_Rect.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_Rect.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_Rect.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_Rect.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_Rect.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_Rect.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_Rect.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_Rect.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_Rect.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_Rect.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_Rect.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_Rect.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_Rect.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_Rect.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_Rect.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_Rect.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Rect.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_Rect.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_Rect.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../IP/IP_Rect.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../IP/IP_Rect.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../IP/IP_Rect.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../IP/IP_Rect.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../IP/IP_Rect.o: ../../IP/IP_Define.h
../../IP/IP_Polygon.o: ../../IP/IP_Polygon.h ../../FB/Intrinsic.h
../../IP/IP_Polygon.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_Polygon.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/Skl.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Polygon.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_Polygon.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_Polygon.o: ../../FB/Monitor.h ../../FB/vec3.h
../../IP/IP_Polygon.o: ../../FB/vec3_func.h ../../FB/vec3f_func.h
../../IP/IP_Polygon.o: ../../FB/Parallel_node.h ../../FB/MonCompo.h
../../IP/IP_Polygon.o: ../../FB/basic_func.h ../../FB/Sampling.h
../../IP/IP_Polygon.o: ../../FB/Interval_Mngr.h ../../FB/TPControl.h
../../IP/IP_Polygon.o: ../../FB/IDtable.h ../../IP/IP_Define.h
../../IP/IP_Cylinder.o: ../../IP/IP_Cylinder.h ../../FB/Intrinsic.h
../../IP/IP_Cylinder.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_Cylinder.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/Skl.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Cylinder.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_Cylinder.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_Cylinder.o: ../../FB/Monitor.h ../../FB/vec3.h
../../IP/IP_Cylinder.o: ../../FB/vec3_func.h ../../FB/vec3f_func.h
../../IP/IP_Cylinder.o: ../../FB/Parallel_node.h ../../FB/MonCompo.h
../../IP/IP_Cylinder.o: ../../FB/basic_func.h ../../FB/Sampling.h
../../IP/IP_Cylinder.o: ../../FB/Interval_Mngr.h ../../FB/TPControl.h
../../IP/IP_Cylinder.o: ../../FB/IDtable.h ../../IP/IP_Define.h
../../IP/IP_Step.o: ../../IP/IP_Step.h ../../FB/Intrinsic.h
../../IP/IP_Step.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_Step.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_Step.o: /usr/local/sphere/include/Skl.h
../../IP/IP_Step.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_Step.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_Step.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_Step.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_Step.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Step.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Step.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_Step.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_Step.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_Step.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_Step.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_Step.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_Step.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_Step.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_Step.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_Step.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_Step.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_Step.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_Step.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_Step.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_Step.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_Step.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_Step.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Step.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_Step.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_Step.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../IP/IP_Step.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../IP/IP_Step.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../IP/IP_Step.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../IP/IP_Step.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../IP/IP_Step.o: ../../IP/IP_Define.h
../../IP/IP_PPLT2D.o: ../../IP/IP_PPLT2D.h ../../FB/Intrinsic.h
../../IP/IP_PPLT2D.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_PPLT2D.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/Skl.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_PPLT2D.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_PPLT2D.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_PPLT2D.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../IP/IP_PPLT2D.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../IP/IP_PPLT2D.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../IP/IP_PPLT2D.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../IP/IP_PPLT2D.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../IP/IP_PPLT2D.o: ../../IP/IP_Define.h
../../IP/IP_Duct.o: ../../IP/IP_Duct.h ../../FB/Intrinsic.h
../../IP/IP_Duct.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_Duct.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_Duct.o: /usr/local/sphere/include/Skl.h
../../IP/IP_Duct.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_Duct.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_Duct.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_Duct.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_Duct.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Duct.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Duct.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_Duct.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_Duct.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_Duct.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_Duct.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_Duct.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_Duct.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_Duct.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_Duct.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_Duct.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_Duct.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_Duct.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_Duct.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_Duct.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_Duct.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_Duct.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_Duct.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Duct.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_Duct.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_Duct.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../IP/IP_Duct.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../IP/IP_Duct.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../IP/IP_Duct.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../IP/IP_Duct.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../IP/IP_Duct.o: ../../IP/IP_Define.h
../../IP/IP_SHC1D.o: ../../IP/IP_SHC1D.h ../../FB/Intrinsic.h
../../IP/IP_SHC1D.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_SHC1D.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/Skl.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_SHC1D.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_SHC1D.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_SHC1D.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../IP/IP_SHC1D.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../IP/IP_SHC1D.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../IP/IP_SHC1D.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../IP/IP_SHC1D.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../IP/IP_SHC1D.o: ../../IP/IP_Define.h
../../IP/IP_Sphere.o: ../../IP/IP_Sphere.h ../../FB/Intrinsic.h
../../IP/IP_Sphere.o: ../../FB/FB_Define.h ../../FB/mydebug.h
../../IP/IP_Sphere.o: ../../FB/BndOuter.h ../../FB/Control.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/Skl.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/endianUtil.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/SklSolverBase.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/SklBase.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/SklDefine.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/SklReservWord.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklSolverConfig.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgBase.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgSteer.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgInFile.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgHistory.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgParameter.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklIBCConfig.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklOBCConfig.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklMBoxConfig.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklMDMConfig.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/base/SklFileIOCntl.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgDatas.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/base/SklOutputDataSet.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklDataManager.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklArayBase.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklScalar.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklAray.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/parallel/SklParaManager.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/parallel/SklParaIF.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/parallel/SklParaComponent.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklVector2D.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklVector3D.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/coord/SklCrdMngr.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/coord/SklCoord.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/config/SklCfgOutFile.h
../../IP/IP_Sphere.o: /usr/local/sphere/include/SklUtil.h ../../FB/Medium.h
../../IP/IP_Sphere.o: ../../FB/Component.h ../../FB/FBUtility.h
../../IP/IP_Sphere.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
../../IP/IP_Sphere.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
../../IP/IP_Sphere.o: ../../FB/MonCompo.h ../../FB/basic_func.h
../../IP/IP_Sphere.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
../../IP/IP_Sphere.o: ../../FB/TPControl.h ../../FB/IDtable.h
../../IP/IP_Sphere.o: ../../IP/IP_Define.h

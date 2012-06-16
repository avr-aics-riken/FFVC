#
# SPHERE - Skeleton for PHysical and Engineering REsearch
#
# Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
#

SPH_APP_SRCS = SklCreateSolver.C SklDeleteSolver.C
SOLV_FACT_SRCS = 
SOLV_FACT_SRCS += SklFactoryCBC.C

SRCS = $(SPH_APP_SRCS) $(SOLV_FACT_SRCS)

.SUFFIXES: .o .C
OBJS = $(SRCS:.C=.o)

$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $(XML2FLAGS) $(SPHERE_CFLAGS) $(REALOPT) -o $(TARGET) $(OBJS) \
	$(LDFLAGS) $(MPICH_LDFLAGS) $(SPHERE_LDFLAGS) $(SPH_SOLV_LDFAGS) \
	$(SPH_USR_DEF_LIBS) $(SPH_SOLV_LIBS) $(SPHERE_LIBS) $(MPICH_LIBS) $(XML2LIBS) $(LIBS)
	mv $(TARGET) ../bin/.

.C.o:
	$(CXX) $(CXXFLAGS) $(XML2FLAGS) $(SPHERE_CFLAGS) $(SPH_SOLV_FLAGS) $(REALOPT) -c $<

clean:
	$(RM) $(OBJS)

depend:
	makedepend -I. -I$(SPHEREDIR)/include $(SOLVER_INC_DIR) -Y *.C

# DO NOT DELETE

SklCreateSolver.o: /usr/local/sphere/include/SklGloval.h
SklCreateSolver.o: /usr/local/sphere/include/SklDefine.h
SklCreateSolver.o: /usr/local/sphere/include/SklReservWord.h
SklCreateSolver.o: /usr/local/sphere/include/SklSolverBase.h
SklCreateSolver.o: /usr/local/sphere/include/Skl.h
SklCreateSolver.o: /usr/local/sphere/include/endianUtil.h
SklCreateSolver.o: /usr/local/sphere/include/SklBase.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklSolverConfig.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgBase.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgBase.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgSteer.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgInFile.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgHistory.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgOutFile.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgParameter.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklIBCConfig.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklOBCConfig.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklMBoxConfig.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklMDMConfig.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklCreateSolver.o: /usr/local/sphere/include/base/SklFileIOCntl.h
SklCreateSolver.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklCreateSolver.o: /usr/local/sphere/include/base/SklOutputDataSet.h
SklCreateSolver.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
SklCreateSolver.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
SklCreateSolver.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklDataManager.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklArayBase.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklScalar.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklAray.h
SklCreateSolver.o: /usr/local/sphere/include/parallel/SklParaManager.h
SklCreateSolver.o: /usr/local/sphere/include/parallel/SklParaIF.h
SklCreateSolver.o: /usr/local/sphere/include/parallel/SklParaComponent.h
SklCreateSolver.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
SklCreateSolver.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
SklCreateSolver.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
SklCreateSolver.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklVector2D.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklVector3D.h
SklCreateSolver.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
SklCreateSolver.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
SklCreateSolver.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
SklCreateSolver.o: /usr/local/sphere/include/coord/SklCrdMngr.h
SklCreateSolver.o: /usr/local/sphere/include/coord/SklCoord.h SklSolverType.h
SklCreateSolver.o: SklFactoryCBC.h
SklCreateSolver.o: /usr/local/sphere/include/SklSolvFactoryBase.h
SklDeleteSolver.o: /usr/local/sphere/include/SklGloval.h
SklDeleteSolver.o: /usr/local/sphere/include/SklDefine.h
SklDeleteSolver.o: /usr/local/sphere/include/SklReservWord.h
SklDeleteSolver.o: /usr/local/sphere/include/SklSolverBase.h
SklDeleteSolver.o: /usr/local/sphere/include/Skl.h
SklDeleteSolver.o: /usr/local/sphere/include/endianUtil.h
SklDeleteSolver.o: /usr/local/sphere/include/SklBase.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklSolverConfig.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgBase.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgBase.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgSteer.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgInFile.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgHistory.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgOutFile.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgParameter.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklIBCConfig.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklOBCConfig.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklMBoxConfig.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklMDMConfig.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklDeleteSolver.o: /usr/local/sphere/include/base/SklFileIOCntl.h
SklDeleteSolver.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklDeleteSolver.o: /usr/local/sphere/include/base/SklOutputDataSet.h
SklDeleteSolver.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
SklDeleteSolver.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
SklDeleteSolver.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklDataManager.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklArayBase.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklScalar.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklAray.h
SklDeleteSolver.o: /usr/local/sphere/include/parallel/SklParaManager.h
SklDeleteSolver.o: /usr/local/sphere/include/parallel/SklParaIF.h
SklDeleteSolver.o: /usr/local/sphere/include/parallel/SklParaComponent.h
SklDeleteSolver.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
SklDeleteSolver.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
SklDeleteSolver.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
SklDeleteSolver.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklVector2D.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklVector3D.h
SklDeleteSolver.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
SklDeleteSolver.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
SklDeleteSolver.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
SklDeleteSolver.o: /usr/local/sphere/include/coord/SklCrdMngr.h
SklDeleteSolver.o: /usr/local/sphere/include/coord/SklCoord.h SklSolverType.h
SklDeleteSolver.o: SklFactoryCBC.h
SklDeleteSolver.o: /usr/local/sphere/include/SklSolvFactoryBase.h
SklFactoryCBC.o: SklFactoryCBC.h
SklFactoryCBC.o: /usr/local/sphere/include/SklSolvFactoryBase.h
SklFactoryCBC.o: /usr/local/sphere/include/SklGloval.h
SklFactoryCBC.o: /usr/local/sphere/include/SklDefine.h
SklFactoryCBC.o: /usr/local/sphere/include/SklReservWord.h
SklFactoryCBC.o: /usr/local/sphere/include/SklSolverBase.h
SklFactoryCBC.o: /usr/local/sphere/include/Skl.h
SklFactoryCBC.o: /usr/local/sphere/include/endianUtil.h
SklFactoryCBC.o: /usr/local/sphere/include/SklBase.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklSolverConfig.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgBase.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgBase.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgSteer.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgInFile.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgHistory.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgOutFile.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgParameter.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklIBCConfig.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklOBCConfig.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklMBoxConfig.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklMDMConfig.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklFactoryCBC.o: /usr/local/sphere/include/base/SklFileIOCntl.h
SklFactoryCBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklFactoryCBC.o: /usr/local/sphere/include/base/SklOutputDataSet.h
SklFactoryCBC.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
SklFactoryCBC.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
SklFactoryCBC.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklDataManager.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklArayBase.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklScalar.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklAray.h
SklFactoryCBC.o: /usr/local/sphere/include/parallel/SklParaManager.h
SklFactoryCBC.o: /usr/local/sphere/include/parallel/SklParaIF.h
SklFactoryCBC.o: /usr/local/sphere/include/parallel/SklParaComponent.h
SklFactoryCBC.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
SklFactoryCBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
SklFactoryCBC.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
SklFactoryCBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklVector2D.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklVector3D.h
SklFactoryCBC.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
SklFactoryCBC.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
SklFactoryCBC.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
SklFactoryCBC.o: /usr/local/sphere/include/coord/SklCrdMngr.h
SklFactoryCBC.o: /usr/local/sphere/include/coord/SklCoord.h
SklFactoryCBC.o: ../CBC/SklSolverCBC.h ../../FB/FB_Define.h
SklFactoryCBC.o: ../../FB/mydebug.h ../CBC/CBC_Define.h
SklFactoryCBC.o: ../CBC/FortranFuncCBC.h /usr/local/sphere/include/SklUtil.h
SklFactoryCBC.o: ../CBC/CBC_SetBC.h ../../IP/IP_Duct.h ../../FB/Intrinsic.h
SklFactoryCBC.o: ../../FB/FB_Define.h ../../FB/BndOuter.h ../../FB/Control.h
SklFactoryCBC.o: ../../FB/Medium.h ../../FB/Component.h ../../FB/FBUtility.h
SklFactoryCBC.o: ../../FB/Monitor.h ../../FB/vec3.h ../../FB/vec3_func.h
SklFactoryCBC.o: ../../FB/vec3f_func.h ../../FB/Parallel_node.h
SklFactoryCBC.o: ../../FB/MonCompo.h ../../FB/basic_func.h
SklFactoryCBC.o: ../../FB/Sampling.h ../../FB/Interval_Mngr.h
SklFactoryCBC.o: ../../FB/TPControl.h ../../FB/IDtable.h ../../IP/IP_Define.h
SklFactoryCBC.o: ../../IP/IP_PPLT2D.h ../../IP/IP_SHC1D.h ../../IP/IP_PMT.h
SklFactoryCBC.o: ../../IP/IP_Rect.h ../../IP/IP_Step.h ../../IP/IP_Cylinder.h
SklFactoryCBC.o: ../../IP/IP_Polygon.h ../../IP/IP_Sphere.h

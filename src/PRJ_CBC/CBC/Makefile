#
# SPHERE - Skeleton for PHysicalSvCube and Engineering REsearch
#
# Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
#

TARGET = libsphsolvCBC.a

CSRCS = \


CXXSRCS = \
       SklSolverCBCInitialize.C \
       SklSolverCBCLoop.C \
       SklSolverCBCPost.C \
       SklSolverCBCUsage.C \
       SklSolverCBC.C \
       CBC_Heat.C \
       CBC_SetBC.C \
       NS_FS_E_CBC.C \
       NS_FS_E_CDS.C \
       PS_E_CBC.C \
       IF_TRP_VOF.C \


FCSRCS = \


F90SRCS = \
       ../../F_CBC/BCvec_cc.f90 \
       ../../F_CBC/BCprs_cc.f90 \
       ../../F_CBC/cbc_3d.f90 \
       ../../F_CBC/cbc_poisson.f90 \
       ../../F_CBC/cbc_utility.f90 \
       ../../F_CBC/cbc_pscalar.f90 \
       ../../F_CBC/cbc_forcing.f90 \
       ../../F_CDS/cds_poisson.f90 \
       ../../F_CDS/cds_vector.f90 \
       ../../F_VOF/c3d_vof.f90 \


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

SklSolverCBCInitialize.o: SklSolverCBC.h CBC_Define.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/Skl.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/endianUtil.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/SklSolverBase.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/SklBase.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/SklDefine.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/SklReservWord.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklSolverConfig.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgSteer.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgInFile.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgHistory.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgOutFile.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgParameter.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklIBCConfig.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklOBCConfig.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklMBoxConfig.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklMDMConfig.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/base/SklFileIOCntl.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/base/SklOutputDataSet.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklDataManager.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklArayBase.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklScalar.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklAray.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/parallel/SklParaManager.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/parallel/SklParaIF.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/parallel/SklParaComponent.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklVector2D.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklVector3D.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/coord/SklCrdMngr.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/coord/SklCoord.h
SklSolverCBCInitialize.o: FortranFuncCBC.h
SklSolverCBCInitialize.o: /usr/local/sphere/include/SklUtil.h CBC_SetBC.h
SklSolverCBCLoop.o: SklSolverCBC.h CBC_Define.h
SklSolverCBCLoop.o: /usr/local/sphere/include/Skl.h
SklSolverCBCLoop.o: /usr/local/sphere/include/endianUtil.h
SklSolverCBCLoop.o: /usr/local/sphere/include/SklSolverBase.h
SklSolverCBCLoop.o: /usr/local/sphere/include/SklBase.h
SklSolverCBCLoop.o: /usr/local/sphere/include/SklDefine.h
SklSolverCBCLoop.o: /usr/local/sphere/include/SklReservWord.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklSolverConfig.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgSteer.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgInFile.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgHistory.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgOutFile.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgParameter.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklIBCConfig.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklOBCConfig.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklMBoxConfig.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklMDMConfig.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBCLoop.o: /usr/local/sphere/include/base/SklFileIOCntl.h
SklSolverCBCLoop.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBCLoop.o: /usr/local/sphere/include/base/SklOutputDataSet.h
SklSolverCBCLoop.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
SklSolverCBCLoop.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
SklSolverCBCLoop.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklDataManager.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklArayBase.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklScalar.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklAray.h
SklSolverCBCLoop.o: /usr/local/sphere/include/parallel/SklParaManager.h
SklSolverCBCLoop.o: /usr/local/sphere/include/parallel/SklParaIF.h
SklSolverCBCLoop.o: /usr/local/sphere/include/parallel/SklParaComponent.h
SklSolverCBCLoop.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
SklSolverCBCLoop.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
SklSolverCBCLoop.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
SklSolverCBCLoop.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklVector2D.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklVector3D.h
SklSolverCBCLoop.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
SklSolverCBCLoop.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
SklSolverCBCLoop.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
SklSolverCBCLoop.o: /usr/local/sphere/include/coord/SklCrdMngr.h
SklSolverCBCLoop.o: /usr/local/sphere/include/coord/SklCoord.h
SklSolverCBCLoop.o: FortranFuncCBC.h /usr/local/sphere/include/SklUtil.h
SklSolverCBCLoop.o: CBC_SetBC.h
SklSolverCBCPost.o: SklSolverCBC.h CBC_Define.h
SklSolverCBCPost.o: /usr/local/sphere/include/Skl.h
SklSolverCBCPost.o: /usr/local/sphere/include/endianUtil.h
SklSolverCBCPost.o: /usr/local/sphere/include/SklSolverBase.h
SklSolverCBCPost.o: /usr/local/sphere/include/SklBase.h
SklSolverCBCPost.o: /usr/local/sphere/include/SklDefine.h
SklSolverCBCPost.o: /usr/local/sphere/include/SklReservWord.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklSolverConfig.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgSteer.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgInFile.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgHistory.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgOutFile.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgParameter.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklIBCConfig.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklOBCConfig.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklMBoxConfig.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklMDMConfig.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBCPost.o: /usr/local/sphere/include/base/SklFileIOCntl.h
SklSolverCBCPost.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBCPost.o: /usr/local/sphere/include/base/SklOutputDataSet.h
SklSolverCBCPost.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
SklSolverCBCPost.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
SklSolverCBCPost.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklDataManager.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklArayBase.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklScalar.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklAray.h
SklSolverCBCPost.o: /usr/local/sphere/include/parallel/SklParaManager.h
SklSolverCBCPost.o: /usr/local/sphere/include/parallel/SklParaIF.h
SklSolverCBCPost.o: /usr/local/sphere/include/parallel/SklParaComponent.h
SklSolverCBCPost.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
SklSolverCBCPost.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
SklSolverCBCPost.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
SklSolverCBCPost.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklVector2D.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklVector3D.h
SklSolverCBCPost.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
SklSolverCBCPost.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
SklSolverCBCPost.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
SklSolverCBCPost.o: /usr/local/sphere/include/coord/SklCrdMngr.h
SklSolverCBCPost.o: /usr/local/sphere/include/coord/SklCoord.h
SklSolverCBCPost.o: FortranFuncCBC.h /usr/local/sphere/include/SklUtil.h
SklSolverCBCPost.o: CBC_SetBC.h
SklSolverCBCUsage.o: SklSolverCBC.h CBC_Define.h
SklSolverCBCUsage.o: /usr/local/sphere/include/Skl.h
SklSolverCBCUsage.o: /usr/local/sphere/include/endianUtil.h
SklSolverCBCUsage.o: /usr/local/sphere/include/SklSolverBase.h
SklSolverCBCUsage.o: /usr/local/sphere/include/SklBase.h
SklSolverCBCUsage.o: /usr/local/sphere/include/SklDefine.h
SklSolverCBCUsage.o: /usr/local/sphere/include/SklReservWord.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklSolverConfig.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgSteer.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgInFile.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgHistory.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgOutFile.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgParameter.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklIBCConfig.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklOBCConfig.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklMBoxConfig.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklMDMConfig.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBCUsage.o: /usr/local/sphere/include/base/SklFileIOCntl.h
SklSolverCBCUsage.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBCUsage.o: /usr/local/sphere/include/base/SklOutputDataSet.h
SklSolverCBCUsage.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
SklSolverCBCUsage.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
SklSolverCBCUsage.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklDataManager.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklArayBase.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklScalar.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklAray.h
SklSolverCBCUsage.o: /usr/local/sphere/include/parallel/SklParaManager.h
SklSolverCBCUsage.o: /usr/local/sphere/include/parallel/SklParaIF.h
SklSolverCBCUsage.o: /usr/local/sphere/include/parallel/SklParaComponent.h
SklSolverCBCUsage.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
SklSolverCBCUsage.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
SklSolverCBCUsage.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
SklSolverCBCUsage.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklVector2D.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklVector3D.h
SklSolverCBCUsage.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
SklSolverCBCUsage.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
SklSolverCBCUsage.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
SklSolverCBCUsage.o: /usr/local/sphere/include/coord/SklCrdMngr.h
SklSolverCBCUsage.o: /usr/local/sphere/include/coord/SklCoord.h
SklSolverCBCUsage.o: FortranFuncCBC.h /usr/local/sphere/include/SklUtil.h
SklSolverCBCUsage.o: CBC_SetBC.h
SklSolverCBC.o: /usr/local/sphere/include/SklGloval.h
SklSolverCBC.o: /usr/local/sphere/include/SklDefine.h
SklSolverCBC.o: /usr/local/sphere/include/SklReservWord.h
SklSolverCBC.o: /usr/local/sphere/include/SklSolverBase.h
SklSolverCBC.o: /usr/local/sphere/include/Skl.h
SklSolverCBC.o: /usr/local/sphere/include/endianUtil.h
SklSolverCBC.o: /usr/local/sphere/include/SklBase.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklSolverConfig.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgBase.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgSteer.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgInFile.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgHistory.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgOutFile.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgParameter.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklIBCConfig.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklOBCConfig.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklMBoxConfig.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklMDMConfig.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
SklSolverCBC.o: /usr/local/sphere/include/base/SklFileIOCntl.h
SklSolverCBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
SklSolverCBC.o: /usr/local/sphere/include/base/SklOutputDataSet.h
SklSolverCBC.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
SklSolverCBC.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
SklSolverCBC.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklDataManager.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklArayBase.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklScalar.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklAray.h
SklSolverCBC.o: /usr/local/sphere/include/parallel/SklParaManager.h
SklSolverCBC.o: /usr/local/sphere/include/parallel/SklParaIF.h
SklSolverCBC.o: /usr/local/sphere/include/parallel/SklParaComponent.h
SklSolverCBC.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
SklSolverCBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
SklSolverCBC.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
SklSolverCBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklVector2D.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklVector3D.h
SklSolverCBC.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
SklSolverCBC.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
SklSolverCBC.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
SklSolverCBC.o: /usr/local/sphere/include/coord/SklCrdMngr.h
SklSolverCBC.o: /usr/local/sphere/include/coord/SklCoord.h SklSolverCBC.h
SklSolverCBC.o: CBC_Define.h FortranFuncCBC.h
SklSolverCBC.o: /usr/local/sphere/include/SklUtil.h CBC_SetBC.h
CBC_Heat.o: SklSolverCBC.h CBC_Define.h /usr/local/sphere/include/Skl.h
CBC_Heat.o: /usr/local/sphere/include/endianUtil.h
CBC_Heat.o: /usr/local/sphere/include/SklSolverBase.h
CBC_Heat.o: /usr/local/sphere/include/SklBase.h
CBC_Heat.o: /usr/local/sphere/include/SklDefine.h
CBC_Heat.o: /usr/local/sphere/include/SklReservWord.h
CBC_Heat.o: /usr/local/sphere/include/config/SklSolverConfig.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgDatas.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgBase.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgBase.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgSteer.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgInFile.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgHistory.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgOutFile.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgParameter.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
CBC_Heat.o: /usr/local/sphere/include/config/SklIBCConfig.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
CBC_Heat.o: /usr/local/sphere/include/config/SklOBCConfig.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
CBC_Heat.o: /usr/local/sphere/include/config/SklMBoxConfig.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
CBC_Heat.o: /usr/local/sphere/include/config/SklMDMConfig.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
CBC_Heat.o: /usr/local/sphere/include/base/SklFileIOCntl.h
CBC_Heat.o: /usr/local/sphere/include/config/SklCfgDatas.h
CBC_Heat.o: /usr/local/sphere/include/base/SklOutputDataSet.h
CBC_Heat.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
CBC_Heat.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
CBC_Heat.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklDataManager.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklArayBase.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklScalar.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklAray.h
CBC_Heat.o: /usr/local/sphere/include/parallel/SklParaManager.h
CBC_Heat.o: /usr/local/sphere/include/parallel/SklParaIF.h
CBC_Heat.o: /usr/local/sphere/include/parallel/SklParaComponent.h
CBC_Heat.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
CBC_Heat.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
CBC_Heat.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
CBC_Heat.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklVector2D.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklVector3D.h
CBC_Heat.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
CBC_Heat.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
CBC_Heat.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
CBC_Heat.o: /usr/local/sphere/include/coord/SklCrdMngr.h
CBC_Heat.o: /usr/local/sphere/include/coord/SklCoord.h FortranFuncCBC.h
CBC_Heat.o: /usr/local/sphere/include/SklUtil.h CBC_SetBC.h
CBC_SetBC.o: CBC_SetBC.h FortranFuncCBC.h CBC_Define.h
NS_FS_E_CBC.o: SklSolverCBC.h CBC_Define.h /usr/local/sphere/include/Skl.h
NS_FS_E_CBC.o: /usr/local/sphere/include/endianUtil.h
NS_FS_E_CBC.o: /usr/local/sphere/include/SklSolverBase.h
NS_FS_E_CBC.o: /usr/local/sphere/include/SklBase.h
NS_FS_E_CBC.o: /usr/local/sphere/include/SklDefine.h
NS_FS_E_CBC.o: /usr/local/sphere/include/SklReservWord.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklSolverConfig.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgBase.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgBase.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgSteer.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgInFile.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgHistory.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgOutFile.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgParameter.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklIBCConfig.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklOBCConfig.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklMBoxConfig.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklMDMConfig.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
NS_FS_E_CBC.o: /usr/local/sphere/include/base/SklFileIOCntl.h
NS_FS_E_CBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
NS_FS_E_CBC.o: /usr/local/sphere/include/base/SklOutputDataSet.h
NS_FS_E_CBC.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
NS_FS_E_CBC.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
NS_FS_E_CBC.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklDataManager.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayBase.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklAray.h
NS_FS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaManager.h
NS_FS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaIF.h
NS_FS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaComponent.h
NS_FS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
NS_FS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
NS_FS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
NS_FS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklVector2D.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklVector3D.h
NS_FS_E_CBC.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
NS_FS_E_CBC.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
NS_FS_E_CBC.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
NS_FS_E_CBC.o: /usr/local/sphere/include/coord/SklCrdMngr.h
NS_FS_E_CBC.o: /usr/local/sphere/include/coord/SklCoord.h FortranFuncCBC.h
NS_FS_E_CBC.o: /usr/local/sphere/include/SklUtil.h CBC_SetBC.h
NS_FS_E_CDS.o: SklSolverCBC.h CBC_Define.h /usr/local/sphere/include/Skl.h
NS_FS_E_CDS.o: /usr/local/sphere/include/endianUtil.h
NS_FS_E_CDS.o: /usr/local/sphere/include/SklSolverBase.h
NS_FS_E_CDS.o: /usr/local/sphere/include/SklBase.h
NS_FS_E_CDS.o: /usr/local/sphere/include/SklDefine.h
NS_FS_E_CDS.o: /usr/local/sphere/include/SklReservWord.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklSolverConfig.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgDatas.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgBase.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgBase.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgSteer.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgInFile.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgHistory.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgOutFile.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgParameter.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklIBCConfig.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklOBCConfig.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklMBoxConfig.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklMDMConfig.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
NS_FS_E_CDS.o: /usr/local/sphere/include/base/SklFileIOCntl.h
NS_FS_E_CDS.o: /usr/local/sphere/include/config/SklCfgDatas.h
NS_FS_E_CDS.o: /usr/local/sphere/include/base/SklOutputDataSet.h
NS_FS_E_CDS.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
NS_FS_E_CDS.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
NS_FS_E_CDS.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklDataManager.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklArayBase.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklScalar.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklAray.h
NS_FS_E_CDS.o: /usr/local/sphere/include/parallel/SklParaManager.h
NS_FS_E_CDS.o: /usr/local/sphere/include/parallel/SklParaIF.h
NS_FS_E_CDS.o: /usr/local/sphere/include/parallel/SklParaComponent.h
NS_FS_E_CDS.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
NS_FS_E_CDS.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
NS_FS_E_CDS.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
NS_FS_E_CDS.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklVector2D.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklVector3D.h
NS_FS_E_CDS.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
NS_FS_E_CDS.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
NS_FS_E_CDS.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
NS_FS_E_CDS.o: /usr/local/sphere/include/coord/SklCrdMngr.h
NS_FS_E_CDS.o: /usr/local/sphere/include/coord/SklCoord.h FortranFuncCBC.h
NS_FS_E_CDS.o: /usr/local/sphere/include/SklUtil.h CBC_SetBC.h
PS_E_CBC.o: SklSolverCBC.h CBC_Define.h /usr/local/sphere/include/Skl.h
PS_E_CBC.o: /usr/local/sphere/include/endianUtil.h
PS_E_CBC.o: /usr/local/sphere/include/SklSolverBase.h
PS_E_CBC.o: /usr/local/sphere/include/SklBase.h
PS_E_CBC.o: /usr/local/sphere/include/SklDefine.h
PS_E_CBC.o: /usr/local/sphere/include/SklReservWord.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklSolverConfig.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgBase.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgBase.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgSteer.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgInFile.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgHistory.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgOutFile.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgParameter.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklIBCConfig.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklOBCConfig.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklMBoxConfig.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklMDMConfig.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
PS_E_CBC.o: /usr/local/sphere/include/base/SklFileIOCntl.h
PS_E_CBC.o: /usr/local/sphere/include/config/SklCfgDatas.h
PS_E_CBC.o: /usr/local/sphere/include/base/SklOutputDataSet.h
PS_E_CBC.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
PS_E_CBC.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
PS_E_CBC.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklDataManager.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayBase.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklAray.h
PS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaManager.h
PS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaIF.h
PS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaComponent.h
PS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
PS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
PS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
PS_E_CBC.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklVector2D.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklVector3D.h
PS_E_CBC.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
PS_E_CBC.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
PS_E_CBC.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
PS_E_CBC.o: /usr/local/sphere/include/coord/SklCrdMngr.h
PS_E_CBC.o: /usr/local/sphere/include/coord/SklCoord.h FortranFuncCBC.h
PS_E_CBC.o: /usr/local/sphere/include/SklUtil.h CBC_SetBC.h
IF_TRP_VOF.o: SklSolverCBC.h CBC_Define.h /usr/local/sphere/include/Skl.h
IF_TRP_VOF.o: /usr/local/sphere/include/endianUtil.h
IF_TRP_VOF.o: /usr/local/sphere/include/SklSolverBase.h
IF_TRP_VOF.o: /usr/local/sphere/include/SklBase.h
IF_TRP_VOF.o: /usr/local/sphere/include/SklDefine.h
IF_TRP_VOF.o: /usr/local/sphere/include/SklReservWord.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklSolverConfig.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgDatas.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgBase.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgErrMsg.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgDomainInfo.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgNodeBase.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgBase.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgUsrDefNode.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgSteer.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgSteerInputData.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgInFile.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgSteerOutputData.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgHistory.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgOutFile.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgParameter.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgUsrDefTopNode.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklIBCConfig.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgInnerBnd.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklOBCConfig.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgOuterBnd.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklMBoxConfig.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgMBoxTable.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklMDMConfig.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgMdmTable.h
IF_TRP_VOF.o: /usr/local/sphere/include/base/SklFileIOCntl.h
IF_TRP_VOF.o: /usr/local/sphere/include/config/SklCfgDatas.h
IF_TRP_VOF.o: /usr/local/sphere/include/base/SklOutputDataSet.h
IF_TRP_VOF.o: /usr/local/sphere/include/base/SklFileIO_SPH.h
IF_TRP_VOF.o: /usr/local/sphere/include/base/SklFileIO_FDV.h
IF_TRP_VOF.o: /usr/local/sphere/include/base/SklFileIO_SVX.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklDataManager.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklArayBase.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklScalar.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklArayIdx1.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklAray.h
IF_TRP_VOF.o: /usr/local/sphere/include/parallel/SklParaManager.h
IF_TRP_VOF.o: /usr/local/sphere/include/parallel/SklParaIF.h
IF_TRP_VOF.o: /usr/local/sphere/include/parallel/SklParaComponent.h
IF_TRP_VOF.o: /usr/local/sphere/include/parallel/SklParaManagerEv.h
IF_TRP_VOF.o: /usr/local/sphere/include/parallel/SklParaNodeInfoEv.h
IF_TRP_VOF.o: /usr/local/sphere/include/parallel/SklParaManagerMb.h
IF_TRP_VOF.o: /usr/local/sphere/include/parallel/SklParaNodeInfoMb.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklScalar2D.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklArayIdx2.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklScalar3D.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklArayIdx3.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklScalar4D.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklArayIdx4.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklScalar4DEx.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklVector2D.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklVector2DEx.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklVector3D.h
IF_TRP_VOF.o: /usr/local/sphere/include/dataclass/SklVector3DEx.h
IF_TRP_VOF.o: /usr/local/sphere/include/fileio/SklVoxDataCntl.h
IF_TRP_VOF.o: /usr/local/sphere/include/fileio/SklVoxDataSet.h
IF_TRP_VOF.o: /usr/local/sphere/include/coord/SklCrdMngr.h
IF_TRP_VOF.o: /usr/local/sphere/include/coord/SklCoord.h FortranFuncCBC.h
IF_TRP_VOF.o: /usr/local/sphere/include/SklUtil.h CBC_SetBC.h

###################################################################
#
# CAERU Library
#
# Copyright (c) All right reserved. 2012
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

TARGET = libFB.a

CSRCS = \


CXXSRCS = \
  Control.C \
  FBUtility.C \
  Intrinsic.C \
  TPControl.C \
  ParseBC.C \
  ParseMat.C \
  FileIO.C \
  Component.C \
  BndOuter.C \
  Interval_Mngr.C \
  Alloc.C \
  VoxInfo.C \
  History.C \
  SetBC.C \
  CompoFraction.C \
  dfi.C \
  DataHolder.C \
  Monitor.C \
  MonCompo.C \
  Sampling.C \
  PLOT3D.C \
  PLOT3D_read.C \
  PLOT3D_write.C \



F90SRCS = \
  FB_util.f90 \
  PLT3D.f90 \
  PLT3D_read.f90 \
  PLT3D_write.f90 \

SRCS  = $(CSRCS) $(CXXSRCS) $(F90SRCS)

.SUFFIXES: .o .c .C .f90
COBJS = $(CSRCS:.c=.o)
CXXOBJS = $(CXXSRCS:.C=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS  = $(COBJS) $(CXXOBJS) $(F90OBJS)

$(TARGET):$(OBJS)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.c.o:
	$(CC) $(CFLAGS) $(UDEF_INC_PATH) -o $@ -c $<

.C.o:
	$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) -o $@ -c $<

.f90.o:
	$(F90) $(F90FLAGS) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -I. -Y $(SRCS)

# DO NOT DELETE

Control.o: Control.h DomainInfo.h FB_Define.h mydebug.h Medium.h Component.h
Control.o: FBUtility.h BndOuter.h Interval_Mngr.h TPControl.h PLOT3D_read.h
Control.o: PLOT3D.h FB_Ffunc.h PLOT3D_write.h
FBUtility.o: FBUtility.h FB_Define.h mydebug.h
Intrinsic.o: Intrinsic.h DomainInfo.h FB_Define.h mydebug.h BndOuter.h
Intrinsic.o: Control.h Medium.h Component.h FBUtility.h Interval_Mngr.h
Intrinsic.o: TPControl.h PLOT3D_read.h PLOT3D.h FB_Ffunc.h PLOT3D_write.h
Intrinsic.o: vec3.h vec3_func.h vec3f_func.h
TPControl.o: TPControl.h FB_Define.h mydebug.h
ParseBC.o: ParseBC.h DomainInfo.h FB_Define.h mydebug.h FBUtility.h
ParseBC.o: BndOuter.h Control.h Medium.h Component.h Interval_Mngr.h
ParseBC.o: TPControl.h PLOT3D_read.h PLOT3D.h FB_Ffunc.h PLOT3D_write.h
ParseBC.o: vec3.h vec3_func.h vec3f_func.h Intrinsic.h
ParseMat.o: ParseMat.h FB_Define.h mydebug.h Medium.h Component.h TPControl.h
ParseMat.o: FBUtility.h
FileIO.o: FileIO.h DomainInfo.h FB_Define.h mydebug.h FB_Ffunc.h
Component.o: Component.h FB_Define.h mydebug.h
BndOuter.o: BndOuter.h FB_Define.h mydebug.h
Interval_Mngr.o: Interval_Mngr.h FB_Define.h mydebug.h
Alloc.o: Alloc.h
VoxInfo.o: VoxInfo.h DomainInfo.h FBUtility.h FB_Define.h mydebug.h
VoxInfo.o: Component.h Medium.h SetBC.h BndOuter.h Control.h Interval_Mngr.h
VoxInfo.o: TPControl.h PLOT3D_read.h PLOT3D.h FB_Ffunc.h PLOT3D_write.h
VoxInfo.o: Intrinsic.h vec3.h vec3_func.h vec3f_func.h CompoFraction.h
VoxInfo.o: FindexS3D.h
History.o: History.h Control.h DomainInfo.h FB_Define.h mydebug.h Medium.h
History.o: Component.h FBUtility.h BndOuter.h Interval_Mngr.h TPControl.h
History.o: PLOT3D_read.h PLOT3D.h FB_Ffunc.h PLOT3D_write.h
SetBC.o: SetBC.h DomainInfo.h FB_Define.h mydebug.h BndOuter.h Control.h
SetBC.o: Medium.h Component.h FBUtility.h Interval_Mngr.h TPControl.h
SetBC.o: PLOT3D_read.h PLOT3D.h FB_Ffunc.h PLOT3D_write.h Intrinsic.h vec3.h
SetBC.o: vec3_func.h vec3f_func.h
CompoFraction.o: CompoFraction.h FB_Define.h mydebug.h vec3.h vec3_func.h
CompoFraction.o: vec3f_func.h
dfi.o: dfi.h FB_Define.h mydebug.h util_Path.h
DataHolder.o: DataHolder.h FB_Define.h mydebug.h TPControl.h DomainInfo.h
Monitor.o: Monitor.h DomainInfo.h FB_Define.h mydebug.h vec3.h vec3_func.h
Monitor.o: vec3f_func.h MonCompo.h basic_func.h Component.h Sampling.h
Monitor.o: FBUtility.h TPControl.h Control.h Medium.h BndOuter.h
Monitor.o: Interval_Mngr.h PLOT3D_read.h PLOT3D.h FB_Ffunc.h PLOT3D_write.h
MonCompo.o: MonCompo.h DomainInfo.h FB_Define.h mydebug.h vec3.h vec3_func.h
MonCompo.o: vec3f_func.h basic_func.h Component.h Sampling.h FBUtility.h
Sampling.o: Sampling.h FB_Define.h mydebug.h FBUtility.h vec3.h vec3_func.h
Sampling.o: vec3f_func.h basic_func.h
PLOT3D.o: PLOT3D.h DomainInfo.h FB_Define.h mydebug.h FB_Ffunc.h
PLOT3D_read.o: PLOT3D_read.h PLOT3D.h DomainInfo.h FB_Define.h mydebug.h
PLOT3D_read.o: FB_Ffunc.h
PLOT3D_write.o: PLOT3D_write.h PLOT3D.h DomainInfo.h FB_Define.h mydebug.h
PLOT3D_write.o: FB_Ffunc.h

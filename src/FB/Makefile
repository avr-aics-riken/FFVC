###################################################################
#
# FFV : Frontflow / violet Cartesian
#
# Copyright (c) 2012 All right reserved.
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
  dfiinfo.C \
  DataHolder.C \
  Monitor.C \
  MonCompo.C \
  Sampling.C \



F90SRCS = \
  FB_util.f90 \


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
Control.o: FBUtility.h BndOuter.h Interval_Mngr.h TPControl.h
FBUtility.o: FBUtility.h FB_Define.h mydebug.h
Intrinsic.o: Intrinsic.h DomainInfo.h FB_Define.h mydebug.h BndOuter.h
Intrinsic.o: Control.h Medium.h Component.h FBUtility.h Interval_Mngr.h
Intrinsic.o: TPControl.h vec3.h vec3_func.h vec3f_func.h
TPControl.o: TPControl.h FB_Define.h mydebug.h
ParseBC.o: ParseBC.h DomainInfo.h FB_Define.h mydebug.h FBUtility.h
ParseBC.o: BndOuter.h Control.h Medium.h Component.h Interval_Mngr.h
ParseBC.o: TPControl.h vec3.h vec3_func.h vec3f_func.h Intrinsic.h
ParseMat.o: ParseMat.h FB_Define.h mydebug.h Medium.h Component.h TPControl.h
ParseMat.o: FBUtility.h
FileIO.o: FileIO.h DomainInfo.h FB_Define.h mydebug.h FB_Ffunc.h FBUtility.h
Component.o: Component.h FB_Define.h mydebug.h
BndOuter.o: BndOuter.h FB_Define.h mydebug.h
Interval_Mngr.o: Interval_Mngr.h FB_Define.h mydebug.h
Alloc.o: Alloc.h FB_Define.h mydebug.h
VoxInfo.o: VoxInfo.h DomainInfo.h FB_Define.h mydebug.h FBUtility.h
VoxInfo.o: Component.h Medium.h SetBC.h BndOuter.h Control.h Interval_Mngr.h
VoxInfo.o: TPControl.h Intrinsic.h vec3.h vec3_func.h vec3f_func.h
VoxInfo.o: CompoFraction.h FindexS3D.h fill_bid.h fill_mid.h
History.o: History.h Control.h DomainInfo.h FB_Define.h mydebug.h Medium.h
History.o: Component.h FBUtility.h BndOuter.h Interval_Mngr.h TPControl.h
SetBC.o: SetBC.h DomainInfo.h FB_Define.h mydebug.h BndOuter.h Control.h
SetBC.o: Medium.h Component.h FBUtility.h Interval_Mngr.h TPControl.h
SetBC.o: Intrinsic.h vec3.h vec3_func.h vec3f_func.h
CompoFraction.o: CompoFraction.h FB_Define.h mydebug.h vec3.h vec3_func.h
CompoFraction.o: vec3f_func.h
dfi.o: dfi.h FB_Define.h mydebug.h util_Path.h
dfiinfo.o: dfiinfo.h DomainInfo.h FB_Define.h mydebug.h TPControl.h
DataHolder.o: DataHolder.h FB_Define.h mydebug.h TPControl.h DomainInfo.h
Monitor.o: Monitor.h DomainInfo.h FB_Define.h mydebug.h vec3.h vec3_func.h
Monitor.o: vec3f_func.h MonCompo.h basic_func.h Component.h Sampling.h
Monitor.o: FBUtility.h TPControl.h Control.h Medium.h BndOuter.h
Monitor.o: Interval_Mngr.h
MonCompo.o: MonCompo.h DomainInfo.h FB_Define.h mydebug.h vec3.h vec3_func.h
MonCompo.o: vec3f_func.h basic_func.h Component.h Sampling.h FBUtility.h
Sampling.o: Sampling.h FB_Define.h mydebug.h vec3.h vec3_func.h vec3f_func.h
Sampling.o: basic_func.h

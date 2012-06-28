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
#       SetBC.C \
#       Monitor.C \
#       MonCompo.C \
#       Sampling.C \
#       DataHolder.C \
#       dfi.C \
#       CompoFraction.C \

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
FileIO.o: FileIO.h DomainInfo.h FB_Define.h mydebug.h FBUtility.h FB_Ffunc.h
Component.o: Component.h FB_Define.h mydebug.h
BndOuter.o: BndOuter.h FB_Define.h mydebug.h
Interval_Mngr.o: Interval_Mngr.h FB_Define.h mydebug.h
Alloc.o: Alloc.h
VoxInfo.o: VoxInfo.h DomainInfo.h FBUtility.h FB_Define.h mydebug.h
VoxInfo.o: Component.h Medium.h SetBC.h BndOuter.h Control.h Interval_Mngr.h
VoxInfo.o: TPControl.h Intrinsic.h vec3.h vec3_func.h vec3f_func.h
History.o: History.h Control.h DomainInfo.h FB_Define.h mydebug.h Medium.h
History.o: Component.h FBUtility.h BndOuter.h Interval_Mngr.h TPControl.h

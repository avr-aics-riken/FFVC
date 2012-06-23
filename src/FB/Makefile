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
       TPControl.C 
#       ParseBC.C \
#       ParseMat.C \
#       SetBC.C \
#       Monitor.C \
#       MonCompo.C \
#       Interval_Mngr.C \
#       Sampling.C \
#       DataHolder.C \
#       dfi.C \
#       FileIO.C \
#       History.C \
#       Alloc.C \
#       VoxInfo.C \
#       BndOuter.C \
#       Component.C \
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

Control.o: Control.h FB_Define.h mydebug.h Medium.h Component.h FBUtility.h
Control.o: BndOuter.h Interval_Mngr.h TPControl.h
FBUtility.o: FBUtility.h FB_Define.h mydebug.h
Intrinsic.o: Intrinsic.h FB_Define.h mydebug.h BndOuter.h Control.h Medium.h
Intrinsic.o: Component.h FBUtility.h Interval_Mngr.h TPControl.h vec3.h
Intrinsic.o: vec3_func.h vec3f_func.h
TPControl.o: TPControl.h FB_Define.h mydebug.h

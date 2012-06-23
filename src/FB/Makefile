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
#       Alloc.C \
#       VoxInfo.C \
#       BndOuter.C \
#       Component.C \
#       CompoFraction.C \
       Control.C \
       FBUtility.C \
#       FileIO.C \
#       History.C \
       Intrinsic.C \
#       ParseBC.C \
#       ParseMat.C \
#       SetBC.C \
#       Monitor.C \
#       MonCompo.C \
#       Interval_Mngr.C \
#       Sampling.C \
#       DataHolder.C \
#       dfi.C \
       TPControl.C \


FCSRCS = \


F90SRCS = \
       FB_util.f90 \


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

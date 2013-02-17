###################################################################
#
# CAERU Library
#
# Copyright (c) 2012-2013  All right reserved.
#
# Institute of Industrial Science, University of Tokyo, Japan. 
#
###################################################################

TARGET = libIP.a

CSRCS = \


CXXSRCS = \
	IP_PMT.C \
	IP_Rect.C \
	IP_Polygon.C \
	IP_Cylinder.C \
	IP_Step.C \
	IP_PPLT2D.C \
	IP_Duct.C \
	IP_SHC1D.C \
	IP_Sphere.C \
  IP_Jet.C \


F90SRCS = \


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

IP_PMT.o: IP_PMT.h IP_Define.h
IP_Rect.o: IP_Rect.h IP_Define.h
IP_Polygon.o: IP_Polygon.h IP_Define.h
IP_Cylinder.o: IP_Cylinder.h IP_Define.h
IP_Step.o: IP_Step.h IP_Define.h
IP_PPLT2D.o: IP_PPLT2D.h IP_Define.h
IP_Duct.o: IP_Duct.h IP_Define.h
IP_SHC1D.o: IP_SHC1D.h IP_Define.h
IP_Sphere.o: IP_Sphere.h IP_Define.h
IP_Jet.o: IP_Jet.h IP_Define.h

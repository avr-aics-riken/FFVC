#
# SPHERE - Skeleton for PHysicalSvCube and Engineering REsearch
#
# Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
#

TARGET = libsphnsolvIP.a

CSRCS = \


CXXSRCS = \
       ../../IP/IP_PPLT2D.C \
       ../../IP/IP_Duct.C \
       ../../IP/IP_SHC1D.C \
       ../../IP/IP_PMT.C \
       ../../IP/IP_Rect.C \
       ../../IP/IP_Cylinder.C \
       ../../IP/IP_Step.C \
       ../../IP/IP_Polygon.C \
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

../../IP/IP_PPLT2D.o: ../../IP/IP_PPLT2D.h ../../IP/IP_Define.h
../../IP/IP_Duct.o: ../../IP/IP_Duct.h ../../IP/IP_Define.h
../../IP/IP_SHC1D.o: ../../IP/IP_SHC1D.h ../../IP/IP_Define.h
../../IP/IP_PMT.o: ../../IP/IP_PMT.h ../../IP/IP_Define.h
../../IP/IP_Rect.o: ../../IP/IP_Rect.h ../../IP/IP_Define.h
../../IP/IP_Cylinder.o: ../../IP/IP_Cylinder.h ../../IP/IP_Define.h
../../IP/IP_Step.o: ../../IP/IP_Step.h ../../IP/IP_Define.h
../../IP/IP_Polygon.o: ../../IP/IP_Polygon.h ../../IP/IP_Define.h
../../IP/IP_Sphere.o: ../../IP/IP_Sphere.h ../../IP/IP_Define.h

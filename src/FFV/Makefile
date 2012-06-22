###################################################################
#
# CAERU Library
#
# Copyright (c) All right reserved. 2012
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

TARGET = ffv

RM = \rm -f

###################################################################
#
# CPM_LIBS;  openmpi, textparser, cpmlibのライブラリパスとライブラリの指定
# CPM_FLAGS; openmpi, textparser, cpmlibのインクルードパスの指定
#
#
#
#
#

OMP_FLAGS=-openmp -openmp-report=2

CC=icc
CFLAGS=-O3
CXX=icpc
CXXFLAGS=-O3 $(OMP_FLAGS) -DTD_USE_NAMESPACE
FC=ifort
FCFLAGS=-O3
F90=ifort
F90FLAGS=-O3 -Warn unused -fpp $(OMP_FLAGS) -D_STATIC
LDFLAGS=-L/opt/intel/composerxe/lib
LIBS=-lifport -lifcore

UDEF_INC_PATH=-I../../Cutlib-2.0.3/include -I../../Polylib_2_0_3/include -I../../PMlib-1.0/include
UDEF_LIB_PATH=-L../../Polylib_2_0_3/lib -L../../Cutlib-2.0.3/lib -L../../PMlib-1.0/lib
UDEF_LIBS=-lMPIPolylib -lcutlib -lpmlib

CPM_LIBS = '/usr/local/cpm/bin/cpm-config --cflags'
CPM_FLAGS = '/usr/local/cpm/bin/cpm-config --libs'

REALOPT=float
SPH_EXTERNAL_HEADER_PATH=-I../../FB -I../../F_CBC -I../../IP 

###################################################################

CSRCS = \


CXXSRCS = \
       main.C \


FCSRCS = \


F90SRCS = \


SRCS  = $(CSRCS) $(CXXSRCS) $(FCSRCS) $(F90SRCS)

.SUFFIXES: .o .c .C .f .f90
COBJS = $(CSRCS:.c=.o)
CXXOBJS = $(CXXSRCS:.C=.o)
FCOBJS = $(FCSRCS:.f=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS  = $(COBJS) $(CXXOBJS) $(FCOBJS) $(F90OBJS)


###################################################################

$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $(CPM_FLAGS) $(UDEF_INC_PATH) $(REALOPT) -o $(TARGET) $(OBJS) \
	$(LDFLAGS) $(UDEF_LIB_PATH) $(CPM_LIBS) $(LIBS) $(UDEF_LIBS)
#	mv $(TARGET) ../bin/.

.c.o:
	$(CC) $(CFLAGS) $(CPM_FLAGS) $(UDEF_INC_PATH) $(REALOPT) -o $@ -c $<

.C.o:
	$(CXX) $(CXXFLAGS) $(CPM_FLAGS) $(UDEF_INC_PATH) $(REALOPT) -o $@ -c $<

.f.o:
	$(FC) $(FCFLAGS) $(CPM_FLAGS) -o $@ -c $<

.f90.o:
	$(F90) $(F90FLAGS) $(CPM_FLAGS) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -I. -Y $(CXXSRCS) $(CSRCS)

# DO NOT DELETE

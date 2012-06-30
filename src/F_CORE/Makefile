###################################################################
#
# CAERU Library
#
# Copyright (c) All right reserved. 2012
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

TARGET = libFFVcore.a

CSRCS = \


CXXSRCS = \

F90SRCS = \
	ffv_forcing.f90 \
	ffv_pbc.f90 \
	ffv_utility.f90 \
	ffv_vbc_inner.f90 \
	ffv_vbc_outer.f90 \
	ffv_velocity_binary.f90 \
  ffv_pscalar.f90 \
#ffv_poisson.f90 \



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
	$(F90) $(F90FLAGS) -I. -I../FB -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -I. -I../FB -I../FFV -Y $(SRCS)

# DO NOT DELETE

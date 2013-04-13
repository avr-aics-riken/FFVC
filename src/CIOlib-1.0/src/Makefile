###################################################################
#
# CIOlib - Cartesian Input / Output library
#
# Copyright (c) RIKEN, All right reserved. 2013
#
###################################################################


TARGET = ../lib/libCIO.a

CSRCS = \

CXXSRCS = \
  cio_DFI.C \
  cio_TextParser.C \
  cio_DFI_SPH.C \
  cio_DFI_BOV.C \
  cio_DFI_Write.C \

F90SRCS = \
  cio_DFI_COARSE.f90 \
  cio_DFI_util.f90 \


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
	$(F90) $(F90FLAGS) $(UDEF_INC_PATH) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -I. -Y $(SRCS)

# DO NOT DELETE

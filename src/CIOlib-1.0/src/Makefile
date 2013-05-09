###################################################################
#
# CIOlib - Cartesian Input / Output library
#
# Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
# All right reserved.
#
###################################################################

TARGET = ../lib/libCIO.a

CSRCS = \

CXXSRCS = \
  cio_DFI.C \
  cio_TextParser.C \
  cio_DFI_SPH.C \
  cio_DFI_BOV.C \
  cio_DFI_Write.C

F90SRCS = \
  cio_DFI_COARSE.f90 \
  cio_DFI_util.f90


SRCS  = $(CSRCS) $(CXXSRCS) $(F90SRCS)


.SUFFIXES: .o .c .C .f90
COBJS = $(CSRCS:.c=.o)
CXXOBJS = $(CXXSRCS:.C=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS  = $(COBJS) $(CXXOBJS) $(F90OBJS)


$(TARGET):$(OBJS)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)


.C.o:
	$(CXX) $(CXXFLAGS) -I../include $(UDEF_INC_PATH) -o $@ -c $<

.f90.o:
	$(F90) $(F90FLAGS) -I. -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -I. -I../include -Y $(SRCS)

# DO NOT DELETE

cio_DFI.o: ../include/cio_DFI.h ../include/cio_TextParser.h
cio_DFI.o: ../include/cio_Define.h ../include/cio_DFI_Ffunc.h
cio_DFI.o: ../include/cio_DFI_SPH.h ../include/cio_DFI_BOV.h
cio_TextParser.o: ../include/cio_TextParser.h ../include/cio_Define.h
cio_DFI_SPH.o: ../include/cio_DFI.h ../include/cio_TextParser.h
cio_DFI_SPH.o: ../include/cio_Define.h ../include/cio_DFI_Ffunc.h
cio_DFI_SPH.o: ../include/cio_DFI_SPH.h
cio_DFI_BOV.o: ../include/cio_DFI.h ../include/cio_TextParser.h
cio_DFI_BOV.o: ../include/cio_Define.h ../include/cio_DFI_Ffunc.h
cio_DFI_BOV.o: ../include/cio_DFI_BOV.h
cio_DFI_Write.o: ../include/cio_DFI.h ../include/cio_TextParser.h
cio_DFI_Write.o: ../include/cio_Define.h ../include/cio_DFI_Ffunc.h

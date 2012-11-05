###################################################################
#
# FFV : Frontflow / violet Cartesian
#
# Copyright (c) 2012 All right reserved.
#
# Institute of Industrial Science, University of Tokyo, Japan. 
#
###################################################################

TARGET = libPLOT3D.a

CSRCS = \


CXXSRCS = \
  PLOT3D.C \
  PLOT3D_read.C \
  PLOT3D_write.C \
  ffv_PLOT3D.C \


F90SRCS = \
  PLT3D.f90 \
  PLT3D_read.f90 \
  PLT3D_write.f90 \

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

PLOT3D.o: PLOT3D.h
PLOT3D_read.o: PLOT3D_read.h PLOT3D.h
PLOT3D_write.o: PLOT3D_write.h PLOT3D.h
ffv_PLOT3D.o: ffv_PLOT3D.h PLOT3D_write.h PLOT3D.h PLOT3D_read.h

SHELL=/bin/sh
TARGET = c++sample

CPMLIBDIR = ../../CPMlib
include $(CPMLIBDIR)/Makefile.inc

CXXSRCS = main.C \


F77SRCS = \


F90SRCS = \


.SUFFIXES: .C .f .f90
CXXOBJS = $(CXXSRCS:.C=.o)
F77OBJS = $(F77SRCS:.f=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS    = $(CXXOBJS) $(F77OBJS) $(F90OBJS)

$(TARGET):	$(OBJS)
	$(CXX) -o $(@) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

.C.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.f.o:
	$(F77) $(F77FLAGS) -c -o $@ $<

.f90.o:
	$(F90) $(F90FLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJS) $(TARGET) *.mod

depend:
	makedepend -Y $(CXXFLAGS) $(CXXSRCS)

rebuild:	clean $(TARGET)


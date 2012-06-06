SHELL=/bin/sh
TARGET = mconvp

F90        = ifort
F90FLAGS   = -O3
LDFLAGS    =
LIBS       =

F90SRCS = mconvp.f90\


.SUFFIXES: .f90
F90OBJS = $(F90SRCS:.f90=.o)
OBJS    = $(F90OBJS)

$(TARGET):	$(OBJS)
	$(F90) -o $(@) $(F90FLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

.f90.o:
	$(F90) $(F90FLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJS) $(TARGET) *.mod

depend:
	makedepend -Y $(CXXFLAGS) $(CXXSRCS)

rebuild:	clean $(TARGET)


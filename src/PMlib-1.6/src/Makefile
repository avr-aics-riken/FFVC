###################################################################
#
# PMlib - Performance Monitor library
#
# Copyright (c) RIKEN. All right reserved. 2010-2013
#
# 2010-2011  VCAD System Research Program, RIKEN
# 2012-2013  Advanced Institute for Computational Science (AICS), RIKEN
# 
###################################################################

PMLIB = ../lib/libpmlib.a

all: $(PMLIB)

OBJS = PerfWatch.o PerfMonitor.o


$(PMLIB): $(OBJS)
	$(AR) $@ $(OBJS) 
	$(RANLIB) $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -I../include $(UDEF_INC_PATH) -c $<

clean:
	rm -rf *.o $(PMLIB)

depend:
	makedepend -I. -I../include -Y PerfWatch.cpp PerfMonitor.cpp
  
# DO NOT DELETE

PerfWatch.o: ../include/PerfWatch.h
PerfMonitor.o: ../include/PerfMonitor.h ../include/PerfWatch.h

#####################################################################################
#
#  環境に合わせて、make_settingマクロを編集後、make
#
#####################################################################################

PMLIB = ../lib/libpmlib.a

all: $(PMLIB)

OBJS = PerfWatch.o PerfMonitor.o


# INTEL
CFLAGS = -O3 -Wall -I../include -I$(MPI_DIR)/include
#CFLAGS = -O3 -I../include -I$(MPI_DIR)/include

# FX10
#CFLAGS = -Kfast -I../include



$(PMLIB): $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

clean:
	rm -rf *.o $(PMLIB)

depend:
	makedepend -I. -I../include -I$(MPI_DIR)/include -Y PerfWatch.cpp PerfMonitor.cpp
# DO NOT DELETE

PerfWatch.o: ../include/PerfWatch.h
PerfMonitor.o: ../include/PerfMonitor.h ../include/PerfWatch.h

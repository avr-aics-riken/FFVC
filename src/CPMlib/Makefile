#####################################################################################
# Makefile for CPM library
#####################################################################################
SHELL=/bin/sh
CPMLIB = ./lib/libcpmlib.a

include Makefile.inc
CXXFLAGS += -I./include

CXXSRCS = src/cpm_DomainInfo.cpp \
          src/cpm_VoxelInfo.cpp \
          src/cpm_ParaManager.cpp \
          src/cpm_ParaManager_MPI.cpp \
          src/cpm_ParaManager_Alloc.cpp \
          src/cpm_ParaManager_frtIF.cpp \
          src/cpm_TextParser.cpp \
          src/cpm_TextParserDomain.cpp \

.SUFFIXES: .cpp
CXXOBJS = $(CXXSRCS:.cpp=.o)
OBJS    = $(CXXOBJS)


$(CPMLIB):	$(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJS) $(CPMLIB)

depend:
	makedepend -Y $(CXXFLAGS) $(CXXSRCS)

rebuild:	clean $(CPMLIB)


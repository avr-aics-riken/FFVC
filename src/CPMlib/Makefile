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

# DO NOT DELETE

src/cpm_DomainInfo.o: ./include/cpm_DomainInfo.h ./include/cpm_Base.h
src/cpm_DomainInfo.o: ./include/cpm_Define.h /opt/openmpi/include/mpi.h
src/cpm_DomainInfo.o: ./include/cpm_Version.h
src/cpm_VoxelInfo.o: ./include/cpm_VoxelInfo.h ./include/cpm_Base.h
src/cpm_VoxelInfo.o: ./include/cpm_Define.h /opt/openmpi/include/mpi.h
src/cpm_VoxelInfo.o: ./include/cpm_Version.h ./include/cpm_DomainInfo.h
src/cpm_ParaManager.o: ./include/cpm_ParaManager.h ./include/cpm_Base.h
src/cpm_ParaManager.o: ./include/cpm_Define.h /opt/openmpi/include/mpi.h
src/cpm_ParaManager.o: ./include/cpm_Version.h ./include/cpm_DomainInfo.h
src/cpm_ParaManager.o: ./include/cpm_VoxelInfo.h ./include/cpm_ObjList.h
src/cpm_ParaManager.o: ./include/inline/cpm_ParaManager_inline.h
src/cpm_ParaManager.o: ./include/inline/cpm_ParaManager_BndComm.h
src/cpm_ParaManager.o: ./include/inline/cpm_ParaManager_BndCommEx.h
src/cpm_ParaManager_MPI.o: ./include/cpm_ParaManager.h ./include/cpm_Base.h
src/cpm_ParaManager_MPI.o: ./include/cpm_Define.h /opt/openmpi/include/mpi.h
src/cpm_ParaManager_MPI.o: ./include/cpm_Version.h ./include/cpm_DomainInfo.h
src/cpm_ParaManager_MPI.o: ./include/cpm_VoxelInfo.h ./include/cpm_ObjList.h
src/cpm_ParaManager_MPI.o: ./include/inline/cpm_ParaManager_inline.h
src/cpm_ParaManager_MPI.o: ./include/inline/cpm_ParaManager_BndComm.h
src/cpm_ParaManager_MPI.o: ./include/inline/cpm_ParaManager_BndCommEx.h
src/cpm_ParaManager_Alloc.o: ./include/cpm_ParaManager.h ./include/cpm_Base.h
src/cpm_ParaManager_Alloc.o: ./include/cpm_Define.h
src/cpm_ParaManager_Alloc.o: /opt/openmpi/include/mpi.h
src/cpm_ParaManager_Alloc.o: ./include/cpm_Version.h
src/cpm_ParaManager_Alloc.o: ./include/cpm_DomainInfo.h
src/cpm_ParaManager_Alloc.o: ./include/cpm_VoxelInfo.h
src/cpm_ParaManager_Alloc.o: ./include/cpm_ObjList.h
src/cpm_ParaManager_Alloc.o: ./include/inline/cpm_ParaManager_inline.h
src/cpm_ParaManager_Alloc.o: ./include/inline/cpm_ParaManager_BndComm.h
src/cpm_ParaManager_Alloc.o: ./include/inline/cpm_ParaManager_BndCommEx.h
src/cpm_ParaManager_frtIF.o: ./include/cpm_ParaManager.h ./include/cpm_Base.h
src/cpm_ParaManager_frtIF.o: ./include/cpm_Define.h
src/cpm_ParaManager_frtIF.o: /opt/openmpi/include/mpi.h
src/cpm_ParaManager_frtIF.o: ./include/cpm_Version.h
src/cpm_ParaManager_frtIF.o: ./include/cpm_DomainInfo.h
src/cpm_ParaManager_frtIF.o: ./include/cpm_VoxelInfo.h
src/cpm_ParaManager_frtIF.o: ./include/cpm_ObjList.h
src/cpm_ParaManager_frtIF.o: ./include/inline/cpm_ParaManager_inline.h
src/cpm_ParaManager_frtIF.o: ./include/inline/cpm_ParaManager_BndComm.h
src/cpm_ParaManager_frtIF.o: ./include/inline/cpm_ParaManager_BndCommEx.h
src/cpm_TextParser.o: ./include/cpm_TextParser.h ./include/cpm_Base.h
src/cpm_TextParser.o: ./include/cpm_Define.h /opt/openmpi/include/mpi.h
src/cpm_TextParser.o: ./include/cpm_Version.h
src/cpm_TextParser.o: /my_apps/cbc_git/src/TP/include/TextParser.h
src/cpm_TextParser.o: /my_apps/cbc_git/src/TP/include/TextParserCommon.h
src/cpm_TextParserDomain.o: ./include/cpm_TextParserDomain.h
src/cpm_TextParserDomain.o: ./include/cpm_TextParser.h ./include/cpm_Base.h
src/cpm_TextParserDomain.o: ./include/cpm_Define.h /opt/openmpi/include/mpi.h
src/cpm_TextParserDomain.o: ./include/cpm_Version.h
src/cpm_TextParserDomain.o: /my_apps/cbc_git/src/TP/include/TextParser.h
src/cpm_TextParserDomain.o: /my_apps/cbc_git/src/TP/include/TextParserCommon.h
src/cpm_TextParserDomain.o: ./include/cpm_DomainInfo.h

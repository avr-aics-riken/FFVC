###
###########################################
# CPM - Cartesian Partition Manager       #
#                                         #
# Copyright (C) 2012 University of Tokyo. #
###########################################
###
lib_LIBRARIES = libcpmlib.a

libcpmlib_a_CXXFLAGS = @MPICH_CFLAGS@ @TP_CFLAGS@ @REALOPT@ -I$(top_builddir)/include

libcpmlib_a_SOURCES = \
    cpm_DomainInfo.cpp \
    cpm_VoxelInfo.cpp \
    cpm_ParaManager.cpp \
    cpm_ParaManager_MPI.cpp \
    cpm_ParaManager_Alloc.cpp \
    cpm_ParaManager_frtIF.cpp \
    cpm_TextParser.cpp \
    cpm_TextParserDomain.cpp

cpmlib_includedir = $(includedir)
cpmlib_include_HEADERS = \
    $(top_builddir)/include/cpm_Base.h \
    $(top_builddir)/include/cpm_Define.h \
    $(top_builddir)/include/cpm_DomainInfo.h \
    $(top_builddir)/include/cpm_ObjList.h \
    $(top_builddir)/include/cpm_ParaManager.h \
    $(top_builddir)/include/cpm_TextParser.h \
    $(top_builddir)/include/cpm_TextParserDomain.h \
    $(top_builddir)/include/cpm_Version.h \
    $(top_builddir)/include/cpm_VoxelInfo.h \
    $(top_builddir)/include/cpm_fparam.fi

cpmlib_inline_includedir = $(includedir)/inline
cpmlib_inline_include_HEADERS = \
    $(top_builddir)/include/inline/cpm_ParaManager_BndComm.h \
    $(top_builddir)/include/inline/cpm_ParaManager_BndCommEx.h \
    $(top_builddir)/include/inline/cpm_ParaManager_inline.h


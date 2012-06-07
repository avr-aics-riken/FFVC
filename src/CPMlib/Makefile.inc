#####################################################################################
# Makefile settings for CPM library
#####################################################################################

# library settings
TEXTPARSER = ../TP
MPICHDIR   = /opt/openmpi
MPICHLIB   = -lmpi -lmpi_f77
FRT_R8     = -r8
REAL_TYPE  = FLOAT
#REAL_TYPE  = DOUBLE

# compiler settings
CXX        = icpc
CXXFLAGS   = -O3
F77        = ifort
F77FLAGS   = -O3
F90        = ifort
F90FLAGS   = -O3
LDFLAGS    =
LIBS       =

AR     = ar cru
RANLIB = ranlib

# other settings
CXXFLAGS += -D_DEBUG
#CXXFLAGS += -D_CPM_NO_INLINE_
#CXXFLAGS += -D_BUFSIZE_DOUBLE_
#CXXFLAGS += -D_BUFSIZE_LONG_DOUBLE_

ifeq ("$(REAL_TYPE)","DOUBLE")
  CXXFLAGS += -D_REAL_IS_DOUBLE_
  F77FLAGS += $(FRT_R8)
  F90FLAGS += $(FRT_R8)
endif

ifneq ("$(MPICHDIR)_","_")
  CXXFLAGS += -I$(MPICHDIR)/include
  F77FLAGS += -I$(MPICHDIR)/include
  F90FLAGS += -I$(MPICHDIR)/include
  LDFLAGS  += -L$(MPICHDIR)/lib
  LIBS     += $(MPICHLIB)
endif

ifneq ("$(CPMLIBDIR)_","_")
  CXXFLAGS += -I$(CPMLIBDIR)/include
  F77FLAGS += -I$(CPMLIBDIR)/include
  F90FLAGS += -I$(CPMLIBDIR)/include
  LDFLAGS  += -L$(CPMLIBDIR)/lib
  LIBS     += -lcpmlib
endif

ifneq ("$(TEXTPARSER)_","_")
  CXXFLAGS += -I$(TEXTPARSER)/include
  LDFLAGS  += -L$(TEXTPARSER)/lib
#  LIBS     += -lTextParser
  LIBS     += $(TEXTPARSER)/lib/libTextParser.a
endif

LIBS += -lstdc++


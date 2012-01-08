#
# SPHERE - Skeleton for PHysical and Engineering REsearch
#
# Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
#

BASEDIR  = $(PWD)
TARGET   = sphere

include project_local_settings

ifeq ("$(SPH_PARA_MODULE)","MPI")
  SPHERE_LIBS += -lsphparampi
else
  MPICH_DIR =
  MPICH_CFLAGS =
  MPICH_LDFLAGS =
  MPICH_LIBS =
  SPHERE_LIBS += -lsphparadmy
endif

ifeq ("$(SPH_CSM_MODULE)","CSM")
  SPHERE_LIBS += -lcsm
  SPHERE_LIBS += $(CSM_LIBS)
  SPHERE_CFLAGS += $(CSM_CFLAGS)
  SPHERE_LDFLAGS += $(CSM_LDFLAGS)
else
  SPHERE_LIBS += -lcsmdmy
endif

# depends on the platform
AR = ar cru
RANLIB = ranlib
RM = \rm -f

# for app
SPH_SOLV_FLAGS  = 
SPH_SOLV_LDFAGS = 
SPH_SOLV_LIBS   = 

SPH_SOLV_FLAGS  += -I../CBC
SPH_SOLV_LDFAGS += -L../CBC
SPH_SOLV_LIBS   += -lsphsolvCBC
SPH_SOLV_FLAGS  += -I../IP
SPH_SOLV_LDFAGS += -L../IP
SPH_SOLV_LIBS   += -lsphnsolvIP
SPH_SOLV_FLAGS  += -I../FB
SPH_SOLV_LDFAGS += -L../FB
SPH_SOLV_LIBS   += -lsphnsolvFB

SPH_SOLV_FLAGS += $(SPH_EXTERNAL_HEADER_PATH)


# for User DEFined settings
SPH_SOLV_FLAGS += $(UDEF_INC_PATH)
TMP_LIBS = $(UDEF_LIB_UPPER)
TMP_LIBS += $(SPHERE_LIBS)
TMP_LIBS += $(UDEF_LIB_LOWER)
SPHERE_LIBS := $(TMP_LIBS)
SPHERE_LDFLAGS += $(UDEF_LIB_PATH)
SPHERE_CFLAGS += $(UDEF_OPT)

all:
	( \
	cd CBC; \
	make \
		SPHEREDIR='$(SPHEREDIR)' \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS) $(SPH_SOLV_FLAGS)' \
		CC='$(CC)' \
		CFLAGS='$(CFLAGS) $(SPH_SOLV_FLAGS)' \
		FC='$(FC)' \
		FCFLAGS='$(FCFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		LDFLAGS='$(LDFLAGS)' \
		MPICH_DIR='$(MPICH_DIR)' \
		MPICH_CFLAGS='$(MPICH_CFLAGS)' \
		MPICH_LDFLAGS='$(MPICH_LDFLAGS)' \
		MPICH_LIBS='$(MPICH_LIBS)' \
		XML2FLAGS='$(XML2FLAGS)' \
		XML2LIBS='$(XML2LIBS)' \
		REALOPT='$(REALOPT)' \
		SPHERE_CFLAGS='$(SPHERE_CFLAGS)' \
		SPH_DEVICE='$(SPH_DEVICE)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
	)
	( \
	cd IP; \
	make \
		SPHEREDIR='$(SPHEREDIR)' \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS) $(SPH_SOLV_FLAGS)' \
		CC='$(CC)' \
		CFLAGS='$(CFLAGS) $(SPH_SOLV_FLAGS)' \
		FC='$(FC)' \
		FCFLAGS='$(FCFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		LDFLAGS='$(LDFLAGS)' \
		MPICH_DIR='$(MPICH_DIR)' \
		MPICH_CFLAGS='$(MPICH_CFLAGS)' \
		MPICH_LDFLAGS='$(MPICH_LDFLAGS)' \
		MPICH_LIBS='$(MPICH_LIBS)' \
		XML2FLAGS='$(XML2FLAGS)' \
		XML2LIBS='$(XML2LIBS)' \
		REALOPT='$(REALOPT)' \
		SPHERE_CFLAGS='$(SPHERE_CFLAGS)' \
		SPH_DEVICE='$(SPH_DEVICE)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
	)
	( \
	cd FB; \
	make \
		SPHEREDIR='$(SPHEREDIR)' \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS) $(SPH_SOLV_FLAGS)' \
		CC='$(CC)' \
		CFLAGS='$(CFLAGS) $(SPH_SOLV_FLAGS)' \
		FC='$(FC)' \
		FCFLAGS='$(FCFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		LDFLAGS='$(LDFLAGS)' \
		MPICH_DIR='$(MPICH_DIR)' \
		MPICH_CFLAGS='$(MPICH_CFLAGS)' \
		MPICH_LDFLAGS='$(MPICH_LDFLAGS)' \
		MPICH_LIBS='$(MPICH_LIBS)' \
		XML2FLAGS='$(XML2FLAGS)' \
		XML2LIBS='$(XML2LIBS)' \
		REALOPT='$(REALOPT)' \
		SPHERE_CFLAGS='$(SPHERE_CFLAGS)' \
		SPH_DEVICE='$(SPH_DEVICE)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
	)
	( \
	cd app; \
	make \
		TARGET='$(TARGET)' \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		CC='$(CC)' \
		CFLAGS='$(CFLAGS)' \
		FC='$(FC)' \
		FCFLAGS='$(FCFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		LDFLAGS='$(LDFLAGS)' \
		LIBS='$(LIBS)' \
		SPH_USR_DEF_LIBS='$(SPH_USR_DEF_LIBS)' \
		MPICH_DIR='$(MPICH_DIR)' \
		MPICH_CFLAGS='$(MPICH_CFLAGS)' \
		MPICH_LDFLAGS='$(MPICH_LDFLAGS)' \
		MPICH_LIBS='$(MPICH_LIBS)' \
		XML2FLAGS='$(XML2FLAGS)' \
		XML2LIBS='$(XML2LIBS)' \
		REALOPT='$(REALOPT)' \
		SPHERE_CFLAGS='$(SPHERE_CFLAGS)' \
		SPHERE_LDFLAGS='$(SPHERE_LDFLAGS)' \
		SPHERE_LIBS='$(SPHERE_LIBS)' \
		SPH_DEVICE='$(SPH_DEVICE)' \
		SPH_SOLV_FLAGS='$(SPH_SOLV_FLAGS)' \
		SPH_SOLV_LDFAGS='$(SPH_SOLV_LDFAGS)' \
		SPH_SOLV_LIBS='$(SPH_SOLV_LIBS)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		SPH_SOLV_FLAGS='$(SPH_SOLV_FLAGS)' \
	)

clean:
	(cd CBC; make clean)
	(cd IP; make clean)
	(cd FB; make clean)
	(cd app; make clean)

allclean: clean
	(cd bin; $(RM) $(TARGET) )

depend:
	(cd CBC; make SPHEREDIR='$(SPHEREDIR)' depend)
	(cd IP; make SPHEREDIR='$(SPHEREDIR)' depend)
	(cd FB; make SPHEREDIR='$(SPHEREDIR)' depend)
	(cd app; make SPHEREDIR='$(SPHEREDIR)' SOLVER_INC_DIR="-I../CBC -I../IP -I../FB " depend)



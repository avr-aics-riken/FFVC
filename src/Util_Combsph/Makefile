###################################################################
#
# CAERU Library
#
# Copyright (c) All right reserved. 2012
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

BASEDIR  = $(PWD)
TARGET   = combsph

include ../FFV/make_setting

all:
	( \
	cd ../PMlib-1.4/src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		MPI_DIR='$(MPI_DIR)' \
	)

	( \
	cd ../FB; \
	make \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		CC='$(CC)' \
		CFLAGS='$(CFLAGS)' \
		FC='$(FC)' \
		FCFLAGS='$(FCFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
	)

	( \
	make -f MakeCOMB \
		TARGET='$(TARGET)' \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		LDFLAGS='$(LDFLAGS)' \
		LIBS='$(LIBS)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
		UDEF_LIB_PATH='$(UDEF_LIB_PATH)' \
		UDEF_LIBS='$(UDEF_LIBS)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
	)

clean: 
	(make -f MakeCOMB clean)
  
mostlyclean: clean
	(cd ../FB; make clean)
	(cd ../PMlib-1.4/src; make clean)

allclean: mostlyclean
	(cd ../../bin; $(RM) $(TARGET) )

depend:
	(cd ../PMlib-1.4/src; make depend)
	(cd ../FB; make depend)
	(make -f MakeCOMB depend)

# DO NOT DELETE

main.o: comb.h FileIO_sph.h endianUtil.h COMB_Define.h
comb.o: comb.h FileIO_sph.h endianUtil.h COMB_Define.h
comb_sph.o: comb.h FileIO_sph.h endianUtil.h COMB_Define.h
comb_plot3d.o: comb.h FileIO_sph.h endianUtil.h COMB_Define.h
FileIO_sph.o: FileIO_sph.h endianUtil.h
FileIO_read_sph.o: FileIO_sph.h endianUtil.h

if MPI
LIBFORMPI=libTextParser_mpi.la
libTextParser_mpi_la_SOURCES=TextParser.cpp TextParserTree.cpp TextParserElement.cpp
libTextParser_mpi_la_CXXFLAGS=-I../include -DBUILD_MPI
else
LIBFORMPI=
endif
lib_LTLIBRARIES = libTextParser.la libTextParser_f90api.la ${LIBFORMPI}

libTextParser_la_SOURCES=TextParser.cpp TextParserTree.cpp TextParserElement.cpp

CXXFLAGS=-I../include 
libTextParser_la_CXXFLAGS=-I../include 


libTextParser_f90api_la_SOURCES=TextParser_fapi.f90
#libTextParser_f90api_la_FCFLAGS=-I../include
FCFLAGS=-I../include

DISTCLEANFILES=*~

TextParser.o : ../include/TextParser.h ../include/TextParserTree.h ../include/TextParserCommon.h ../include/TextParserElement.h
TextParserTree.o: ../include/TextParserTree.h ../include/TextParserCommon.h ../include/TextParserElement.h
TextParserElement.o: ../include/TextParserCommon.h ../include/TextParserElement.h
TextParser_fapi.o:../include/TextParser.inc ../include/TextParser.h
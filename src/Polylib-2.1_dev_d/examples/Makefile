TARGET_NOMPI = test
TARGET_MPI = test_mpi test_mpi2 test_mpi3

CXX=mpicxx
CXXFLAGS=-I/usr/local/TextParser/include -I../include/  -std=c++0x 
CXXLDADD= -L/usr/local/TextParser/lib -L../lib/ -lTextParser_mpi -lPolylib -lMPIPolylib
CXXLDADD_NOMPI= -L/usr/local/TextParser/lib -L../lib/ -lTextParser -lPolylib


.PHONY:all clean

all: $(TARGET_NOMPI) $(TARGET_MPI)

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $<
test.o: test.cxx
	g++ -c $(CXXFLAGS) $<
test_mpi.o: test_mpi.cxx
	$(CXX) -c $(CXXFLAGS) $<
test_mpi2.o: test_mpi2.cxx
	$(CXX) -c $(CXXFLAGS) $<
test_mpi3.o: test_mpi3.cxx
	$(CXX) -c $(CXXFLAGS) $<
CarGroup.o: CarGroup.cxx
	$(CXX) -c $(CXXFLAGS) $<
MyGroupFactory.o: MyGroupFactory.cxx
	$(CXX) -c $(CXXFLAGS) $<

test: test.o
	g++ -o test test.o $(CXXLDADD_NOMPI)

test_mpi:test_mpi.o
	mpicxx -o test_mpi test_mpi.o $(CXXLDADD)
test_mpi2:test_mpi2.o
	mpicxx -o test_mpi2 test_mpi2.o $(CXXLDADD)

test_mpi3:test_mpi3.o CarGroup.o MyGroupFactory.o
	mpicxx -o test_mpi3 test_mpi3.o CarGroup.o MyGroupFactory.o $(CXXLDADD)

clean :
	-rm *.o $(TARGET_NOMPI) $(TARGET_MPI) polylib_config_*.tpp Polylib_*
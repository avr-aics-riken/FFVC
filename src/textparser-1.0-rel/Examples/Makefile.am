if MPI
MPIPROGRAMS=Example3_cpp_mpi
Example3_cpp_mpi_SOURCES=Example3_cpp_mpi.cpp
Example3_cpp_mpi_LDADD= -lstdc++ -L../src -lTextParser_mpi
else
MPIPROGRAMS=
endif

noinst_PROGRAMS=Example1_cpp Example2_cpp Example3_cpp Example4_cpp \
Example1_c Example2_c Example3_c Example4_c \
Example1_f90 Example2_f90 Example3_f90 Example4_f90 \
${MPIPROGRAMS}

# Example3_cpp_mpi_newAPI

#Example_f90_tmp Example_tmp Example_tmp2 
# Example_c_tmp

DISTCLEANFILES=*~
CLEANFILES=Output*.txt nodq/Output*.txt

LDADD=-L../src -lTextParser -lstdc++
CXXFLAGS=-I../include
CFLAGS=-I../include
FCFLAGS=-I../include

Example1_cpp_SOURCES=Example1_cpp.cpp
Example1_c_SOURCES=Example1_c.c
Example1_f90_SOURCES=Example1_f90.f90
Example2_cpp_SOURCES=Example2_cpp.cpp
Example2_c_SOURCES=Example2_c.c
Example2_f90_SOURCES=Example2_f90.f90
Example3_cpp_SOURCES=Example3_cpp.cpp
Example3_c_SOURCES=Example3_c.c
Example3_f90_SOURCES=Example3_f90.f90
Example4_cpp_SOURCES=Example4_cpp.cpp
Example4_c_SOURCES=Example4_c.c
Example4_f90_SOURCES=Example4_f90.f90

#Example3_cpp_mpi_newAPI_SOURCES=Example3_cpp_mpi_newAPI.cpp

#Example_f90_tmp_SOURCES=Example_f90_tmp.f90

Example1_f90_LDADD=  -lstdc++ -L../src -lTextParser -lTextParser_f90api
Example2_f90_LDADD=  -lstdc++ -L../src -lTextParser -lTextParser_f90api
Example3_f90_LDADD=  -lstdc++ -L../src -lTextParser -lTextParser_f90api
Example4_f90_LDADD=  -lstdc++ -L../src -lTextParser -lTextParser_f90api
#Example_f90_tmp_LDADD=-L../src -lstdc++ -lTextParser -lTextParser_f90api

#Example3_cpp_mpi_newAPI_LDADD= -lstdc++ -L../src -lTextParser -lTextParser_mpi

Example1_f90.o:../include/TextParser.inc ../src/libTextParser.la ../src/libTextParser_f90api.la
Example2_f90.o:../include/TextParser.inc ../src/libTextParser.la ../src/libTextParser_f90api.la
Example3_f90.o:../include/TextParser.inc ../src/libTextParser.la ../src/libTextParser_f90api.la
Example4_f90.o:../include/TextParser.inc ../src/libTextParser.la ../src/libTextParser_f90api.la
Example1_cpp.o:../include/TextParser.h ../src/libTextParser.la
Example2_cpp.o:../include/TextParser.h ../src/libTextParser.la
Example3_cpp.o:../include/TextParser.h ../src/libTextParser.la
Example4_cpp.o:../include/TextParser.h ../src/libTextParser.la
Example1_c.o:../include/TextParser.h ../src/libTextParser.la
Example2_c.o:../include/TextParser.h ../src/libTextParser.la
Example3_c.o:../include/TextParser.h ../src/libTextParser.la
Example4_c.o:../include/TextParser.h ../src/libTextParser.la


dist_noinst_DATA=\
incorrect_basic_19.txt \
incorrect_basic_18.txt \
incorrect_basic_17.txt \
incorrect_basic_16.txt \
incorrect_basic_15.txt \
incorrect_basic_14.txt \
incorrect_basic_13.txt \
incorrect_basic_12.txt \
incorrect_basic_11.txt \
incorrect_basic_10.txt \
incorrect_basic_9.txt \
incorrect_basic_8.txt \
incorrect_basic_7.txt \
incorrect_basic_6.txt \
incorrect_basic_5.txt \
incorrect_basic_4.txt \
incorrect_basic_3.txt \
incorrect_basic_2.txt \
incorrect_basic_1.txt \
incorrect_label_16.txt \
incorrect_label_15.txt \
incorrect_label_14.txt \
incorrect_label_13.txt \
incorrect_label_12.txt \
incorrect_label_11.txt \
incorrect_label_10.txt \
incorrect_label_9.txt \
incorrect_label_8.txt \
incorrect_label_7.txt \
incorrect_label_6.txt \
incorrect_label_2.txt \
incorrect_label_1.txt \
incorrect_labelarray_3.txt \
incorrect_labelarray_2.txt \
incorrect_labelarray_1.txt \
correct_labelarray_4.txt \
correct_labelarray_3.txt \
correct_labelarray_2.txt \
correct_labelarray_1.txt \
correct_cond_20.txt \
correct_cond_19.txt \
correct_cond_18.txt \
correct_cond_17.txt \
correct_cond_16.txt \
correct_cond_15.txt \
correct_cond_14.txt \
correct_cond_13.txt \
correct_cond_12.txt \
correct_cond_11.txt \
correct_cond_10.txt \
correct_cond_9.txt \
correct_cond_8.txt \
correct_cond_7.txt \
correct_cond_6.txt \
correct_cond_5.txt \
correct_cond_4.txt \
correct_cond_3.txt \
correct_cond_2.txt \
correct_cond_1.txt \
incorrect_cond_5.txt \
incorrect_cond_4.txt \
incorrect_cond_3.txt \
incorrect_cond_2.txt \
incorrect_cond_1.txt \
correct_label_8.txt \
correct_label_7.txt \
correct_label_6.txt \
correct_label_5.txt \
correct_label_4.txt \
correct_label_3.txt \
correct_label_2.txt \
correct_label_1.txt \
correct_string_4.txt \
correct_string_3.txt \
correct_string_2.txt \
correct_string_1.txt \
correct_basic_13.txt \
correct_basic_12.txt \
correct_basic_11.txt \
correct_basic_10.txt \
correct_basic_9.txt \
correct_basic_8.txt \
correct_basic_7.txt \
correct_basic_6.txt \
correct_basic_5.txt \
correct_basic_4.txt \
correct_basic_3.txt \
correct_basic_2.txt \
correct_basic_1.txt
#incorrect_label_5.txt 
#incorrect_label_4.txt 
#incorrect_label_3.txt 
#incorrect_label_17.txt


.PHONY:$(dist_noinst_DATA)
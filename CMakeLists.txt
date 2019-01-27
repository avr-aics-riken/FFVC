###################################################################################
#
# FFV-C : Frontflow / violet Cartesian
#
# Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
# All rights reserved.
#
# Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology(RIIT),
# Kyushu University.
# All rights reserved.
#
###################################################################################
#
#-*- CMAKE -*-#
# Process this file with cmake to produce a configure script.

# Following parameters must be given at cmake command line
# INSTALL_DIR, CXX, CC, F90

# Build options, see Readme.md for details.
#
# -D INSTALL_DIR=install_directory
#
# -D with_MPI=(yes|no)
#
# -D real_type=(float|double)
#
# -D enable_OPENMP=(yes|no)
#
# -D with_TP=installed_directory
#
# -D with_PM=installed_directory
#
# -D with_CPM=installed_directory
#
# -D with_CDM=installed_directory
#
# -D with_PL=installed_directory
#
# -D with_PAPI=installed_directory or OFF
#

cmake_minimum_required(VERSION 2.6)


#######
# Build options
#######

option(with_MPI "Enable MPI" "ON")
option(enable_OPENMP "Enable OpenMP" "ON")
option(real_type "Type of floating point" "OFF")
option(with_PAPI "Enable PAPI" "OFF")


#######
# Project setting
#######

project(FFVC C CXX Fortran)

if(INSTALL_DIR STREQUAL "")
  set(INSTALL_DIR "/usr/local/FFVC")
endif()

# project()でCMAKE_INSTALL_PREFIX=/usr/localとなるので、再定義
set (CMAKE_INSTALL_PREFIX "${INSTALL_DIR}")

set(CMAKE_BUILD_TYPE Release)


#######
# Policy setting
#######

cmake_policy(SET CMP0012 NEW)


#######
set(PROJECT_VERSION "2.8.0")
set(LIB_REVISION "20190127_1630")
#######



#######
# Cmake module
#######

set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileOptionSelector.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/ffvc_package.cmake)



#######
# Options
#######

# Fortran free form
FreeForm()

AddOptimizeOption()

# OpennMP
checkOpenMP()


# Real type
precision()


# Serial flag for CPMlib header
if(NOT with_MPI)
  SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -DDISABLE_MPI")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_MPI")
endif()


#######
# Display options
#######

message(" ")
message("${PROJECT_NAME} Version : ${PROJECT_VERSION}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "CMAKE_CXX_FLAGS       : "      ${CMAKE_CXX_FLAGS})
message( STATUS "CMAKE_Fortran_COMPILER: "      ${CMAKE_Fortran_COMPILER})
message( STATUS "CMAKE_Fortran_FLAGS   : "      ${CMAKE_Fortran_FLAGS})

message( STATUS "MPI support           : "      ${with_MPI})
message( STATUS "OpenMP support        : "      ${enable_OPENMP})
message( STATUS "Type of floating point: "      ${real_type})
message( STATUS "TextParser support    : "      ${with_TP})
message( STATUS "PMlib support         : "      ${with_PM})
message( STATUS "PAPI support          : "      ${with_PAPI})
message( STATUS "CPMlib support        : "      ${with_CPM})
message( STATUS "Polylib support       : "      ${with_PL})
message( STATUS "CDMlib support        : "      ${with_CDM})

message(" ")


if(CMAKE_CXX_COMPILER MATCHES ".*FCCpx$")
else()
  if(with_MPI)
    include(FindMPI)

    if(MPI_CXX_FOUND)
      message("MPI-CXX-CMP  : ${MPI_CXX_COMPILER}")
      message("MPI-CXX-FLAGS: ${MPI_CXX_COMPILE_FLAGS}")
      message("MPI-CXX-INC  : ${MPI_CXX_INCLUDE_PATH}")
      message("MPI-LNK-FLAGS: ${MPI_CXX_LINK_FLAGS}")
      message("MPI-LIBS     : ${MPI_CXX_LIBRARIES}")
    endif(MPI_CXX_FOUND)
  endif()
endif()

if(USE_F_TCS STREQUAL "YES")
else()
  if(with_MPI)
    include(FindMPI)

    if(MPI_CXX_FOUND)
      message("MPI-CXX-CMP  : ${MPI_CXX_COMPILER}")
      message("MPI-CXX-FLAGS: ${MPI_CXX_COMPILE_FLAGS}")
      message("MPI-CXX-INC  : ${MPI_CXX_INCLUDE_PATH}")
      message("MPI-LNK-FLAGS: ${MPI_CXX_LINK_FLAGS}")
      message("MPI-LIBS     : ${MPI_CXX_LIBRARIES}")
    endif(MPI_CXX_FOUND)
  endif()
endif()


message(" ")
message("PROJECT : FFVC " ${PROJECT_VERSION} " : " ${LIB_REVISION})
message(" ")



#######
# find Libraries to depend on
#######

# TextParser
IF(with_TP)
  SET(TP_DIR "${with_TP}")
  SET(TP_INC "${TP_DIR}/include")
  SET(TP_LIB "${TP_DIR}/lib")
ELSE()
  MESSAGE("Error: can not find TPlib.")
ENDIF()

# PMlib
IF(with_PM)
  SET(PM_DIR "${with_PM}")
  SET(PM_INC "${PM_DIR}/include")
  SET(PM_LIB "${PM_DIR}/lib")
ELSE()
  MESSAGE("Error: can not find PMlib.")
ENDIF()

# PAPI
if(with_PAPI STREQUAL "OFF")
elseif(with_PAPI)
  SET(PAPI_DIR "${with_PAPI}")
  SET(PAPI_INC "${PAPI_DIR}/include")
  SET(PAPI_LIB "${PAPI_DIR}/lib")
else()
  MESSAGE("Error: can not find PAPI")
endif()


# Polylib
IF(with_PL)
  SET(PL_DIR "${with_PL}")
  SET(PL_INC "${PL_DIR}/include")
  SET(PL_LIB "${PL_DIR}/lib")
ELSE()
  MESSAGE("Error: can not find Polylib.")
ENDIF()

# CPMlib
IF(with_CPM)
  SET(CPM_DIR "${with_CPM}")
  SET(CPM_INC "${CPM_DIR}/include")
  SET(CPM_LIB "${CPM_DIR}/lib")
ELSE()
  MESSAGE("Error: can not find CPMlib.")
ENDIF()

# CDMlib
IF(with_CDM)
  SET(CDM_DIR "${with_CDM}")
  SET(CDM_INC "${CDM_DIR}/include")
  SET(CDM_LIB "${CDM_DIR}/lib")
  add_definitions(-D_CDM_OUTPUT)
ELSE()
  MESSAGE("Error: can not find CDMlib.")
ENDIF()


add_definitions(-DHAVE_CONFIG_H)



#######
# Check header files
#######

include(CheckIncludeFiles)

CHECK_INCLUDE_FILES(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)



#######
# Sub Directories
#######

add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(examples)
add_subdirectory(src_util)
add_subdirectory(pt_tool)


#######
# configure files
#######

configure_file(src/FFV/ffv_Version.h.in ${PROJECT_BINARY_DIR}/src/FFV/ffv_Version.h @ONLY)



#######
# Install files
#######

install(FILES  ChangeLog.md
               Readme.md
               License.txt
               ${PROJECT_BINARY_DIR}/CMakeCache.txt
        DESTINATION share)



#######
# Installer
#######

ffvc_package()

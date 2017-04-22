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

add_subdirectory(ASD)
add_subdirectory(F_CORE)
add_subdirectory(F_LS)
add_subdirectory(F_VOF)
add_subdirectory(FB)
add_subdirectory(FFV)
add_subdirectory(FILE_IO)
add_subdirectory(Geometry)
add_subdirectory(IP)


include_directories(
       ${PROJECT_BINARY_DIR}/src/FFV       # ffv_Version.h
       ${PROJECT_SOURCE_DIR}/src/ASD
       ${PROJECT_SOURCE_DIR}/src/F_CORE
       ${PROJECT_SOURCE_DIR}/src/F_LS
       ${PROJECT_SOURCE_DIR}/src/F_VOF
       ${PROJECT_SOURCE_DIR}/src/FB
       ${PROJECT_SOURCE_DIR}/src/FFV
       ${PROJECT_SOURCE_DIR}/src/FILE_IO
       ${PROJECT_SOURCE_DIR}/src/Geometry
       ${PROJECT_SOURCE_DIR}/src/IP
       ${TP_INC}
       ${PM_INC}
       ${PL_INC}
       ${CPM_INC}
       ${CDM_INC}
)

link_directories(
       ${PROJECT_BINARY_DIR}/src/F_CORE
       ${PROJECT_BINARY_DIR}/src/F_LS
       ${PROJECT_BINARY_DIR}/src/F_VOF
       ${PROJECT_BINARY_DIR}/src/FB
       ${PROJECT_BINARY_DIR}/src/FFV
       ${PROJECT_BINARY_DIR}/src/FILE_IO
       ${PROJECT_BINARY_DIR}/src/Geometry
       ${PROJECT_BINARY_DIR}/src/IP
       ${TP_LIB}
       ${PM_LIB}
       ${PL_LIB}
       ${CPM_LIB}
       ${CDM_LIB}
)

set(ffvc_src main.C)

if(with_MPI)
  set(ffvc_target "ffvc-mpi")
else()
  set(ffvc_target "ffvc")
endif()

add_executable(${ffvc_target} ${ffvc_src})
target_link_libraries(${ffvc_target} FFV FIO FCORE FLS ASD GEOM IP FB)

install(TARGETS ${ffvc_target} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


###

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
# Copyright (c) 2016-2020 Research Institute for Information Technology(RIIT),
# Kyushu University.
# All rights reserved.
#
###################################################################################

include_directories(
       ${CPM_INC}
       ${PROJECT_SOURCE_DIR}/src/FB
       ${PROJECT_SOURCE_DIR}/src/F_LS
)

set(ffvc_files
       ffv_blas.f90
       ffv_SOR.f90
)

add_library(FLS STATIC ${ffvc_files})

SET_TARGET_PROPERTIES(FLS PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(FLS ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})

###

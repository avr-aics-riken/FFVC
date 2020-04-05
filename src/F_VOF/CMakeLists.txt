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
       ${PROJECT_SOURCE_DIR}/src/F_VOF
)

set(ffvc_files
       c3d_vof.f90
)

add_library(FVOF STATIC ${ffvc_files})

SET_TARGET_PROPERTIES(FVOF PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(FVOF ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})

###

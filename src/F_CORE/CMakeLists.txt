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
       ${PROJECT_SOURCE_DIR}/src/F_CORE
)

set(ffvc_files
       ffv_Ffunc.h
       ffv_forcing.f90
       ffv_pbc.f90
       ffv_pscalar.f90
       ffv_utility.f90
       ffv_vbc_inner.f90
       ffv_vbc_outer.f90
       ffv_vbc_outer_face.f90
       ffv_vbc_outer_flux.f90
       ffv_velocity_binary.f90
       force.h
       ffv_velocity_cds.f90
       FB_util.f90
       io_sph.f90
)

add_library(FCORE STATIC ${ffvc_files})

SET_TARGET_PROPERTIES(FCORE PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(FCORE ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})

###

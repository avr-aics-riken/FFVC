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


include_directories(
       ${PROJECT_BINARY_DIR}/src/FFV  # ffv_Version.h
       ${PROJECT_SOURCE_DIR}/src/FB
       ${PROJECT_SOURCE_DIR}/src/Geometry
       ${PROJECT_SOURCE_DIR}/src/ASD
       ${TP_INC}
       ${PM_INC}
       ${PL_INC}
       ${CPM_INC}
)

link_directories(
      ${TP_LIB}
      ${PL_LIB}
)

set(ffvc_files
       SubDomain.C
       ASDmodule.C
)

add_library(ASD STATIC ${ffvc_files})


if(with_MPI)
  SET(ffvc_libs "-lPOLYmpi -lTPmpi")
else()
  SET(ffvc_libs "-lPOLY -lTP")
endif()

target_link_libraries(ASD "${ffvc_libs}")

###

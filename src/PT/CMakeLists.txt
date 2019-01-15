###################################################################################
#
# FFV-C : Frontflow / violet Cartesian
#
# Copyright (c) 2016-2019 Research Institute for Information Technology(RIIT),
# Kyushu University.
# All rights reserved.
#
###################################################################################


include_directories(
       ${PROJECT_BINARY_DIR}/src/FFV     # ffv_Version.h
       ${PROJECT_SOURCE_DIR}/src/FB
       ${TP_INC}
       ${PL_INC}
       ${PM_INC}
       ${CPM_INC}
)

link_directories(
      ${TP_LIB}
      ${PM_LIB}
      ${CPM_LIB}
)

set(ffvc_files
       Chunk.C
       Tracking.C
       Cloud.C
       PtComm.C
)

add_library(PT STATIC ${ffvc_files})

if(with_MPI)
  SET(ffvc_libs "-lCPMmpi -lPMmpi -lTPmpi")
else()
  SET(ffvc_libs "-lCPM -lPM -lTP")
endif()

target_link_libraries(PT "${ffvc_libs}")

###

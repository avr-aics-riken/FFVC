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
       ${PROJECT_BINARY_DIR}/src/FFV    # ffv_Version.h
       ${PROJECT_SOURCE_DIR}/src/FB
       ${PROJECT_SOURCE_DIR}/src/Geometry
       ${TP_INC}
       ${PM_INC}
       ${PL_INC}
       ${CPM_INC}
       ${CDM_INC}
)

link_directories(
      ${TP_LIB}
      ${PM_LIB}
      ${CPM_LIB}
      ${CDM_LIB}
)

IF(with_PAPI)
  include_directories(${PAPI_INC})
  link_directories(${PAPI_LIB})
  # link static library
  set(PAPI_STATIC_LIB ${PAPI_LIB}/libpapi.a
                      ${PAPI_LIB}/libpfm.a
                      )
endif()

set(ffvc_files
       Alloc.C
       BndOuter.C
       Component.C
       Control.C
       DataHolder.C
       FBUtility.C
       History.C
       Intrinsic.C
       IterationControl.C
       MonCompo.C
       Monitor.C
       ParseBC.C
       ParseMat.C
       Sampling.C
       SetBC.C
       VoxInfo.C
)

add_library(FB STATIC ${ffvc_files})

if(with_MPI)
  SET(ffvc_libs "-lCDMmpi -lCPMmpi -lPMmpi -lTPmpi")
else()
  SET(ffvc_libs "-lCDM -lCPM -lPM -lTP")
endif()

if(with_PAPI)
  if(TARGET_ARCH STREQUAL "FX100")
  # on AICS FX100 static PAPI libs does not match the other system libs.
  # so the following dynamic libs must be used. bonehead Fujitsu...
    set(ffvc_libs "${ffvc_libs} -lpapi_ext -lpapi -lpfm")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  ## GNU linker ld options are better passed using -Wl, rather than -B compiler phase option
    set(ffvc_libs "${ffvc_libs} -lpapi_ext ${PAPI_STATIC_LIB}")
  else()
    ## Intel compiler & others
    set(ffvc_libs "${ffvc_libs} -lpapi_ext ${PAPI_STATIC_LIB}")
  endif()
endif()

target_link_libraries(FB "${ffvc_libs}")

###

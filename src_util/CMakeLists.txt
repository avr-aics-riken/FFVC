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
       ${PROJECT_BINARY_DIR}/src_util
       ${PROJECT_SOURCE_DIR}/src/FB
       ${TP_INC}
       ${CPM_INC}
       ${CDM_INC}
       ${PM_INC}
)


link_directories(
      ${PROJECT_BINARY_DIR}/src/FB
      ${CDM_LIB}
      ${CPM_LIB}
      ${PM_LIB}
      ${TP_LIB}
)


set(comb_files
    main.C
    comb.C
    comb_sph.C
    comb_avs.C
    FileIO_sph.C
    FileIO_read_sph.C
)


if(with_MPI)
  SET(comb_libs "-lCDMmpi -lCPMmpi -lPMmpi -lTPmpi")
else()
  SET(comb_libs "-lCDM -lCPM -lPM -lTP")
endif()

add_executable(combsph ${comb_files})
target_link_libraries(combsph "${comb_libs}")

install(TARGETS combsph DESTINATION bin)

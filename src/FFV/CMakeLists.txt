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
       ${PROJECT_BINARY_DIR}/src/FFV    # ffv_Version.h
       ${PROJECT_SOURCE_DIR}/src/FB
       ${PROJECT_SOURCE_DIR}/src/F_CORE
       ${PROJECT_SOURCE_DIR}/src/F_LS
       ${PROJECT_SOURCE_DIR}/src/IP
       ${PROJECT_SOURCE_DIR}/src/Geometry
       ${PROJECT_SOURCE_DIR}/src/FILE_IO
       ${PROJECT_SOURCE_DIR}/src/FFV
       ${TP_INC}
       ${PM_INC}
       ${PL_INC}
       ${CPM_INC}
       ${CDM_INC}
)

link_directories(
#      ${TP_LIB}
#      ${PM_LIB}
#      ${CPM_LIB}
#      ${CDM_LIB}
)

set(ffvc_files
       NS_FS_E_Binary.C
       NS_FS_E_CDS.C
       PS_Binary.C
       ffv.C
       ffv_Alloc.C
       ffv_Filter.C
       ffv_Heat.C
       ffv_Initialize.C
       ffv_LS.C
       ffv_Loop.C
       ffv_Post.C
       ffv_SetBC.C
       ffv_TerminateCtrl.C
       dryrunBC.C
)

add_library(FFV STATIC ${ffvc_files})

###

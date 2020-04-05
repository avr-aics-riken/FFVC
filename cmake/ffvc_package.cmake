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

##
## Packaging
##

macro (ffvc_package)
set (CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS TRUE)
include (InstallRequiredSystemLibraries)

set (CPACK_PACKAGE_NAME "${PROJECT}")
set (CPACK_PACKAGE_VENDOR "FFVC")

set (CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set (CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")

set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Frontflow violet / Cartesian")
set (CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/License.txt")
set (CPACK_RESOURCE_FILE_README  "${PROJECT_SOURCE_DIR}/Readme.md")

set (CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
set (CPACK_SOURCE_GENERATOR "TGZ")

# config for source installer
set (CPACK_SOURCE_OUTPUT_CONFIG_FILE "${PROJECT_BINARY_DIR}/CPackSourceConfig.cmake")

set (CPACK_SOURCE_IGNORE_FILES "/BUILD/;/build/;/.git/;.#;/#;.DS_Store" CACHE STRING "CPACK will ignore this files" )

# config for binary installer
#set (CPACK_OUTPUT_CONFIG_FILE "${PROJECT_BINARY_DIR}/CPackConfig.cmake")

include (CPack)
endmacro()

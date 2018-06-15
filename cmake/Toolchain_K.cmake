set(CMAKE_SYSTEM_NAME Linux)

include(CMakeForceCompiler)

# CMAKE_FORCE_Fortran_Compiler is not supported ver. 2.6
# A : cmake version >= 2.7
# B : cmake version < 2.7

if (CMAKE_MAJOR_VERSION GREATER 2)
  set(build_rule "A")
else()
  if(CMAKE_MINOR_VERSION GREATER 6)
    set(build_rule "A")
  else()
    set(build_rule "B")
  endif()
endif()


if(with_MPI)
  CMAKE_FORCE_C_COMPILER(mpifccpx GNU)
  CMAKE_FORCE_CXX_COMPILER(mpiFCCpx GNU)
  if (build_rule STREQUAL "A")
    CMAKE_FORCE_Fortran_COMPILER(mpifrtpx GNU)
  else()
    set(CMAKE_Fortran_COMPILER mpifrtpx)
    set(CMAKE_Fortran_COMPILER_WORKS true)
  endif()
else()
  CMAKE_FORCE_C_COMPILER(fccpx GNU)
  CMAKE_FORCE_CXX_COMPILER(FCCpx GNU)
  if (build_rule STREQUAL "A")
    CMAKE_FORCE_Fortran_COMPILER(frtpx GNU)
  else()
    set(CMAKE_Fortran_COMPILER frtpx)
    set(CMAKE_Fortran_COMPILER_WORKS true)
  endif()
endif()

set(CMAKE_FIND_ROOT_PATH /opt/FJSVtclang/1.2.0)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

set(TARGET_ARCH "K")
set(USE_F_TCS "YES")

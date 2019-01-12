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

##
## Compile option selector
##


macro (AddOptimizeOption)
  if (TARGET_ARCH STREQUAL "INTEL_F_TCS")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast,parallel,optmsg=2 -V -Nsrc -Xg -Nfjcex -std=c++11")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast,parallel,optmsg=2 -V -Nsrc -Xg")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Cpp -Kfast,parallel,optmsg=2 -V -Qt")
    # -Xg   : gcc compatible flag
    # -fPIC : PIC flag
    # -Nfjcex : to link PMlib

   elseif (USE_F_TCS STREQUAL "YES")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Nsrc -x0 -Xg -Nfjcex")
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Nsrc -x0 -Xg")
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Cpp -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Qt")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall")
    set(CMAKE_Fortran_FLAGS "-O3 -cpp -Wall --free-line-length-none")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -O3 -qopt-report=5 -DMPICH_IGNORE_CXX_SEEK -std=c++11")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHOST -O3 -qopt-report=5")
    set(CMAKE_Fortran_FLAGS "-xHOST -O3 -Warn unused -fpp -qopt-report=5")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -Mipa=fast,inline -O4 -Minfo=intensity,vect -cpp=mm")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fast -Mipa=fast,inline -O4 -Minfo=intensity,vect -cpp=mm")
    set(CMAKE_Fortran_FLAGS "-fast -Mipa=fast,inline -O4 -Minfo=intensity,vect -cpp=mm")

  else()
    message("using default option")
  endif()
endmacro()


macro (AddSSE)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fastsse")
  endif()
endmacro()


macro (FreeForm)
  if(CMAKE_Fortran_COMPILER MATCHES ".*frtpx$")
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

  elseif(TARGET_ARCH STREQUAL "INTEL_F_TCS")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Free")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mfree")

  endif()
endmacro()


macro(C99)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endmacro()


macro(CPP11)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  endif()
endmacro()


macro(checkOpenMP)
  if(enable_OPENMP)
    if(USE_F_TCS STREQUAL "YES")
      set(OpenMP_C_FLAGS "-Kopenmp")
      set(OpenMP_CXX_FLAGS "-Kopenmp")
      set(OpenMP_Fortran_FLAGS "-Kopenmp")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(OpenMP_C_FLAGS "-fopenmp")
      set(OpenMP_CXX_FLAGS "-fopenmp")
      set(OpenMP_Fortran_FLAGS "-fopenmp")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      set(OpenMP_C_FLAGS "-qopenmp")
      set(OpenMP_CXX_FLAGS "-qopenmp")
      set(OpenMP_Fortran_FLAGS "-qopenmp")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
      set(OpenMP_C_FLAGS "-mp")
      set(OpenMP_CXX_FLAGS "-mp")
      set(OpenMP_Fortran_FLAGS "-mp")
    else()
      find_package(OpenMP REQUIRED)
    endif()

    # OpenMP_*_FLAGSにはfind_package(OpenMP REQUIRED)でオプションフラグが設定される
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endmacro()


macro(precision)
  if(real_type STREQUAL "OFF")
  # nothing, default is float
  set(real_type "float")

  elseif(real_type STREQUAL "float")
  # nothing

  elseif(real_type STREQUAL "double")
    ADD_DEFINITIONS(-D_REAL_IS_DOUBLE_)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_REAL_IS_DOUBLE_")

    if(CMAKE_Fortran_COMPILER MATCHES ".*frtpx$")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -CcdRR8")
    elseif(CMAKE_Fortran_COMPILER MATCHES ".*frt$")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -CcdRR8")

    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")

    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")

    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
    endif()

  else() # neither 'float' nor 'double'
    message("@@@@@@@@@@@")
    message("FATAL ERROR : Invalid floating type : ${real_type}")
    message("@@@@@@@@@@@")
  ENDIF()
endmacro()

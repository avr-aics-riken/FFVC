#ifndef _COMB_DEFINE_H_
#define _COMB_DEFINE_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 @file   COMB_Define.h
 @brief  FlowBase Definition Header
 @author aics
 */


#define LOG_OUT_   if(lflag)
#define LOG_OUTV_  if(lflagv)
#define STD_OUT_   if(pflag) 
#define STD_OUTV_  if(pflagv) 

#define OUTFORMAT_IS_SPH    0
#define OUTFORMAT_IS_PLOT3D 1
#define OUTFORMAT_IS_AVS    2

#define OUTPUT_FORMAT_TYPE_UNKNOWN  -1
#define OUTPUT_FORMAT_TYPE_BINARY    0
#define OUTPUT_FORMAT_TYPE_ASCII     1

#define OUTPUT_DIV_FUNC_ON  0
#define OUTPUT_DIV_FUNC_OFF 1

// 定常or非定常解析
#define STEADY   0
#define UNSTEADY 1

#define FILE_PATH_LENGTH 64

#define DFI_LINE_LENGTH 256

#define TT_OTHER_ENDIAN		1
#define TT_LITTLE_ENDIAN	2
#define TT_BIG_ENDIAN		3

#define REAL_UNKNOWN 0
#define SPH_FLOAT    1
#define SPH_DOUBLE   2

#define SPH_DATA_UNKNOWN 0
#define SPH_SCALAR       1
#define SPH_VECTOR       2

#endif // _COMB_DEFINE_H_

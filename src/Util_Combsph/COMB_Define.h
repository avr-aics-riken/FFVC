#ifndef _COMB_DEFINE_H_
#define _COMB_DEFINE_H_

// #################################################################
//
// Combine sph files and output
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan.
//
// #################################################################

/**
 @file   COMB_Define.h
 @brief  FlowBase Definition Header
 @author kero
 */

#include "mydebug.h"

#define LOG_OUT_   if(lflag)
#define LOG_OUTV_  if(lflagv)
#define STD_OUT_   if(pflag) 
#define STD_OUTV_  if(pflagv) 

#define OUTFORMAT_IS_SPH    0
#define OUTFORMAT_IS_PLOT3D 1

#define DFI_LINE_LENGTH 256

#define TT_OTHER_ENDIAN		1
#define TT_LITTLE_ENDIAN	2
#define TT_BIG_ENDIAN		3

//#define OUTPUT_REAL_UNKNOWN 0
//#define OUTPUT_FLOAT        1
//#define OUTPUT_DOUBLE       2

#define REAL_UNKNOWN 0
#define SPH_FLOAT    1
#define SPH_DOUBLE   2

#define SPH_DATA_UNKNOWN 0
#define SPH_SCALAR       1
#define SPH_VECTOR       2

#endif // _COMB_DEFINE_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################
//
///
/// @file  type.h
/// @brief クロスプラットホームのデータ型宣言
///

#ifndef __FFV_FILEIO_TYPE_H__
#define __FFV_FILEIO_TYPE_H__


#ifdef _WIN32

typedef short          int16_t;
typedef unsigned short uint16_t;
typedef int            int32_t;
typedef unsigned int   uint32_t;
typedef long           int64_t;
typedef unsigned long  uint64_t;

#else  // _WIN32

#include <stdint.h>

#endif // _WIN32

typedef bool          b8;   ///< 論理型
typedef char          s8;   ///< 符号付き 8bit整数型
typedef unsigned char u8;   ///< 符号なし 8bit整数型
typedef int16_t       s16;  ///< 符号付き16bit整数型
typedef uint16_t      u16;  ///< 符号なし16bit整数型
typedef int32_t       s32;  ///< 符号付き32bit整数型
typedef uint32_t      u32;  ///< 符号なし32bit整数型
typedef int64_t       s64;  ///< 符号付き64bit整数型
typedef uint64_t      u64;  ///< 符号なし64bit整数型
typedef float         f32;  ///< 32bit浮動小数点 (単精度浮動小数点)
typedef double        f64;  ///< 64bit浮動小数点 (倍精度浮動小数点)

#endif // __FFV_FILEIO_TYPE_H__

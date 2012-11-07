#ifndef _LAYOUT_DEFINE_H_
#define _LAYOUT_DEFINE_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/**
 @file   LAYOUT_Define.h
 @brief  FlowBase Definition Header
 @author kero
 */

#include "mydebug.h"

//#ifdef _REAL_IS_DOUBLE_
//  #define REAL_TYPE double
//#else
//  /** 実数型の指定
//   * - デフォルトでは、REAL_TYPE=float
//   * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
//   *   REAL_TYPE=doubleになる
//   */
//  #define REAL_TYPE float
//#endif

#define OUTFORMAT_IS_SPH    0
#define OUTFORMAT_IS_PLOT3D 1

// general
//#define FB_FILE_PATH_LENGTH 64
#define FB_BUFF_LENGTH      256

#endif // _LAYOUT_DEFINE_H_

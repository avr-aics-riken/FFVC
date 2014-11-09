#ifndef _FB_ALLOC_H_
#define _FB_ALLOC_H_

//##################################################################################
//
// Flow Base class
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

/** 
 * @file   Alloc.h
 * @brief  FlowBase Files class Header
 * @author aics
 */

#include <string.h>
#include "FB_Define.h"


class Alloc {
  
public:
  /** コンストラクタ */
  Alloc() {}
  
  /**　デストラクタ */
  ~Alloc() {}
  
public:

  /**
   @brief データ領域をアロケートする
   @retval エラーコード
   @param[in]      sz       計算内部領域のサイズ
   @param[in]      gc       ガイドセルサイズ
   */
  
  
  static double* Double_S3D(const int* sz, const int gc);
  
  static float* Float_S3D(const int* sz, const int gc);
  
  static float* Float_S4D(const int* sz, const int gc, const int dnum);
  
  static int* Int_S3D(const int* sz, const int gc);
  
	static REAL_TYPE* Real_S3D(const int* sz, const int gc);
  
  static REAL_TYPE* Real_S4D(const int* sz, const int gc, const int dnum);
  
  static REAL_TYPE* Real_V3D(const int* sz, const int gc);

  static unsigned* Uint_S3D(const int* sz, const int gc);

};

#endif // _FB_ALLOC_H_

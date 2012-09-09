#ifndef _FB_ALLOC_H_
#define _FB_ALLOC_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   Alloc.h
 * @brief  FlowBase Files class Header
 * @author kero
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

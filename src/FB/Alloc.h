#ifndef _FB_ALLOC_H_
#define _FB_ALLOC_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file Alloc.h
 * @brief FlowBase Files class Header
 * @author kero
 */

#include "cpm_Define.h"
#include "limits.h"


class Alloc {
  
public:
  /** コンストラクタ */
  Alloc() {}
  
  /**　デストラクタ */
  ~Alloc() {}
  
public:

  /**
   @brief データ領域をアロケートする（Scalar:float）
   @retval エラーコード
   @param[in]      sz       計算内部領域のサイズ
   @param[in]      gc       ガイドセルサイズ
   @param[in/out]  mc       メモリ使用量のカウンタ
   */
  static float* Float_S3D(const int* sz, const int gc, float &mc);
  
  
  /**
   @brief データ領域をアロケートする（Scalar4:REAL_TYPE）
   @retval エラーコード
   @param[in]      sz       計算内部領域のサイズ
   @param[in]      gc       ガイドセルサイズ
   @param[in]      dnum     4つめのサイズ
   @param[in/out]  mc       メモリ使用量のカウンタ
   */
  static float* Float_S4D(const int* sz, const int gc, const int dnum, float &mc);
  
  
  /**
   @brief データ領域をアロケートする（Scalar:int）
   @retval エラーコード
   @param[in]      sz       計算内部領域のサイズ
   @param[in]      gc       ガイドセルサイズ
   @param[in/out]  mc       メモリ使用量のカウンタ
   */
  static int* Int_S3D(const int* sz, const int gc, float &mc);
  
  
  /**
   @brief データ領域をアロケートする（Scalar:REAL_TYPE）
   @retval エラーコード
   @param[in]      sz       計算内部領域のサイズ
   @param[in]      gc       ガイドセルサイズ
   @param[in/out]  mc       メモリ使用量のカウンタ
   */
	static REAL_TYPE* Real_S3D(const int* sz, const int gc, float &mc);
  
  
  /**
   @brief データ領域をアロケートする（Vector:REAL_TYPE）
   @retval エラーコード
   @param[in]      sz       計算内部領域のサイズ
   @param[in]      gc       ガイドセルサイズ
   @param[in/out]  mc       メモリ使用量のカウンタ
   */
  static REAL_TYPE* Real_V3D(const int* sz, const int gc, float &mc);

  
  /**
   @brief データ領域をアロケートする（Scalar:unsigned）
   @retval エラーコード
   @param[in]      sz       計算内部領域のサイズ
   @param[in]      gc       ガイドセルサイズ
   @param[in/out]  mc       メモリ使用量のカウンタ
   */
  static unsigned* Uint_S3D(const int* sz, const int gc, float &mc);

};

#endif // _FB_ALLOC_H_

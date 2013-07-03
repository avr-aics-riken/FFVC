//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 @file   IP_Polygon.C
 @brief  IP_Polygon class
 @author kero
 */

#include "IP_Polygon.h"


// #################################################################
/* @brief 矩形の計算領域のセルIDを設定する
 * @param [in,out] mid   媒質情報の配列
 * @param [in]     R     Controlクラスのポインタ
 * @param [in]     G_org グローバルな原点（無次元）
 * @param [in]     Nmax  Controlクラスのポインタ
 * @param [in]     mat   MediumListクラスのポインタ
 */
void IP_Polygon::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int ref = R->Mode.Base_Medium;

  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ref) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = ref;
      }
    }
  }
}

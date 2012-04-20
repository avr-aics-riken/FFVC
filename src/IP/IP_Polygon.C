/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Polygon.C
//@brief IP_Polygon class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_Polygon.h"

/**
 @fn bool IP_Polygon::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
 @brief 領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 @note 基点，領域サイズ，ピッチが有効，分割数はそれらから計算
 */
void IP_Polygon::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
{
  
  // 等分割のチェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  
  // 分割数を計算
  //sz[0] = (unsigned)ceil(wth[0]/pch[0]);
  //sz[1] = (unsigned)ceil(wth[1]/pch[1]);
  //sz[2] = (unsigned)ceil(wth[2]/pch[2]);
}

/**
 @fn void IP_Polygon::setup(int* mid, Control* R, REAL_TYPE* G_org)
 @brief 計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_Polygon::setup(int* mid, Control* R, REAL_TYPE* G_org)
{
  int m = (imax+2*guide)*(jmax+2*guide)*(kmax+2*guide);
  int ref = R->Mode.Base_Medium;

#pragma omp parallel for firstprivate(m, ref) schedule(static)
  for (int i=0; i<m; i++) mid[i] = ref;
}

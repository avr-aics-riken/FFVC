/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file IP_Polygon.C
//@brief IP_Polygon class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_Polygon.h"

/**
 @fn bool IP_Polygon::setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3])
 @brief 領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
void IP_Polygon::setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3])
{
  wth[0] = pch[0]*(SKL_REAL)sz[0];
  wth[1] = pch[1]*(SKL_REAL)sz[1];
  wth[2] = pch[2]*(SKL_REAL)sz[2];
  
  // チェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    assert(0);
  }
  if ( ((unsigned)(wth[0]/pch[0]) != sz[0]) ||
       ((unsigned)(wth[1]/pch[1]) != sz[1]) ||
       ((unsigned)(wth[2]/pch[2]) != sz[2]) ) {
    Hostonly_ printf("Error : Invalid parameters among 'VoxelSize', 'VoxelPitch', and 'VoxelWidth' in DomainInfo section.\n");
    assert(0);
  }
}

/**
 @fn void IP_Polygon::setup(int* mid, Control* R, SKL_REAL* G_org)
 @brief 計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_Polygon::setup(int* mid, Control* R, SKL_REAL* G_org)
{
  int i,m;
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  
  // Initialize  全領域をfluidにしておく
  m = (imax+2*guide)*(jmax+2*guide)*(kmax+2*guide);
  for (i=0; i<m; i++) mid[i]=mid_fluid;

}
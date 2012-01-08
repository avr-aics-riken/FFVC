/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file IP_PMT.C
//@brief IP_PMT class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_PMT.h"

/**
 @fn bool IP_PMT::setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3])
 @brief 領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
void IP_PMT::setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3])
{
  // forced
  if (R->Unit.Param != NONDIMENSIONAL) {
    Hostonly_ printf("\tError : PMT class is designed for only non-dimensional parameter\n");
    assert(0);
  }
  
  // 性能測定モードをOnにする
  R->Hide.PM_Test = ON;

  // 等ピッチのチェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    assert(0);
  }
  if ( pch[0] <= 0.0 ) {
    Hostonly_ printf("Error : 'VoxelPitch' must be grater than zero.\n");
    assert(0);
  }
  
  wth[0] = pch[0]*(SKL_REAL)sz[0];
  wth[1] = pch[1]*(SKL_REAL)sz[1];
  wth[2] = pch[2]*(SKL_REAL)sz[2];
  
  org[0] = -0.5*wth[0];
  org[1] = -0.5*wth[1];
  org[2] = -0.5*wth[2];
  
  // Setting depends on Example,  INTRINSIC
  if ( (sz[0]/2*2 != sz[0]) || (sz[1]/2*2 != sz[1]) || (sz[2]/2*2 != sz[2]) ) {
    printf("\tDimension size must be even for all direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
    assert(0);
  }
}

/**
 @fn void IP_PMT::setup(int* mid, Control* R, SKL_REAL* G_org)
 @brief 計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_PMT::setup(int* mid, Control* R, SKL_REAL* G_org)
{
  int i,j,k;
  unsigned m;

  // Inner
  for (k=1; k<=(int)kmax; k++) {
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
        m = SklUtil::getFindexS3D(size, guide, i, j, k);
        mid[m] = 1;
      }
    }
  }
}
/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file IP_RSP.C
//@brief IP_RSP class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_RSP.h"

/**
 @fn bool IP_RSP::setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3])
 @brief RSPの領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
void IP_RSP::setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3])
{
  // forced
  if (R->Unit.Param != DIMENSIONAL) {
    Hostonly_ printf("\tError : RSP class is designed for only dimensional parameter\n");
    assert(0);
  }

  // 分割数を取得する
  if ( !SklUtil::getCellInfo(R->NoDimension, sz, org, pch, wth) ) assert(0);
  
  // Z方向のチェック
  if ( sz[2] != 3 ) {
    Hostonly_ printf("\tError : The size of Z-direction must be 3.\n");
    assert(0);
  }
  
  // 領域アスペクト比のチェック
  if ( sz[0]*5 != sz[1] ) {
    Hostonly_ printf("\tError : The number of division must be 1:5 (=X:Y)\n");
    assert(0);
  }

  // 二次元の領域設定
  wth[0] =  0.02;
  wth[1] =  0.1;
  org[0] = -0.01;
  org[1] =  0.0;
  
  pch[0] = wth[0] / (SKL_REAL)sz[0];
  pch[1] = wth[1] / (SKL_REAL)sz[1];
  
  if ( (pch[0]-pch[1])/pch[0] > 1.0e-3 ) { // 桁落ち防止
    Hostonly_ printf("\tVoxel width must be same between X(%e) and Y(%e) direction.\n", pch[0], pch[1]);
    assert(0);
  }
  
  // Z方向は3層に合わせて調整
  pch[2] =  pch[0];
  wth[2] =  pch[2]*3.0;
  org[2] = -pch[2]*1.5;
  
  // 加速時間の警告
  Hostonly_ printf("\t##### ATTENTION : Rayleigh's Problem requires NO acceleration #####\n");
}

/**
 @fn void IP_RSP::setup(int* mid, Control* R, SKL_REAL* G_org)
 @brief RSPの計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_RSP::setup(int* mid, Control* R, SKL_REAL* G_org)
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
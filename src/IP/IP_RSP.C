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

//@file IP_RSP.C
//@brief IP_RSP class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_RSP.h"

/**
 @fn bool IP_RSP::setDomainParameter(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE reg[3], REAL_TYPE pch[3])
 @brief RSPの領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
void IP_RSP::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
  // forced
  if (R->Unit.Param != DIMENSIONAL) {
    Hostonly_ printf("\tError : RSP class is designed for only dimensional parameter\n");
    Exit(0);
  }

  // 分割数を取得する
  if ( !FBUtility::getCellInfo(3, sz, org, pch, reg) ) Exit(0);
  
  // Z方向のチェック
  if ( sz[2] != 3 ) {
    Hostonly_ printf("\tError : The size of Z-direction must be 3.\n");
    Exit(0);
  }
  
  // 領域アスペクト比のチェック
  if ( sz[0]*5 != sz[1] ) {
    Hostonly_ printf("\tError : The number of division must be 1:5 (=X:Y)\n");
    Exit(0);
  }

  // 二次元の領域設定
  reg[0] =  0.02;
  reg[1] =  0.1;
  org[0] = -0.01;
  org[1] =  0.0;
  
  pch[0] = reg[0] / (REAL_TYPE)sz[0];
  pch[1] = reg[1] / (REAL_TYPE)sz[1];
  
  if ( (pch[0]-pch[1])/pch[0] > 1.0e-3 ) { // 桁落ち防止
    Hostonly_ printf("\tVoxel width must be same between X(%e) and Y(%e) direction.\n", pch[0], pch[1]);
    Exit(0);
  }
  
  // Z方向は3層に合わせて調整
  pch[2] =  pch[0];
  reg[2] =  pch[2]*3.0;
  org[2] = -pch[2]*1.5;
  
  // 加速時間の警告
  Hostonly_ printf("\t##### ATTENTION : Rayleigh's Problem requires NO acceleration #####\n");
}

/**
 @fn void IP_RSP::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
 @brief RSPの計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 @param mat
 */
void IP_RSP::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  size_t m;
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // Inner
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd); //FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = 1;
      }
    }
  }
}

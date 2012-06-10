/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_PMT.C
//@brief IP_PMT class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_PMT.h"

/**
 @fn bool IP_Rect::getTP(Control* R, TPControl* tpCntl)
 @brief パラメータを取得する
 */
bool IP_Rect::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 媒質指定
  label="/Parameter/Intrinsic_Example/Fluid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Fluid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_fluid = str;
  
  label="/Parameter/Intrinsic_Example/Solid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Solid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_solid = str;
  
  return true;
}

/**
 @fn bool IP_PMT::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
 @brief 領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
void IP_PMT::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
{
  // forced
  if (R->Unit.Param != NONDIMENSIONAL) {
    Hostonly_ printf("\tError : PMT class is designed for only non-dimensional parameter\n");
    Exit(0);
  }
  
  // 性能測定モードをOnにする
  R->Hide.PM_Test = ON;

  // 等ピッチのチェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  if ( pch[0] <= 0.0 ) {
    Hostonly_ printf("Error : 'VoxelPitch' must be grater than zero.\n");
    Exit(0);
  }
  
  wth[0] = pch[0]*(REAL_TYPE)sz[0];
  wth[1] = pch[1]*(REAL_TYPE)sz[1];
  wth[2] = pch[2]*(REAL_TYPE)sz[2];
  
  org[0] = -0.5*wth[0];
  org[1] = -0.5*wth[1];
  org[2] = -0.5*wth[2];
  
  // Setting depends on Example,  INTRINSIC
  if ( (sz[0]/2*2 != sz[0]) || (sz[1]/2*2 != sz[1]) || (sz[2]/2*2 != sz[2]) ) {
    printf("\tDimension size must be even for all direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
    Exit(0);
  }
}

/**
 @fn void IP_PMT::setup(int* mid, Control* R, REAL_TYPE* G_org)
 @brief 計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_PMT::setup(int* mid, Control* R, REAL_TYPE* G_org)
{
  unsigned m, m_sz[3], gd;
  m_sz[0] = size[0];
  m_sz[1] = size[1];
  m_sz[2] = size[2];
  gd = guide;
  
#pragma omp parallel for firstprivate(m_sz, gd) schedule(static)
  for (int k=1; k<=(int)m_sz[2]; k++) {
    for (int j=1; j<=(int)m_sz[1]; j++) {
      for (int i=1; i<=(int)m_sz[0]; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = 1;
      }
    }
  }
}

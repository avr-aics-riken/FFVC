// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file IP_PMT.C
 * @brief IP_PMT class
 * @author kero
 */

#include "IP_PMT.h"


// パラメータを取得する
bool IP_PMT::getTP(Control* R, TPControl* tpCntl)
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


// 領域情報を設定する
void IP_PMT::setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
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
  
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = pch[1]*(REAL_TYPE)sz[1];
  reg[2] = pch[2]*(REAL_TYPE)sz[2];
  
  org[0] = -0.5*reg[0];
  org[1] = -0.5*reg[1];
  org[2] = -0.5*reg[2];
  
  // Setting depends on Example,  INTRINSIC
  if ( (sz[0]/2*2 != sz[0]) || (sz[1]/2*2 != sz[1]) || (sz[2]/2*2 != sz[2]) ) {
    printf("\tDimension size must be even for all direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
    Exit(0);
  }
}


// 計算領域のセルIDを設定する
void IP_PMT::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  size_t m;
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd); //FBUtility::getFindexS3D(sz, gd, i, j, k);
        mid[m] = 1;
      }
    }
  }
}

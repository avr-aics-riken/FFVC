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
 @file IP_Polygon.C
 @brief IP_Polygon class
 @author kero
 */

#include "IP_Polygon.h"


// 領域情報を設定する
void IP_Polygon::setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  
  // 等分割のチェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  
  // 分割数を計算
  //sz[0] = (unsigned)ceil(reg[0]/pch[0]);
  //sz[1] = (unsigned)ceil(reg[1]/pch[1]);
  //sz[2] = (unsigned)ceil(reg[2]/pch[2]);
}



// 計算領域のセルIDを設定する
void IP_Polygon::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  // ローカルにコピー
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  size_t m = (imax+2*guide)*(jmax+2*guide)*(kmax+2*guide);
  int ref = R->Mode.Base_Medium;

#pragma omp parallel for firstprivate(m, ref) schedule(static)
  for (size_t i=0; i<m; i++) mid[i] = ref;
}

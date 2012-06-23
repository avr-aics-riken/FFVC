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
 * @file IP_Rect.C
 * @brief IP_Rect class
 * @author kero
 */

#include "IP_Rect.h"


// パラメータを取得する
bool IP_Rect::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 分割数の偶数チェックオプション
  label="/Parameter/Intrinsic_Example/Check_Even";

  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Check_Even' in 'Intrinsic_Example'\n");
    return false;
  }

  if     ( !strcasecmp(str.c_str(), "yes") ) {
    even = ON;
  }
  else if( !strcasecmp(str.c_str(), "no") ) {
    even = OFF;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid 'Check_Even' in 'Intrinsic_Example'\n");
    return false;
  }
  
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


// Rectの領域情報を設定する
void IP_Rect::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
{
  wth[0] = pch[0]*(REAL_TYPE)sz[0];
  wth[1] = pch[1]*(REAL_TYPE)sz[1];
  wth[2] = pch[2]*(REAL_TYPE)sz[2];
  
  // 偶数のチェック
  if ( even == ON ) {
    if ( sz[0]/2*2 != sz[0] ) {
      printf("\tDimension size must be even for x direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
      Exit(0);
    }
    if ( sz[1]/2*2 != sz[1] ) {
      printf("\tDimension size must be even for y direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
      Exit(0);
    }
    if ( sz[2]/2*2 != sz[2] ) {
      printf("\tDimension size must be even for z direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
      Exit(0);
    }
  }
}


// 計算領域のセルIDを設定する
void IP_Rect::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  unsigned m, m_sz[3], gd;
  m_sz[0] = size[0];
  m_sz[1] = size[1];
  m_sz[2] = size[2];
  gd = guide;
  
  int id_fluid;
  
  if ( (id_fluid = find_ID_from_Label(mat, Nmax, m_fluid)) == 0 ) {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }

#pragma omp parallel for firstprivate(m_sz, gd) schedule(static)
  for (int k=1; k<=(int)m_sz[2]; k++) {
    for (int j=1; j<=(int)m_sz[1]; j++) {
      for (int i=1; i<=(int)m_sz[0]; i++) {
        m = FBUtility::getFindexS3D(m_sz, gd, i, j, k);
        mid[m] = id_fluid;
      }
    }
  }
}

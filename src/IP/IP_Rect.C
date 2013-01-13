// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   IP_Rect.C
 * @brief  IP_Rect class
 * @author kero
 */

#include "IP_Rect.h"


// #################################################################
// パラメータを取得する
bool IP_Rect::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 分割数の偶数チェックオプション
  label="/Parameter/IntrinsicExample/CheckEven";

  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }

  if     ( !strcasecmp(str.c_str(), "yes") )
  {
    even = ON;
  }
  else if( !strcasecmp(str.c_str(), "no") )
  {
    even = OFF;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  // 媒質指定
  label="/Parameter/IntrinsicExample/FluidMedium";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label="/Parameter/IntrinsicExample/SolidMedium";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  return true;
}


// #################################################################
// Rectの領域情報を設定する
void IP_Rect::setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = pch[1]*(REAL_TYPE)sz[1];
  reg[2] = pch[2]*(REAL_TYPE)sz[2];
  
  // 偶数のチェック
  if ( even == ON )
  {
    if ( sz[0]/2*2 != sz[0] )
    {
      printf("\tDimension size must be even for x direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
      Exit(0);
    }
    if ( sz[1]/2*2 != sz[1] )
    {
      printf("\tDimension size must be even for y direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
      Exit(0);
    }
    if ( sz[2]/2*2 != sz[2] )
    {
      printf("\tDimension size must be even for z direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
      Exit(0);
    }
  }
}


// #################################################################
// 計算領域のセルIDを設定する
void IP_Rect::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int id_fluid;
  
  if ( (id_fluid = R->find_ID_from_Label(mat, Nmax, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }

#pragma omp parallel for firstprivate(ix, jx, kx, gd, id_fluid) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = id_fluid;
      }
    }
  }
}

// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   IP_PPLT2D.C
 * @brief  IP_PPLT2D class
 * @author kero
 */

#include "IP_PPLT2D.h"


// #################################################################
//@brief パラメータを取得する
bool IP_PPLT2D::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 媒質指定
  label = "/Parameter/IntrinsicExample/Fluid_medium";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label = "/Parameter/IntrinsicExample/SolidMedium";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  return true;
}


// #################################################################
// PPLT2Dの領域情報を設定する
void IP_PPLT2D::setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
  // forced
  if (R->Unit.Param != NONDIMENSIONAL)
  {
    Hostonly_ printf("\tError : PPLT2D class is designed for only non-dimensional parameter\n");
    Exit(0);
  }
  
  // Z方向のチェック
  if ( sz[2] != 3 )
  {
    Hostonly_ printf("\tError : The size of Z-direction must be 3.\n");
    Exit(0);
  }
  
  // 領域アスペクト比のチェック
  if ( sz[0] != sz[1]*2 )
  {
    Hostonly_ printf("\tError : The number of division must be 2:1 (=X:Y)\n");
    Exit(0);
  }

  // 二次元の領域設定
  reg[0] =  2.0;
  reg[1] =  1.0;
  org[0] = -1.0;
  org[1] =  0.0;
  
  pch[0] = reg[0] / (REAL_TYPE)sz[0];
  pch[1] = reg[1] / (REAL_TYPE)sz[1];
  
  if ( pch[0] != pch[1] )
  {
    Hostonly_ printf("\tVoxel width must be same between X and Y direction.\n");
    Exit(0);
  }
  
  // Z方向は3層に合わせて調整
  pch[2] =  pch[0];
  reg[2] =  pch[0]*3.0;
  org[2] = -pch[2]*1.5;
}


// #################################################################
// PPLT2Dの計算領域のセルIDを設定する
void IP_PPLT2D::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
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
  
  // Inner
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

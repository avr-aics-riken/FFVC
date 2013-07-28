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

/** 
 * @file   IP_PPLT2D.C
 * @brief  IP_PPLT2D class
 * @author kero
 */

#include "IP_PPLT2D.h"


// #################################################################
/*
 * @brief パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_PPLT2D::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 媒質指定
  label = "/Parameter/IntrinsicExample/FluidMedium";
  
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
/* @brief 領域パラメータを設定する
 * @param [in]     R   Controlクラスのポインタ
 * @param [in]     sz  分割数
 * @param [in,out] org 計算領域の基点
 * @param [in,out] reg 計算領域のbounding boxサイズ
 * @param [in,out] pch セル幅
 */
void IP_PPLT2D::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  // forced
  if (R->Unit.Param != NONDIMENSIONAL)
  {
    Hostonly_ printf("\tError : PPLT2D class is designed for only non-dimensional parameter\n");
    Exit(0);
  }
  
  // 2D
  mode = dim_2d;
  
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
  
  // Z方向は3層に合わせて調整
  pch[2] =  pch[0];
  reg[2] =  pch[0]*3.0;
  org[2] = -pch[2]*1.5;
}


// #################################################################
/*
 * @brief 並行平板の計算領域のセルIDを設定する
 * @param [in,out] mid      媒質情報の配列
 * @param [in]     R        Controlクラスのポインタ
 * @param [in]     G_org    グローバルな原点（無次元）
 * @param [in]     NoMedium 媒質数
 * @param [in]     mat      MediumListクラスのポインタ
 */
void IP_PPLT2D::setup(int* mid, Control* R, REAL_TYPE* G_org, const int NoMedium, const MediumList* mat)
{
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int id_fluid;
  
  if ( (id_fluid = R->findIDfromLabel(mat, NoMedium, m_fluid)) == 0 )
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

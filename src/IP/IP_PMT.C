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
 * @file   IP_PMT.C
 * @brief  IP_PMT class
 * @author kero
 */

#include "IP_PMT.h"


// #################################################################
/*
 * @brief パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_PMT::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
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
/*
 * @brief 領域パラメータを設定する
 * @param [in]     R   Controlクラスのポインタ
 * @param [in]     sz  分割数
 * @param [in,out] org 計算領域の基点
 * @param [in,out] reg 計算領域のbounding boxサイズ
 * @param [in,out] pch セル幅
 */
void IP_PMT::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{  
  // forced
  if (R->Unit.Param != NONDIMENSIONAL)
  {
    Hostonly_ printf("\tError : PMT class is designed for only non-dimensional parameter\n");
    Exit(0);
  }
  
  // 性能測定モードをOnにする
  R->Hide.PM_Test = ON;
  
  // 偶数チェック
  even = ON;
  
  org[0] = -0.5*reg[0];
  org[1] = -0.5*reg[1];
  org[2] = -0.5*reg[2];
}


// #################################################################
/*
 * @brief Cavityの計算領域のセルIDを設定する
 * @param [in,out] mid      媒質情報の配列
 * @param [in]     R        Controlクラスのポインタ
 * @param [in]     G_org    グローバルな原点（無次元）
 * @param [in]     NoMedium 媒質数
 * @param [in]     mat      MediumListクラスのポインタ
 */
void IP_PMT::setup(int* mid, Control* R, REAL_TYPE* G_org, const int NoMedium, const MediumList* mat)
{
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int id_fluid;
  
  if ( (id_fluid = R->find_ID_from_Label(mat, NoMedium, m_fluid)) == 0 )
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

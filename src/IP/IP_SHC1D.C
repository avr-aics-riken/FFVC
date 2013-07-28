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
 @file   IP_SHC1D.C
 @brief  IP_SHC1D class
 @author kero
 */

#include "IP_SHC1D.h"


// #################################################################
/*
 * @brief パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_SHC1D::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 媒質指定
  label="/Parameter/IntrinsicExample/FluidMedium";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  
  label="/Parameter/IntrinsicExample/SolidMedium";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  
  label="/Parameter/IntrinsicExample/InactiveSolidMedium";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_inactive = str;
  
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
void IP_SHC1D::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
  // forced
  if (R->Unit.Param != DIMENSIONAL) {
    Hostonly_ printf("\tError : SHC1D class is designed for only DIMENSIONAL parameter\n");
    Exit(0);
  }
	
  int i,j,k;
  i = sz[0] -2;
  j = sz[1];
  k = sz[2];
  
  if ( (i != j) || (j != k) || (k != i) )
  {
    printf("\tDimension size error (%d %d %d)\n", sz[0], sz[1], sz[2]);
    Exit(0);
  }
  
	pch[0] = 1.0 / (REAL_TYPE)i;
  pch[1] = pch[0];
  pch[2] = pch[0];
	
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = (REAL_TYPE)j * pch[1];
  reg[2] = (REAL_TYPE)k * pch[2];
  org[0] = -1.0*pch[0];
  org[1] = -0.5*reg[1];
  org[2] = -0.5*reg[2];

}


// #################################################################
/*
 * @brief 矩形の計算領域のセルIDを設定する
 * @param [in,out] mid      媒質情報の配列
 * @param [in]     R        Controlクラスのポインタ
 * @param [in]     G_org    グローバルな原点（無次元）
 * @param [in]     NoMedium 媒質数
 * @param [in]     mat      MediumListクラスのポインタ
 */
void IP_SHC1D::setup(int* mid, Control* R, REAL_TYPE* G_org, const int NoMedium, const MediumList* mat)
{
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  int id_fluid, id_solid, id_solid_inactive;
  
  if ( (id_fluid = R->findIDfromLabel(mat, NoMedium, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  if ( (id_solid = R->findIDfromLabel(mat, NoMedium, m_solid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_solid.c_str());
    Exit(0);
  }
  
  if ( (id_solid_inactive = R->findIDfromLabel(mat, NoMedium, m_inactive)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_inactive.c_str());
    Exit(0);
  }
  
  
  // SOLID_CONDUCTIONモードでは，fluidセルはInactive
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        mid[_F_IDX_S3D(i, j, k, ix, jx, kx, gd)] = id_fluid;
      }
    }
  }
  
  
  // fin
  int mj = jx/2 + 1;
  int mk = kx/2 + 1;
  
  for (int i=2; i<=ix-1; i++) {
    mid[_F_IDX_S3D(i, mj, mk, ix, jx, kx, gd)] = id_solid;
  }
  

  // vertical wall >> Inactive
  int mi = 1;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      mid[_F_IDX_S3D(mi, j, k, ix, jx, kx, gd)] = id_solid_inactive;
    }
  } 
  
}

// #################################################################
/*
 * @brief 矩形の計算領域のセルIDを設定する
 * @param [in,out] bid  カット点の境界条件ID
 */
void IP_SHC1D::setup_bc(int* bid)
{
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // カット情報
  int id_isothermal   = 4;
  int id_adiabatic    = 5;
  int id_heattransfer = 6;
  
  // Direction
  // 0 - w
  // 1 - e
  // 2 - s
  // 3 - n
  // 4 - b
  // 5 - t
  
  // ISOTHERMAL on vertical wall
  int mi = 2;
  int mj = jx/2 + 1;
  int mk = kx/2 + 1;
  
  // i=2のw側にBC
  bid[_F_IDX_S3D(mi, mj, mk, ix, jx, kx, gd)] = id_isothermal << 0;
  
  // i=ix-1のe側に断熱
  mi = ix-1;
  bid[_F_IDX_S3D(mi, mj, mk, ix, jx, kx, gd)] = id_adiabatic << 5;
  
  for (int i=2; i<=ix-1; i++) {
    int b = 0;
    b |= id_heattransfer << 10; // s
    b |= id_heattransfer << 15; // n
    b |= id_heattransfer << 20; // b
    b |= id_heattransfer << 25; // t
    bid[_F_IDX_S3D(i, mj, mk, ix, jx, kx, gd)] = b;
  }
  
}

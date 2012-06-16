/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

/**
 @file IP_SHC1D.C
 @brief IP_SHC1D class
 @author keno, FSI Team, RIKEN
 */

#include "IP_SHC1D.h"

//@brief パラメータを取得する
bool IP_SHC1D::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 媒質指定
  label="/Parameter/Intrinsic_Example/Inactive_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Inactive_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_inactive = str;
  
  label="/Parameter/Intrinsic_Example/Fluid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Fluid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_fluid = str;
  
  label="/Parameter/Intrinsic_Example/Fin_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Fin_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_fin = str;
  
  label="/Parameter/Intrinsic_Example/isothermal_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'driver_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_isothermal = str;
  
  label="/Parameter/Intrinsic_Example/adiabatic_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'adiabatic_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_adiabatic = str;
  
  return true;
}

/**
 @fn void IP_SHC1D::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
 @brief Cavity3Dの領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
void IP_SHC1D::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
{
  // forced
  if (R->Unit.Param != DIMENSIONAL) {
    Hostonly_ printf("\tError : SHC1D class is designed for only dimensional parameter\n");
    Exit(0);
  }
	
	pch[0] = 1.0 / (REAL_TYPE)(sz[0]-2);
  pch[1] = pch[0];
  pch[2] = pch[0];
	
  wth[0] = pch[0]*(REAL_TYPE)sz[0];
  wth[1] = 5.0*pch[1];
  wth[2] = 5.0*pch[2];
  org[0] = -1.0*pch[0];
  org[1] = -2.5*pch[1];
  org[2] = -2.5*pch[2];

  // Setting depends on Example,  INTRINSIC
  if ( (sz[1] != 5) || (sz[2] != 5) ) {
    printf("\tSHC1D case requires jmax = kmax = 5\n");
    Exit(0);
  }
}

/**
 @fn void IP_SHC1D::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, IDtable* itbl)
 @brief モデルIDのセットアップ
 @param mid ボクセル配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 @param itbl IDtable
 */
void IP_SHC1D::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, IDtable* itbl)
{
  unsigned i, j, k, m;
  int gd;
  const int mid_inactive=600;
  const int mid_fluid=1;
  const int mid_fin=610;
  const int mid_isothermal=500;
  const int mid_adiabatic=520;
  
  gd = (int)guide;

  // Initialize and for outer solid
  for (k=1-gd; k<=(int)kmax+gd; k++) {
    for (j=1-gd; j<=(int)jmax+gd; j++) {
      for (i=1-gd; i<=(int)imax+gd; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = mid_inactive;
      }
    }
  }

  // inner fluid
  for (k=2; k<=kmax-1; k++) {
    for (j=2; j<=jmax-1; j++) {
      for (i=2; i<=imax-1; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = mid_fluid;
      }
    }
  }

  // fin
  j = k = 3;
  for (i=3; i<=imax-2; i++) {
    m = FBUtility::getFindexS3D(size, guide, i, j, k);
    mid[m] = mid_fin;
  }

  // iso-thermal fin
  i = 2;
  j = k = 3;
  m = FBUtility::getFindexS3D(size, guide, i, j, k);
  mid[m] = mid_isothermal;
  
  // adiabatic fin
  i = imax-1;
  j = k = 3;
  m = FBUtility::getFindexS3D(size, guide, i, j, k);
  mid[m] = mid_adiabatic;  
  
}

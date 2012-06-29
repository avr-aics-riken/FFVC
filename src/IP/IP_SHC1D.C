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
 @file IP_SHC1D.C
 @brief IP_SHC1D class
 @author kero
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


// 領域情報を設定する
void IP_SHC1D::setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
  // forced
  if (R->Unit.Param != DIMENSIONAL) {
    Hostonly_ printf("\tError : SHC1D class is designed for only dimensional parameter\n");
    Exit(0);
  }
	
	pch[0] = 1.0 / (REAL_TYPE)(sz[0]-2);
  pch[1] = pch[0];
  pch[2] = pch[0];
	
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = 5.0*pch[1];
  reg[2] = 5.0*pch[2];
  org[0] = -1.0*pch[0];
  org[1] = -2.5*pch[1];
  org[2] = -2.5*pch[2];

  // Setting depends on Example,  INTRINSIC
  if ( (sz[1] != 5) || (sz[2] != 5) ) {
    printf("\tSHC1D case requires jmax = kmax = 5\n");
    Exit(0);
  }
}


// モデルIDのセットアップ
void IP_SHC1D::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  const int mid_inactive=600;
  const int mid_fluid=1;
  const int mid_fin=610;
  const int mid_isothermal=500;
  const int mid_adiabatic=520;
  
  size_t m;
  
  // ローカルにコピー
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gd = guide;

  // Initialize and for outer solid
  for (int k=1-gd; k<=kmax+gd; k++) {
    for (int j=1-gd; j<=jmax+gd; j++) {
      for (int i=1-gd; i<=imax+gd; i++) {
        m = FBUtility::getFindexS3D(size, gd, i, j, k);
        mid[m] = mid_inactive;
      }
    }
  }

  // inner fluid
  for (int k=2; k<=kmax-1; k++) {
    for (int j=2; j<=jmax-1; j++) {
      for (int i=2; i<=imax-1; i++) {
        m = FBUtility::getFindexS3D(size, gd, i, j, k);
        mid[m] = mid_fluid;
      }
    }
  }

  // fin
  int i, j, k;
  
  j = k = 3;
  for (i=3; i<=imax-2; i++) {
    m = FBUtility::getFindexS3D(size, gd, i, j, k);
    mid[m] = mid_fin;
  }

  // iso-thermal fin
  i = 2;
  j = k = 3;
  m = FBUtility::getFindexS3D(size, gd, i, j, k);
  mid[m] = mid_isothermal;
  
  // adiabatic fin
  i = imax-1;
  j = k = 3;
  m = FBUtility::getFindexS3D(size, gd, i, j, k);
  mid[m] = mid_adiabatic;  
  
}

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Rect.C
//@brief IP_Rect class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_Rect.h"

/**
 @fn bool IP_Rect::getXML(SklSolverConfig* CF, Control* R)
 @brief パラメータを取得する
 @param CF コンフィギュレーションツリー
 @param R Controlクラスのポインタ
 */
bool IP_Rect::getXML(SklSolverConfig* CF, Control* R)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  const char *str=NULL;
  
  if ( !(elemTop = CF->GetTop(PARAMETER)) ) return false;
  
  if( !(elmL1 = elemTop->GetElemFirst("Intrinsic_Example")) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of 'Intrinsic_Example'\n");
    return false;
  }
  
  // 分割数の偶数チェックオプション
  if ( !elmL1->GetValue("Check_Even", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Check_Even' in 'Intrinsic_Example'\n");
    return false;
  }
  if ( !strcasecmp(str, "yes") ) {
    even = ON;
  }
  else if ( !strcasecmp(str, "no") ) {
    even = OFF;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid 'Check_Even' in 'Intrinsic_Example'\n");
    return false;
  }
  
  return true;
}

/**
 @fn bool IP_Rect::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
 @brief Rectの領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
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

/**
 @fn void IP_Rect::setup(int* mid, Control* R, REAL_TYPE* G_org)
 @brief 計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_Rect::setup(int* mid, Control* R, REAL_TYPE* G_org)
{
  unsigned m, m_sz[3], gd;
  m_sz[0] = size[0];
  m_sz[1] = size[1];
  m_sz[2] = size[2];
  gd = guide;

#pragma omp parallel for firstprivate(m_sz, gd) schedule(static)
  for (int k=1; k<=(int)m_sz[2]; k++) {
    for (int j=1; j<=(int)m_sz[1]; j++) {
      for (int i=1; i<=(int)m_sz[0]; i++) {
        m = FBUtility::getFindexS3D(m_sz, gd, i, j, k);
        mid[m] = 1;
      }
    }
  }
}

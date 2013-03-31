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
/*
 * @brief パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_Rect::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 2D or 3D mode
  label="/Parameter/IntrinsicExample/Dimension";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if     ( !strcasecmp(str.c_str(), "2d") ) {
    mode = dim_2d;
  }
  else if( !strcasecmp(str.c_str(), "3d") ) {
    mode = dim_3d;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  
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
/* @brief パラメータの表示
 * @param [in] fp ファイルポインタ
 * @param [in] R  コントロールクラスのポインタ
 */
void IP_Rect::printPara(FILE* fp, const Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Rectangular Class Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tCheck even number                  :  %s\n", (even == ON)?"Yes":"No");

}


// #################################################################
/* @brief 領域パラメータを設定する
 * @param [in]     R   Controlクラスのポインタ
 * @param [in]     sz  分割数
 * @param [in,out] org 計算領域の基点
 * @param [in,out] reg 計算領域のbounding boxサイズ
 * @param [in,out] pch セル幅
 */
void IP_Rect::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
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
    if ( (mode == dim_3d) && (sz[2]/2*2 != sz[2]) )
    {
      printf("\tDimension size must be even for z direction (%d %d %d)\n", sz[0], sz[1], sz[2]);
      Exit(0);
    }
  }
  
  
}


// #################################################################
/* @brief 矩形の計算領域のセルIDを設定する
 * @param [in,out] mid   媒質情報の配列
 * @param [in]     R     Controlクラスのポインタ
 * @param [in]     G_org グローバルな原点（無次元）
 * @param [in]     Nmax  MediumList配列のサイズ
 * @param [in]     mat   MediumListクラスのポインタ
 */
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

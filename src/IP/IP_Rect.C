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
 * @file   IP_Rect.C
 * @brief  IP_Rect class
 * @author kero
 */

#include "IP_Rect.h"


// #################################################################
/* @brief パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_Rect::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  
  // 2D or 3D mode
  label="/IntrinsicExample/Dimension";
  
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
  label="/IntrinsicExample/CheckEven";

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
  label="/IntrinsicExample/FluidMedium";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label="/IntrinsicExample/SolidMedium";
  
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
/* @brief 矩形の計算領域のセルIDを設定する
 * @param [in,out] mid      媒質情報の配列
 * @param [in]     R        Controlクラスのポインタ
 * @param [in]     G_org    グローバルな原点（無次元）
 * @param [in]     NoMedium 媒質数
 * @param [in]     mat      MediumListクラスのポインタ
 */
void IP_Rect::setup(int* mid, Control* R, REAL_TYPE* G_org, const int NoMedium, const MediumList* mat)
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

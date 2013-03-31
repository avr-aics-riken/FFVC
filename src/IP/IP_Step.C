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
 * @file   IP_Step.C
 * @brief  IP_Step class
 * @author kero
 */

#include "IP_Step.h"


// #################################################################
/*
 * @brief パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_Step::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
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
  
  // x-dir step
  label="/Parameter/IntrinsicExample/StepLength";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  width = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // z-dir step
  label="/Parameter/IntrinsicExample/StepHeight";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  height = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label="/Parameter/IntrinsicExample/DriverLength";
  if ( tpCntl->GetValue(label, &ct ) ) {
    drv_length = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if ( drv_length < 0.0 ) {
    Hostonly_ stamped_printf("\tError : Value of '%s' must be positive.\n", label.c_str());
    return false;
  }
  
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
  
  // Only driver is specified
  if ( drv_length > 0.0 )
  {
    label="/Parameter/IntrinsicExample/DriverMedium";
    if ( !(tpCntl->GetValue(label, &str )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver = str;
    
    label="/Parameter/IntrinsicExample/DriverFaceMedium";
    if ( !(tpCntl->GetValue(label, &str )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver_face = str;
  }
  
  return true;
}


// #################################################################
/**
 * @brief パラメータの表示
 * @param [in] fp ファイルポインタ
 * @param [in] R  コントロールクラスのポインタ
 */
void IP_Step::printPara(FILE* fp, const Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Backstep Class Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tStep Width (x-dir.)    [m] / [-]   : %12.5e / %12.5e\n", width, width/RefL);
  fprintf(fp,"\tStep Height(z-dir.)    [m] / [-]   : %12.5e / %12.5e\n", height, height/RefL);
  if ( drv_length > 0.0 ) {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
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
void IP_Step::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = pch[1]*(REAL_TYPE)sz[1];
  reg[2] = pch[2]*(REAL_TYPE)sz[2];
  
  // チェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  
  if ( (reg[0] != (REAL_TYPE)sz[0]*pch[0]) ||
       (reg[1] != (REAL_TYPE)sz[1]*pch[1]) ||
       (reg[2] != (REAL_TYPE)sz[2]*pch[2]) ) {
    Hostonly_ printf("Error : Invalid parameters among 'GlobalRegion', 'GlobalPitch', and 'GlobalVoxel' in DomainInfo section.\n");
    Exit(0);
  }

  // 次元とサイズ
  if (mode == dim_2d) {
    if (size[2] != 3) {
      Hostonly_ printf("Error : VoxelSize kmax must be 3 if 2-dimensional.\n");
    }
  }
}


// #################################################################
/*
 * @brief 計算領域のセルIDを設定する
 * @param [in,out] mid   媒質情報の配列
 * @param [in]     R     Controlクラスのポインタ
 * @param [in]     G_org グローバルな原点（無次元）
 * @param [in]     Nmax  Controlクラスのポインタ
 * @param [in]     mat   MediumListクラスのポインタ
 */
void IP_Step::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int mid_fluid, mid_solid, mid_driver, mid_driver_face;
  
  // 流体
  if ( (mid_fluid = R->find_ID_from_Label(mat, Nmax, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  // 固体
  if ( (mid_solid = R->find_ID_from_Label(mat, Nmax, m_solid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_solid.c_str());
    Exit(0);
  }
  
  if ( drv_length > 0.0 )
  {
    // ドライバ部
    if ( (mid_driver = R->find_ID_from_Label(mat, Nmax, m_driver)) == 0 )
    {
      Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_driver.c_str());
      Exit(0);
    }
    
    // ドライバ流出面
    if ( (mid_driver_face = R->find_ID_from_Label(mat, Nmax, m_driver_face)) == 0 )
    {
      Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_driver_face.c_str());
      Exit(0);
    }
  }


  REAL_TYPE ox, oy, oz, dh;
  
  // ノードローカルの無次元値
  ox = origin[0];
  oy = origin[1];
  oz = origin[2];
  dh = deltaX;
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // length, widthなどは有次元値
  REAL_TYPE len = G_origin[0] + (drv_length+width)/R->RefLength; // グローバルな無次元位置
  REAL_TYPE ht  = G_origin[1] + height/R->RefLength;
  
  // Initialize  内部領域をfluidにしておく
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = mid_fluid;
      }
    }
  }
  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 ) {
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_driver, ox, dh, len) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          REAL_TYPE x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) mid[m] = mid_driver;
        }
      }
    }
    
  }
  
  // ドライバの下流面にIDを設定
  if ( drv_length > 0.0 )
  {
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_driver, mid_fluid, mid_driver_face) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i,   j, k, ix, jx, kx, gd);
          size_t m1= _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
          if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
            mid[m] = mid_driver_face;
          }
        }
      }
    }
    
  }

  // ステップ部分を上書き
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_solid, ox, oy, dh, len, ht) \
schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*dh + dh*(i-1);
        REAL_TYPE y = oy + 0.5*dh + dh*(j-1);
        if ( (x < len) && (y < ht) )
        {
          mid[m] = mid_solid;
        }
      }
    }
  }
  
}

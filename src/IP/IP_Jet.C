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
 * @file   IP_Jet.C
 * @brief  IP_Jet class
 * @author kero
 */

#include "IP_Jet.h"


// #################################################################
/** パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_Jet::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;

  
  // 2D or 3D mode
  label="/Parameter/IntrinsicExample/Mode";
  
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
  
  
  // Core
  label="/Parameter/IntrinsicExample/RadiusOfCore";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  r0 = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  label="/Parameter/IntrinsicExample/AngularVelocityOfCore";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  omg0 = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL / RefV;
  }
  
  
  // Ring1
  label="/Parameter/IntrinsicExample/InnerRadiusOfRing1";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  r1i = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  label="/Parameter/IntrinsicExample/OuterRadiusOfRing1";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  r1o = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  label="/Parameter/IntrinsicExample/AngularVelocityOfRing1";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  omg1 = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL / RefV;
  }
  
  // Ring2
  label="/Parameter/IntrinsicExample/InnerRadiusOfRing2";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  r2i = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  label="/Parameter/IntrinsicExample/OuterRadiusOfRing2";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  r2o = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  label="/Parameter/IntrinsicExample/AngularVelocityOfRing2";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  omg2 = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL / RefV;
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
/** 領域を設定する
 * @param [in]     R   Controlクラスのポインタ
 * @param [in]     sz  分割数
 * @param [in,out] org 計算領域の基点
 * @param [in,out] reg 計算領域のbounding boxサイズ
 * @param [in,out] pch セル幅
 */
void IP_Jet::setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  RefV = R->RefVelocity;
  
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = pch[1]*(REAL_TYPE)sz[1];
  reg[2] = pch[2]*(REAL_TYPE)sz[2];
  
  // チェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) )
  {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  
  if ( ((int)(reg[0]/pch[0]) != sz[0]) ||
       ((int)(reg[1]/pch[1]) != sz[1]) ||
       ((int)(reg[2]/pch[2]) != sz[2]) )
  {
    Hostonly_ printf("Error : Invalid parameters among 'GlobalRegion', 'GlobalPitch', and 'GlobalVoxel' in DomainInfo section.\n");
    Exit(0);
  }
  
  // 次元とサイズ
  if (mode == dim_2d)
  {
    if (size[2] != 3)
    {
      Hostonly_ printf("Error : VoxelSize kmax must be 3 if 2-dimensional.\n");
    }
  }
  
}



// #################################################################
/**
 * @brief パラメータの表示
 * @param [in] fp ファイルポインタ
 * @param [in] R  コントロールクラスのポインタ
 */
void IP_Jet::printPara(FILE* fp, const Control* R)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  // パターンのチェック
  
  pat_0 = ( r0 == 0.0 ) ? 0 : 1;
  
  pat_1 = ( (r1i == 0.0) && (r1o == 0.0) ) ? 0 : 1;
  
  pat_2 = ( (r2i == 0.0) && (r2o == 0.0) ) ? 0 : 1;
  
  if ( pat_1 == 1 )
  {
    if ( r1i < r0 )
    {
      stamped_printf("\tInner Radius of Ring1 must be greater than Radius of Core\n");
      Exit(0);
    }
    
    if ( r1o < r1i )
    {
      stamped_printf("\tInner Radius of Ring1 must be smaller than Outer Radius\n");
      Exit(0);
    }
  }

  if ( pat_2 == 1)
  {
    if ( r2i < r1o )
    {
      stamped_printf("\tInner Radius of Ring2 must be greater than Outer Radius of Ring1\n");
      Exit(0);
    }
    
    if ( r2o < r2i )
    {
      stamped_printf("\tInner Radius of Ring2 must be smaller than Outer Radius\n");
      Exit(0);
    }
  }
  
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Backstep Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tCore  Radius           [m]   / [-] : %12.5e / %12.5e\n", r0, r0/RefL);
  fprintf(fp,"\t      Angular Velocity [m/s] / [-] : %12.5e / %12.5e\n", omg0, omg0*RefL/RefV);
  
  if ( pat_1 == 1 )
  {
    fprintf(fp,"\tRing1 Inner Radius     [m]   / [-] : %12.5e / %12.5e\n", r1i, r1i/RefL);
    fprintf(fp,"\t      Outer Radius     [m]   / [-] : %12.5e / %12.5e\n", r1o, r1o/RefL);
    fprintf(fp,"\t      Angular Velocity [m/s] / [-] : %12.5e / %12.5e\n", omg1, omg1*RefL/RefV);
  }
  
  if ( pat_2 == 1)
  {
    fprintf(fp,"\tRing2 Inner Radius     [m]   / [-] : %12.5e / %12.5e\n", r2i, r2i/RefL);
    fprintf(fp,"\t      Outer Radius     [m]   / [-] : %12.5e / %12.5e\n", r2o, r2o/RefL);
    fprintf(fp,"\t      Angular Velocity [m/s] / [-] : %12.5e / %12.5e\n", omg2, omg2*RefL/RefV);
  }
}



// #################################################################
/** 計算領域のセルIDを設定する
 * @param [in,out] mid   媒質情報の配列
 * @param [in]     R     Controlクラスのポインタ
 * @param [in]     G_org グローバルな原点（無次元）
 * @param [in]     Nmax  Controlクラスのポインタ
 * @param [in]     mat   MediumListクラスのポインタ
 */
void IP_Jet::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
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

  
  // ノードローカルの無次元値
  REAL_TYPE ox = origin[0];
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  REAL_TYPE dh = deltaX;
  
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

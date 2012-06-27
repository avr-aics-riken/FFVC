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
 * @file IP_Step.C
 * @brief IP_Step class
 * @author kero
 */

#include "IP_Step.h"


//@brief パラメータを取得する
bool IP_Step::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  // 2D or 3D mode
  label="/Parameter/Intrinsic_Example/mode";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'mode' in 'Intrinsic_Example'\n");
    return false;
  }
  
  if     ( !strcasecmp(str.c_str(), "2d") ) {
    mode = dim_2d;
  }
  else if( !strcasecmp(str.c_str(), "3d") ) {
    mode = dim_3d;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid 'mode' in 'Intrinsic_Example'\n");
    return false;
  }
  
  // x-dir step
  label="/Parameter/Intrinsic_Example/width";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'width' in 'Intrinsic_Example'\n");
    return false;
  }
  else{
	  width = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // z-dir step
  label="/Parameter/Intrinsic_Example/height";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'height' in 'Intrinsic_Example'\n");
    return false;
  }
  else{
	  height = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label="/Parameter/Intrinsic_Example/Driver";
  if ( tpCntl->GetValue(label, &ct ) ) {
    drv_length = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Driver' in 'Intrinsic_Example'\n");
    return false;
  }
  
  if ( drv_length < 0.0 ) {
    Hostonly_ stamped_printf("\tError : Value of 'Driver' in 'Intrinsic_Example' must be positive.\n");
    return false;
  }
  
  // 媒質指定
  label="/Parameter/Intrinsic_Example/Fluid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Fluid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_fluid = str;
  
  label="/Parameter/Intrinsic_Example/Solid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Solid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_solid = str;
  
  label="/Parameter/Intrinsic_Example/driver_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'driver_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_driver = str;
  
  label="/Parameter/Intrinsic_Example/driver_face_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'driver_face_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_driver_face = str;
  
  return true;
}


// パラメータの表示
void IP_Step::printPara(FILE* fp, Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Backstep Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tStep Width (x-dir.)    [m] / [-]   : %12.5e / %12.5e\n", width, width/RefL);
  fprintf(fp,"\tStep Height(z-dir.)    [m] / [-]   : %12.5e / %12.5e\n", height, height/RefL);
  if ( drv_length > 0.0 ) {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
}


// 領域情報を設定する
void IP_Step::setDomain(Control* R, const int* sz, const REAL_TYPE* org, REAL_TYPE* reg, const REAL_TYPE* pch)
{
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = pch[1]*(REAL_TYPE)sz[1];
  reg[2] = pch[2]*(REAL_TYPE)sz[2];
  
  // チェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  if ( ((int)(reg[0]/pch[0]) != sz[0]) ||
       ((int)(reg[1]/pch[1]) != sz[1]) ||
       ((int)(reg[2]/pch[2]) != sz[2]) ) {
    Hostonly_ printf("Error : Invalid parameters among 'VoxelSize', 'VoxelPitch', and 'VoxelWidth' in DomainInfo section.\n");
    Exit(0);
  }

  // 次元とサイズ
  if (mode == dim_2d) {
    if (size[2] != 3) {
      Hostonly_ printf("Error : VoxelSize kmax must be 3 if 2-dimensional.\n");
    }
  }
}


// 計算領域のセルIDを設定する
void IP_Step::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  int mid_driver=3;       /// ドライバ部
  int mid_driver_face=4;  /// ドライバ流出面

  REAL_TYPE x, y, z, dh, len, ht;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  REAL_TYPE ox_g, oy_g, oz_g;
  
  // ノードローカルの無次元値
  ox = R->org[0];
  oy = R->org[1];
  oz = R->org[2];
  Lx = R->Lbx[0];
  Ly = R->Lbx[1];
  Lz = R->Lbx[2];
  dh = R->dh;

  ox_g = G_org[0];
  oy_g = G_org[1];
  oz_g = G_org[2];
  
  size_t m;
  
  // ローカルにコピー
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gd = guide;

  // length, widthなどは有次元値
  len = ox_g + (drv_length+width)/R->RefLength; // グローバルな無次元位置
  ht  = oy_g + height/R->RefLength;
  
  // Initialize  内部領域をfluidにしておく
  for (int k=1; k<=kmax; k++) {
    for (int j=1; j<=jmax; j++) {
      for (int i=1; i<=imax; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = mid_fluid;
      }
    }
  }
  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 ) {
    for (int k=1; k<=kmax; k++) {
      for (int j=1; j<=jmax; j++) {
        for (int i=1; i<=imax; i++) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) mid[m] = mid_driver;
        }
      }
    }  
  }
  
  // ドライバの下流面にIDを設定
  if ( drv_length > 0.0 ) {
    
    size_t m1;
    
    for (int k=1; k<=kmax; k++) {
      for (int j=1; j<=jmax; j++) {
        for (int i=1; i<=imax; i++) {
          m = FBUtility::getFindexS3D(size, guide, i,   j, k);
          m1= FBUtility::getFindexS3D(size, guide, i+1, j, k);
          if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
            mid[m] = mid_driver_face;
          }
        }
      }
    }    
  }

  // ステップ部分を上書き
  for (int k=1; k<=kmax; k++) {
    for (int j=1; j<=jmax; j++) {
      for (int i=1; i<=imax; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        x = ox + 0.5*dh + dh*(i-1);
        y = oy + 0.5*dh + dh*(j-1);
        if ( (x < len) && (y < ht) ) {
          mid[m] = mid_solid;
        }
      }
    }
  }
}

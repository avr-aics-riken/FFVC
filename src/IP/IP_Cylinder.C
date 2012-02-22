/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Cylinder.C
//@brief IP_Cylinder class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_Cylinder.h"

/**
 @fn bool IP_Cylinder::getXML(SklSolverConfig* CF, Control* R)
 @brief Cylinderのパラメータを取得する
 @param CF コンフィギュレーションツリー
 @param R Controlクラスのポインタ
 */
bool IP_Cylinder::getXML(SklSolverConfig* CF, Control* R)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  const char *str=NULL;
  REAL_TYPE ct=0.0;
  
  if ( !(elemTop = CF->GetTop(PARAMETER)) ) return false;
  
  if( !(elmL1 = elemTop->GetElemFirst("Intrinsic_Example")) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of 'Intrinsic_Example'\n");
    return false;
  }
  
  // 次元
  if ( !elmL1->GetValue("mode", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Mode' in 'Intrinsic_Example'\n");
    return false;
  }
  if ( !strcasecmp(str, "2d") ) {
    mode = dim_2d;
  }
  else if ( !strcasecmp(str, "3d") ) {
    mode = dim_3d;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid mode in 'Intrinsic_Example'\n");
    return false;
  }
  
  // x-dir step
  if ( elmL1->GetValue(CfgIdt("width"), &ct) ) {
    width = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Depth' in 'Intrinsic_Example'\n");
    return false;
  }
  
  // z-dir step
  if ( elmL1->GetValue(CfgIdt("height"), &ct) ) {
    height = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Height' in 'Intrinsic_Example'\n");
    return false;
  }
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  if ( elmL1->GetValue(CfgIdt("Driver"), &ct) ) {
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
  
  return true;
}

/**
 @fn void IP_Cylinder::printPara(FILE* fp, Control* R)
 @brief パラメータの表示
 @param fp ファイルポインタ
 @param R コントロールクラスのポインタ
 */
void IP_Cylinder::printPara(FILE* fp, Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    assert(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Duct Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tStep Width (x-dir.)    [m] / [-]   : %12.5e / %12.5e\n", width, width/RefL);
  fprintf(fp,"\tStep Height(z-dir.)    [m] / [-]   : %12.5e / %12.5e\n", height, height/RefL);
  if ( drv_length > 0.0 ) {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
}

/**
 @fn bool IP_Cylinder::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
 @brief Cylinderの領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
void IP_Cylinder::setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3])
{
  wth[0] = pch[0]*(REAL_TYPE)sz[0];
  wth[1] = pch[1]*(REAL_TYPE)sz[1];
  wth[2] = pch[2]*(REAL_TYPE)sz[2];
  
  // チェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    assert(0);
  }
  if ( ((unsigned)(wth[0]/pch[0]) != sz[0]) ||
       ((unsigned)(wth[1]/pch[1]) != sz[1]) ||
       ((unsigned)(wth[2]/pch[2]) != sz[2]) ) {
    Hostonly_ printf("Error : Invalid parameters among 'VoxelSize', 'VoxelPitch', and 'VoxelWidth' in DomainInfo section.\n");
    assert(0);
  }

  // 次元とサイズ
  if (mode == dim_2d) {
    if (kmax != 3) {
      Hostonly_ printf("Error : VoxelSize kmax must be 3 if 2-dimensional.\n");
    }
  }
}

/**
 @fn void IP_Cylinder::setup(int* mid, Control* R, REAL_TYPE* G_org)
 @brief Cylinderの計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_Cylinder::setup(int* mid, Control* R, REAL_TYPE* G_org)
{
  int i,j,k, gd;
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  int mid_driver=3;       /// ドライバ部
  int mid_driver_face=4;  /// ドライバ流出面
  unsigned m;
  REAL_TYPE x, y, z, dh, len, ht;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  
  ox = R->org[0];
  oy = R->org[1];
  oz = R->org[2];
  Lx = R->Lbx[0];
  Ly = R->Lbx[1];
  Lz = R->Lbx[2];
  dh = R->dh;
  gd = (int)guide;

  len= drv_length/R->RefLength;
  ht = height/R->RefLength;
  
  // Initialize  内部領域をfluidにしておく
  for (k=1; k<=(int)kmax; k++) { 
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
        m = SklUtil::getFindexS3D(size, guide, i, j, k);
        mid[m] = mid_fluid;
      }
    }
  }
  
  // ドライバ部分
  if ( drv_length > 0.0 ) {
    if ( pn.nID[X_MINUS] < 0 ) {
      for (k=1; k<=(int)kmax; k++) {
        for (j=1; j<=(int)jmax; j++) {
          for (i=1; i<=(int)imax; i++) {
            m = SklUtil::getFindexS3D(size, guide, i, j, k);
            x = ox + 0.5*dh + dh*(i-1);
            if ( x < ox+len ) mid[m] = mid_driver;
          }
        }
      }
    }     
  }
  
  // ドライバの下流面にIDを設定
  if ( drv_length > 0.0 ) {
    
    unsigned m1;
    for (k=1; k<=(int)kmax; k++) {
      for (j=1; j<=(int)jmax; j++) {
        for (i=1; i<=(int)imax; i++) {
          m = SklUtil::getFindexS3D(size, guide, i,   j, k);
          m1= SklUtil::getFindexS3D(size, guide, i+1, j, k);
          if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
            mid[m] = mid_driver_face;
          }
        }
      }
    }    
  }

  // ステップ部分を上書き
  for (k=1; k<=(int)kmax; k++) {
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
        m = SklUtil::getFindexS3D(size, guide, i, j, k);
        x = ox + 0.5*dh + dh*(i-1);
        y = oy + 0.5*dh + dh*(j-1);
        if ( (x < ox+len) && (y < oy+ht) ) {
          mid[m] = mid_solid;
        }
      }
    }
  }
}

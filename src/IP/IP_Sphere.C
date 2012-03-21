/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Shere.C
//@brief IP_Shere class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_Sphere.h"

/**
 @fn bool IP_Sphere::getXML(SklSolverConfig* CF, Control* R)
 @brief パラメータを取得する
 @param CF コンフィギュレーションツリー
 @param R Controlクラスのポインタ
 */
bool IP_Sphere::getXML(SklSolverConfig* CF, Control* R)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  const char *str=NULL;
  REAL_TYPE ct=0.0;
  
  if ( !(elemTop = CF->GetTop(PARAMETER)) ) return false;
  
  if( !(elmL1 = elemTop->GetElemFirst("Intrinsic_Example")) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of 'Intrinsic_Example'\n");
    return false;
  }
  
  // offset
  if ( elmL1->GetValue(CfgIdt("offset"), &ct) ) {
    offset = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Offset' in 'Intrinsic_Example'\n");
    return false;
  }
  if ( offset < 0.0 ) {
    Hostonly_ stamped_printf("\tParsing error : offset must be positive.\n");
    return false;
  }
  
  // radius
  if ( elmL1->GetValue(CfgIdt("radius"), &ct) ) {
    radius = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Radius' in 'Intrinsic_Example'\n");
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
 @fn void IP_Sphere::printPara(FILE* fp, Control* R)
 @brief パラメータの表示
 @param fp ファイルポインタ
 @param R コントロールクラスのポインタ
 */
void IP_Sphere::printPara(FILE* fp, Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Sphere Parameters\n\n");
  
  fprintf(fp,"\tOffset                 [m] / [-]   : %12.5e / %12.5e\n", offset, offset/RefL);
  fprintf(fp,"\tRadius of Sphere       [m] / [-]   : %12.5e / %12.5e\n", radius, radius/RefL);
  if ( drv_length > 0.0 ) {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
}

/**
 @fn bool IP_Sphere::setDomain(Control* R, unsigned m_sz[3], REAL_TYPE m_org[3], REAL_TYPE m_wth[3], REAL_TYPE m_pch[3])
 @brief 領域情報を設定する
 @param R Controlクラスのポインタ
 @param[in]  m_sz グローバル計算領域のセルサイズ
 @param[out] m_org グローバル計算領域の基点
 @param[out] m_wth グローバル計算領域のbounding boxサイズ
 @param[in]  m_pch セルピッチ
 */
void IP_Sphere::setDomain(Control* R, unsigned m_sz[3], REAL_TYPE m_org[3], REAL_TYPE m_wth[3], REAL_TYPE m_pch[3])
{
  pch = m_pch;

  sz  = (int)m_sz;
  
  // チェック
  if ( (pch.x != pch.y) || (pch.y != pch.z) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  
  // 領域サイズ
  wth = pch * (float)sz;
  m_wth[0] = wth.x;
  m_wth[1] = wth.y;
  m_wth[2] = wth.z;
  
  // 領域中心座標
  ctr.x = offset;
  ctr.y = 0.0;
  ctr.z = 0.0;
  
  // 球の中心座標
  os.x = 0.0;
  os.y = 0.0;
  os.z = 0.0;
  
  // 領域の基点
  org = ctr - 0.5*wth;
  m_org[0] = org.x;
  m_org[1] = org.y;
  m_org[2] = org.z;
  
  // 球のbbox
  box_min.x = os.x - radius;

}

//@fn void IP_Sphere::bbox_index(int* st, int* ed)
//@brief コンポーネントの属するセルインデクスを求める
void IP_Sphere::bbox_index(int* st, int* ed)
{
  find_index(st, box_min);
  find_index(ed, box_max);
}

//@fn void IP_Sphere::find_index(int* w, const FB::Vec3f p)
//@brief 点pの属するセルインデクスを求める
//@note Fortran index
void IP_Sphere::find_index(int* w, const FB::Vec3f p)
{
  FB::Vec3f q = (p-org)/pch;
  
  w[0] = (int)ceil(q.x);
  w[1] = (int)ceil(q.y);
  w[2] = (int)ceil(q.z);
}

/**
 @fn void IP_Sphere::setup(int* mid, Control* R, REAL_TYPE* G_org)
 @brief 計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_Sphere::setup(int* mid, Control* R, REAL_TYPE* G_org)
{
  int i,j,k, gd;
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  int mid_driver=3;       /// ドライバ部
  int mid_driver_face=4;  /// ドライバ流出面
  unsigned m;
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
  gd = (int)guide;
  ox_g = G_org[0];
  oy_g = G_org[1];
  oz_g = G_org[2];

  // length, widthなどは有次元値
  len = ox_g + (drv_length+width)/R->RefLength; // グローバルな無次元位置
  ht  = oy_g + height/R->RefLength;
  
  // Initialize  内部領域をfluidにしておく
  for (k=1; k<=(int)kmax; k++) { 
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = mid_fluid;
      }
    }
  }
  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 ) {
    for (k=1; k<=(int)kmax; k++) {
      for (j=1; j<=(int)jmax; j++) {
        for (i=1; i<=(int)imax; i++) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) mid[m] = mid_driver;
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
  for (k=1; k<=(int)kmax; k++) {
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
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

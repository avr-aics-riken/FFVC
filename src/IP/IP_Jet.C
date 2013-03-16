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
/* @brief Jetの流入強化条件
 * @param [in,out] v     速度
 * @param [out]    sum   速度の積算 \sum{v}
 * @param [in,out] flop  flop count
 */
void IP_Jet::vobc_jet_inflow(REAL_TYPE* v)
{
  
  // グローバル
  REAL_TYPE dh = deltaX;
  REAL_TYPE ox_g = G_origin[0];
  REAL_TYPE oy_g = G_origin[1];
  REAL_TYPE oz_g = G_origin[2];
  
  // ノードローカル
  REAL_TYPE ox = origin[0];
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE r, ri, ro, x, y, z;
  REAL_TYPE u1_in = (q1 / a1) / RefV;
  REAL_TYPE u2_in = (q2 / a2) / RefV;
  
  // X-側のJet吹き出し部設定
  if ( nID[X_MINUS] < 0 )
  {
    int i=0;
    
    // Ring1
    ri = r1i;
    ro = r1o;
#pragma omp parallel for firstprivate(i, ix, jx, kx, gd, ri, ro, omg1, dh, u1_in) \
private(x, y, z, r) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        x = ox + ( (REAL_TYPE)i-0.5 ) * dh;
        y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
        z = oz + ( (REAL_TYPE)k-0.5 ) * dh;

        r = sqrt(y*y + z*z);
        
        if ( (ri < r) && (r < ro) )
        {
          v[_F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd)] = u1_in;
          v[_F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd)] = -omg1 * z;
          v[_F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd)] =  omg1 * y;
        }
      }
    }
    
    // Ring2
    ri = r2i;
    ro = r2o;
#pragma omp parallel for firstprivate(i, ix, jx, kx, gd, ri, ro, omg2, dh, u2_in) \
private(x, y, z, r) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        x = ox + ( (REAL_TYPE)i-0.5 ) * dh;
        y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
        z = oz + ( (REAL_TYPE)k-0.5 ) * dh;
        
        r = sqrt(y*y + z*z);
        
        if ( (ri < r) && (r < ro) )
        {
          v[_F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd)] = u2_in;
          v[_F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd)] = -omg2 * z;
          v[_F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd)] =  omg2 * y;
        }
      }
    }
    
    
  } // X_MINUS面の処理
  
  
}
                             
                     
// #################################################################
/* @brief パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_Jet::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  const double pai = (double)(2.0*asin(1.0));

  
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
  
  
  // Ring1
  label="/Parameter/IntrinsicExample/Ring1/UseRing";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  if     ( !strcasecmp(str.c_str(), "yes") ) {
    pat_1 = ON;
  }
  else if( !strcasecmp(str.c_str(), "no") ) {
    pat_1 = OFF;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  if ( pat_1 == ON )
  {
    label="/Parameter/IntrinsicExample/Ring1/InnerRadius";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r1i = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
    }
    
    label="/Parameter/IntrinsicExample/RIng1/OuterRadius";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r1o = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
    }
    
    label="/Parameter/IntrinsicExample/RIng1/RotationFrequency";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      n1 = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL / RefV;
    }
    
    label="/Parameter/IntrinsicExample/Ring1/InletMassFlow";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      q1 = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL * RefL * RefV;
    }
    
    omg1 = 2.0 * pai * n1;
    a1 = pai * (r1o*r1o - r1i*r1i);
  }
  
  
  
  // Ring2
  label="/Parameter/IntrinsicExample/Ring2/UseRing";
  
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  if     ( !strcasecmp(str.c_str(), "yes") ) {
    pat_2 = ON;
  }
  else if( !strcasecmp(str.c_str(), "no") ) {
    pat_2 = OFF;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  if ( pat_2 == ON )
  {
    label="/Parameter/IntrinsicExample/Ring2/InnerRadius";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r2i = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
    }
    
    label="/Parameter/IntrinsicExample/Ring2/OuterRadius";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r2o = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
    }
    
    label="/Parameter/IntrinsicExample/Ring2/RotationFrequency";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      n2 = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL / RefV;
    }
    
    label="/Parameter/IntrinsicExample/Ring2/InletMassFlow";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      q2 = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL * RefL * RefV;
    }
    
    omg2 = 2.0 * pai * n2;
    a2 = pai * (r2o*r2o - r2i*r2i);
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
/* @brief 領域を設定する
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
  
  if ((reg[0] != (REAL_TYPE)sz[0]*pch[0]) ||
      (reg[1] != (REAL_TYPE)sz[1]*pch[1]) ||
      (reg[2] != (REAL_TYPE)sz[2]*pch[2]) ) {
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
/*
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
  
  // REAL_TYPE r1i, r1o;        ///< Ring1の内外径
  // REAL_TYPE r2i, r2o;        ///< Ring2の内外径
  
  /*
   
   0           r1i   r1o  r2i   r2o
   |------------+-----+----+-----+-------> Radius
   |    Solid   |  F  |  S |  F  |  Solid
   
   */
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Backstep Parameters\n\n");
  
  
  // Ring1のチェック
  if ( pat_1 == ON )
  {
    if ( r1i < 0.0 )
    {
      stamped_printf("\tInner Radius of Ring1 must be positive\n");
      Exit(0);
    }
    
    if ( r1o < r1i )
    {
      stamped_printf("\tInner Radius of Ring1 must be smaller than Outer Radius\n");
      Exit(0);
    }
  }

  // Ring2のチェック
  if ( pat_2 == ON)
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
  
  
  // 次元モード
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  
  // Ring1
  if ( pat_1 == ON )
  {
    fprintf(fp,"\tRing1 Inner Radius     [m]   / [-] : %12.5e / %12.5e\n", r1i, r1i/RefL);
    fprintf(fp,"\t      Outer Radius     [m]   / [-] : %12.5e / %12.5e\n", r1o, r1o/RefL);
    fprintf(fp,"\t      Massflow         [m^3/s]/[-] : %12.5e / %12.5e\n", q1, q1/(RefL*RefL*RefV));
    fprintf(fp,"\t      Area             [m^2] / [-] : %12.5e / %12.5e\n", a1, a1/(RefL*RefL));
    fprintf(fp,"\t      Inlet U(x-dir)   [m/s] / [-] : %12.5e / %12.5e\n", q1/a1, (q1/a1)/RefV);
    fprintf(fp,"\t      Angular Velocity [rad/s]/[-] : %12.5e / %12.5e\n", omg1, omg1*RefL/RefV);
    fprintf(fp,"\t      Rotation Freq.   [1/s] / [-] : %12.5e / %12.5e\n", n1, n1);
  }
  
  // Ring2
  if ( pat_2 == ON )
  {
    fprintf(fp,"\tRing2 Inner Radius     [m]   / [-] : %12.5e / %12.5e\n", r2i, r2i/RefL);
    fprintf(fp,"\t      Outer Radius     [m]   / [-] : %12.5e / %12.5e\n", r2o, r2o/RefL);
    fprintf(fp,"\t      Massflow         [m^3/s]/[-] : %12.5e / %12.5e\n", q2, q2/(RefL*RefL*RefV));
    fprintf(fp,"\t      Area             [m^2] / [-] : %12.5e / %12.5e\n", a2, a2/(RefL*RefL));
    fprintf(fp,"\t      Inlet U(x-dir)   [m/s] / [-] : %12.5e / %12.5e\n", q2/a2, (q2/a2)/RefV);
    fprintf(fp,"\t      Angular Velocity [rad/s]/[-] : %12.5e / %12.5e\n", omg2, omg2*RefL/RefV);
    fprintf(fp,"\t      Rotation Freq.   [1/s] / [-] : %12.5e / %12.5e\n", n2, n2);
  }
}



// #################################################################
/* @brief 計算領域のセルIDを設定する
 * @param [in,out] mid   媒質情報の配列
 * @param [in]     R     Controlクラスのポインタ
 * @param [in]     G_org グローバルな原点（無次元）
 * @param [in]     Nmax  Controlクラスのポインタ
 * @param [in]     mat   MediumListクラスのポインタ
 */
void IP_Jet::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 隣接ランクのIDを取得 nID[6]
  const int* nID = paraMngr->GetNeighborRankID();
  
  
  // グローバルな値
  REAL_TYPE dh = deltaX;
  REAL_TYPE ox_g = G_origin[0];
  REAL_TYPE oy_g = G_origin[1];
  REAL_TYPE oz_g = G_origin[2];

  // ノードローカルの無次元値
  FB::Vec3f org_l;
  org_l.x = (float)origin[0];
  org_l.y = (float)origin[1];
  org_l.z = (float)origin[2];
  REAL_TYPE Lx = region[0];
  REAL_TYPE Ly = region[1];
  REAL_TYPE Lz = region[2];
  
  
  // Initialize  全領域をfluidにしておく
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid) \
schedule(static)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = mid_fluid;
      }
    }
  }
  
  
  FB::Vec3f base, b, t;
  float ph = (float)dh;
  float r, ri, ro;
  
  // X-側のJet吹き出し部設定
  if ( nID[X_MINUS] < 0 )
  {
    int i=0;
    
    // デフォルトでガイドセルをSolidにする
#pragma omp parallel for firstprivate(i, ix, jx, kx, gd, mid_solid) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = mid_solid;
      }
    }
    
    // Ring1
    ri = r1i;
    ro = r1o;
#pragma omp parallel for firstprivate(i, ix, jx, kx, gd, mid_fluid, ri, ro, org_l, ph) \
private(b, r, base) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        base.assign((float)i-0.5, (float)j-0.5, (float)k-0.5);
        b = org_l + base*ph;
        r = b.length();
        
        if ( (ri < r) && (r < ro) )
        {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          mid[m] = mid_fluid;
        }
        
      }
    }
   
    // Ring2
    ri = r2i;
    ro = r2o;
#pragma omp parallel for firstprivate(i, ix, jx, kx, gd, mid_fluid, ri, ro, org_l, ph) \
private(b, r, base) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        base.assign((float)i-0.5, (float)j-0.5, (float)k-0.5);
        b = org_l + base*ph;
        r = b.length();
        
        if ( (ri < r) && (r < ro) )
        {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          mid[m] = mid_fluid;
        }
        
      }
    }
    
    
    
  } // X_MINUS面の処理
  
}

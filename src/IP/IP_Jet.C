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
/* @brief Jetの流入境界条件による発散値の修正
 * @param [in,out] div   発散値
 * @param [in] bv        BCindex V
 * @param [in,out] vf    セルフェイス速度    
 * @param [in,out] flop  flop count
 * @retval         sum   流入量（無次元）
 */
REAL_TYPE IP_Jet::divJetInflow(REAL_TYPE* div, const int* bv, REAL_TYPE* vf, double& flop)
{
  
  // グローバル
  REAL_TYPE dh = deltaX;
  
  // ノードローカル
  //REAL_TYPE ox = origin[0];
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE sum = 0.0;
  
  // Ring1
  if ( pat_1 == ON )
  {
    REAL_TYPE ri = r1i;
    REAL_TYPE ro = r1o;
    REAL_TYPE u_in = q1 / a1;
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ri, ro, dh, u_in) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dh;
        
        REAL_TYPE r = sqrt(y*y + z*z);
        
        if ( (ri < r) && (r < ro) )
        {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          div[m] -= u_in;
          vf[_F_IDX_V3D(0, j, k, 0, ix, jx, kx, gd)] = u_in;
        }
      }
    }
    
    sum += q1;
    flop += (double)jx * (double)kx * 21.0; // DP 31.0
  }
  
  
  // Ring2
  if ( pat_2 == ON )
  {
    REAL_TYPE ri = r2i;
    REAL_TYPE ro = r2o;
    REAL_TYPE u_in = q2 / a2;
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ri, ro, dh, u_in) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dh;
        
        REAL_TYPE r = sqrt(y*y + z*z);
        
        if ( (ri < r) && (r < ro) )
        {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          div[m] -= u_in;
          vf[_F_IDX_V3D(0, j, k, 0, ix, jx, kx, gd)] = u_in;
        }
      }
    }
    
    sum += q2;
    flop += (double)jx * (double)kx * 21.0; // DP 31.0
  }
  
  return sum;
}


// #################################################################
/* @brief Jetの流入境界条件　Xマイナス方向のみ
 * @param [in,out] wv    疑似速度
 * @param [in]     rei   レイノルズ数の逆数
 * @param [in]     v0    速度ベクトル（n-step）
 * @param [in]     bv    BCindex V
 * @param [in,out] flop  flop count
 */
void IP_Jet::vobc_pv_JetInflow(REAL_TYPE* wv,
                               const REAL_TYPE rei,
                               const REAL_TYPE* v0,
                               const int* bv,
                               double& flop)
{
  
  // グローバル
  REAL_TYPE dh = deltaX;
  REAL_TYPE dh2= rei / (dh * dh);
  
  // ノードローカル
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  
  // 壁面速度
  REAL_TYPE uy = 0.0;
  REAL_TYPE uz = 0.0;

  
  int i = 1;
  REAL_TYPE rin_1  = r1i;
  REAL_TYPE rout_1 = r1o;
  REAL_TYPE rin_2  = r2i;
  REAL_TYPE rout_2 = r2o;
  REAL_TYPE uin_1  = q1 / a1;
  REAL_TYPE uin_2  = q2 / a2;
  REAL_TYPE omg_1  = omg1;
  REAL_TYPE omg_2  = omg2;
  
#pragma omp parallel for firstprivate(i, ix, jx, kx, gd, dh, dh2, uy, uz, \
rin_1, rout_1, omg_1, uin_1, rin_2, rout_2, omg_2, uin_2) \
schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      
      REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
      REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dh;
      
      REAL_TYPE r = sqrt(y*y + z*z);
      
      size_t m0 = _F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd);
      size_t m1 = _F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd);
      size_t m2 = _F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd);
      
      if ( (rin_1 < r) && (r < rout_1) ) // Ring1
      {
        REAL_TYPE ur = uin_1;
        REAL_TYPE vr = -omg_1 * z;
        REAL_TYPE wr =  omg_1 * y;
        
        REAL_TYPE c  = uin_1;
        REAL_TYPE ac = fabs(c);
        
        REAL_TYPE up = v0[m0];
        REAL_TYPE vp = v0[m1];
        REAL_TYPE wp = v0[m2];
        
        REAL_TYPE ex = up - ur;
        REAL_TYPE ey = vp - vr;
        REAL_TYPE ez = wp - wr;
        
        REAL_TYPE fu = 0.5*(c*(up+ur) - ac*ex);
        REAL_TYPE fv = 0.5*(c*(vp+vr) - ac*ey);
        REAL_TYPE fw = 0.5*(c*(wp+wr) - ac*ez);

        wv[m0] += (fu*dh - ex*dh2);
        wv[m1] += (fv*dh - ey*dh2);
        wv[m2] += (fw*dh - ez*dh2);
      }
      else if ( (rin_2 < r) && (r < rout_2) ) // Ring2
      {
        REAL_TYPE ur = uin_2;
        REAL_TYPE vr = -omg_2 * z;
        REAL_TYPE wr =  omg_2 * y;
        
        REAL_TYPE c  = uin_2;
        REAL_TYPE ac = fabs(c);
        
        REAL_TYPE up = v0[m0];
        REAL_TYPE vp = v0[m1];
        REAL_TYPE wp = v0[m2];
        
        REAL_TYPE ex = up - ur;
        REAL_TYPE ey = vp - vr;
        REAL_TYPE ez = wp - wr;
        
        REAL_TYPE fu = 0.5*(c*(up+ur) - ac*ex);
        REAL_TYPE fv = 0.5*(c*(vp+vr) - ac*ey);
        REAL_TYPE fw = 0.5*(c*(wp+wr) - ac*ez);
        
        wv[m0] += (fu*dh - ex*dh2);
        wv[m1] += (fv*dh - ey*dh2);
        wv[m2] += (fw*dh - ez*dh2);
      }
      else // Wall
      {
        wv[m1] += (uy - v0[m1]) * dh2 * 2.0;
        wv[m2] += (uz - v0[m2]) * dh2 * 2.0;
      }
    }
  }
  flop += (double)jx * (double)kx * 27.0; // DP 37.0
}



// #################################################################
/* @brief Jetの流入境界条件をガイドセルに代入
 * @param [in,out] v     セルセンター速度
 */
void IP_Jet::vobcJetInflowGC(REAL_TYPE* v)
{
  
  // グローバル
  REAL_TYPE dh = deltaX;
  
  // ノードローカル
  //REAL_TYPE ox = origin[0];
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;


  // Ring1
  if ( pat_1 == ON )
  {
    REAL_TYPE ri = r1i;
    REAL_TYPE ro = r1o;
    REAL_TYPE u_in = q1 / a1;
    REAL_TYPE omg = omg1;
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ri, ro, omg, dh, u_in) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dh;
        
        REAL_TYPE r = sqrt(y*y + z*z);
        
        if ( (ri < r) && (r < ro) )
        {
          v[_F_IDX_V3D(0, j, k, 0, ix, jx, kx, gd)] = u_in;
          v[_F_IDX_V3D(0, j, k, 1, ix, jx, kx, gd)] = -omg * z;
          v[_F_IDX_V3D(0, j, k, 2, ix, jx, kx, gd)] =  omg * y;
        }
      }
    }
  }

  
  // Ring2
  if ( pat_2 == ON )
  {
    REAL_TYPE ri = r2i;
    REAL_TYPE ro = r2o;
    REAL_TYPE u_in = q2 / a2;
    REAL_TYPE omg = omg2;
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ri, ro, omg, dh, u_in) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dh;
        
        REAL_TYPE r = sqrt(y*y + z*z);
        
        if ( (ri < r) && (r < ro) )
        {
          v[_F_IDX_V3D(0, j, k, 0, ix, jx, kx, gd)] = u_in;
          v[_F_IDX_V3D(0, j, k, 1, ix, jx, kx, gd)] = -omg * z;
          v[_F_IDX_V3D(0, j, k, 2, ix, jx, kx, gd)] =  omg * y;
        }
      }
    }
  }

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
  
  // パラメータは無次元で保持
  
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
      r1i = ( R->Unit.Param == DIMENSIONAL ) ? ct/RefL : ct;
    }
    
    label="/Parameter/IntrinsicExample/RIng1/OuterRadius";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r1o = ( R->Unit.Param == DIMENSIONAL ) ? ct/RefL : ct;
    }
    
    label="/Parameter/IntrinsicExample/RIng1/RotationFrequency";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      n1 = ( R->Unit.Param == DIMENSIONAL ) ? ct * RefL / RefV : ct;
    }
    
    label="/Parameter/IntrinsicExample/Ring1/InletMassFlow";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      q1 = ( R->Unit.Param == DIMENSIONAL ) ? ct / (RefL * RefL * RefV) : ct;
    }
    
    omg1 = 2.0 * pai * n1;
    a1 = pai * (r1o*r1o - r1i*r1i) / (RefL*RefL);
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
      r2i = ( R->Unit.Param == DIMENSIONAL ) ? ct/RefL : ct;
    }
    
    label="/Parameter/IntrinsicExample/Ring2/OuterRadius";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r2o = ( R->Unit.Param == DIMENSIONAL ) ? ct/RefL : ct;
    }
    
    label="/Parameter/IntrinsicExample/Ring2/RotationFrequency";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      n2 = ( R->Unit.Param == DIMENSIONAL ) ? ct * RefL / RefV : ct;
    }
    
    label="/Parameter/IntrinsicExample/Ring2/InletMassFlow";
    if ( !(tpCntl->GetValue(label, &ct )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      q2 = ( R->Unit.Param == DIMENSIONAL ) ? ct / (RefL * RefL * RefV) : ct;
    }
    
    omg2 = 2.0 * pai * n2;
    a2 = pai * (r2o*r2o - r2i*r2i) / (RefL*RefL);
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
  fprintf(fp,"\n\t>> Intrinsic Jet Parameters\n\n");
  
  
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
  fprintf(fp,"\tDimension Mode                     :  %s\n\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  
  // Ring1
  if ( pat_1 == ON )
  {
    fprintf(fp,"\tRing1 Inner Radius     [m]   / [-] : %12.5e / %12.5e\n", r1i*RefL, r1i);
    fprintf(fp,"\t      Outer Radius     [m]   / [-] : %12.5e / %12.5e\n", r1o*RefL, r1o);
    fprintf(fp,"\t      Massflow         [m^3/s]/[-] : %12.5e / %12.5e\n", q1*(RefL*RefL*RefV), q1);
    fprintf(fp,"\t      Area             [m^2] / [-] : %12.5e / %12.5e\n", a1*(RefL*RefL), a1);
    fprintf(fp,"\t      Inlet U(x-dir)   [m/s] / [-] : %12.5e / %12.5e\n", q1/a1*RefV, q1/a1);
    fprintf(fp,"\t      Angular Velocity [rad/s]/[-] : %12.5e / %12.5e\n", omg1*RefV/RefL, omg1);
    fprintf(fp,"\t      Rotation Freq.   [1/s] / [-] : %12.5e / %12.5e\n", n1*RefV/RefL, n1);
    fprintf(fp,"\n");
  }
  
  // Ring2
  if ( pat_2 == ON )
  {
    fprintf(fp,"\tRing2 Inner Radius     [m]   / [-] : %12.5e / %12.5e\n", r2i*RefL, r2i);
    fprintf(fp,"\t      Outer Radius     [m]   / [-] : %12.5e / %12.5e\n", r2o*RefL, r2o);
    fprintf(fp,"\t      Massflow         [m^3/s]/[-] : %12.5e / %12.5e\n", q2*(RefL*RefL*RefV), q2);
    fprintf(fp,"\t      Area             [m^2] / [-] : %12.5e / %12.5e\n", a2*(RefL*RefL), a2);
    fprintf(fp,"\t      Inlet U(x-dir)   [m/s] / [-] : %12.5e / %12.5e\n", q2/a2*RefV, q2/a2);
    fprintf(fp,"\t      Angular Velocity [rad/s]/[-] : %12.5e / %12.5e\n", omg2*RefV/RefL, omg2);
    fprintf(fp,"\t      Rotation Freq.   [1/s] / [-] : %12.5e / %12.5e\n", n2*RefV/RefL, n2);
  }
}



// #################################################################
/* @brief 計算領域のセルIDを設定する
 * @param [in,out] mid   媒質情報の配列
 * @param [in]     R     Controlクラスのポインタ
 * @param [in]     G_org グローバルな原点（無次元）
 * @param [in]     Nmax  媒質数
 * @param [in]     mat   MediumListクラスのポインタ
 */
void IP_Jet::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int mid_fluid;        /// 流体
  int mid_solid;        /// 固体
  
  // 媒質設定のチェック
  if ( (mid_fluid = R->find_ID_from_Label(mat, Nmax, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  if ( (mid_solid = R->find_ID_from_Label(mat, Nmax, m_solid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_solid.c_str());
    Exit(0);
  }
  
  
  REAL_TYPE dh = deltaX;
  
  // ノードローカル
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  
  // 隣接ランクのIDを取得 nID[6]
  const int* nID = paraMngr->GetNeighborRankID();
  
  
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
  

  
  // X-側の場合に，Jet吹き出し部設定
  if ( nID[X_MINUS] < 0 )
  {
    // デフォルトでガイドセルをSolidにする
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_solid) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        size_t m = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
        mid[m] = mid_solid;
      }
    }
    
    
    // Ring1
    if ( pat_1 == ON )
    {
      REAL_TYPE ri = r1i;
      REAL_TYPE ro = r1o;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid, ri, ro, oy, oz, dh) \
schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          
          REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
          REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dh;
          
          REAL_TYPE r = sqrt(y*y + z*z);
          
          if ( (ri < r) && (r < ro) )
          {
            mid[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = mid_fluid;
          }
          
        }
      }
    }

   
    // Ring2
    if ( pat_2 == ON )
    {
      REAL_TYPE ri = r2i;
      REAL_TYPE ro = r2o;
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid, ri, ro, oy, oz, dh) \
schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          
          REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dh;
          REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dh;
          
          REAL_TYPE r = sqrt(y*y + z*z);
          
          if ( (ri < r) && (r < ro) )
          {
            mid[_F_IDX_S3D(0, j, k, ix, jx, kx, gd)] = mid_fluid;
          }
          
        }
      }
    }
    
  } // X_MINUS面の処理
  
}

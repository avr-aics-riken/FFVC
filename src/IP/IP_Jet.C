//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
 * @file   IP_Jet.C
 * @brief  IP_Jet class
 * @author aics
 */

#include "IP_Jet.h"


// #################################################################
// Jetの流入境界条件による発散値の修正
void IP_Jet::divJetInflow(REAL_TYPE* div, const int face, double& flop)
{
  // X_MINUS面の外部境界面のみ
  if ( nID[face] >= 0) return;

  REAL_TYPE dx = pitch[0];
  REAL_TYPE dy = pitch[1];
  REAL_TYPE dz = pitch[2];
  
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
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ri, ro, u_in, dx, dy, dz) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dy;
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dz;
        
        REAL_TYPE r = sqrt(y*y + z*z);
        
        if ( (ri < r) && (r < ro) )
        {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          div[m] -= (u_in / dx);
          //vf[_F_IDX_V3D(0, j, k, 0, ix, jx, kx, gd)] = u_in;
        }
      }
    }
    
    flop += (double)jx * (double)kx * 37.0;
  }
  
  
  // Ring2
  if ( pat_2 == ON )
  {
    REAL_TYPE ri = r2i;
    REAL_TYPE ro = r2o;
    REAL_TYPE u_in = q2 / a2;
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ri, ro, u_in, dx, dy, dz) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dy;
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dz;
        
        REAL_TYPE r = sqrt(y*y + z*z);
        
        if ( (ri < r) && (r < ro) )
        {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          div[m] -= (u_in / dx);
          //vf[_F_IDX_V3D(0, j, k, 0, ix, jx, kx, gd)] = u_in;
        }
      }
    }
    
    flop += (double)jx * (double)kx * 37.0;
  }
  
}


// #################################################################
// 流束型流入境界条件
void IP_Jet::vobc_pv_JetInflow(REAL_TYPE* wv,
                               const int face, 
                               const REAL_TYPE rei,
                               const REAL_TYPE* v0,
                               double& flop)
{
  // X_MINUS面の外部境界面のみ
  if ( nID[face] >= 0) return;
  

  REAL_TYPE dx = pitch[0];
  REAL_TYPE dy = pitch[1];
  REAL_TYPE dz = pitch[2];
  
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
  
#pragma omp parallel for firstprivate(i, ix, jx, kx, gd, uy, uz, \
rin_1, rout_1, omg_1, uin_1, rin_2, rout_2, omg_2, uin_2, dx, dy, dz) \
schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      
      REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dy;
      REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dz;
      
      REAL_TYPE r = sqrt(y*y + z*z);
      
      size_t m0 = _F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd);
      size_t m1 = _F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd);
      size_t m2 = _F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd);
      
      flop += 28.0;
      
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

        wv[m0] += (fu*dx - ex * rei / (dx * dx));
        wv[m1] += (fv*dy - ey * rei / (dy * dy));
        wv[m2] += (fw*dz - ez * rei / (dz * dz));
        
        flop += 61.0;
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
        
        wv[m0] += (fu*dx - ex * rei / (dx * dx));
        wv[m1] += (fv*dy - ey * rei / (dy * dy));
        wv[m2] += (fw*dz - ez * rei / (dz * dz));
        
        flop += 61.0;
      }
      else // Wall
      {
        wv[m1] += (uy - v0[m1]) * 2.0 * rei / (dy * dy);
        wv[m2] += (uz - v0[m2]) * 2.0 * rei / (dz * dz);
        
        flop += 26.0;
      }
    }
  }

}



// #################################################################
// Jetの流入境界条件をガイドセルに代入
void IP_Jet::vobcJetInflowGC(REAL_TYPE* v, const int face)
{
  // X_MINUS面の外部境界面のみ
  if ( nID[face] >= 0) return;
  
  //REAL_TYPE dx = pitch[0];
  REAL_TYPE dy = pitch[1];
  REAL_TYPE dz = pitch[2];
  
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
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ri, ro, omg, u_in, dy, dz) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dy;
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dz;
        
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
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, ri, ro, omg, u_in, dy, dz) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        
        REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dy;
        REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dz;
        
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
//  パラメータをロード
bool IP_Jet::getTP(Control* R, TextParser* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  const double pai = (double)(2.0*asin(1.0));

  
  // 2D or 3D mode
  label = "/IntrinsicExample/Dimension";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if     ( !strcasecmp(str.c_str(), "2d") )
  {
    mode = dim_2d;
  }
  else if( !strcasecmp(str.c_str(), "3d") )
  {
    mode = dim_3d;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  // パラメータは無次元で保持
  
  // Ring1
  label = "/IntrinsicExample/Ring1/UseRing";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  if     ( !strcasecmp(str.c_str(), "yes") )
  {
    pat_1 = ON;
  }
  else if( !strcasecmp(str.c_str(), "no") )
  {
    pat_1 = OFF;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  if ( pat_1 == ON )
  {
    label = "/IntrinsicExample/Ring1/InnerRadius";
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r1i = ( R->Unit.Param == DIMENSIONAL ) ? ct/RefL : ct;
    }
    
    label = "/IntrinsicExample/RIng1/OuterRadius";
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r1o = ( R->Unit.Param == DIMENSIONAL ) ? ct/RefL : ct;
    }
    
    label = "/IntrinsicExample/Ring1/RotationFrequency";
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      n1 = ( R->Unit.Param == DIMENSIONAL ) ? ct * RefL / RefV : ct;
    }
    
    label = "/IntrinsicExample/Ring1/InletMassFlow";
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      q1 = ( R->Unit.Param == DIMENSIONAL ) ? ct / (RefL * RefL * RefV) : ct;
    }
    
    omg1 = 2.0 * pai * n1;
    a1 = pai * (r1o*r1o - r1i*r1i);
  }
  
  
  
  // Ring2
  label = "/IntrinsicExample/Ring2/UseRing";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  if     ( !strcasecmp(str.c_str(), "yes") )
  {
    pat_2 = ON;
  }
  else if( !strcasecmp(str.c_str(), "no") )
  {
    pat_2 = OFF;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  if ( pat_2 == ON )
  {
    label = "/IntrinsicExample/Ring2/InnerRadius";
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r2i = ( R->Unit.Param == DIMENSIONAL ) ? ct/RefL : ct;
    }
    
    label = "/IntrinsicExample/Ring2/OuterRadius";
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      r2o = ( R->Unit.Param == DIMENSIONAL ) ? ct/RefL : ct;
    }
    
    label = "/IntrinsicExample/Ring2/RotationFrequency";
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      n2 = ( R->Unit.Param == DIMENSIONAL ) ? ct * RefL / RefV : ct;
    }
    
    label = "/IntrinsicExample/Ring2/InletMassFlow";
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else{
      q2 = ( R->Unit.Param == DIMENSIONAL ) ? ct / (RefL * RefL * RefV) : ct;
    }
    
    omg2 = 2.0 * pai * n2;
    a2 = pai * (r2o*r2o - r2i*r2i);
  }
  
  
  
  // 媒質指定
  label = "/IntrinsicExample/FluidMedium";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label = "/IntrinsicExample/SolidMedium";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  return true;
}



// #################################################################
// パラメータの表示
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
  
  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Jet Class Parameters\n\n");
  
  
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
// 外部境界の設定
void IP_Jet::setOBC(const int face, int* bcd, Control* R, REAL_TYPE* G_org, const int NoMedium, const MediumList* mat, int* cutL, int* cutU, int*bid)
{
  int mid_fluid;        /// 流体
  int mid_solid;        /// 固体  
  
  // 媒質設定のチェック
  if ( (mid_fluid = FBUtility::findIDfromLabel(mat, NoMedium, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  if ( (mid_solid = FBUtility::findIDfromLabel(mat, NoMedium, m_solid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_solid.c_str());
    Exit(0);
  }
  
  
  //REAL_TYPE dx = pitch[0];
  REAL_TYPE dy = pitch[1];
  REAL_TYPE dz = pitch[2];
  
  // ノードローカル
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
  // 隣接ランクのIDを取得 nID[6]
  const int* nID = paraMngr->GetNeighborRankID(procGrp);
  
  
  
  // X-側の場合に，Jet吹き出し部設定
  if ( face == X_minus )
  {
    if ( nID[X_minus] < 0 )
    {
      // デフォルトでガイドセルをSolidにする
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_solid) \
schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          
          // 媒質エントリ
          size_t m = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
          setMediumID(bcd[m], mid_solid);
          
          // 交点
          size_t l = _F_IDX_S3D(1  , j  , k  , ix, jx, kx, gd);
          int r = quantize9(0.5); /// 壁面までの距離
          setCutL9(cutL[l], r, X_minus);
          
          // 境界ID
          setBit5(bid[l], mid_solid, X_minus);
        }
      }
      
      
      // Ring1
      if ( pat_1 == ON )
      {
        REAL_TYPE ri = r1i;
        REAL_TYPE ro = r1o;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid, ri, ro, oy, oz, dy, dz) \
schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            
            REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dy;
            REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dz;
            
            REAL_TYPE r = sqrt(y*y + z*z);
            
            if ( (ri < r) && (r < ro) )
            {
              size_t m = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
              setMediumID(bcd[m], mid_fluid);
              
              size_t l = _F_IDX_S3D(1  , j  , k  , ix, jx, kx, gd);
              int r = quantize9(1.0);
              setCutL9(cutL[l], r, X_minus);
              setBit5(bid[l], 0, X_minus);
            }
            
          }
        }
      } // ring1
      
      
      // Ring2
      if ( pat_2 == ON )
      {
        REAL_TYPE ri = r2i;
        REAL_TYPE ro = r2o;
        
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid, ri, ro, oy, oz, dy, dz) \
schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            
            REAL_TYPE y = oy + ( (REAL_TYPE)j-0.5 ) * dy;
            REAL_TYPE z = oz + ( (REAL_TYPE)k-0.5 ) * dz;
            
            REAL_TYPE r = sqrt(y*y + z*z);
            
            if ( (ri < r) && (r < ro) )
            {
              size_t m = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
              setMediumID(bcd[m], mid_fluid);
              
              size_t l = _F_IDX_S3D(1  , j  , k  , ix, jx, kx, gd);
              int r = quantize9(1.0);
              setCutL9(cutL[l], r, X_minus);
              setBit5(bid[l], 0, X_minus);
            }
            
          }
        }
      }// ring2
    }
    
  } // face
  
}

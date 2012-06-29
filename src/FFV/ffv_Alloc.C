// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################


/** 
 * @file   ffv_Alloc.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"

// Adams-Bashforth法に用いる配列のアロケーション
void FFV::allocArray_AB2 (double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_abf
  if ( !(dc_abf = paraMngr->AllocRealV3DEx(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
}



// 平均値処理に用いる配列のアロケーション
void FFV::allocArray_Average (double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_ap
  if ( !(dc_ap = paraMngr->AllocFloatS3D(guide)) ) Exit(0);
  total += mc * (double)sizeof(float);
  
  // dc_av
  if ( !(dc_av = paraMngr->AllocRealV3DEx(guide)) ) Exit(0);
  total += mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  if ( C.isHeatProblem() ) {
    // dc_at
    if ( !(dc_at = paraMngr->AllocFloatS3D(guide)) ) Exit(0);
    total += mc * (double)sizeof(float);
  }
}



// 粗格子読み込みに用いる配列のアロケーション
void FFV::allocArray_CoarseMesh(const int* r_size, double &prep)
{
  double mc = (double)(r_size[0] * r_size[1] * r_size[2]);
  
  // dc_r_p
  if ( !(dc_r_p = Alloc::Real_S3D(r_size, guide)) ) Exit(0);
  prep += mc * (double)sizeof(REAL_TYPE);
  
  // dc_r_v
  if ( !(dc_r_v = Alloc::Real_V3D(r_size, guide)) ) Exit(0);
  prep += mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  
  if ( C.isHeatProblem() ) {
    
    // dc_r_t
    if ( !(dc_r_t = Alloc::Real_S3D(r_size, guide)) ) Exit(0);
    prep += mc * (double)sizeof(REAL_TYPE);
    
  }
}



// コンポーネント体積率の配列のアロケーション
void FFV::allocArray_CompoVF(double &prep, double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_cvf
  if ( !(dc_cvf = paraMngr->AllocFloatS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(float);
  total+= mc * (double)sizeof(float);
}



// カット情報の配列
void FFV::allocArray_Cut(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_cut
  if ( !(dc_cut = paraMngr->AllocFloatS4DEx(6, guide)) ) Exit(0);
  total+= mc * (double)sizeof(float) * 6.0;
  
  // dc_bid
  if ( !(dc_bid = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  total += mc * (double)sizeof(int);
}



// 熱の主計算部分に用いる配列のアロケーション
void FFV::allocArray_Heat(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_t
  if ( !(dc_t = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // dc_t0
  if ( !(dc_t0 = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // dc_qbc
  if ( !(dc_qbc = paraMngr->AllocRealV3DEx(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
}



// 体積率の配列のアロケーション
void FFV::allocArray_Interface(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_vof
  if ( !(dc_vof = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
}



// LES計算に用いる配列のアロケーション
void FFV::allocArray_LES(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_vt
  if ( !(dc_vt = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
}



// 主計算部分に用いる配列のアロケーション
void FFV::allocArray_Main(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_v
  if ( !(dc_v = paraMngr->AllocRealV3DEx(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // dc_vc
  if ( !(dc_vc = paraMngr->AllocRealV3DEx(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // dc_v0
  if ( !(dc_v0 = paraMngr->AllocRealV3DEx(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // dc_wv
  if ( !(dc_wv = paraMngr->AllocRealV3DEx(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // dc_wvex
  if ( !(dc_wvex = paraMngr->AllocRealV3DEx(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // dc_p
  if ( !(dc_p = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // dc_p0
  if ( !(dc_p0 = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // dc_wk2
  if ( !(dc_wk2 = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
}



// 前処理に用いる配列のアロケーション
void FFV::allocArray_Prep(double &prep, double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_ws
  if ( !(dc_ws = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(REAL_TYPE);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // dc_mid
  if ( !(dc_mid = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  
  // dc_bcd
  if ( !(dc_bcd = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
  
  // dc_bcp
  if ( !(dc_bcp = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
  
  // dc_bcv
  if ( !(dc_bcv = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
  
  if ( C.isHeatProblem() ) {
    
    // dc_bh1
    if ( !(dc_bh1 = paraMngr->AllocIntS3D(guide)) ) Exit(0);
    prep += mc * (double)sizeof(int);
    total+= mc * (double)sizeof(int);
    
    // dc_bh2
    if ( !(dc_bh2 = paraMngr->AllocIntS3D(guide)) ) Exit(0);
    prep += mc * (double)sizeof(int);
    total+= mc * (double)sizeof(int);
  }
}



// Runge-Kutta法に用いる配列のアロケーション
void FFV::allocArray_RK(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // dc_dp
  if ( !(dc_dp = paraMngr->AllocRealS3D(guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
}

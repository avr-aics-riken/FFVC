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
 * @file   ffv_Heat.C
 * @brief  FFV BC Class
 * @author kero
 */

#include <math.h>
#include <time.h>
#include <fcntl.h>

#ifndef _WIN32
#include <unistd.h>
#include <strings.h>
#else
#include "sph_win32_util.h"
#endif
#include <sys/types.h>

#include "ffv.h"



// #################################################################
// 移流項のEuler陽解法による時間積分
void FFV::ps_ConvectionEE(REAL_TYPE* tc, const REAL_TYPE delta_t, const int* bd, const REAL_TYPE* t0, double& flop)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  REAL_TYPE dt = delta_t;
  
  flop += (double)(ix*jx*kx)* 3.0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dt) schedule(static)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        tc[m] = t0[m] + dt * tc[m] * GET_SHIFT_F(bd[m], STATE_BIT); // 対流項の評価は流動なので，状態を参照
      }
    }
  }
}



// #################################################################
// Boussinesq浮力項の計算
void FFV::Buoyancy(REAL_TYPE* v, const REAL_TYPE dgr, const REAL_TYPE* t, const int* bd, double& flop)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE d = dgr;
  
  flop += (double)(ix*jx*kx)* 3.0;

#pragma omp parallel for firstprivate(ix, jx, kx, gd, d) schedule(static)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        size_t l = _F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd); // k方向が重力方向
        v[l] += d * t[m]* GET_SHIFT_F(bd[m], STATE_BIT); // 対流項の計算なので，状態を参照
      }
    }
  }
}



// #################################################################
// 単媒質に対する熱伝導方程式を陰解法で解く
void FFV::ps_LS(IterationCtl* IC, const double rhs_nrm, const double r0)
{
  double flop = 0.0;      /// 浮動小数点演算数
  double res=0.0;         /// 残差
  double b2=0.0;          /// 反復式のソースベクトルのノルム
  double nrm = 0.0;       ///
  REAL_TYPE dt = deltaT;  /// 時間積分幅
  
  // d_t   温度 n+1 step
  // d_t0  温度 n step
  // d_qbc 境界条件の熱流束
  // d_ws  対流項のみの部分段階
  // d_bh2 BCindex H2
  
  unsigned int wait_num=0;
  int req[12];
  
	
  switch (IC->getLS()) 
  {
    case SOR:
      
      // >>> Passive scalar Diffusion subsection 3
      TIMING_start(tm_heat_diff_sct_3);
      
      // 反復処理
      TIMING_start(tm_heat_diff_PSOR);
      
      flop = 0.0;
      res = ps_Diff_SM_PSOR(d_t, b2, dt, d_qbc, d_bh2, d_ws, IC, flop);
      
      TIMING_stop(tm_heat_diff_PSOR, flop);
      
      // 外部周期境界条件
      TIMING_start(tm_heat_diff_OBC);
      BC.OuterTBCperiodic(d_t);
      TIMING_stop(tm_heat_diff_OBC, 0.0);
      
      // 温度の同期
      if ( numProc > 1 )
      {
        TIMING_start(tm_heat_update_comm);
        if ( paraMngr->BndCommS3D(d_t, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
        TIMING_stop(tm_heat_update_comm, face_comm_size*guide);
      }
      
      // 残差の集約
      if ( numProc > 1 )
      {
        TIMING_start(tm_heat_diff_res_comm);
        REAL_TYPE tmp_wk[2], m_tmp[2];
        tmp_wk[0] = m_tmp[0] = res;
        tmp_wk[1] = m_tmp[1] = b2;
        if ( paraMngr->Allreduce(tmp_wk, m_tmp, 2, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop(tm_heat_diff_res_comm, 2.0*numProc*2.0*sizeof(REAL_TYPE) );
        
        res = sqrt( m_tmp[0]/(REAL_TYPE)G_Acell ); // 残差のRMS
        b2  = sqrt( m_tmp[1]/(REAL_TYPE)G_Acell ); // ソースベクトルのRMS
      }
      
      TIMING_stop(tm_heat_diff_sct_3, 0.0);
      // <<< Passive scalar Diffusion subsection 3
      break;

    default:
      printf("\tInvalid Linear Solver for Heat\n");
      Exit(0);
      break;
  }
  
  // Residual resを上書き
  switch ( IC->getNormType() )
  {
    case r_b:
    case r_r0:
      
      TIMING_start(tm_poi_src_nrm);
      res = 0.0;
      flop = 0.0;
      //poi_residual_(&res, size, &guide, d_p, d_ws, d_bcp, &flop);
      TIMING_stop(tm_poi_src_nrm, flop);
      
      if ( numProc > 1 )
      {
        TIMING_start(tm_poi_src_comm);
        double m_tmp = res;
        if ( paraMngr->Allreduce(&m_tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop(tm_poi_src_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
      }
      
      res = sqrt(res);
      break;
  }
  
  // 残差の保存
  switch ( IC->getNormType() )
  {
      
    case dx_b:
      IC->setNormValue( res/rhs_nrm );
      break;
      
    case r_b:
      IC->setNormValue( res/rhs_nrm );
      break;
      
    case r_r0:
      IC->setNormValue( res/r0 );
      break;
  }
  
}


// #################################################################
// 単媒質に対する熱伝導方程式をEuler陽解法で解く
REAL_TYPE FFV::ps_Diff_SM_EE(REAL_TYPE* t, const REAL_TYPE dt, const REAL_TYPE* qbc, const int* bh2, const REAL_TYPE* ws, double& flop)
{
  REAL_TYPE g_p, g_w, g_e, g_s, g_n, g_b, g_t;
  REAL_TYPE t_p, t_w, t_e, t_s, t_n, t_b, t_t;
  REAL_TYPE      a_w, a_e, a_s, a_n, a_b, a_t;
  REAL_TYPE dth1, dth2, delta, res;
  REAL_TYPE dh = (REAL_TYPE)deltaX;    /// 空間格子幅
  int s;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  dth1 = dt/dh;
  dth2 = dth1*C.getRcpPeclet()/dh;
  res  = 0.0;
  flop += (double)(ix*jx*kx)* 50.0;

#pragma omp parallel for firstprivate(ix, jx, kx, gd, dth1, dth2) \
            private(t_p, t_w, t_e, t_s, t_n, t_b, t_t) \
            private(g_p, g_w, g_e, g_s, g_n, g_b, g_t) \
            private(a_w, a_e, a_s, a_n, a_b, a_t) \
            private(s, delta) \
            schedule(static) reduction(+:res)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
#include "FindexS3D.h"
        
        t_p = t[m_p];
        t_w = t[m_w];
        t_e = t[m_e];
        t_s = t[m_s];
        t_n = t[m_n];
        t_b = t[m_b];
        t_t = t[m_t];

        s   = bh2[m_p];
        
        g_w = GET_SHIFT_F(s, GMA_W); // gamma(BC)
        g_e = GET_SHIFT_F(s, GMA_E);
        g_s = GET_SHIFT_F(s, GMA_S);
        g_n = GET_SHIFT_F(s, GMA_N);
        g_b = GET_SHIFT_F(s, GMA_B);
        g_t = GET_SHIFT_F(s, GMA_T);
        
        a_w = GET_SHIFT_F(s, ADIABATIC_W); // gamma(A)
        a_e = GET_SHIFT_F(s, ADIABATIC_E);
        a_s = GET_SHIFT_F(s, ADIABATIC_S);
        a_n = GET_SHIFT_F(s, ADIABATIC_N);
        a_b = GET_SHIFT_F(s, ADIABATIC_B);
        a_t = GET_SHIFT_F(s, ADIABATIC_T);
        
        g_p = (REAL_TYPE)( (s>>H_DIAG) & 0x7 ); // 3bitを取り出す

        delta = (dth2*( g_w * a_w * t_w  // west  
                      + g_e * a_e * t_e  // east  
                      + g_s * a_s * t_s  // south 
                      + g_n * a_n * t_n  // north 
                      + g_b * a_b * t_b  // bottom
                      + g_t * a_t * t_t  // top   
                      - g_p *       t_p
                      )
                - dth1*(-(1.0-g_w)*a_w * qbc[_F_IDX_V3DEX(0, i-1, j  , k  , ix, jx, kx, gd)]
                        +(1.0-g_e)*a_e * qbc[_F_IDX_V3DEX(0, i  , j  , k  , ix, jx, kx, gd)]
                        -(1.0-g_s)*a_s * qbc[_F_IDX_V3DEX(1, i  , j-1, k  , ix, jx, kx, gd)]
                        +(1.0-g_n)*a_n * qbc[_F_IDX_V3DEX(1, i  , j  , k  , ix, jx, kx, gd)]
                        -(1.0-g_b)*a_b * qbc[_F_IDX_V3DEX(2, i  , j  , k-1, ix, jx, kx, gd)]
                        +(1.0-g_t)*a_t * qbc[_F_IDX_V3DEX(2, i  , j  , k  , ix, jx, kx, gd)]
                        ) )* GET_SHIFT_F(s, ACTIVE_BIT);
        t[m_p] = ws[m_p] + delta;
        res += delta*delta;
      }
    }
  }
  
  return res;
}


// #################################################################
// 単媒質に対する熱伝導方程式をEuler陰解法で解く
double FFV::ps_Diff_SM_PSOR(REAL_TYPE* t, double& b2, const REAL_TYPE dt, const REAL_TYPE* qbc, const int* bh2, const REAL_TYPE* ws, IterationCtl* IC, double& flop)
{
  REAL_TYPE g_p, g_w, g_e, g_s, g_n, g_b, g_t;
  REAL_TYPE t_p, t_w, t_e, t_s, t_n, t_b, t_t;
  REAL_TYPE a_p, a_w, a_e, a_s, a_n, a_b, a_t;
  REAL_TYPE dth1, dth2;
  REAL_TYPE s0, dd, delta, sb;
  double bb=0.0;
  REAL_TYPE omg;
  double res; // 残差の自乗和
  REAL_TYPE dh = (REAL_TYPE)deltaX;    /// 空間格子幅
  int s;

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  dth1 = dt/dh;
  dth2 = dth1*C.getRcpPeclet()/dh;
  omg = IC->getOmega();
  res = b2 = 0.0;
  flop += (double)(ix*jx*kx)* 58.0;

#pragma omp parallel for firstprivate(ix, jx, kx, gd, dth1, dth2, omg) \
            private(t_p, t_w, t_e, t_s, t_n, t_b, t_t) \
            private(g_p, g_w, g_e, g_s, g_n, g_b, g_t) \
            private(a_p, a_w, a_e, a_s, a_n, a_b, a_t) \
            private(s, dd, sb, s0, delta) \
            schedule(static) reduction(+:res) reduction(+:bb)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
#include "FindexS3D.h"
        
        t_p = t[m_p];
        t_w = t[m_w];
        t_e = t[m_e];
        t_s = t[m_s];
        t_n = t[m_n];
        t_b = t[m_b];
        t_t = t[m_t];

        s   = bh2[m_p];
        a_p = GET_SHIFT_F(s, ACTIVE_BIT);
        g_w = GET_SHIFT_F(s, GMA_W); // gamma(BC)
        g_e = GET_SHIFT_F(s, GMA_E);
        g_s = GET_SHIFT_F(s, GMA_S);
        g_n = GET_SHIFT_F(s, GMA_N);
        g_b = GET_SHIFT_F(s, GMA_B);
        g_t = GET_SHIFT_F(s, GMA_T);
        
        a_w = GET_SHIFT_F(s, ADIABATIC_W); // gamma(A)
        a_e = GET_SHIFT_F(s, ADIABATIC_E);
        a_s = GET_SHIFT_F(s, ADIABATIC_S);
        a_n = GET_SHIFT_F(s, ADIABATIC_N);
        a_b = GET_SHIFT_F(s, ADIABATIC_B);
        a_t = GET_SHIFT_F(s, ADIABATIC_T);

        g_p = (REAL_TYPE)( (s>>H_DIAG) & 0x7 ); // 3bitを取り出す
        dd = 1.0 / (1.0 + dth2*g_p);
        
        sb = -dth1*(-(1.0-g_w)*a_w * qbc[_F_IDX_V3DEX(0, i-1, j  , k  , ix, jx, kx, gd)]
                    +(1.0-g_e)*a_e * qbc[_F_IDX_V3DEX(0, i  , j  , k  , ix, jx, kx, gd)]
                    -(1.0-g_s)*a_s * qbc[_F_IDX_V3DEX(1, i  , j-1, k  , ix, jx, kx, gd)]
                    +(1.0-g_n)*a_n * qbc[_F_IDX_V3DEX(1, i  , j  , k  , ix, jx, kx, gd)]
                    -(1.0-g_b)*a_b * qbc[_F_IDX_V3DEX(2, i  , j  , k-1, ix, jx, kx, gd)]
                    +(1.0-g_t)*a_t * qbc[_F_IDX_V3DEX(2, i  , j  , k  , ix, jx, kx, gd)]
                    );
        bb += (double)(sb*sb * a_p);
        s0 = sb + ws[m_p]
          + dth2*( g_w * a_w * t_w  // west  
                 + g_e * a_e * t_e  // east  
                 + g_s * a_s * t_s  // south 
                 + g_n * a_n * t_n  // north 
                 + g_b * a_b * t_b  // bottom
                 + g_t * a_t * t_t  // top  
                 );
        delta = (dd * s0 - t_p) * a_p;
        t[m_p] += omg*delta;
        res += (double)(delta*delta);
      }
    }
  }
  
  b2 = bb;

	return res;
}

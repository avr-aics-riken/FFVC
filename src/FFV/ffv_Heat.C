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
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   ffv_Heat.C
 * @brief  FFV BC Class
 * @author aics
 */

#include "ffv.h"
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


// #################################################################
// 単媒質に対する熱伝導方程式を陰解法で解く
void FFV::ps_LS(LinearSolver* IC, const double b_l2, const double r0_l2)
{
  double flop = 0.0;      /// 浮動小数点演算数
  double res=0.0;         /// 残差
  double b2=0.0;          /// 反復式のソースベクトルのノルム
  double nrm = 0.0;       ///
  REAL_TYPE dt = deltaT;  /// 時間積分幅
  double var[3];          /// 誤差、残差、解
  double x_l2;            /// 解ベクトルのL2ノルム
  
  // d_ie  内部エネルギー n+1 step
  // d_ie0 内部エネルギー n step
  // d_qbc 境界条件の熱流束
  // d_ws  対流項のみの部分段階
  // d_bcd BCindex B
  
  unsigned int wait_num=0;
  int req[12];
  
	
  switch (IC->getLS()) 
  {
    case SOR:
      
      // 反復処理
      TIMING_start("Thermal_Diff_PSOR");
      
      flop = 0.0;
      res = ps_Diff_SM_PSOR(d_ie, b2, dt, d_qbc, d_bcd, d_ws, IC, flop);
      
      TIMING_stop("Thermal_Diff_PSOR", flop);
      
      // 外部周期境界条件
      TIMING_start("Thermal_Diff_OBC_Face");
      BC.OuterTBCperiodic(d_ie, ensPeriodic);
      TIMING_stop("Thermal_Diff_OBC_Face", 0.0);
      
      // 温度の同期
      if ( numProc > 1 )
      {
        TIMING_start("Sync_Thermal_Update");
        if ( paraMngr->BndCommS3D(d_ie, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
        TIMING_stop("Sync_Thermal_Update", face_comm_size*guide*sizeof(REAL_TYPE));
      }
      
      // 残差の集約
      if ( numProc > 1 )
      {
        TIMING_start("A_R_Thermal_Diff_Res");
        REAL_TYPE tmp_wk[2], m_tmp[2];
        tmp_wk[0] = m_tmp[0] = res;
        tmp_wk[1] = m_tmp[1] = b2;
        if ( paraMngr->Allreduce(tmp_wk, m_tmp, 2, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop("A_R_Thermal_Diff_Res", 2.0*numProc*2.0*sizeof(REAL_TYPE) );
        
        res = sqrt( m_tmp[0] ); // 残差のRMS
        b2  = sqrt( m_tmp[1] ); // ソースベクトルのRMS
      }
      break;

    default:
      printf("\tInvalid Linear Solver for Heat\n");
      Exit(0);
      break;
  }
  
  // Residual resを上書き
  //TIMING_start(tm_poi_src_nrm);
  //res = 0.0;
  //flop = 0.0;
  //poi_residual_(&res, size, &guide, d_p, d_ws, d_bcp, &flop);
  //TIMING_stop(tm_poi_src_nrm, flop);
  
  if ( numProc > 1 )
  {
    TIMING_start("All_Reduce");
    double m_tmp = res;
    if ( paraMngr->Allreduce(&m_tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("All_Reduce", 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
  }
  
  res = sqrt(res);
  
  
  // 残差の保存
  double ErrEPS = IC->getErrCriterion();
  
  switch ( IC->getResType() )
  {
      case nrm_r_b:
      IC->setResidual( (b_l2<ErrEPS) ? res/ErrEPS : res/b_l2 );
      break;
      
      case nrm_r_x:
      IC->setResidual( (x_l2<ErrEPS) ? res/ErrEPS : res/x_l2 );
      break;
      
      case nrm_r_r0:
      IC->setResidual( (r0_l2<ErrEPS) ? res/ErrEPS : res/r0_l2 );
      break;
      
      default:
      printf("\tInvalid Residual Norm for Pressure\n");
      Exit(0);
      break;
  }
  
}


// #################################################################
// 単媒質に対する熱伝導方程式をEuler陰解法で解く
double FFV::ps_Diff_SM_PSOR(REAL_TYPE* t, double& b_l2, const REAL_TYPE dt, const REAL_TYPE* qbc, const int* bh, const REAL_TYPE* ws, IterationCtl* IC, double& flop)
{
  REAL_TYPE g_p, g_w, g_e, g_s, g_n, g_b, g_t;
  REAL_TYPE t_p, t_w, t_e, t_s, t_n, t_b, t_t;
  REAL_TYPE a_p, a_w, a_e, a_s, a_n, a_b, a_t;
  REAL_TYPE dth1, dth2;
  REAL_TYPE s0, dd, delta, sb;
  double bb=0.0;
  REAL_TYPE omg;
  double res; // 残差の自乗和
  REAL_TYPE dh = (REAL_TYPE)pitch[0];    /// 空間格子幅
  int s;

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  dth1 = dt/dh;
  dth2 = dth1*C.getRcpPeclet()/dh;
  omg = IC->getOmega();
  res = b_l2 = 0.0;
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
        
#include "../FB/FindexS3D.h"
        
        t_p = t[m_p];
        t_w = t[m_w];
        t_e = t[m_e];
        t_s = t[m_s];
        t_n = t[m_n];
        t_b = t[m_b];
        t_t = t[m_t];

        s   = bh[m_p];
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
        
        sb = -dth1*(-(1.0-g_w)*a_w * qbc[_F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd)]
                    +(1.0-g_e)*a_e * qbc[_F_IDX_S4DEX(1, i, j, k, 6, ix, jx, kx, gd)]
                    -(1.0-g_s)*a_s * qbc[_F_IDX_S4DEX(2, i, j, k, 6, ix, jx, kx, gd)]
                    +(1.0-g_n)*a_n * qbc[_F_IDX_S4DEX(3, i, j, k, 6, ix, jx, kx, gd)]
                    -(1.0-g_b)*a_b * qbc[_F_IDX_S4DEX(4, i, j, k, 6, ix, jx, kx, gd)]
                    +(1.0-g_t)*a_t * qbc[_F_IDX_S4DEX(5, i, j, k, 6, ix, jx, kx, gd)]
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
  
  b_l2 = bb;

	return res;
}

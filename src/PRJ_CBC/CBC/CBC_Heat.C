/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file 3D_Heat.C
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

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

#include "SklSolverCBC.h"

/**
 @fn void SklSolverCBC::ps_ConvectionEE(SKL_REAL* tc, SKL_REAL dt, unsigned* bd, SKL_REAL* t0, SKL_REAL& flop)
 @brief 移流項のEuler陽解法による時間積分
 @param tc[in/out] 対流項の流束の和/部分段階の温度
 @param dt 時間積分幅
 @param bd BCindex ID
 @param t0 nステップの温度
 @param flop 浮動小数演算数
 @note
    - tc = t0 + dt/dh*sum_flux(n)
 */
void SklSolverCBC::ps_ConvectionEE(SKL_REAL* tc, SKL_REAL dt, unsigned* bd, SKL_REAL* t0, SKL_REAL& flop)
{
  int i,j,k;
  unsigned register m;
  
  flop += (SKL_REAL)(ixc*jxc*kxc*3);
  
  for (k=1; k<=kxc; k++) {
    for (j=1; j<=jxc; j++) {
      for (i=1; i<=ixc; i++) {
        m = SklUtil::getFindexS3D(size, guide, i  , j  , k  );
        tc[m] = t0[m] + dt * tc[m] * GET_SHIFT_F(bd[m], STATE_BIT); // 対流項の評価は流動なので，状態を参照
      }
    }
  }
}

/**
 @fn void SklSolverCBC::Buoyancy(SKL_REAL* v, SKL_REAL dgr, SKL_REAL* t, unsigned* bd, SKL_REAL& flop)
 @brief Boussinesq浮力項の計算
 @retval 浮動小数演算数
 @param v 速度
 @param dgr 係数 
 @param t 温度
 @param bd BCindex ID
 @param flop 浮動小数点演算数
 */
void SklSolverCBC::Buoyancy(SKL_REAL* v, SKL_REAL dgr, SKL_REAL* t, unsigned* bd, SKL_REAL& flop)
{
  int i,j,k;
  unsigned register m, l;
  
  flop += (SKL_REAL)(ixc*jxc*kxc*2);
  
  for (k=1; k<=kxc; k++) {
    for (j=1; j<=jxc; j++) {
      for (i=1; i<=ixc; i++) {
        m = SklUtil::getFindexS3D(size, guide, i  , j  , k  );
        l = SklUtil::getFindexV3DEx(size, guide, 2, i  , j  , k  ); // k方向が重力方向
        v[l] += dgr*t[m]* GET_SHIFT_F(bd[m], STATE_BIT); // 対流項の計算なので，状態を参照
      }
    }
  }
}

/**
 @fn void SklSolverCBC::ps_LS(ItrCtl* IC)
 @brief 単媒質に対する熱伝導方程式を陰解法で解く
 @param IC IterationCtlクラス
 */
void SklSolverCBC::ps_LS(ItrCtl* IC)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  SKL_REAL *t=NULL;             /// 温度 n+1 step
  SKL_REAL *t0=NULL;             /// 温度 n step
  SKL_REAL *qbc=NULL;           /// 境界条件の熱流束
  SKL_REAL *ws=NULL;            /// 対流項のみの部分段階
  unsigned *bh2=NULL;           /// BCindex H2
	SKL_REAL res=0.0;             /// 残差
  SKL_REAL b2=0.0;              /// 反復式のソースベクトルのノルム
  SKL_REAL flop_count=0.0;
  SKL_REAL np_f = (SKL_REAL)para_mng->GetNodeNum(pn.procGrp); /// 全ノード数
  SKL_REAL comm_size = count_comm_size(size, guide);       /// 通信面1面あたりの通信量
  SKL_REAL dt = SklGetDeltaT(); /// 時間積分幅
  SKL_REAL res_r =0.0;
  SKL_REAL nrm = 0.0;
  
  unsigned int wait_num=0;
  int req[12];
  
  if( !(t   = dc_t->GetData()) )    assert(0);
  if( !(t0  = dc_t0->GetData()) )   assert(0);
  if( !(qbc = dc_qbc->GetData()) )  assert(0);
  if( !(bh2 = dc_bh2->GetData()) )  assert(0);
  if( !(ws  = dc_ws->GetData()) )   assert(0);
	
  switch (IC->get_LS()) {
    case SOR:
      
      // >>> Passive scalar Diffusion subsection 3
      TIMING_start(tm_heat_diff_sct_3);
      
      // 反復処理
      TIMING_start(tm_heat_diff_PSOR);
      
      flop_count = 0.0;
      res = ps_Diff_SM_PSOR(t, b2, dt, qbc, bh2, ws, IC, flop_count);
      
      TIMING_stop(tm_heat_diff_PSOR, flop_count);
      
      // 外部境界条件
      TIMING_start(tm_heat_diff_OBC);
      BC.OuterTBC(dc_t);
      TIMING_stop(tm_heat_diff_OBC, 0.0);
      
      // 温度の同期
      if ( para_mng->IsParallel() ) {
        TIMING_start(tm_heat_update_comm);
        dc_t->CommBndCell(guide);
        TIMING_stop(tm_heat_update_comm, comm_size*(SKL_REAL)guide);
      }
      
      // 残差の集約
      if ( para_mng->IsParallel() ) {
        TIMING_start(tm_heat_diff_res_comm);
        SKL_REAL tmp_wk[2], m_tmp[2];
        tmp_wk[0] = m_tmp[0] = res;
        tmp_wk[1] = m_tmp[1] = b2;
        para_mng->Allreduce(tmp_wk, m_tmp, 2, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
        TIMING_stop( tm_heat_diff_res_comm, 2.0*np_f*2.0*(SKL_REAL)sizeof(SKL_REAL) );
        
        res = sqrt( m_tmp[0]/(SKL_REAL)G_Acell ); // 残差のRMS
        b2  = sqrt( m_tmp[1]/(SKL_REAL)G_Acell ); // ソースベクトルのRMS
      }
      
      // 残差の保存
      res_r = (b2 == 0.0) ? res : res/b2;
      nrm = ( IC->get_normType() == ItrCtl::t_res_l2_r ) ? res_r : res;
      IC->set_normValue( nrm );
      
      TIMING_stop(tm_heat_diff_sct_3, 0.0);
      // <<< Passive scalar Diffusion subsection 3
      break;

    default:
      printf("\tInvalid Linear Solver for Heat\n");
      assert(0);
      break;
  }
}

/**
 @fn SKL_REAL SklSolverCBC::ps_Diff_SM_EE(SKL_REAL* t, SKL_REAL dt, SKL_REAL* qbc, unsigned* bx, SKL_REAL* ws, SKL_REAL& flop)
 @brief 単媒質に対する熱伝導方程式をEuler陽解法で解く
 @retval 拡散項の変化量Δθの絶対値
 @param t n+1時刻の温度場
 @param dt 時間積分幅
 @param qbc 境界条件熱流束
 @param bh2 BCindex H2
 @param ws 部分段階の温度
 @param flop 浮動小数点演算数
 */
SKL_REAL SklSolverCBC::ps_Diff_SM_EE(SKL_REAL* t, SKL_REAL dt, SKL_REAL* qbc, unsigned* bh2, SKL_REAL* ws, SKL_REAL& flop)
{
  int i, j, k;
  unsigned m_p, m_w, m_e, m_s, m_n, m_b, m_t;
  SKL_REAL g_p, g_w, g_e, g_s, g_n, g_b, g_t;
  SKL_REAL t_p, t_w, t_e, t_s, t_n, t_b, t_t;
  SKL_REAL      a_w, a_e, a_s, a_n, a_b, a_t;
  SKL_REAL dth1, dth2, delta, res;
  unsigned register s;

  dth1 = dt/C.dh;
  dth2 = dth1*C.getRcpPeclet()/C.dh;
  res  = 0.0;
  flop += (SKL_REAL)(ixc*jxc*kxc)* 50.0;

  for (k=1; k<=kxc; k++) {
    for (j=1; j<=jxc; j++) {
      for (i=1; i<=ixc; i++) {
        m_p = SklUtil::getFindexS3D(size, guide, i  , j  , k  );
        m_w = SklUtil::getFindexS3D(size, guide, i-1, j  , k  );
        m_e = SklUtil::getFindexS3D(size, guide, i+1, j  , k  );
        m_s = SklUtil::getFindexS3D(size, guide, i  , j-1, k  );
        m_n = SklUtil::getFindexS3D(size, guide, i  , j+1, k  );
        m_b = SklUtil::getFindexS3D(size, guide, i  , j  , k-1);
        m_t = SklUtil::getFindexS3D(size, guide, i  , j  , k+1);
        
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
        
        g_p = (SKL_REAL)( (s>>H_DIAG) & 0x7 ); // 3bitを取り出す

        delta = (dth2*( g_w * a_w * t_w  // west  
                      + g_e * a_e * t_e  // east  
                      + g_s * a_s * t_s  // south 
                      + g_n * a_n * t_n  // north 
                      + g_b * a_b * t_b  // bottom
                      + g_t * a_t * t_t  // top   
                      - g_p *       t_p
                      )
                - dth1*(-(1.0-g_w)*a_w * qbc[SklUtil::getFindexV3DEx(size, guide, 0, i-1, j  , k  )]
                        +(1.0-g_e)*a_e * qbc[SklUtil::getFindexV3DEx(size, guide, 0, i  , j  , k  )]
                        -(1.0-g_s)*a_s * qbc[SklUtil::getFindexV3DEx(size, guide, 1, i  , j-1, k  )]
                        +(1.0-g_n)*a_n * qbc[SklUtil::getFindexV3DEx(size, guide, 1, i  , j  , k  )]
                        -(1.0-g_b)*a_b * qbc[SklUtil::getFindexV3DEx(size, guide, 2, i  , j  , k-1)]
                        +(1.0-g_t)*a_t * qbc[SklUtil::getFindexV3DEx(size, guide, 2, i  , j  , k  )]
                        ) )* GET_SHIFT_F(s, ACTIVE_BIT);
        t[m_p] = ws[m_p] + delta;
        res += delta*delta;
      }
    }
  }
  
  return (res);
}

/**
 @fn SKL_REAL SklSolverCBC::ps_Diff_SM_PSOR(SKL_REAL* t, SKL_REAL& b2, SKL_REAL dt, SKL_REAL* qbc, unsigned* bh2, SKL_REAL* ws, ItrCtl* IC, SKL_REAL& flop)
 @brief 単媒質に対する熱伝導方程式をEuler陰解法で解く
 @retval ローカルノードの変化量の自乗和
 @param t n+1時刻の温度場
 @param b2[out] ソースベクトルの自乗和
 @param dt 時間積分幅
 @param qbc 境界条件熱流束
 @param bh2 BCindex H2
 @param ws 部分段階の温度
 @param IC ItrCtlクラス
 @param flop[out] 浮動小数点演算数
 @todo 
    - 残差確認，並列時と逐次の比較
 */
SKL_REAL SklSolverCBC::ps_Diff_SM_PSOR(SKL_REAL* t, SKL_REAL& b2, SKL_REAL dt, SKL_REAL* qbc, unsigned* bh2, SKL_REAL* ws, ItrCtl* IC, SKL_REAL& flop)
{
  int i,j,k;
  unsigned m_p, m_w, m_e, m_s, m_n, m_b, m_t;
  SKL_REAL g_p, g_w, g_e, g_s, g_n, g_b, g_t;
  SKL_REAL t_p, t_w, t_e, t_s, t_n, t_b, t_t;
  SKL_REAL a_p, a_w, a_e, a_s, a_n, a_b, a_t;
  SKL_REAL dth1, dth2;
  SKL_REAL s0, dd, delta, sb;
  SKL_REAL omg;
  SKL_REAL res; // 残差の自乗和
  unsigned register s;

  dth1 = dt/C.dh;
  dth2 = dth1*C.getRcpPeclet()/C.dh;
  omg = IC->get_omg();
  res = b2 = 0.0;
  flop += (SKL_REAL)(ixc*jxc*kxc)* 58.0;

  for (k=1; k<=kxc; k++) {
    for (j=1; j<=jxc; j++) {
      for (i=1; i<=ixc; i++) {
        m_p = SklUtil::getFindexS3D(size, guide, i  , j  , k  );
        m_w = SklUtil::getFindexS3D(size, guide, i-1, j  , k  );
        m_e = SklUtil::getFindexS3D(size, guide, i+1, j  , k  );
        m_s = SklUtil::getFindexS3D(size, guide, i  , j-1, k  );
        m_n = SklUtil::getFindexS3D(size, guide, i  , j+1, k  );
        m_b = SklUtil::getFindexS3D(size, guide, i  , j  , k-1);
        m_t = SklUtil::getFindexS3D(size, guide, i  , j  , k+1);
        
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

        g_p = (SKL_REAL)( (s>>H_DIAG) & 0x7 ); // 3bitを取り出す
        dd = 1.0 / (1.0 + dth2*g_p);
        
        sb = -dth1*(-(1.0-g_w)*a_w * qbc[SklUtil::getFindexV3DEx(size, guide, 0, i-1, j  , k  )]
                    +(1.0-g_e)*a_e * qbc[SklUtil::getFindexV3DEx(size, guide, 0, i  , j  , k  )]
                    -(1.0-g_s)*a_s * qbc[SklUtil::getFindexV3DEx(size, guide, 1, i  , j-1, k  )]
                    +(1.0-g_n)*a_n * qbc[SklUtil::getFindexV3DEx(size, guide, 1, i  , j  , k  )]
                    -(1.0-g_b)*a_b * qbc[SklUtil::getFindexV3DEx(size, guide, 2, i  , j  , k-1)]
                    +(1.0-g_t)*a_t * qbc[SklUtil::getFindexV3DEx(size, guide, 2, i  , j  , k  )]
                    );
        b2 += sb*sb * a_p;
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
        res += delta*delta;
      }
    }
  }

	return res;
}

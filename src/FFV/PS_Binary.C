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
 * @file   PS_Binary.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"

// 温度の移流拡散方程式をEuler陽解法/Adams-Bashforth法で解く
void FFV::PS_Binary()
{
  // local variables
  double flop;                         /// 浮動小数演算数
  double rhs_nrm = 0.0;                /// 反復解法での定数項ベクトルのL2ノルム
  double res_init = 0.0;               /// 反復解法での初期残差ベクトルのL2ノルム
  double convergence=0.0;              /// 定常収束モニター量
  double res=0.0;                      /// 残差
  
  REAL_TYPE dt = deltaT;               /// 時間積分幅
  REAL_TYPE dh = (REAL_TYPE)deltaX;    /// 空間格子幅
  REAL_TYPE pei=C.getRcpPeclet();      /// ペクレ数の逆数
  
  REAL_TYPE half = 0.5;                /// 定数
  REAL_TYPE one = 1.0;                 /// 定数
  REAL_TYPE zero = 0.0;                /// 定数
  int cnv_scheme = C.CnvScheme;        /// 対流項スキーム
  
  IterationCtl* ICt = &IC[ic_tmp1];  /// 拡散項の反復

  // point Data
  // d_v   セルセンタ速度 v^{n+1}
  // d_t   温度 t^n -> t^{n+1}
  // d_t0  温度 t^n
  // d_qbc 温度のソース項
  // d_ws  ワーク
  // d_bh1 温度のビットフラグ
  // d_bh2 温度のビットフラグ
  // d_bcv 速度のビットフラグ
  
  
  // >>> Passive scalar Convection section
  TIMING_start(tm_heat_convection_sct);
  
  // n stepの値をdc_t0に保持，dc_tはn+1レベルの値として利用
  TIMING_start(tm_copy_array);
  U.xcopy(d_t0, size, guide, d_t, one, kind_scalar);
  TIMING_stop(tm_copy_array, 0.0);
  
  
  // 指定境界条件の参照値を代入する
  TIMING_start(tm_heat_spec_temp);
  //BC.assignTemp(d_t0, d_bh1, CurrentTime, &C);
  TIMING_stop(tm_heat_spec_temp, 0.0);
  
  
  // 対流項の寄与
  if ( C.KindOfSolver != SOLID_CONDUCTION) // 流れの場合，対流項の寄与分のみを積分しdc_tで保持
  {
    TIMING_start(tm_heat_cnv);
    flop = 0.0;
    int swt = 0; // 断熱壁
    ps_muscl_(d_ws, size, &guide, &dh, &cnv_scheme, v00, d_v, d_t0, d_bcv, d_bcp, d_bh1, d_bh2, &swt, &flop);
    TIMING_stop(tm_heat_cnv, flop);

		// 対流フェイズの流束型境界条件
    TIMING_start(tm_heat_cnv_BC);
    flop=0.0;
		BC.OuterTBCconvection(d_ws, d_bh1, d_vf, d_t0, CurrentTime, &C, v00);
    TIMING_stop(tm_heat_cnv_BC, flop);
		
    // 時間積分
    TIMING_start(tm_heat_cnv_EE);
    flop = 0.0;
    ps_ConvectionEE(d_ws, dt, d_bh2, d_t0, flop);
    TIMING_stop(tm_heat_cnv_EE, flop);
  }
  else // 熱伝導の場合，対流項の寄与分はないので前ステップの値
  {
    TIMING_start(tm_copy_array);
    U.xcopy(d_ws, size, guide, d_t0, one, kind_scalar);
    TIMING_stop(tm_copy_array, 0.0);
  }

  
  // 外部周期境界条件
  TIMING_start(tm_heat_diff_OBC);
  BC.OuterTBCperiodic(d_ws);
  TIMING_stop(tm_heat_diff_OBC, 0.0);
  
  
  TIMING_stop(tm_heat_convection_sct, 0.0);
  // <<< Passive scalar Convection section
  
  
  
  // 拡散項の計算
  // >>> Passive scalar Diffusion section
  TIMING_start(tm_heat_diffusion_sct);
  
  
  // >>> Passive scalar Diffusion subsection 1
  TIMING_start(tm_heat_diff_sct_1);
  
  
  // 内部境界条件 体積要素
  TIMING_start(tm_heat_diff_IBC_vol);
  BC.InnerTBCvol(d_ws, d_bh2, dt);
  TIMING_stop(tm_heat_diff_IBC_vol, 0.0);
  
  
  // 部分段階の温度の同期
  if ( numProc > 1 )
  {
    TIMING_start(tm_heat_diff_comm);
    if ( paraMngr->BndCommS3D(d_ws, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    TIMING_stop(tm_heat_diff_comm, face_comm_size*guide);
  }
  
  TIMING_stop(tm_heat_diff_sct_1, 0.0);
  // <<< Passive scalar Diffusion subsection 1
  
  
  
  // >>> Passive scalar Diffusion subsection 2
  TIMING_start(tm_heat_diff_sct_2);
  
  
  // 熱流束境界条件のクリア qbcは積算するため
  TIMING_start(tm_assign_const);
  U.xset(d_qbc, size, guide, zero, kind_vector);
  TIMING_stop(tm_assign_const, 0.0);
  
  
  // 内部境界条件　熱流束型の境界条件は時間進行の前
  TIMING_start(tm_heat_diff_IBC_face);
  BC.InnerTBCface(d_qbc, d_bh1, d_ws, d_t0); // 境界値はt^{n}から計算すること
  TIMING_stop(tm_heat_diff_IBC_face, 0.0);
  
  
  TIMING_start(tm_heat_diff_OBC_face);
  BC.OuterTBCdiffusion(d_qbc, d_ws, d_t0, &C);
  TIMING_stop(tm_heat_diff_OBC_face, 0.0);
  
  
  // 境界条件の熱流束の同期 >>　不要？
  if ( numProc > 1 )
  {
    TIMING_start(tm_heat_diff_QBC_comm);
    if ( paraMngr->BndCommV3DEx(d_qbc, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    TIMING_stop(tm_heat_diff_QBC_comm, face_comm_size*3.0); // 3成分
  }
  
  TIMING_stop(tm_heat_diff_sct_2, 0.0);
  // <<< Passive scalar Diffusion subsection 2
  
  
  if ( C.AlgorithmH == Heat_EE_EE ) // 陽的時間進行
  {
    // >>> Passive scalar Diffusion subsection 3
    TIMING_start(tm_heat_diff_sct_3);
    
    TIMING_start(tm_heat_diff_EE);
    flop = 0.0;
    //res = ps_Diff_SM_EE(t, dt, qbc, bh2, ws, flop); // resは拡散項のみの絶対残差
    ps_diff_ee_(d_t, size, &guide, &res, &dh, &dt, &pei, d_qbc, d_bh2, d_ws, &flop);
    TIMING_stop(tm_heat_diff_EE, flop);
    

    
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
      double tmp = res;
      if ( paraMngr->Allreduce(&tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      TIMING_stop( tm_heat_diff_res_comm, 2.0*numProc*sizeof(REAL_TYPE) ); // 双方向 x ノード数
    }
    ICt->setNormValue( sqrt(res/(double)G_Acell) ); // RMS
    
    
    TIMING_stop(tm_heat_diff_sct_3, 0.0);
    // <<< Passive scalar Diffusion subsection 3
  }
  else // 陰解法
  {
    // 反復初期値
    TIMING_start(tm_copy_array);
    U.xcopy(d_t, size, guide, d_ws, one, kind_scalar);
    TIMING_stop(tm_copy_array, 0.0);
    
    for (ICt->setLoopCount(0); ICt->getLoopCount()< ICt->getMaxIteration(); ICt->incLoopCount())
    {

      // 線形ソルバー
      ps_LS(ICt, rhs_nrm, res_init);
      
      switch (ICt->getNormType())
      {
        case dx_b:
        case r_b:
        case r_r0:
          convergence = ICt->getNormValue();
          break;
          
        default:
          stamped_printf("\tInvalid convergence type\n");
          Exit(0);
      }
      
      if ( convergence < ICt->getCriterion() ) break;
    }
    
  }
  
  TIMING_stop(tm_heat_diffusion_sct, 0.0);
  // <<< Passive scalar Diffusion section
  
  
  // >>> Passive scalar Post
  TIMING_start(tm_heat_loop_post_sct);
  
  // 変数のカットオフオプション
  TIMING_start(tm_heat_range);
  if ( C.Hide.Range_Limit == Control::Range_Cutoff )
  {
      fb_limit_scalar_(d_t, size, &guide);
  }
  TIMING_stop(tm_heat_range, 0.0);
  
  TIMING_stop(tm_heat_loop_post_sct, 0.0);
  // <<< Passive scalar Post
  
  
  // ノルムの増加率が規定値をこえたら，終了
  if (CM_H.previous != 0.0 )
  {
    CM_H.rate = convergence / CM_H.previous;
  }
  else
  {
    CM_H.rate = 1.0;
  }
  CM_H.previous = convergence;
  
}

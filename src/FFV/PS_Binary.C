//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   PS_Binary.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"

// 温度の移流拡散方程式をEuler陽解法/Adams-Bashforth法で解く
void FFV::PS_Binary()
{
  // local variables
  double flop;                         /// 浮動小数演算数
  double b_l2 = 0.0;                   /// 反復解法での定数項ベクトルのL2ノルム
  double res0_l2 = 0.0;                /// 反復解法での初期残差ベクトルのL2ノルム
  double res=0.0;                      /// 残差
  
  REAL_TYPE dt = deltaT;               /// 時間積分幅
  REAL_TYPE dh = (REAL_TYPE)deltaX;    /// 空間格子幅
  REAL_TYPE pei=C.getRcpPeclet();      /// ペクレ数の逆数
  REAL_TYPE coef = C.RefDensity * C.RefSpecificHeat * C.RefVelocity * C.RefLength;
  
  REAL_TYPE half = 0.5;                /// 定数
  REAL_TYPE one = 1.0;                 /// 定数
  REAL_TYPE zero = 0.0;                /// 定数
  int cnv_scheme = C.CnvScheme;        /// 対流項スキーム
  
  LinearSolver* LSt = &LS[ic_tmp1];    /// 拡散項の反復

  // point Data
  // d_v   セルセンタ速度 v^{n+1}
  // d_ie  内部エネルギー ie^n -> ie^{n+1}
  // d_ie0 内部エネルギー ie^n
  // d_qbc 熱流束のソース項
  // d_ws  ワーク
  // d_bcd BCindex B
  // d_cdf Component Directional BC Flag
  
  
  // n stepの値をd_ie0に保持，d_ieはn+1レベルの値として利用
  TIMING_start("Copy_Array");
  U.copyS3D(d_ie0, size, guide, d_ie, one);
  TIMING_stop("Copy_Array", 0.0);
  
  
  // 積算量のクリア，移動熱量は対流部と拡散部に分けて積算するため
  for (int i=0; i<NOFACE; i++) {
    C.H_Dface[i] = 0.0;
  }
  
  
  // 対流項の寄与
  if ( C.KindOfSolver != SOLID_CONDUCTION) // 流れの場合，対流項の寄与分のみを積分しdc_tで保持
  {
    TIMING_start("Thermal_Convection");
    flop = 0.0;
    int swt = 0; // 断熱壁
    ps_muscl_(d_ws, size, &guide, &dh, &cnv_scheme, v00, d_v, d_ie0, d_bcp, d_cdf, d_bcd, &swt, &flop);
    TIMING_stop("Thermal_Convection", flop);

		// 対流フェイズの流束型境界条件
    TIMING_start("Thermal_Convection_BC");
    flop=0.0;
		BC.TBCconvection(d_ws, d_cdf, d_vf, d_ie0, CurrentTime, &C, v00);
    TIMING_stop("Thermal_Convection_BC", flop);
		
    // 時間積分
    TIMING_start("Thermal_Convection_EE");
    flop = 0.0;
    ps_ConvectionEE(d_ws, dt, d_bcd, d_ie0, flop);
    TIMING_stop("Thermal_Convection_EE", flop);
  }
  else // 熱伝導の場合，対流項の寄与分はないので前ステップの値
  {
    TIMING_start("Copy_Array");
    U.copyS3D(d_ws, size, guide, d_ie0, one);
    TIMING_stop("Copy_Array", 0.0);
  }

  
  // 外部周期境界条件
  TIMING_start("Thermal_Diff_Outer_BC");
  BC.OuterTBCperiodic(d_ws, ensPeriodic);
  TIMING_stop("Thermal_Diff_Outer_BC", 0.0);
  
  
  
  // 拡散項の計算
  
  // 内部境界条件 体積要素
  TIMING_start("Thermal_Diff_IBC_Vol");
  BC.InnerTBCvol(d_ws, d_bcd, dt);
  TIMING_stop("Thermal_Diff_IBC_Vol", 0.0);
  
  
  // 部分段階の温度の同期
  if ( numProc > 1 )
  {
    TIMING_start("Sync_Thermal");
    if ( paraMngr->BndCommS3D(d_ws, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("Sync_Thermal", face_comm_size*guide*sizeof(REAL_TYPE));
  }
  
  
  // 熱流束境界条件のクリア qbcは積算するため
  TIMING_start("assign_Const_to_Array");
  U.initS4DEX(d_qbc, size, guide, zero);
  TIMING_stop("assign_Const_to_Array", 0.0);
  
  
  // 内部境界条件　熱流束型の境界条件は時間進行の前
  TIMING_start("Thermal_Diff_IBC_Face");
  BC.InnerTBCface(d_qbc, d_cdf, d_bcd, d_ws, d_ie0); // 境界値はt^{n}から計算すること
  TIMING_stop("Thermal_Diff_IBC_Face", 0.0);
  
  
  TIMING_start("Thermal_Diff_OBC_Face");
  BC.OuterTBCdiffusion(d_qbc, d_ws, d_ie0, d_bcd, &C);
  TIMING_stop("Thermal_Diff_OBC_Face", 0.0);
  
  
  // 境界条件の熱流束の同期 >>　不要？
  if ( numProc > 1 )
  {
    TIMING_start("Sync_Thermal_QBC");
    if ( paraMngr->BndCommS4DEx(d_qbc, 6, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("Sync_Thermal_QBC", face_comm_size*6.0*guide*sizeof(REAL_TYPE)); // 6成分
  }

  
  
  if ( C.AlgorithmH == Heat_EE_EE ) // 陽的時間進行
  {
    
    TIMING_start("Thermal_Diff_EE");
    flop = 0.0;
    //res = ps_Diff_SM_EE(t, dt, qbc, bh, ws, flop); // resは拡散項のみの修正ベクトルの自乗和
    
    int h_mode;
    if ( (C.KindOfSolver == CONJUGATE_HT) || (C.KindOfSolver == CONJUGATE_HT_NATURAL) )
    {
      h_mode = 0;
    }
    else
    {
      h_mode = 1;
    }
    ps_diff_ee_(d_ie, size, &guide, &res, &dh, &dt, d_qbc, d_bcd, d_ws, &C.NoCompo, mat_tbl, &h_mode, &flop);
    
    TIMING_stop("Thermal_Diff_EE", flop);
    

    
    // 外部周期境界条件
    TIMING_start("Thermal_Diff_Outer_BC");
    BC.OuterTBCperiodic(d_ie, ensPeriodic);
    TIMING_stop("Thermal_Diff_Outer_BC", 0.0);
    
    
    // 温度の同期
    if ( numProc > 1 )
    {
      TIMING_start("Sync_Thermal");
      if ( paraMngr->BndCommS3D(d_ie, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      TIMING_stop("Sync_Thermal", face_comm_size*guide*sizeof(REAL_TYPE));
    }
    
    
    // 修正ベクトルの自乗和の集約
    if ( numProc > 1 )
    {
      TIMING_start("A_R_Thermal_Diff_Res");
      double tmp = res;
      if ( paraMngr->Allreduce(&tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      TIMING_stop("A_R_Thermal_Diff_Res", 2.0*numProc*sizeof(REAL_TYPE) ); // 双方向 x ノード数
    }
    LSt->setError( sqrt(res) ); // RMS
    
    
  }
  else // 陰解法
  {
    // 反復初期値
    TIMING_start("Copy_Array");
    U.copyS3D(d_ie, size, guide, d_ws, one);
    TIMING_stop("Copy_Array", 0.0);
    
    for (LSt->setLoopCount(0); LSt->getLoopCount()< LSt->getMaxIteration(); LSt->incLoopCount())
    {

      // 線形ソルバー
      ps_LS(LSt, b_l2, res0_l2);
      
      if ( LSt->isErrConverged() || LSt->isResConverged() ) break;
    }
    
  }
  
  
  // 変数のカットオフオプション
  TIMING_start("Thermal_Range_Cut");
  if ( C.Hide.Range_Limit == Control::Range_Cutoff )
  {
      fb_limit_scalar_(d_ie, size, &guide);
  }
  TIMING_stop("Thermal_Range_Cut", 0.0);

  
}

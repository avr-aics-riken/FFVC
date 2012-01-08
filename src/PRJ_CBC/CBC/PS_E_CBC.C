/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file PS_E_CBS.C
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

#include "SklSolverCBC.h"

/**
 @fn void SklSolverCBC::PS_E_CBC(void)
 @brief 温度の移流拡散方程式をEuler陽解法/Adams-Bashforth法で解く
 */
void SklSolverCBC::PS_E_CBC(void)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  // Dimensions
  SKL_REAL *v=NULL;   /// 速度
  SKL_REAL *t=NULL;   /// 温度 n+1 step
  SKL_REAL *t0=NULL;  /// 温度 n step
  SKL_REAL *qbc=NULL; /// 境界条件の熱流束
  SKL_REAL *ws=NULL;  /// 対流項の寄与分，熱伝導計算時には前ステップの値
  unsigned *bh1=NULL; /// BCindex H1
  unsigned *bh2=NULL; /// BCindex H2
  unsigned *bcv=NULL; /// 速度のビットフラグ
  
  // local variables
  SKL_REAL tm = SklGetTotalTime();  /// 計算開始からの積算時刻
  SKL_REAL dt = SklGetDeltaT();     /// 時間積分幅
  SKL_REAL flop_count=0.0;          /// 浮動小数演算数
  SKL_REAL pei=C.getRcpPeclet();    /// ペクレ数の逆数　
  SKL_REAL res=0.0;                 /// 残差
  SKL_REAL convergence=0.0;         /// 定常収束モニター量
  SKL_REAL comm_size = 0.0;         /// 通信面1面あたりの通信量を全ノードで積算した通信量(Byte)
  SKL_REAL np_f = (SKL_REAL)para_mng->GetNodeNum(pn.procGrp);
  ItrCtl* ICt = &IC[ItrCtl::ic_tdf_ei];  /// 拡散項の反復
  
  comm_size = CU.count_comm_size(size, guide);
  
  // >>> Passive scalar Convection section
  TIMING__ PM.start(tm_heat_convection_sct);

  // point Data
  if( !(v   = dc_v->GetData()) )    assert(0);
  if( !(t   = dc_t->GetData()) )    assert(0);
  if( !(t0  = dc_t0->GetData()) )   assert(0);
  if( !(qbc = dc_qbc->GetData()) )  assert(0);
  if( !(ws  = dc_ws->GetData()) )   assert(0);
  if( !(bh1 = dc_bh1->GetData()) )  assert(0);
  if( !(bh2 = dc_bh2->GetData()) )  assert(0);
  if( !(bcv = dc_bcv->GetData()) )  assert(0);
  
  // n stepの値をdc_t0に保持，dc_tはn+1レベルの値として利用
  TIMING__ PM.start(tm_copy_array);
  CU.copy_SKL_REAL(t0, t, dc_t0->GetArrayLength());
  TIMING__ PM.stop(tm_copy_array, 0.0);
  
  // 速度指定境界条件の参照値を代入する
  TIMING__ PM.start(tm_heat_spec_temp);
  BC.assign_Temp(t0, bh1, tm, &C);
  TIMING__ PM.stop(tm_heat_spec_temp, 0.0);
  
  // 対流項の寄与
  if ( C.KindOfSolver != SOLID_CONDUCTION) { // 流れの場合，対流項の寄与分のみを積分しdc_tで保持
    
    TIMING__ PM.start(tm_heat_cnv);
    flop_count = 0.0;
    int swt = 0; // 断熱壁
    cbc_ps_muscl_(ws, sz, gc, dh, (int*)&C.CnvScheme, v00, v, t0, (int*)bcv, (int*)bh1, (int*)bh2, &swt, &flop_count);
    TIMING__ PM.stop(tm_heat_cnv, flop_count);

		// 対流フェイズの流束型境界条件
    TIMING__ PM.start(tm_heat_cnv_BC);
    flop_count=0.0;
		BC.ps_BC_Convection(ws, bh1, v, t0, tm, &C, v00, flop_count);
    TIMING__ PM.stop(tm_heat_cnv_BC, flop_count);
		
    // 時間積分
    TIMING__ PM.start(tm_heat_cnv_EE);
    flop_count = 0.0;
    ps_ConvectionEE(ws, dt, bh2, t0, flop_count);
    TIMING__ PM.stop(tm_heat_cnv_EE, flop_count);
  }
  else { // 熱伝導の場合，対流項の寄与分はないので前ステップの値
    TIMING__ PM.start(tm_copy_array);
    CU.copy_SKL_REAL(ws, t0, dc_ws->GetArrayLength());
    TIMING__ PM.stop(tm_copy_array, 0.0);
  }

  TIMING__ PM.stop(tm_heat_convection_sct, 0.0);
  // <<< Passive scalar Convection section
  
  
  
  // 拡散項の計算
  // >>> Passive scalar Diffusion section
  TIMING__ PM.start(tm_heat_diffusion_sct);
  
  // 外部境界条件
  TIMING__ PM.start(tm_heat_diff_OBC);
  BC.OuterTBC(dc_ws);
  TIMING__ PM.stop(tm_heat_diff_OBC, 0.0);
  
  
  
  // >>> Passive scalar Diffusion subsection 1
  TIMING__ PM.start(tm_heat_diff_sct_1);
  
  // 内部境界条件 体積要素
  TIMING__ PM.start(tm_heat_diff_IBC_vol);
  flop_count = 0.0;
  BC.InnerTBCvol(ws, bh2, dt, flop_count);
  TIMING__ PM.stop(tm_heat_diff_IBC_vol, flop_count);
  
  // 部分段階の温度の同期
  if ( para_mng->IsParallel() ) {
    TIMING__ PM.start(tm_heat_diff_comm);
    dc_ws->CommBndCell(guide);
    TIMING__ PM.stop(tm_heat_diff_comm, comm_size*(SKL_REAL)guide);
  }
  
  TIMING__ PM.stop(tm_heat_diff_sct_1, 0.0);
  // <<< Passive scalar Diffusion subsection 1
  
  
  
  // >>> Passive scalar Diffusion subsection 2
  TIMING__ PM.start(tm_heat_diff_sct_2);
  
  // 熱流束境界条件のクリア qbcは積算するため
  TIMING__ PM.start(tm_assign_const);
  SklInitializeSKL_REAL(dc_qbc->GetData(), 0.0, dc_qbc->GetArrayLength());
  TIMING__ PM.stop(tm_assign_const, 0.0);
  
  // 内部境界条件　熱流束型の境界条件は時間進行の前
  TIMING__ PM.start(tm_heat_diff_IBC_face);
  flop_count = 0.0;
  BC.InnerTBCface(qbc, bh1, ws, t0, flop_count); // 境界値はt^{n}から計算すること
  TIMING__ PM.stop(tm_heat_diff_IBC_face, flop_count);
  
  TIMING__ PM.start(tm_heat_diff_OBC_face);
  BC.OuterTBCface(qbc, bh1, ws, t0, &C, flop_count);
  TIMING__ PM.stop(tm_heat_diff_OBC_face, flop_count);
  
  // 境界条件の熱流束の同期 >>　不要？
  if ( para_mng->IsParallel() ) {
    TIMING__ PM.start(tm_heat_diff_QBC_comm);
    dc_qbc->CommBndCell(1);
    TIMING__ PM.stop(tm_heat_diff_QBC_comm, comm_size*3.0); // 3成分
  }
  
  TIMING__ PM.stop(tm_heat_diff_sct_2, 0.0);
  // <<< Passive scalar Diffusion subsection 2
  
  
  
  if ( C.AlgorithmH == Control::Heat_EE_EE ) { // 陽的時間進行

    // >>> Passive scalar Diffusion subsection 3
    TIMING__ PM.start(tm_heat_diff_sct_3);
    
    TIMING__ PM.start(tm_heat_diff_EE);
    flop_count = 0.0;
    //res = ps_Diff_SM_EE(t, dt, qbc, bh2, ws, flop_count); // resは拡散項のみの絶対残差
    cbc_ps_diff_ee_(t, sz, gc, &res, dh, &dt, &pei, qbc, (int*)bh2, ws, &flop_count);
    TIMING__ PM.stop(tm_heat_diff_EE, flop_count);
    
    // 外部境界条件
    TIMING__ PM.start(tm_heat_diff_OBC);
    BC.OuterTBC(dc_t);
    TIMING__ PM.stop(tm_heat_diff_OBC, 0.0);
    
    // 温度の同期
    if ( para_mng->IsParallel() ) {
      TIMING__ PM.start(tm_heat_update_comm);
      dc_t->CommBndCell(guide);
      TIMING__ PM.stop(tm_heat_update_comm, comm_size*(SKL_REAL)guide);
    }
    
    // 残差の集約
    if ( para_mng->IsParallel() ) {
      TIMING__ PM.start(tm_heat_diff_res_comm);
      SKL_REAL tmp = res;
      para_mng->Allreduce(&tmp, &res, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
      TIMING__ PM.stop( tm_heat_diff_res_comm, 2.0*np_f*(SKL_REAL)sizeof(SKL_REAL) ); // 双方向 x ノード数
    }
    ICt->set_normValue( sqrt(res/(SKL_REAL)G_Acell) ); // RMS
    
    TIMING__ PM.stop(tm_heat_diff_sct_3, 0.0);
    // <<< Passive scalar Diffusion subsection 3
  }
  else { // 陰解法
    // 反復初期値
    TIMING__ PM.start(tm_copy_array);
    CU.copy_SKL_REAL(t, ws, dc_t->GetArrayLength());
    TIMING__ PM.stop(tm_copy_array, 0.0);
    
    for (ICt->LoopCount=0; ICt->LoopCount< ICt->get_ItrMax(); ICt->LoopCount++) {

      // 線形ソルバー
      ps_LS(ICt);
      
      switch (ICt->get_normType()) {
        case ItrCtl::t_res_l2_a:
          convergence = ICt->get_normValue();
          break;
          
        case ItrCtl::t_res_l2_r:
          convergence = ICt->get_normValue();
          break;
          
        default:
          stamped_printf("\tInvalid convergence type\n");
          assert(0);
      }
      
      if ( convergence < ICt->get_eps() ) break;
    }
    
  }
  
  TIMING__ PM.stop(tm_heat_diffusion_sct, 0.0);
  // <<< Passive scalar Diffusion section
  
  
  // >>> Passive scalar Post
  TIMING__ PM.start(tm_heat_loop_post_sct);
  
  // 変数のカットオフオプション
  TIMING__ PM.start(tm_heat_range);
  CU.CutOffRange (t, cmp, bh1, &C);
  TIMING__ PM.stop(tm_heat_range, 0.0);
  
  TIMING__ PM.stop(tm_heat_loop_post_sct, 0.0);
  // <<< Passive scalar Post
  
}
/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 20010-2012
 *
 */

//@file NS_FS_E_CDS.C
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

#include "SklSolverCBC.h"

/**
 @fn void SklSolverCBC::NS_FS_E_CDS(void)
 @brief Fractional Step法でNavier-Stokes方程式を解く．距離情報近似．
 */
void SklSolverCBC::NS_FS_E_CDS(void)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  // Dimensions
  REAL_TYPE *v =NULL;    /// セルセンタ速度 v^n -> v^{n+1}
  REAL_TYPE *v0=NULL;    /// セルセンタ速度 v^nの保持
  REAL_TYPE *vc=NULL;    /// 疑似速度ベクトル
  REAL_TYPE *wv=NULL;    /// ワーク　陰解法の時の疑似速度ベクトル，射影ステップの境界条件
  REAL_TYPE *p =NULL;    /// 圧力 p^n -> p^{n+1}
  REAL_TYPE *p0=NULL;    /// 圧力 p^nの保持
  REAL_TYPE *src0=NULL;  /// Poissonのソース項0　速度境界を考慮
  REAL_TYPE *src1=NULL;  /// Poissonのソース項1　反復毎に変化するソース項，摩擦速度，発散値
  unsigned *bcd=NULL;   /// IDのビットフラグ
  unsigned *bcp=NULL;   /// 圧力のビットフラグ
  unsigned *bcv=NULL;   /// 速度のビットフラグ
  REAL_TYPE *wss=NULL;   /// ワーク　壁関数利用時のWSS，ベクトル出力時のテンポラリ
  REAL_TYPE *t0=NULL;    /// 温度 t^n 
  REAL_TYPE *vt=NULL;    /// 渦粘性係数
  REAL_TYPE *abf=NULL;   /// Adams-Bashforth用のワーク
  float* cvf=NULL;       /// コンポーネントの体積率
  
  // local variables
  REAL_TYPE tm = SklGetTotalTime();    /// 計算開始からの積算時刻
  REAL_TYPE dt = SklGetDeltaT();       /// 時間積分幅
  REAL_TYPE flop_count;                /// 浮動小数演算数
  REAL_TYPE convergence=0.0;           /// 定常収束モニター量
  REAL_TYPE coef = C.dh / dt;          /// Poissonソース項の係数
  REAL_TYPE Re = C.Reynolds;           /// レイノルズ数
  REAL_TYPE rei = C.getRcpReynolds();  /// レイノルズ数の逆数
  REAL_TYPE b2 = 0.0;                  /// 反復解法での定数項ベクトルのノルム
  REAL_TYPE half = 0.5;                /// 定数                          
  REAL_TYPE np_f = (REAL_TYPE)para_mng->GetNodeNum(pn.procGrp); /// 全ノード数
  REAL_TYPE comm_size;                 /// 通信面1面あたりの通信量
  int wall_prof = (int)C.Mode.Wall_profile; /// 壁面条件（slip/noslip）
  int cnv_scheme = (int)C.CnvScheme;  /// 対流項スキーム
  REAL_TYPE clear_value = 0.0;
  
  comm_size = count_comm_size(size, 1);
  
  int v_mode=0;
  
  ItrCtl* ICp = &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復
  ItrCtl* ICv = &IC[ItrCtl::ic_vis_cn];  /// 粘性項のCrank-Nicolson反復
  
  // >>> Fractional step section
  TIMING_start(tm_frctnl_stp_sct);
  
  // >>> Fractional step sub-section 1
  TIMING_start(tm_frctnl_stp_sct_1); 
  
  // point Data
  if( !(v   = dc_v->GetData()) )   Exit(0);
  if( !(v0  = dc_v0->GetData()) )  Exit(0);
  if( !(vc  = dc_vc->GetData()) )  Exit(0);
  if( !(wv  = dc_wv->GetData()) )  Exit(0);
  if( !(p   = dc_p->GetData()) )   Exit(0);
  if( !(p0  = dc_p0->GetData()) )  Exit(0);
  if( !(src0= dc_ws->GetData()) )  Exit(0);
  if( !(src1= dc_wk2->GetData()) ) Exit(0); // div(u)の値を保持，出力のところまでは再利用しないこと
  if( !(bcd = dc_bcd->GetData()) ) Exit(0);
  if( !(bcp = dc_bcp->GetData()) ) Exit(0);
  if( !(bcv = dc_bcv->GetData()) ) Exit(0);
  if( !(wss = dc_wvex->GetData())) Exit(0);
  
  // 温度計算
  if ( C.isHeatProblem() ) {
    if( !(t0 = dc_t->GetData()) ) Exit(0);
  }
  
  // LES計算
  if ( C.LES.Calc ==ON ) {
    if( !(vt  = dc_vt->GetData()) )  Exit(0); //
  }
  
  // AB法
  if ( (C.AlgorithmF == Control::Flow_FS_AB2) || (C.AlgorithmF == Control::Flow_FS_AB_CN) ) {
    if( !(abf = dc_abf->GetData()) )  Exit(0);
  }
  mark();
  // コンポーネントの体積率
  if ( C.isVfraction() ) {
    if( !(cvf = dc_cvf->GetData()) )  Exit(0);
  }
  
  // IN_OUT境界条件のときのフラグ処理
  if ( BC.has_InOut() ) {
    TIMING_start(tm_flip_bf);
    BC.flipDir_OBC(bcv, &C);
    TIMING_stop(tm_flip_bf, 0.0);
  }
  
  // n stepの値を保持 >> In use (dc_v0, dc_p0)
  TIMING_start(tm_copy_array);
  fb_copy_real_v_(v0, v, sz, gc);
  fb_copy_real_s_(p0, p, sz, gc);
  TIMING_stop(tm_copy_array, 0.0, 2);
  
  // 壁関数指定時の摩擦速度の計算 src0をテンポラリのワークとして利用
  if ( C.Mode.Wall_profile == Control::Log_Law ) {
    TIMING_start(tm_WallFunc);
    flop_count = 0.0;
    cbc_friction_velocity_(src0, sz, gc, dh, &Re, v0, (int*)bcp, range_Yp, range_Ut, v00, &flop_count);
    TIMING_stop(tm_WallFunc, flop_count);
  }
  
  TIMING_stop(tm_frctnl_stp_sct_1, 0.0);
  // <<< Fractional step subsection 1
  
  mark();
  
  // >>> Fractional step sub-section 2
  TIMING_start(tm_frctnl_stp_sct_2);
  
  // 対流項と粘性項の評価 >> In use (dc_vc, dc_wv)
  switch (C.AlgorithmF) {
    case Control::Flow_FS_EE_EE:
    case Control::Flow_FS_AB2:
      TIMING_start(tm_pseudo_vec);
      flop_count = 0.0;
      v_mode = (C.Mode.Wall_profile == Control::Log_Law) ? 2 : 1;
      
      if ( C.LES.Calc == ON ) {
        Hostonly_ printf("not inplemented yet. sorry:-)\n");
        exit(0);
      }
      else {
        cds_pvec_muscl_(vc, sz, gc, dh, &cnv_scheme, v00, &rei, v0, (int*)bcv, (int*)bcp, &v_mode, cut, &flop_count); 
      }
      TIMING_stop(tm_pseudo_vec, flop_count);
      
      TIMING_start(tm_pvec_flux);
      flop_count = 0.0;
      BC.mod_Pvec_Flux(vc, v0, bcv, tm, &C, v_mode, v00, flop_count);
      TIMING_stop(tm_pvec_flux, flop_count);
      break;
      
    case Control::Flow_FS_AB_CN:
      TIMING_start(tm_pseudo_vec);
      flop_count = 0.0;
      v_mode = 0;
      if ( C.LES.Calc == ON ) {
        Hostonly_ printf("not inplemented yet. sorry:-)\n");
      }
      else {
        cds_pvec_muscl_(wv, sz, gc, dh, &cnv_scheme, v00, &rei, v0, (int*)bcv, (int*)bcp, &v_mode, cut, &flop_count); 
      }
      TIMING_stop(tm_pseudo_vec, flop_count);
      
      TIMING_start(tm_pvec_flux);
      flop_count = 0.0;
      BC.mod_Pvec_Flux(wv, v0, bcv, tm, &C, v_mode, v00, flop_count, true);
      TIMING_stop(tm_pvec_flux, flop_count);
      break;
      
    default:
      Exit(0);
  }
  
  // 時間積分
  switch (C.AlgorithmF) {
    case Control::Flow_FS_EE_EE:
      TIMING_start(tm_pvec_ee);
      flop_count = 0.0;
      cbc_ee_ (vc, sz, gc, &dt, v0, (int*)bcd, &flop_count);
      TIMING_stop(tm_pvec_ee, flop_count);
      break;
      
    case Control::Flow_FS_AB2:
      TIMING_start(tm_pvec_ab);
      flop_count = 0.0;
      if ( SklGetCurrentStep() == 1 ) { // 初期とリスタート後，1ステップめ
        cbc_ee_ (vc, sz, gc, &dt, v0, (int*)bcd, &flop_count);
      }
      else {
        cbc_ab2_(vc, sz, gc, &dt, v0, abf, (int*)bcd, v00, &flop_count);
      }
      TIMING_stop(tm_pvec_ab, flop_count);
      break;
      
    case Control::Flow_FS_AB_CN: // 未対応20110918
      TIMING_start(tm_pvec_abcn);
      flop_count = 0.0;
      if ( SklGetCurrentStep() == 1 ) {
        cbc_ee_ (wv, sz, gc, &dt, v0, (int*)bcd, &flop_count);
      }
      else {
        cbc_ab2_(wv, sz, gc, &dt, v0, abf, (int*)bcd, v00, &flop_count);
      }
      TIMING_stop(tm_pvec_abcn, flop_count);
      
      TIMING_start(tm_pvec_abcn_df_ee);
      flop_count = 0.0;
      cbc_vis_ee_(vc, sz, gc, dh, &dt, v00, &rei, wv, v0, (int*)bcv, &half, &flop_count);
      TIMING_stop(tm_pvec_abcn_df_ee, flop_count);
      
      TIMING_start(tm_pvec_abcn_df_ee_BC);
      flop_count = 0.0;
      BC.mod_Vis_EE(vc, v0, half, bcv, tm, dt, v00, flop_count);
      TIMING_stop(tm_pvec_abcn_df_ee_BC, flop_count);
      break;
      
    default:
      Exit(0);
  }
  
  TIMING_stop(tm_frctnl_stp_sct_2, 0.0);
  // <<< Fractional step subsection 2
  
  mark();
  
  // >>> Fractional step sub-section 3
  TIMING_start(tm_frctnl_stp_sct_3);
  
  // FORCINGコンポーネントの疑似速度ベクトルの方向修正
  if ( C.isForcing() == ON ) {
    TIMING_start(tm_forcing);
    flop_count = 0.0;
    BC.mod_Pvec_Forcing(vc, v, bcd, cvf, v00, dt, flop_count);
    TIMING_stop(tm_forcing, flop_count);
  }
  
  // 浮力項
  if ( C.isHeatProblem() && (C.Mode.Buoyancy == BOUSSINESQ) ) {
    TIMING_start(tm_buoyancy);
    REAL_TYPE dgr = dt*C.Grashof*rei*rei;
    flop_count = 3.0;
    Buoyancy(vc, dgr, t0, bcd, flop_count);
    TIMING_stop(tm_buoyancy, flop_count);
  }
  
  // 疑似ベクトルの境界条件
  TIMING_start(tm_pvec_BC);
  flop_count = 0.0;
  BC.OuterVBC_Pseudo(vc, v0, bcv, tm, dt, &C, v00, flop_count);
  BC.OuterVBC_Periodic(dc_vc);
  BC.InnerVBC_Periodic(dc_vc, dc_bcd);
  TIMING_stop(tm_pvec_BC, flop_count);
  
  // 疑似ベクトルの同期
  if ( para_mng->IsParallel() ) {
    TIMING_start(tm_pvec_comm);
    if ( !dc_vc->CommBndCell(1) ) Exit(0);
    TIMING_stop(tm_pvec_comm, comm_size*1.0*3.0); // ガイドセル数 x ベクトル
  }
  
  TIMING_stop(tm_frctnl_stp_sct_3, 0.0);
  // <<< Fractional step subsection 3
  
  mark();
  
  // >>> Fractional step sub-section 4
  TIMING_start(tm_frctnl_stp_sct_4);
  
  // Crank-Nicolson Iteration
  if ( C.AlgorithmF == Control::Flow_FS_AB_CN ) {
    
    TIMING_start(tm_copy_array);
    fb_copy_real_v_(wv, vc, sz, gc);
    TIMING_stop(tm_copy_array, 0.0);
    
    for (ICv->LoopCount=0; ICv->LoopCount< ICv->get_ItrMax(); ICv->LoopCount++) {
      //CN_Itr(ICv);
      if (  ICv->get_normValue() < ICv->get_eps() ) break;
    }
  }
  
  TIMING_stop(tm_frctnl_stp_sct_4, 0.0);
  // <<< Fractional step subsection 4
  
  
  TIMING_stop(tm_frctnl_stp_sct, 0.0);
  // <<< Fractional step section
  
  
  
  // Poissonのソース部分
  // >>> Poisson Source section
  TIMING_start(tm_poi_src_sct);
  
  // vの初期値をvcにしておく
  TIMING_start(tm_copy_array);
  fb_copy_real_v_(v, vc, sz, gc);
  TIMING_stop(tm_copy_array, 0.0);
  
  // 非反復ソース項のゼロクリア src0
  TIMING_start(tm_assign_const);
  fb_set_real_s_(src0, sz, gc, &clear_value);
  TIMING_stop(tm_assign_const, 0.0);
  
  // 非VBC面に対してのみ，セルセンターの値から発散量を計算
  TIMING_start(tm_div_pvec);
  flop_count = 0.0;
  cds_div_(src0, sz, gc, &coef, vc, (int*)bcv, cut, v00, &flop_count);
  TIMING_stop(tm_div_pvec, flop_count);
  
  // Poissonソース項の速度境界条件（VBC）面による修正
  TIMING_start(tm_poi_src_vbc);
  flop_count = 0.0;
  BC.mod_Psrc_VBC(src0, vc, v0, coef, bcv, tm, dt, &C, v00, flop_count, true);
  TIMING_stop(tm_poi_src_vbc, flop_count);
  
  // (Neumann_BCType_of_Pressure_on_solid_wall == grad_NS)　のとき，\gamma^{N2}の処理
  //hogehoge
  
  // 連立一次方程式の定数項の計算は圧力相対残差の場合のみ >> @todo 発散項以外の外力の影響なども含める
  if ( ICp->get_normType() == ItrCtl::p_res_l2_r) {
    TIMING_start(tm_poi_src_nrm);
    b2 = 0.0;
    cbc_div_cnst_(src0, sz, gc, &b2, (int*)bcp, &flop_count);
    b2 = sqrt(b2);
    TIMING_stop(tm_poi_src_nrm, flop_count);
    
    if ( para_mng->IsParallel() ) {
      TIMING_start(tm_poi_src_comm);
      REAL_TYPE m_tmp = b2;
      para_mng->Allreduce(&m_tmp, &b2, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
      TIMING_stop(tm_poi_src_comm, 2.0*np_f*(REAL_TYPE)sizeof(REAL_TYPE) ); // 双方向 x ノード数
    }
  }
  
  TIMING_stop(tm_poi_src_sct, 0.0);
  // <<< Poisson Source section
  
  
  // VP-Iteration
  // >>> Poisson Iteration section
  TIMING_start(tm_poi_itr_sct);
  
  if ( C.Mode.Log_Itr == ON ) {
    TIMING_start(tm_hstry_itr);
    Hostonly_ H->printHistoryItrTitle(fp_i);
    TIMING_stop(tm_hstry_itr, 0.0);
  }
  
  for (ICp->LoopCount=0; ICp->LoopCount< ICp->get_ItrMax(); ICp->LoopCount++) {
    
    
    // >>> Poisson Iteration subsection 1
    TIMING_start(tm_poi_itr_sct_1);
    
    // 反復ソース項のゼロクリア => src1
    TIMING_start(tm_assign_const);
    fb_set_real_s_(src1, sz, gc, &clear_value);
    TIMING_stop(tm_assign_const, 0.0);
    
    // Forcingコンポーネントによるソース項の寄与分
    if ( C.isForcing() == ON ) {
      TIMING_start(tm_force_src);
      flop_count=0.0;
      BC.mod_Psrc_Forcing(src1, v, bcd, cvf, coef, v00, component_array, flop_count);
      TIMING_stop(tm_force_src, flop_count);
    }
    
    // 内部周期境界部分のディリクレソース項
    //TIMING_start(tm_prdc_src);
    //BC.InnerPrdc_Src(dc_wk2, dc_p, dc_bcd);
    //TIMING_stop(tm_prdc_src, flop_count);
    
    TIMING_stop(tm_poi_itr_sct_1, 0.0);
    // <<< Poisson Iteration subsection 1
    
    // 線形ソルバー
    LS_Binary(ICp, b2);
    
    
    // >>> Poisson Iteration subsection 4
    TIMING_start(tm_poi_itr_sct_4);
    
    // 速度のスカラポテンシャルによる射影と発散値
    TIMING_start(tm_prj_vec);
    flop_count = 0.0;
    cds_update_vec_ (v, src1, sz, gc, &dt, dh, vc, p, (int*)bcp, (int*)bcv, cut, v00, &coef, &flop_count); // src1は，反復毎のソース項をワークとして利用
    TIMING_stop(tm_prj_vec, flop_count);
    
    // セルフェイス速度の境界条件による修正
    REAL_TYPE m_av[C.NoBC][2];
    TIMING_start(tm_prj_vec_bc);
    flop_count=0.0;
    BC.mod_div(src1, bcv, coef, tm, v00, m_av[C.NoBC], flop_count);
    TIMING_stop(tm_prj_vec_bc, flop_count);
    
    // セルフェイス速度の境界条件の通信部分
    if ( C.isOutflow() == ON ) {
      if ( !C.isCDS() ) { // Binary
        if ( para_mng->IsParallel() ) {
          REAL_TYPE tmp[C.NoBC][2];
          
          TIMING_start(tm_prj_vec_bc_comm);
          for (int n=1; n<=C.NoBC; n++) {
            tmp[n][0] = m_av[n][0];
            tmp[n][1] = m_av[n][1];
          }
          para_mng->Allreduce(tmp, m_av, 2*C.NoBC, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
          TIMING_stop(tm_prj_vec_bc_comm, 2.0*(REAL_TYPE)C.NoBC*np_f*(REAL_TYPE)sizeof(REAL_TYPE)*2.0 ); // 双方向 x ノード数 x 変数
        }
        
        for (int n=1; n<=C.NoBC; n++) {
          if ( cmp[n].getType() == OUTFLOW ) {
            cmp[n].val[var_Velocity] = m_av[n][0]/m_av[n][1]; // 無次元平均流速
          }
        }
      }
      else { // Cut-Distance
        ;
      }
    }
    
    // Forcingコンポーネントによる速度と発散値の修正
    if ( C.isForcing() == ON ) {
      REAL_TYPE vm[C.NoBC][2]; // モニター用
      REAL_TYPE tmp[C.NoBC][2];
      
      TIMING_start(tm_prj_frc_mod);
      flop_count=0.0;
      BC.mod_Vdiv_Forcing(v, bcd, cvf, src1, dt, C.dh, v00, vm[C.NoBC], component_array, flop_count);
      TIMING_stop(tm_prj_frc_mod, flop_count);
      
      // 通信部分
      TIMING_start(tm_prj_frc_mod_comm);
      if ( para_mng->IsParallel() ) {
        for (int n=1; n<=C.NoBC; n++) {
          tmp[n][0] = vm[n][0];
          tmp[n][1] = vm[n][1];
        }
        para_mng->Allreduce(tmp, vm, 2*C.NoBC, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
      }
      for (int n=1; n<=C.NoBC; n++) {
        if ( cmp[n].isFORCING() ) {
          vm[n][0] /= (REAL_TYPE)cmp[n].getElement();
          vm[n][1] /= (REAL_TYPE)cmp[n].getElement();
          cmp[n].val[var_Velocity] = vm[n][0]; // 平均速度
          cmp[n].val[var_Pressure] = vm[n][1]; // 平均圧力損失量
        }
      }
      TIMING_stop(tm_prj_frc_mod_comm, 2.0*(REAL_TYPE)C.NoBC*(REAL_TYPE)sizeof(REAL_TYPE)*2.0);
    }

    // 周期型の速度境界条件
    TIMING_start(tm_vec_BC);
    flop_count=0.0;
    BC.OuterVBC_Periodic(dc_v);
    BC.InnerVBC_Periodic(dc_v, dc_bcd);
    TIMING_stop(tm_vec_BC, flop_count);
    
    TIMING_stop(tm_poi_itr_sct_4, 0.0);
    // <<< Poisson Iteration subsection 4
 
    
    // ノルムの計算
    convergence = Norm_Poisson(ICp);
    
    /* Forcingコンポーネントによる速度の方向修正(収束判定から除外)  >> TEST
     TIMING_start(tm_prj_frc_dir);
     flop_count=0.0;
     BC.mod_Dir_Forcing(v, bcd, cvf, v00, flop_count);
     TIMING_stop(tm_prj_frc_dir, flop_count);
     */
    
    // 収束判定　性能測定モードのときは収束判定を行わない
    if ( (C.Hide.PM_Test == OFF) && (convergence < ICp->get_eps()) ) break;
  } // end of iteration
  
  TIMING_stop(tm_poi_itr_sct, 0.0);
  // <<< Poisson Iteration section
  
  
  
  /// >>> NS Loop post section
  TIMING_start(tm_NS_loop_post_sct);
  
  // 同期
  if ( para_mng->IsParallel() ) {
    TIMING_start(tm_vectors_comm);
    dc_v->CommBndCell(guide);
    TIMING_stop(tm_vectors_comm, 2*comm_size*(REAL_TYPE)guide*3.0);
  }
  
  // 外部領域境界面での速度や流量を計算 > 外部流出境界条件の移流速度に利用
  TIMING_start(tm_domain_monitor);
  flop_count=0.0;
  DomainMonitor( BC.get_OBC_Ptr(), &C, flop_count);
  TIMING_stop(tm_domain_monitor, flop_count);
  
  // 流出境界のガイドセル値の更新
  TIMING_start(tm_VBC_update);
  flop_count = 0.0;
  BC.InnerVBC(v, bcv, tm, v00, flop_count, true);
  BC.OuterVBC(v, vc, bcv, tm, dt, &C, v00, flop_count);
  TIMING_stop(tm_VBC_update, flop_count);
  
  // 非同期にして隠す
  if (C.LES.Calc==ON) {
    TIMING_start(tm_LES_eddy);
    flop_count = 0.0;
    cbc_eddy_viscosity_(vt, sz, gc, dh, &C.Reynolds, &C.LES.Cs, v, (int*)bcv, range_Ut, range_Yp, v00);
    TIMING_stop(tm_LES_eddy, flop_count);
    
    if ( para_mng->IsParallel() ) {
      TIMING_start(tm_LES_eddy_comm);
      if( !dc_vt->CommBndCell(guide) ) Exit(0);
      TIMING_stop(tm_LES_eddy_comm, comm_size*(REAL_TYPE)guide);
    }
  }
  
  
  // ノルムの増加率が規定値をこえたら，終了
  if (convergence_prev != 0.0 ) {
    convergence_rate = convergence / convergence_prev;
  }
  else {
    convergence_rate = 1.0;
  }
  convergence_prev = convergence;
  
  TIMING_stop(tm_NS_loop_post_sct, 0.0);
  // >>> NS loop post section
}

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file SklSolverCBCLoop.C
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

#include "SklSolverCBC.h"

//@fn int SklSolverCBC::SklSolverLoop(const unsigned int step)
//@brief タイムステップループの処理
//@param step タイムステップ
int
SklSolverCBC::SklSolverLoop(const unsigned int step) {

  SklParaComponent* para_cmp = SklGetParaComponent();
  const SklParaManager* para_mng = para_cmp->GetParaManager();
  
  REAL_TYPE flop_count=0.0;      /// 浮動小数演算数
  REAL_TYPE avrms[6];            /// 平均値 [0]; vel, [1]; prs, [2]; temp, 変動値 [3]; vel, [4]; prs, [5]; temp
  REAL_TYPE vMax=0.0;            /// 最大速度成分
  REAL_TYPE *v=NULL;             /// 速度（セルセンタ）
  REAL_TYPE *p=NULL;             /// 圧力
  REAL_TYPE *t=NULL;             /// 温度
  REAL_TYPE *tp=NULL;            /// 全圧
  unsigned *bcd=NULL;            /// BCindex ID
  REAL_TYPE np_f = (REAL_TYPE)para_mng->GetNodeNum(pn.procGrp);
  
  // point Data
  if( !(v = dc_v->GetData()) )     Exit(0);
  if( !(p = dc_p->GetData()) )     Exit(0);
  if( !(bcd = dc_bcd->GetData()) ) Exit(0);
  
  if ( C.isHeatProblem() ) {
    if( !(t = dc_t->GetData()) )   Exit(0);
  }
  if (C.Mode.TP == ON ) {
    if( !(tp = dc_p0->GetData()) ) Exit(0);
  }
  
  // トリガーのリセット
  for (int i=0; i<Interval_Manager::tg_END; i++) {
    C.Interval[i].resetTrigger();
  }
  
  // Loop section
  TIMING_start(tm_loop_sct);
  
  // 時間進行
  SklIncrementTime();
  unsigned loop_step = SklGetTotalStep();
  double   loop_time = SklGetTotalTime();

  // 参照座標速度をv00に保持する
  copyV00fromRF(loop_time);

  // モニタークラスに参照速度を渡す
  if (C.Sampling.log == ON) MO.set_V00(v00);

  // 速度成分の最大値
  TIMING_start(tm_vmax);
  flop_count = 0.0;
  cbc_vmax_(&vMax, sz, gc, v00, v, &flop_count);
  TIMING_stop(tm_vmax, flop_count);

  if( para_mng->IsParallel() ){
    TIMING_start(tm_vmax_comm);
    REAL_TYPE vMax_tmp = vMax;
    if( !para_mng->Allreduce(&vMax_tmp, &vMax, 1, SKL_ARRAY_DTYPE_REAL, SKL_MAX, pn.procGrp) ) Exit(0);
    TIMING_stop( tm_vmax_comm, 2.0*np_f*(REAL_TYPE)sizeof(REAL_TYPE) ); // 双方向 x ノード数
  }


  // Flow
  if ( C.KindOfSolver != SOLID_CONDUCTION ) {
    TIMING_start(tm_flow_sct);
    
    switch (C.AlgorithmF) {
      case Control::Flow_FS_EE_EE:
      case Control::Flow_FS_AB2:
      case Control::Flow_FS_AB_CN:
        if (C.Mode.ShapeAprx == BINARY) {
          NS_FS_E_CBC();
        }
        else if (C.Mode.ShapeAprx == CUT_INFO) {
          NS_FS_E_CDS();
        }          
        break;
        
      case Control::Flow_FS_RK_CN:
        break;
        
      default:
        break;
    }
    TIMING_stop(tm_flow_sct, 0.0);
  }

  // Heat
  if ( C.isHeatProblem() ) {
    TIMING_start(tm_heat_sct);
    PS_E_CBC();
    TIMING_stop(tm_heat_sct, 0.0);
  }
  
  // Interface Equation
  if ( C.BasicEqs == INCMP_2PHASE ) {
    TIMING_start(tm_vof_sct);
    IF_TRP_VOF();
    TIMING_stop(tm_vof_sct, 0.0);
  }
  


  // >>> ステップループのユーティリティ
  TIMING_start(tm_loop_uty_sct);
  
  //  >>> ステップループのユーティリティ 1
  TIMING_start(tm_loop_uty_sct_1);

  // 時間平均値操作
  if ( (C.Mode.Average == ON) && SklUtil::IsStartAverage(this, C.Interval[Interval_Manager::tg_avstart].getIntervalTime()) ) {
    TIMING_start(tm_average_time);
    flop_count=0.0;
    Averaging_Time(flop_count);
    TIMING_stop(tm_average_time, flop_count);
  }

  // 空間平均値操作と変動量
  TIMING_start(tm_stat_space);
  flop_count=0.0;
  for (int i=0; i<6; i++) avrms[i] = 0.0;
  Variation_Space(avrms, flop_count);
  TIMING_stop(tm_stat_space, flop_count);

  if ( para_mng->IsParallel() ) {
    REAL_TYPE tmp[6];
    TIMING_start(tm_stat_space_comm);
    for (int n=0; n<6; n++) tmp[n] = avrms[n];
    para_mng->Allreduce(tmp, avrms, 6, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp); // 速度，圧力，温度の3変数 x (平均値+変動値)
    TIMING_stop(tm_stat_space_comm, 2.0*np_f*6.0*(REAL_TYPE)sizeof(REAL_TYPE) ); // 双方向 x ノード数 x 6変数
  }

  avrms[0] = sqrt(avrms[0]/(REAL_TYPE)G_Acell);  // 速度の変動量のRMS
  avrms[1] = sqrt(avrms[1]/(REAL_TYPE)G_Acell);  // 圧力の変動量のRMS
  avrms[2] = sqrt(avrms[2]/(REAL_TYPE)G_Acell);  // 温度の変動量のRMS
  avrms[3] = sqrt(avrms[3]/(REAL_TYPE)G_Acell);  // 速度の空間RMS
  avrms[4] = avrms[4]/(REAL_TYPE)G_Acell;        // 圧力の空間平均
  avrms[5] = avrms[5]/(REAL_TYPE)G_Acell;        // 温度の空間平均
  
  //  <<< ステップループのユーティリティ 1
  TIMING_stop(tm_loop_uty_sct_1, 0.0);
  

  
  // 1ステップ後のモニタ処理 -------------------------------
  
  //  >>> ステップループのユーティリティ 2
  TIMING_start(tm_loop_uty_sct_2);
  
  // Historyクラスのタイムスタンプを更新
  H->updateTimeStamp(loop_step, (REAL_TYPE)loop_time, vMax);
  
  // 基本履歴情報をコンソールに出力
  if ( C.Mode.Log_Base == ON) {
    if ( C.Interval[Interval_Manager::tg_console].isTriggered(loop_step, loop_time) ) {
      TIMING_start(tm_hstry_stdout);
      Hostonly_ H->printHistory(mp, avrms, IC, &C);
      TIMING_stop(tm_hstry_stdout, 0.0);
    }
  }
  
  // 瞬時値のデータ出力
  TIMING_start(tm_file_out);
  if ( C.Interval[Interval_Manager::tg_instant].isTriggered(loop_step, loop_time) ) {
    
    flop_count=0.0;
    
    // 通常
    FileOutput(flop_count);     
  }
  
  // 最終ステップ
  if ( m_currentStep == C.Interval[Interval_Manager::tg_compute].getIntervalStep() ) {
    
    // 指定間隔の出力がない場合のみ（重複を避ける）
    if ( !C.Interval[Interval_Manager::tg_instant].isTriggered(loop_step, loop_time) ) {
      if ( C.Hide.PM_Test != ON ) FileOutput(flop_count);
    }
  }
  
  TIMING_stop(tm_file_out, flop_count); 
  
  
  // 平均値のデータ出力 >　アルゴいまいち
  if (C.Mode.Average == ON) {
    // 開始時刻を過ぎているか
    bool j_flag = false;
    if ( C.Interval[Interval_Manager::tg_avstart].isStep() ) {
      if (loop_step >= C.Interval[Interval_Manager::tg_avstart].getIntervalStep()) j_flag=true;
    }
    else {
      if (loop_time >= C.Interval[Interval_Manager::tg_avstart].getIntervalTime()) j_flag=true;
    }

    if ( j_flag ) {
      // 初期化は1回だけ
      if ( !C.Interval[Interval_Manager::tg_average].initTrigger(loop_step, loop_time, (double)SklGetDeltaT(), Interval_Manager::tg_average) ) Exit(0);
      if ( C.Interval[Interval_Manager::tg_average].isTriggered(loop_step, loop_time) ) {

        TIMING_start(tm_file_out);
        
        flop_count=0.0;
        
        // 通常
        AverageOutput(flop_count);
        
        // 最終ステップ
        if ( m_currentStep == C.Interval[Interval_Manager::tg_compute].getIntervalStep() ) { 
          AverageOutput(flop_count);
        }
        
        TIMING_stop(tm_file_out, flop_count);
      }
    }
  }
  
  
  
  // 履歴のファイル出力
  if ( C.Interval[Interval_Manager::tg_history].isTriggered(loop_step, loop_time) ) {
    
    // 基本履歴情報
    if ( C.Mode.Log_Base == ON ) {
      TIMING_start(tm_hstry_base);
      Hostonly_ H->printHistory(fp_b, avrms, IC, &C);
      TIMING_stop(tm_hstry_base, 0.0);
    }
    
    // 壁面履歴情報
    if ( C.Mode.Log_Wall == ON ) {
      TIMING_start(tm_hstry_wall);
      Hostonly_ H->printHistoryWall(fp_w, range_Yp, range_Ut);
      TIMING_stop(tm_hstry_wall, 0.0);
    }
    
    // 流量収支履歴
    if ( C.Mode.Log_Base == ON ) {
      TIMING_start(tm_hstry_dmfx);
      Hostonly_ H->printHistoryDomfx(fp_d, &C);
      TIMING_stop(tm_hstry_dmfx, 0.0);
    }

    // 力の履歴
    if ( C.Hide.PM_Test != ON ) {
      REAL_TYPE force[3];
      
      TIMING_start(tm_cal_force);
      flop_count=0.0;
      // 性能測定モードのときには出力しない
      
      if ( C.isCDS() ) {
        cds_force_(force, sz, gc, p, (int*)bcd, dc_bid->GetData(), &id_of_solid, dh, &flop_count);
      }
      else {
        cbc_force_(force, sz, gc, p, (int*)bcd, dh, &flop_count);
      }
      
      TIMING_stop(tm_cal_force, 0.0);
      
      
      REAL_TYPE tmp_f[3];
      tmp_f[0] = force[0];
      tmp_f[1] = force[1];
      tmp_f[2] = force[2];
      if( para_cmp->IsParallel() ) {
        para_cmp->Allreduce(tmp_f, force, 3, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
      }
      
      TIMING_start(tm_hstry_force);
      Hostonly_ H->printHistoryForce(fp_f, &C, force);
      TIMING_stop(tm_hstry_force, 0.0);
    }
  }

  if (C.Mode.TP == ON ) {
    TIMING_start(tm_total_prs);
    flop_count=0.0;
    fb_totalp_ (tp, sz, gc, v, p, v00, &flop_count);
    TIMING_stop(tm_total_prs, flop_count);
  }

  // コンポーネント履歴
  if ( C.Sampling.log == ON ) {
    if ( C.Interval[Interval_Manager::tg_sampled].isTriggered(loop_step, loop_time) ) {
      TIMING_start(tm_compo_monitor);
      flop_count=0.0;
      MO.samplingInnerBoundary();
      TIMING_stop(tm_compo_monitor, flop_count);

      TIMING_start(tm_hstry_compo);
      Hostonly_ H->printHistoryCompo(fp_c, cmp, &C);
      TIMING_stop(tm_hstry_compo, 0.0);
    }
  }

  // サンプリング履歴
  if ( C.Sampling.log == ON ) {
    if ( C.Interval[Interval_Manager::tg_sampled].isTriggered(loop_step, loop_time) ) {
      TIMING_start(tm_sampling);
      MO.sampling();
      TIMING_stop(tm_sampling, 0.0);
      
      TIMING_start(tm_hstry_sampling);
      MO.print(loop_step, (REAL_TYPE)loop_time);
      TIMING_stop(tm_hstry_sampling, 0.0);
    }
  }
  
  TIMING_stop(tm_loop_uty_sct_2, 0.0);
  //  <<< ステップループのユーティリティ 2


  
  // 発散時の打ち切り
  if ( m_currentStep > 1 ) {
    if ( (convergence_rate > 100.0) ) {
      Hostonly_ {
        printf      ("\tForced termination : converegence rate >> 100.0\n");
        fprintf(fp_b,"\tForced termination : converegence rate >> 100.0\n");
      }
      return -1;
    }
  }
  
  // 計算時間がtimeにより指定されている場合の終了判断
  if ( !C.Interval[Interval_Manager::tg_compute].isStep() ) {
    if ( C.Interval[Interval_Manager::tg_compute].getIntervalTime() < loop_time ) {
      Hostonly_ {
        printf      ("\tFinish : Time = %e\n", loop_time);
        fprintf(fp_b,"\tFinish : Time = %e\n", loop_time);
      }
      return (0);
    }
  }
  
  TIMING_stop(tm_loop_uty_sct, 0.0);
  //  <<< ステップループのユーティリティ
  
  TIMING_stop(tm_loop_sct, 0.0);
  
  return 1;
}

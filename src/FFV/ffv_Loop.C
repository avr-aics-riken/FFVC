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
 * @file   ffv_Loop.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"


// タイムステップループの処理
int FFV::Loop(const unsigned step) 
{
  // 1 step elapse (sec)
  double step_start = cpm_Base::GetWTime();
  double step_end;
  
  double flop_count=0.0;   /// 浮動小数演算数
  double rms_Var[3];       /// 変動値（速度、圧力、温度）
  double avr_Var[3];       /// 平均値
  REAL_TYPE vMax=0.0;      /// 最大速度成分
  
  bool isNormal=true;      /// 発散チェックフラグ

  
  // Loop section
  TIMING_start(tm_loop_sct);
  
  // 時間進行
  CurrentTime += DT.get_DT(); // 戻り値はdouble
  CurrentStep++;
  
  
  // 参照座標速度をv00に保持する
  RF.setV00(CurrentTime);
  RF.copyV00(v00);
  
  // モニタークラスに参照速度を渡す
  if (C.SamplingMode == ON) MO.setV00(v00);
  
  // 速度成分の最大値
  TIMING_start(tm_vmax);
  flop_count = 0.0;
  find_vmax_(&vMax, size, &guide, v00, d_v, &flop_count);
  TIMING_stop(tm_vmax, flop_count);
  
  if ( numProc > 1 ) 
  {
    TIMING_start(tm_vmax_comm);
    REAL_TYPE vMax_tmp = vMax;
    if ( paraMngr->Allreduce(&vMax_tmp, &vMax, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    TIMING_stop( tm_vmax_comm, 2.0*numProc*sizeof(REAL_TYPE) ); // 双方向 x ノード数
  }
  
  
  // Flow
  if ( C.KindOfSolver != SOLID_CONDUCTION )
  {
    TIMING_start(tm_flow_sct);
    
    switch (C.AlgorithmF) 
    {
      case Flow_FS_EE_EE:
      case Flow_FS_AB2:
      case Flow_FS_AB_CN:
        if (C.Mode.ShapeAprx == BINARY)
        {
          NS_FS_E_Binary();
        }
        else if (C.Mode.ShapeAprx == CUT_INFO)
        {
          NS_FS_E_CDS();
        }          
        break;
        
      case Flow_FS_RK_CN:
        break;
        
      default:
        break;
    }
    TIMING_stop(tm_flow_sct, 0.0);
  }
  
  // Heat
  if ( C.isHeatProblem() ) 
  {
    TIMING_start(tm_heat_sct);
    PS_Binary();
    TIMING_stop(tm_heat_sct, 0.0);
  }
  
  
  
  // Interface Equation
  if ( C.BasicEqs == INCMP_2PHASE ) 
  {
    TIMING_start(tm_vof_sct);
    //IF_TRP_VOF();
    TIMING_stop(tm_vof_sct, 0.0);
  }
  
  
  
  // >>> ステップループのユーティリティ
  TIMING_start(tm_loop_uty_sct);
  
  //  >>> ステップループのユーティリティ 1
  TIMING_start(tm_loop_uty_sct_1);
  
  
  // 統計処理操作 >> 毎ステップ
  if ( (C.Mode.Statistic == ON) && C.Interval[Control::tg_statistic].isStarted(CurrentStep, CurrentTime))
  {
    TIMING_start(tm_average_time);
    flop_count=0.0;
    Averaging(flop_count);
    TIMING_stop(tm_average_time, flop_count);
    
    // 乱流統計量
    if ( C.LES.Calc == ON )
    {
      flop_count = 0.0;
      REAL_TYPE accum = (REAL_TYPE)CurrentStepStat;
      calc_rms_v_(d_rms, d_rmsmean, size, &guide, d_v, d_av, &accum, &flop_count);
    }
  }
  
  
  // 空間平均値操作と変動量
  TIMING_start(tm_stat_space);
  flop_count=0.0;
  for (int i=0; i<3; i++) 
  {
    rms_Var[i] = 0.0;
    avr_Var[i] = 0.0;
  }
  VariationSpace(rms_Var, avr_Var, flop_count);
  TIMING_stop(tm_stat_space, flop_count);
  
  
  if ( numProc > 1 )
  {
    /// var_Velocity=0,  > FB_Define.h
    /// var_Pressure,
    /// var_Temperature,
    double src[6], dst[6]; // Vel, Prs, Tempで3*2
    TIMING_start(tm_stat_space_comm);
    
    for (int n=0; n<3; n++) {
      src[n]   = rms_Var[n];
      src[n+3] = avr_Var[n];
    }
    
    if ( paraMngr->Allreduce(src, dst, 6, MPI_SUM) != CPM_SUCCESS) Exit(0); // 変数 x (平均値+変動値)
    
    for (int n=0; n<3; n++) {
      rms_Var[n] = dst[n];
      avr_Var[n] = dst[n+3];
    }
    
    TIMING_stop(tm_stat_space_comm, 2.0*numProc*6.0*2.0*sizeof(double) ); // 双方向 x ノード数 x 変数
  }

  rms_Var[var_Velocity] = sqrt(rms_Var[var_Velocity]); // 速度の変動量の空間総和
  rms_Var[var_Pressure] = sqrt(rms_Var[var_Pressure]); // 圧力の変動量の空間総和
  
  //avr_Var[var_Velocity] /= (double)G_Acell;  // 速度の空間平均は計算しない
  avr_Var[var_Pressure] /= (double)G_Acell;  // 圧力の空間平均
  
  
  
  if ( C.isHeatProblem() ) 
  {
    rms_Var[var_Temperature] = sqrt(rms_Var[var_Temperature]); // 温度の変動量の空間総和
    avr_Var[var_Temperature] /= (double)G_Acell;   // 温度の空間平均
  }
  
  // 発散チェック
  if ( isnan(rms_Var[var_Velocity])
    || isnan(rms_Var[var_Pressure])
    || isnan(DivC.divergence)
      )
  {
    isNormal = false;
  }
  
  if ( C.isHeatProblem() )
  {
    if ( isnan(rms_Var[var_Temperature]) )
    {
      isNormal = false;
    }
  }
  
  
  
  //  <<< ステップループのユーティリティ 1
  TIMING_stop(tm_loop_uty_sct_1, 0.0);
  
  
  
  // 1ステップ後のモニタ処理 -------------------------------
  
  
  
  //  >>> ステップループのユーティリティ 2
  TIMING_start(tm_loop_uty_sct_2);
  
  // Historyクラスのタイムスタンプを更新
  H->updateTimeStamp(CurrentStep, (REAL_TYPE)CurrentTime, vMax);
  
  
  
  // 瞬時値のデータ出力

  if ( C.Hide.PM_Test == OFF )
  {
    // 通常
    if ( C.Interval[Control::tg_basic].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start(tm_file_out);
      flop_count=0.0;
      F->OutputBasicVariables(CurrentStep, CurrentTime, flop_count);
      TIMING_stop(tm_file_out, flop_count);
      
      if ( F->isVtk() )
      {
        output_vtk_((int*)&CurrentStep, G_origin, G_division, G_size, &myRank, size, &pitch[0], &guide, d_v, d_p);
      }
    }

    
    // 最終ステップ
    if ( C.Interval[Control::tg_compute].isLast(CurrentStep, CurrentTime) )
    {
      // 指定間隔の出力がない場合のみ（重複を避ける）
      if ( !C.Interval[Control::tg_basic].isTriggered(CurrentStep, CurrentTime) )
      {
        TIMING_start(tm_file_out);
        flop_count=0.0;
        F->OutputBasicVariables(CurrentStep, CurrentTime, flop_count);
        TIMING_stop(tm_file_out, flop_count);
      }
    }
    
  }

  
  // 統計値のデータ出力 
  if (C.Mode.Statistic == ON) 
  {
    
    // 開始時刻を過ぎているか
    if ( C.Interval[Control::tg_statistic].isStarted(CurrentStep, CurrentTime) )
    {
      // 通常
      if ( C.Interval[Control::tg_statistic].isTriggered(CurrentStep, CurrentTime) ) 
      {
        TIMING_start(tm_file_out);
        flop_count=0.0;
        F->OutputStatisticalVarables(CurrentStep, CurrentTime, CurrentStepStat, CurrentTimeStat, flop_count);
        TIMING_stop(tm_file_out, flop_count);
        
        // special
        output_mean_((int*)&CurrentStep, G_origin, G_region, G_division, G_size, &myRank, size, &pitch[0], &guide, d_av, d_rms, d_rmsmean);
        
      }
      
      // 最終ステップ
      if ( C.Interval[Control::tg_compute].isLast(CurrentStep, CurrentTime) )
      {
        // 指定間隔の出力がない場合のみ（重複を避ける）
        if ( !C.Interval[Control::tg_statistic].isTriggered(CurrentStep, CurrentTime) )
        {
          TIMING_start(tm_file_out);
          flop_count=0.0;
          F->OutputStatisticalVarables(CurrentStep, CurrentTime, CurrentStepStat, CurrentTimeStat, flop_count);
          TIMING_stop(tm_file_out, flop_count);
        }
      }
    }
  }
  
  if (C.varState[var_TotalP] == ON )
  {
    TIMING_start(tm_total_prs);
    flop_count=0.0;
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop_count);
    TIMING_stop(tm_total_prs, flop_count);
  }
  

  
  // サンプリング履歴
  if ( (C.SamplingMode == ON) && C.Interval[Control::tg_sampled].isTriggered(CurrentStep, CurrentTime) )
  {
    TIMING_start(tm_sampling);
    MO.sampling();
    TIMING_stop(tm_sampling, 0.0);
    
    TIMING_start(tm_hstry_sampling);
    MO.print(CurrentStep, (REAL_TYPE)CurrentTime);
    TIMING_stop(tm_hstry_sampling, 0.0);
  }

  

  // 1 step elapse
  step_end = cpm_Base::GetWTime() - step_start;
  
  
  // dynamic_cast<IterationCtl*>(LS)で生じるコンパイラエラー（gnu）対策
  double container[3*ic_END];
  
  for (int i=0; i<ic_END; i++)
  {
    container[3*i+0] = (double)LS[ic_prs1].getLoopCount();
    container[3*i+1] = LS[ic_prs1].getResidual();
    container[3*i+1] = LS[ic_prs1].getError();
  }
  
  
  // 基本履歴情報をコンソールに出力
  if ( C.Mode.Log_Base == ON)
  {
    if ( C.Interval[Control::tg_console].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start(tm_hstry_stdout);
      Hostonly_
      {
        H->printHistory(stdout, rms_Var, avr_Var, container, &C, &DivC, step_end, true);
        
        if ( C.Mode.CCNV == ON )
        {
          H->printCCNV(rms_Var, avr_Var, container, &C, DivC.divergence, step_end);
        }
      }
      TIMING_stop(tm_hstry_stdout, 0.0);
    }
  }
  
  
  // 履歴のファイル出力
  if ( C.Interval[Control::tg_history].isTriggered(CurrentStep, CurrentTime) ) 
  {
    
    if ( C.Mode.Log_Base == ON ) 
    {
      TIMING_start(tm_hstry_base);
      
      // 基本履歴情報
      Hostonly_ H->printHistory(fp_b, rms_Var, avr_Var, container, &C, &DivC, step_end, true);
      
      // コンポーネント
      if ( C.EnsCompo.monitor )
      {
        Hostonly_ H->printHistoryCompo(fp_c, cmp, &C, deltaT);
      }
      
      // 物体の力の計算
      if ( C.EnsCompo.obstacle )
      {
        TIMING_start(tm_cal_force);
        flop_count=0.0;
        
        calcForce(flop_count);
        
        gatherForce(buffer_force);
        
        Hostonly_ H->printHistoryForce(fp_f, cmp, &C, global_force);
        TIMING_stop(tm_cal_force, 0.0);
      }
      
      // 流量収支履歴
      Hostonly_ H->printHistoryDomfx(fp_d, &C, deltaT);
      
      TIMING_stop(tm_hstry_base, 0.0);
    }
    
    // 壁面履歴情報
    if ( C.Mode.Log_Wall == ON ) 
    {
      TIMING_start(tm_hstry_wall);
      //Hostonly_ H->printHistoryWall(fp_w, range_Yp, range_Ut);
      TIMING_stop(tm_hstry_wall, 0.0);
    }
    
  }
  
  TIMING_stop(tm_loop_uty_sct_2, 0.0);
  //  <<< ステップループのユーティリティ 2
  
  
  
  // 発散時の打ち切り
  if ( !isNormal )
  {
    Hostonly_ {
      printf      ("\tForced termination : floating point exception\n");
      if ( C.Mode.Log_Base == ON) fprintf(fp_b,"\tForced termination : floating point exception\n");
    }
    return -1;
  }
  
  
  TIMING_stop(tm_loop_uty_sct, 0.0);
  //  <<< ステップループのユーティリティ
  
  TIMING_stop(tm_loop_sct, 0.0);
  
  
  
  // 終了判断
  if ( C.Interval[Control::tg_compute].isLast(CurrentStep, CurrentTime) )
  {
    Hostonly_
    {
      printf      ("\tFinish : Time = %e\n", CurrentTime);
      if ( C.Mode.Log_Base == ON) fprintf(fp_b,"\tFinish : Time = %e\n", CurrentTime);
    }
    return 0;
  }
  
  return 1;
}

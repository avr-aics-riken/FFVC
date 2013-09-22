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
 * @file   ffv_Loop.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"


// タイムステップループの処理
int FFV::Loop(const unsigned step) 
{
  // 1 step elapse
  step_start = cpm_Base::GetWTime();
  double step_end;
  
  double flop_count=0.0;   /// 浮動小数演算数
  double avr_Var[3];       /// 平均値（速度、圧力、温度）
  double rms_Var[3];       /// 変動値
  REAL_TYPE vMax=0.0;      /// 最大速度成分


  // トリガーのリセット
  for (int i=0; i<Interval_Manager::tg_END; i++) 
  {
    C.Interval[i].resetTrigger();
  }
  
  // Loop section
  TIMING_start(tm_loop_sct);
  
  // 時間進行
  CurrentTime += DT.get_DT(); // 戻り値はdouble
  CurrentStep++;
  
  
  // 参照座標速度をv00に保持する
  copyV00fromRF(CurrentTime);
  
  // モニタークラスに参照速度を渡す
  if (C.Sampling.log == ON) MO.set_V00(v00);
  
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
  
  // 時間平均値操作
  if ( (C.Mode.Average == ON) && C.Interval[Interval_Manager::tg_average].isStarted(CurrentStep, CurrentTime))
  {
    TIMING_start(tm_average_time);
    flop_count=0.0;
    Averaging(flop_count);
    TIMING_stop(tm_average_time, flop_count);
  }
  
  // 空間平均値操作と変動量
  TIMING_start(tm_stat_space);
  flop_count=0.0;
  for (int i=0; i<3; i++) 
  {
    avr_Var[i] = 0.0;
    rms_Var[i] = 0.0;
  }
  VariationSpace(avr_Var, rms_Var, flop_count);
  TIMING_stop(tm_stat_space, flop_count);

  
  
  if ( numProc > 1 ) 
  {
    /// var_Velocity=0,  > FB_Define.h
    /// var_Pressure,
    /// var_Temperature,
    double src[6], dst[6]; // Vel, Prs, Tempで3*2
    TIMING_start(tm_stat_space_comm);
    
    for (int n=0; n<3; n++) {
      src[n]   = avr_Var[n];
      src[n+3] = rms_Var[n];
    }
    
    if ( paraMngr->Allreduce(src, dst, 6, MPI_SUM) != CPM_SUCCESS) Exit(0); // 変数 x (平均値+変動値)
    
    for (int n=0; n<3; n++) {
      avr_Var[n] = dst[n];
      rms_Var[n] = dst[n+3];
    }
    
    TIMING_stop(tm_stat_space_comm, 2.0*numProc*6.0*2.0*sizeof(double) ); // 双方向 x ノード数 x 変数
  }

  avr_Var[var_Velocity] /= (double)G_Acell;  // 速度の空間平均
  avr_Var[var_Pressure] /= (double)G_Acell;  // 圧力の空間平均
  
  rms_Var[var_Velocity] /= (double)G_Acell;  // 速度の変動量
  rms_Var[var_Pressure] /= (double)G_Acell;  // 圧力の変動量
  rms_Var[var_Velocity] = sqrt(rms_Var[var_Velocity]);
  rms_Var[var_Pressure] = sqrt(rms_Var[var_Pressure]);
  
  if ( C.isHeatProblem() ) 
  {
    avr_Var[var_Temperature] /= (double)G_Acell;   // 温度の空間平均
    rms_Var[var_Temperature] /= (double)G_Acell;   // 温度の変動量
    rms_Var[var_Temperature] = sqrt(rms_Var[var_Temperature]);
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
    if ( C.Interval[Interval_Manager::tg_basic].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start(tm_file_out);
      flop_count=0.0;
      OutputBasicVariables(flop_count);
      TIMING_stop(tm_file_out, flop_count);
    }
    
    if ( C.Interval[Interval_Manager::tg_derived].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start(tm_file_out);
      flop_count=0.0;
      OutputDerivedVariables(flop_count);
      TIMING_stop(tm_file_out, flop_count);
    }
    
    // 最終ステップ
    if ( CurrentStep == C.Interval[Interval_Manager::tg_compute].getIntervalStep() )
    {
      // 指定間隔の出力がない場合のみ（重複を避ける）
      if ( !C.Interval[Interval_Manager::tg_basic].isTriggered(CurrentStep, CurrentTime) )
      {
        TIMING_start(tm_file_out);
        flop_count=0.0;
        OutputBasicVariables(flop_count);
        TIMING_stop(tm_file_out, flop_count);
      }
      
      if ( !C.Interval[Interval_Manager::tg_derived].isTriggered(CurrentStep, CurrentTime) )
      {
        TIMING_start(tm_file_out);
        flop_count=0.0;
        OutputDerivedVariables(flop_count);
        TIMING_stop(tm_file_out, flop_count);
      }
    }
    
  }


  
  //  PLOT3D output
  if (C.FIO.Format == plt3d_fmt)
  {
    // 通常
    if ( C.Interval[Interval_Manager::tg_basic].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start(tm_file_out);
      flop_count=0.0;
      PLT3D.post(CurrentStep, CurrentTime, v00, origin, pitch, dfi_mng_Plot3D, flop_count);
      TIMING_stop(tm_file_out, flop_count);
    }
     
    
    // 最終ステップ
    if ( CurrentStep == C.Interval[Interval_Manager::tg_compute].getIntervalStep() )
    {
      // 指定間隔の出力がない場合のみ（重複を避ける）
      if ( !C.Interval[Interval_Manager::tg_basic].isTriggered(CurrentStep, CurrentTime) )
      {
        if ( C.Hide.PM_Test != ON )
        {
          TIMING_start(tm_file_out);
          flop_count=0.0;
          PLT3D.post(CurrentStep, CurrentTime, v00, origin, pitch, dfi_mng_Plot3D, flop_count);
          TIMING_stop(tm_file_out, flop_count);
        }
      }
    }
  }
  
  // 平均値のデータ出力 
  if (C.Mode.Average == ON) 
  {
    
    // 開始時刻を過ぎているか
    if ( C.Interval[Interval_Manager::tg_average].isStarted(CurrentStep, CurrentTime) )
    {
      // 通常
      if ( C.Interval[Interval_Manager::tg_average].isTriggered(CurrentStep, CurrentTime) ) 
      {
        TIMING_start(tm_file_out);
        flop_count=0.0;
        OutputAveragedVarables(flop_count);
        TIMING_stop(tm_file_out, flop_count);
      }
      
      // 最終ステップ
      if ( CurrentStep == C.Interval[Interval_Manager::tg_compute].getIntervalStep() )
      {
        // 指定間隔の出力がない場合のみ（重複を避ける）
        if ( !C.Interval[Interval_Manager::tg_average].isTriggered(CurrentStep, CurrentTime) )
        {
          TIMING_start(tm_file_out);
          flop_count=0.0;
          OutputAveragedVarables(flop_count);
          TIMING_stop(tm_file_out, flop_count);
        }
      }
    }
  }
  
  
  // 1 step elapse
  step_end = cpm_Base::GetWTime() - step_start;
  
  
  // 基本履歴情報をコンソールに出力
  if ( C.Mode.Log_Base == ON)
  {
    if ( C.Interval[Interval_Manager::tg_console].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start(tm_hstry_stdout);
      Hostonly_
      {
        H->printHistory(stdout, avr_Var, rms_Var, IC, &C, step_end, true);
        if ( C.Mode.CCNV == ON )
        {
          H->printCCNV(avr_Var, rms_Var, IC, &C, step_end);
        }
      }
      TIMING_stop(tm_hstry_stdout, 0.0);
    }
  }
  
  
  // 履歴のファイル出力
  if ( C.Interval[Interval_Manager::tg_history].isTriggered(CurrentStep, CurrentTime) ) 
  {
    
    // 基本履歴情報
    if ( C.Mode.Log_Base == ON ) 
    {
      TIMING_start(tm_hstry_base);
      Hostonly_ H->printHistory(fp_b, avr_Var, rms_Var, IC, &C, step_end, true);
      TIMING_stop(tm_hstry_base, 0.0);
    }
    
    // 壁面履歴情報
    if ( C.Mode.Log_Wall == ON ) 
    {
      TIMING_start(tm_hstry_wall);
      //Hostonly_ H->printHistoryWall(fp_w, range_Yp, range_Ut);
      TIMING_stop(tm_hstry_wall, 0.0);
    }
    
    // 流量収支履歴
    if ( C.Mode.Log_Base == ON ) 
    {
      TIMING_start(tm_hstry_dmfx);
      Hostonly_ H->printHistoryDomfx(fp_d, &C, deltaT);
      TIMING_stop(tm_hstry_dmfx, 0.0);
    }
    
    // 力の履歴
    if ( C.Hide.PM_Test != ON ) 
    {
      REAL_TYPE frc[3];
      
      TIMING_start(tm_cal_force);
      flop_count=0.0;
      // 性能測定モードのときには出力しない
      
      if ( C.isBinary() )
      {
        force_(frc, size, &guide, d_p, d_bcd, &deltaX, &flop_count);
      }
      else 
      {
        //cds_force_(frc, size, &guide, d_p, d_bcd, d_bid, &id_of_solid, &deltaX, &flop_count);
      }
      
      TIMING_stop(tm_cal_force, 0.0);
      
      
      REAL_TYPE tmp_f[3];
      tmp_f[0] = frc[0];
      tmp_f[1] = frc[1];
      tmp_f[2] = frc[2];
      if ( numProc > 1 ) 
      {
        if ( paraMngr->Allreduce(tmp_f, frc, 3, MPI_SUM) != CPM_SUCCESS) Exit(0);
      }
      
      TIMING_start(tm_hstry_force);
      Hostonly_ H->printHistoryForce(fp_f, frc);
      TIMING_stop(tm_hstry_force, 0.0);
    }
  }
  
  if (C.Mode.TP == ON ) 
  {
    TIMING_start(tm_total_prs);
    flop_count=0.0;
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop_count);
    TIMING_stop(tm_total_prs, flop_count);
  }
  
  // セルモニターとコンポーネントの履歴
  if ( C.Sampling.log == ON ) 
  {
    if ( C.Interval[Interval_Manager::tg_history].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start(tm_compo_monitor);
      flop_count=0.0;
      MO.samplingInnerBoundary();
      TIMING_stop(tm_compo_monitor, flop_count);
      
      TIMING_start(tm_hstry_compo);
      Hostonly_ H->printHistoryCompo(fp_c, cmp, &C, deltaT);
      TIMING_stop(tm_hstry_compo, 0.0);
    }
  }
  
  // サンプリング履歴
  if ( C.Sampling.log == ON ) 
  {
    if ( C.Interval[Interval_Manager::tg_sampled].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start(tm_sampling);
      MO.sampling();
      TIMING_stop(tm_sampling, 0.0);
      
      TIMING_start(tm_hstry_sampling);
      MO.print(CurrentStep, (REAL_TYPE)CurrentTime);
      TIMING_stop(tm_hstry_sampling, 0.0);
    }
  }
  
  TIMING_stop(tm_loop_uty_sct_2, 0.0);
  //  <<< ステップループのユーティリティ 2
  
  
  
  // 発散時の打ち切り
  if ( CurrentStep > 1 ) 
  {
    
    switch ( C.KindOfSolver )
    {
      case FLOW_ONLY:
        if ( (CM_F.rate > 100.0) )
        {
          Hostonly_ {
            printf      ("\tForced termination : converegence rate >> 100.0\n");
            fprintf(fp_b,"\tForced termination : converegence rate >> 100.0\n");
          }
          return -1;
        }
        break;
        
      case THERMAL_FLOW:
      case THERMAL_FLOW_NATURAL:
      case CONJUGATE_HEAT_TRANSFER:
        if ( (CM_F.rate > 100.0) || (CM_H.rate > 100.0) )
        {
          Hostonly_ {
            printf      ("\tForced termination : converegence rate >> 100.0\n");
            fprintf(fp_b,"\tForced termination : converegence rate >> 100.0\n");
          }
          return -1;
        }
        break;
        
      case SOLID_CONDUCTION:
        if ( (CM_H.rate > 100.0) )
        {
          Hostonly_ {
            printf      ("\tForced termination : converegence rate >> 100.0\n");
            fprintf(fp_b,"\tForced termination : converegence rate >> 100.0\n");
          }
          return -1;
        }
        break;
    }
    
  }
  
  // 計算時間がtimeにより指定されている場合の終了判断
  if ( !C.Interval[Interval_Manager::tg_compute].isStep() )
  {
    if ( C.Interval[Interval_Manager::tg_compute].getIntervalTime() < CurrentTime )
    {
      Hostonly_ 
      {
        printf      ("\tFinish : Time = %e\n", CurrentTime);
        fprintf(fp_b,"\tFinish : Time = %e\n", CurrentTime);
      }
      return 0;
    }
  }
  
  TIMING_stop(tm_loop_uty_sct, 0.0);
  //  <<< ステップループのユーティリティ
  
  TIMING_stop(tm_loop_sct, 0.0);
  
  return 1;
}

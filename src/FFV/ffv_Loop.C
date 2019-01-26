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
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   ffv_Loop.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"
#include "FileSystemUtil.h"


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

  int isNormal = 0;        /// 発散チェックフラグ


  // Loop section
  TIMING_start("Time_Step_Loop_Section");

  // 時間進行
  CurrentTime += DT.get_DT(); // 戻り値はdouble
  CurrentStep++;


  // 参照座標速度をv00に保持する
  RF.setV00(CurrentTime);
  RF.copyV00(v00);

  // モニタークラスに参照速度を渡す
  if (C.SamplingMode == ON) MO.setV00(v00);

  // 速度成分の最大値
  TIMING_start("Search_Vmax");
  flop_count = 0.0;
  find_vmax_(&vMax, size, &guide, v00, d_v, &flop_count);
  TIMING_stop("Search_Vmax", flop_count);

  if ( numProc > 1 )
  {
    TIMING_start("All_Reduce");
    REAL_TYPE vMax_tmp = vMax;
    if ( paraMngr->Allreduce(&vMax_tmp, &vMax, 1, MPI_MAX, procGrp) != CPM_SUCCESS ) Exit(0);
    TIMING_stop( "All_Reduce", 2.0*numProc*sizeof(REAL_TYPE) ); // 双方向 x ノード数
  }


  // Flow
  if ( C.KindOfSolver != SOLID_CONDUCTION )
  {
    TIMING_start("Flow_Section");

    if ( C.Hide.DryRun == OFF )
    {
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
    }
    else
    {
      dryrunBC();
    }

    TIMING_stop("Flow_Section", 0.0);
  }

  // Heat
  if ( C.isHeatProblem() )
  {
    TIMING_start("Heat_Section");
    PS_Binary();
    TIMING_stop("Heat_Section", 0.0);
  }



  // Interface Equation
  if ( C.BasicEqs == INCMP_2PHASE )
  {
    TIMING_start("VOF_Section");
    //IF_TRP_VOF();
    TIMING_stop("VOF_Section", 0.0);
  }



  // >>> ステップループのユーティリティ
  TIMING_start("Loop_Utility_Section");



  // 統計処理操作 >> 毎ステップ
  if ( (C.Mode.Statistic == ON) && C.Interval[Control::tg_statistic].isStarted(CurrentStep, CurrentTime))
  {
    TIMING_start("Averaging");
    flop_count=0.0;
    Averaging(flop_count);
    TIMING_stop("Averaging", flop_count);
    REAL_TYPE accum = (REAL_TYPE)CurrentStepStat;

    // 乱流統計量
    TIMING_start("Turbulence Statistic");

    if ( C.Mode.StatVelocity == ON )
    {
      flop_count = 0.0;
      calc_rms_v_(d_rms_v, d_rms_mean_v, size, &guide, d_v, d_av, &accum, &flop_count);
    }

    if ( C.Mode.StatPressure == ON )
    {
      flop_count = 0.0;
      calc_rms_s_(d_rms_p, d_rms_mean_p, size, &guide, d_p, d_ap, &accum, &flop_count);
    }

    if ( C.Mode.StatTemperature == ON )
    {
      flop_count = 0.0;
      U.convArrayIE2Tmp(d_ws, size, guide, d_ie, d_bcd, mat_tbl, C.BaseTemp, C.DiffTemp, C.Unit.File, flop_count);
      U.convArrayIE2Tmp(d_ie0, size, guide, d_ae, d_bcd, mat_tbl, C.BaseTemp, C.DiffTemp, C.Unit.File, flop_count); // d_ie0はワーク
      calc_rms_s_(d_rms_t, d_rms_mean_t, size, &guide, d_ws, d_ie0, &accum, &flop_count);
    }


    if ( C.Mode.ReynoldsStress == ON )
    {
      flop_count = 0.0;

      // レイノルズ応力テンソル (瞬時値, 時間平均値)
      calc_reynolds_stress_(d_R, d_aR, size, &guide, d_v, d_av, &accum, &flop_count);
      /*
      // (1) 生成項 (時間平均値)
      calc_production_rate_(d_aP, size, pitch, &guide, d_av, d_R, d_bid, &accum, &flop_count);

      // (2) 散逸項 (時間平均値)
      calc_dissipation_rate_(d_aE, size, pitch, &guide, &C.RefKviscosity, d_v, d_av, d_bid, &accum, &flop_count);

      // (3) 乱流拡散項 (時間平均値)
      calc_turb_transport_rate_(d_aT, size, pitch, &guide, d_v, d_av, d_R, d_bid, &accum, &flop_count);

      // (4) 速度圧力勾配相関項 (時間平均値)
      calc_vel_pregrad_term_(d_aPI, size, pitch, &guide, d_v, d_av, d_p, d_ap, d_bcp, &accum, &flop_count);

      // (5) チャネル乱流統計量 (時間・主流方向・スパン方向平均値) の出力
      if ( C.Mode.ChannelOutputMean == ON )
      {
        if ( (CurrentStepStat % C.Mode.ChannelOutputIter == 0) || (CurrentStep == 1) )
        {
          F->OutputMean(d_av, d_rms_mean_v, d_aR, d_aP, d_aE, d_aT, d_aPI, myRank, size, CurrentStepStat, pitch, &guide, flop_count);
        }
      }*/
    }

    TIMING_stop("Turbulence Statistic", flop_count);

  }


  // 空間平均値操作と変動量
  TIMING_start("Variation_Space");
  flop_count=0.0;
  for (int i=0; i<3; i++)
  {
    rms_Var[i] = 0.0;
    avr_Var[i] = 0.0;
  }
  VariationSpace(rms_Var, avr_Var, flop_count);
  TIMING_stop("Variation_Space", flop_count);


  if ( numProc > 1 )
  {
    /// var_Velocity=0,  > FB_Define.h
    /// var_Pressure,
    /// var_Temperature,
    double src[6], dst[6]; // Vel, Prs, Tempで3*2
    TIMING_start("A_R_variation_space");

    for (int n=0; n<3; n++) {
      src[n]   = rms_Var[n];
      src[n+3] = avr_Var[n];
    }

    if ( paraMngr->Allreduce(src, dst, 6, MPI_SUM, procGrp) != CPM_SUCCESS) Exit(0); // 変数 x (平均値+変動値)

    for (int n=0; n<3; n++) {
      rms_Var[n] = dst[n];
      avr_Var[n] = dst[n+3];
    }

    TIMING_stop("A_R_variation_space", 2.0*numProc*6.0*2.0*sizeof(double) ); // 双方向 x ノード数 x 変数
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
	if ( std::isnan(rms_Var[var_Velocity])
			|| std::isnan(rms_Var[var_Pressure])
			|| std::isnan(DivC.divergence)
      )
  {
    isNormal = 1;
  }

  if ( C.isHeatProblem() )
  {
		if ( std::isnan(rms_Var[var_Temperature]) )
    {
      isNormal = 1;
    }
  }



  // 1ステップ後のモニタ処理 -------------------------------



  // 瞬時値のデータ出力

  if ( C.Hide.PM_Test == OFF )
  {
    // 通常
    if ( C.Interval[Control::tg_basic].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start("File_Output");
      flop_count=0.0;
      F->OutputBasicVariables(CurrentStep, CurrentTime, flop_count);
      TIMING_stop("File_Output", flop_count);

      if ( F->isVtk() )
      {
        // Uzawa
        //int cs = CurrentStep;
        //output_vtk_(&cs, G_origin, G_division, G_size, &myRank, size, pitch, &guide, d_v, d_p);
      }
    }


    // 最終ステップ
    if ( C.Interval[Control::tg_compute].isLast(CurrentStep, CurrentTime) )
    {
      // 指定間隔の出力がない場合のみ（重複を避ける）
      if ( !C.Interval[Control::tg_basic].isTriggered(CurrentStep, CurrentTime) )
      {
        TIMING_start("File_Output");
        flop_count=0.0;
        F->OutputBasicVariables(CurrentStep, CurrentTime, flop_count);
        TIMING_stop("File_Output", flop_count);
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
        TIMING_start("File_Output");
        flop_count=0.0;
        F->OutputStatisticalVarables(CurrentStep, CurrentTime, CurrentStepStat, CurrentTimeStat, flop_count);
        TIMING_stop("File_Output", flop_count);

        // レイノルズ応力
        if (C.Mode.ReynoldsStress == ON)
        {
          char tmp_fname[100];
          int stp  = (int)CurrentStep;
          REAL_TYPE tm = (REAL_TYPE)CurrentTime;
          int stp2 = (int)CurrentStepStat;
          REAL_TYPE tm2 = (REAL_TYPE)CurrentTimeStat;
          int rtype = (8 == sizeof(REAL_TYPE)) ? 2:1;
          REAL_TYPE m_org[3]={(REAL_TYPE)origin[0], (REAL_TYPE)origin[1], (REAL_TYPE)origin[2]};
          REAL_TYPE m_pch[3]={(REAL_TYPE)pitch[0], (REAL_TYPE)pitch[1], (REAL_TYPE)pitch[0]}; // cartesian

          if ( !BVX_IO::CreateDirectory("Reynolds_Stress") ) Exit(0);

          // 瞬時値
          int c = 1;
          sprintf( tmp_fname, "./Reynolds_Stress/tau_uu_%08d_id%06d.sph", stp, myRank );
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_R, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_uv_%08d_id%06d.sph", stp, myRank );
          c= 2;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_R, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_uw_%08d_id%06d.sph", stp, myRank );
          c= 3;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_R, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_vv_%08d_id%06d.sph", stp, myRank );
          c= 4;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_R, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_vw_%08d_id%06d.sph", stp, myRank );
          c= 5;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_R, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_ww_%08d_id%06d.sph", stp, myRank );
          c= 6;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_R, &c, tmp_fname, &rtype);

          // 平均値
          c = 1;
          sprintf( tmp_fname, "./Reynolds_Stress/tau_uu_a_%08d_id%06d.sph", stp, myRank );
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_aR, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_uv_a_%08d_id%06d.sph", stp, myRank );
          c= 2;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_aR, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_uw_a_%08d_id%06d.sph", stp, myRank );
          c= 3;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_aR, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_vv_a_%08d_id%06d.sph", stp, myRank );
          c= 4;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_aR, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_vw_a_%08d_id%06d.sph", stp, myRank );
          c= 5;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_aR, &c, tmp_fname, &rtype);

          sprintf( tmp_fname, "./Reynolds_Stress/tau_ww_a_%08d_id%06d.sph", stp, myRank );
          c= 6;
          sph_write_tensor_(size, &guide, m_org, m_pch, &stp2, &tm2, d_aR, &c, tmp_fname, &rtype);
        }
      }

      // 最終ステップ
      if ( C.Interval[Control::tg_compute].isLast(CurrentStep, CurrentTime) )
      {
        // 指定間隔の出力がない場合のみ（重複を避ける）
        if ( !C.Interval[Control::tg_statistic].isTriggered(CurrentStep, CurrentTime) )
        {
          TIMING_start("File_Output");
          flop_count=0.0;
          F->OutputStatisticalVarables(CurrentStep, CurrentTime, CurrentStepStat, CurrentTimeStat, flop_count);
          TIMING_stop("File_Output", flop_count);
        }
      }

    }
  }



  if (C.varState[var_TotalP] == ON )
  {
    TIMING_start("Total_Pressure");
    flop_count=0.0;
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop_count);
    TIMING_stop("Total_Pressure", flop_count);
  }



  // サンプリング履歴
  if ( (C.SamplingMode == ON) && C.Interval[Control::tg_sampled].isTriggered(CurrentStep, CurrentTime) )
  {
    TIMING_start("Sampling");
    MO.sampling();
    MO.print(CurrentStep, (REAL_TYPE)CurrentTime);
    TIMING_stop("Sampling", 0.0);
  }


  //#################################################
  // 粒子追跡
  TIMING_start("ParticleTracking_Section");
  TR->tracking(CurrentStep, CurrentTime);
  TIMING_stop("ParticleTracking_Section", 0.0);
  //#################################################

  // 1 step elapse
  step_end = cpm_Base::GetWTime() - step_start;


  // dynamic_cast<IterationCtl*>(LS)で生じるコンパイラエラー（gnu）対策
  double container[3*ic_END];

  for (int i=0; i<ic_END; i++)
  {
    container[3*i+0] = (double)LS[ic_prs1].getLoopCount();
    container[3*i+1] = LS[ic_prs1].getResidual();
    container[3*i+2] = LS[ic_prs1].getError();
  }
    
    
    // Historyクラスのタイムスタンプを更新
    H->updateTimeStamp(CurrentStep, (REAL_TYPE)CurrentTime, vMax);


  // 基本履歴情報をコンソールに出力
  if ( C.Mode.Log_Base == ON)
  {
    if ( C.Interval[Control::tg_console].isTriggered(CurrentStep, CurrentTime) )
    {
      TIMING_start("History_out");
      Hostonly_
      {
        H->printHistory(stdout, rms_Var, avr_Var, container, &C, &DivC, step_end, true);

        if ( C.Mode.CCNV == ON )
        {
          H->printCCNV(rms_Var, avr_Var, container, &C, DivC.divergence, step_end);
        }
      }
      TIMING_stop("History_out", 0.0);
    }
  }


  // 履歴のファイル出力
  if ( C.Interval[Control::tg_history].isTriggered(CurrentStep, CurrentTime) )
  {

    if ( C.Mode.Log_Base == ON )
    {
      // 基本履歴情報
      TIMING_start("History_out");
      Hostonly_ H->printHistory(fp_b, rms_Var, avr_Var, container, &C, &DivC, step_end, true);
      TIMING_stop("History_out", 0.0);

      // コンポーネント
      if ( C.EnsCompo.monitor )
      {
        TIMING_start("History_out");
        Hostonly_ H->printHistoryCompo(fp_c, cmp, &C, deltaT);
        TIMING_stop("History_out", 0.0);
      }

      // 物体の力の計算
      if ( C.EnsCompo.obstacle )
      {
        TIMING_start("Force_Calculation");
        flop_count=0.0;
        calcForce(flop_count);
        TIMING_stop("Force_Calculation", flop_count);


        TIMING_start("History_out");
        Hostonly_ H->printHistoryForce(cmp, cmp_force_global);
        TIMING_stop("History_out", 0.0);
      }

      // 物体に作用する力の平均値
      if ( C.Mode.Statistic == ON && C.Interval[Control::tg_statistic].isStarted(CurrentStep, CurrentTime) )
      {
        TIMING_start("History_out");
        Hostonly_ H->printForceAvr(cmp, cmp_force_avr);
        TIMING_stop("History_out", 0.0);
      }


      // 流量収支履歴
      TIMING_start("History_out");
      Hostonly_ H->printHistoryDomfx(fp_d, &C, deltaT);
      TIMING_stop("History_out", 0.0);

    }

    // 壁面履歴情報
    if ( C.Mode.Log_Wall == ON )
    {
      //Hostonly_ H->printHistoryWall(fp_w, range_Yp, range_Ut);
    }

  }

  TIMING_stop("Loop_Utility_Section", 0.0);
  TIMING_stop("Time_Step_Loop_Section", 0.0);


  // 発散時の打ち切り
  if ( numProc > 1 )
  {
    int tmp = isNormal;
    if ( paraMngr->Allreduce(&tmp, &isNormal, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  if ( isNormal > 0 )
  {
    Hostonly_ {
      printf      ("\tForced termination : floating point exception\n");
      if ( C.Mode.Log_Base == ON) fprintf(fp_b,"\tForced termination : floating point exception\n");
    }
    return -1;
  }



  // 終了判断
  if ( C.Interval[Control::tg_compute].isLast(CurrentStep, CurrentTime) )
  {
    Hostonly_
    {
      printf      ("\tFinish : Time = %e\n", CurrentTime);
      //if ( C.Mode.Log_Base == ON) fprintf(fp_b,"\tFinish : Time = %e\n", CurrentTime); >> suppress for history plot
    }
    return 0;
  }

  return 1;
}

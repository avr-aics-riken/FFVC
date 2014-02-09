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
 * @file   ffv_Restart.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"



// #################################################################
/**
 * @brief リスタートプロセス
 * @param [in]     fp     ファイルポインタ
 */
void FFV::Restart(FILE* fp)
{
  double flop_task;
  double g[4];
  
  // 初期スタートのステップ，時間を設定する
  if ( (C.Start == initial_start) || (C.Hide.PM_Test == ON)  )
  {
    CurrentStep = 0;
    CurrentTime = 0.0;
    
    // V00の値のセット．モードがONの場合はV00[0]=1.0に設定，そうでなければtmに応じた値
    if ( C.CheckParam == ON ) RF.setV00(CurrentTime, true);
    else                      RF.setV00(CurrentTime);
    
    RF.copyV00(g);
    for (int i=0; i<4; i++) v00[i]=(REAL_TYPE)g[i];
    
    return;
  }
  
  

  switch (C.Start)
  {
      
    // 同一解像度・同一分割数のリスタート
    case restart_sameDiv_sameRes: 
      Hostonly_ fprintf(stdout, "\t>> Restart with same resolution and same num. of division\n\n");
      Hostonly_ fprintf(fp, "\t>> Restart with same resolution and same num. of division\n\n");
      break;
      
    case restart_sameDiv_refinement: // 同一分割数・リファインメント
      Hostonly_ fprintf(stdout, "\t>> Restart with refinemnt and same num. of division\n\n");
      Hostonly_ fprintf(fp, "\t>> Restart with refinemnt and same num. of division\n\n");
      break;
      
    case restart_diffDiv_sameRes:    // 異なる分割数・同一解像度
      Hostonly_ fprintf(stdout, "\t>> Restart with same resolution and different division\n\n");
      Hostonly_ fprintf(fp, "\t>> Restart with same resolution and different division\n\n");
      break;
      
    case restart_diffDiv_refinement: // 異なる分割数・リファインメント
      Hostonly_ fprintf(stdout, "\t>> Restart with refinement and different division\n\n");
      Hostonly_ fprintf(fp, "\t>> Restart with refinement and different division\n\n");
      break;
      
    default:
      Exit(0);
      break;
  }
  
  flop_task = 0.0;
  RestartInstantaneous(fp, flop_task);
  
}


// #################################################################
/**
 * @brief リスタート時の平均値ファイル読み込み
 * @param [in]  fp   ファイルポインタ
 * @param [out] flop 浮動小数点演算数
 */
void FFV::RestartAvrerage (FILE* fp, double& flop)
{
  std::string fname;
  std::string fmt(C.file_fmt_ext);
  
  unsigned m_Session_step = C.Interval[Control::tg_compute].getStartStep(); ///< セッションの開始ステップ
  double   m_Session_time = C.Interval[Control::tg_compute].getStartTime(); ///< セッションの開始時刻
  
  // ガイド出力
  int gs = C.GuideOut;
  
  
  // まだ平均値開始時刻になっていなければ，何もしない
  if ( C.Interval[Control::tg_average].getMode() == IntervalManager::By_step )
  {
    if ( m_Session_step >= C.Interval[Control::tg_average].getStartStep() )
    {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tStep : base=%u current=%u\n", m_Session_step, CurrentStep);
      Hostonly_ fprintf(fp, "\tStep : base=%u current=%u\n", m_Session_step, CurrentStep);
    }
    else
    {
      return;
    }
  }
  else if ( C.Interval[Control::tg_average].getMode() == IntervalManager::By_time )
  {
    if ( m_Session_time >= C.Interval[Control::tg_average].getStartTime() )
    {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tTime : base=%e[sec.]/%e[-] current=%e[-]\n", m_Session_time*C.Tscale, m_Session_time, CurrentTime);
      Hostonly_ fprintf(fp, "\tTime : base=%e[sec.]/%e[-] current=%e[-]\n", m_Session_time*C.Tscale, m_Session_time, CurrentTime);
    }
    else
    {
      return;
    }
  }
  else
  {
    Exit(0);
  }
  
  
  // 現在のセッションの領域分割数の取得
  int gdiv[3] = {1, 1, 1};
  
  if ( numProc > 1)
  {
    const int* m_div = paraMngr->GetDivNum();
    for (int i=0; i<3; i++ ) gdiv[i]=m_div[i];
  }
  
  
  // エラーコード
  CIO::E_CIO_ERRORCODE cio_error;
  
  
  // Averaged dataの初期化
  if ( C.Mode.Average == ON && C.Interval[Control::tg_average].isStarted(CurrentStep, CurrentTime) )
  {
    DFI_IN_PRSA = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_prsa, G_size, gdiv, cio_error);
    if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
    
    DFI_IN_VELA = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_vela, G_size, gdiv, cio_error);
    if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
    
    if ( DFI_IN_PRSA == NULL || DFI_IN_VELA == NULL ) Exit(0);
    
    
    if ( C.isHeatProblem() )
    {
      DFI_IN_TEMPA = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_tempa, G_size, gdiv, cio_error);
      if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
      if ( DFI_IN_TEMPA == NULL ) Exit(0);
    }
  }
  
  
  
  
  unsigned step_avr = 0;
  double time_avr = 0.0;
  
  
  // Pressure
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  
  //自身の領域終点インデックス
  int tail[3];
  for (int i=0; i<3; i++) tail[i] = head[i]+size[i]-1;
  
  double r_time;
  if ( DFI_IN_PRSA->ReadData(d_ap,
                             C.Restart_step,
                             guide,
                             G_size,
                             gdiv,
                             head,
                             tail,
                             r_time,
                             false,
                             step_avr,
                             time_avr) != CIO::E_CIO_SUCCESS ) Exit(0);
  
  if( d_ap == NULL ) Exit(0);
  
  
  CurrentStep_Avr = step_avr;
  CurrentTime_Avr = time_avr;
  
  if ( DFI_IN_VELA->ReadData(d_wo,
                             C.Restart_step,
                             guide,
                             G_size,
                             gdiv,
                             head,
                             tail,
                             r_time,
                             false,
                             step_avr,
                             time_avr) != CIO::E_CIO_SUCCESS ) Exit(0);
  
  if( d_wo == NULL ) Exit(0);
  
  REAL_TYPE refv = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  REAL_TYPE scale = (REAL_TYPE)step_avr;
  REAL_TYPE u0[4];
  u0[0] = v00[0];
  u0[1] = v00[1];
  u0[2] = v00[2];
  u0[3] = v00[3];
  
  fb_vin_nijk_(d_av, size, &guide, d_wo, u0, &refv, &flop);
  
  if ( (step_avr != CurrentStep_Avr) || (time_avr != CurrentTime_Avr) ) // 圧力とちがう場合
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  
  // Temperature
  if ( C.isHeatProblem() )
  {
    if ( DFI_IN_TEMPA->ReadData(d_ae,
                                C.Restart_step,
                                guide,
                                G_size,
                                gdiv,
                                head,
                                tail,
                                r_time,
                                false,
                                step_avr,
                                time_avr) != CIO::E_CIO_SUCCESS ) Exit(0);
    
    if ( d_ae == NULL ) Exit(0);
    
    if ( (step_avr != CurrentStep_Avr) || (time_avr != CurrentTime_Avr) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
}



// #################################################################
/**
 * @brief リスタートの最大値と最小値の表示
 * @param [in]  fp   ファイルポインタ
 * @param [out] flop 浮動小数点演算数
 */
void FFV::RestartDisplayMinmax(FILE* fp, double& flop)
{
  Hostonly_ fprintf(stdout, "\n\tNon-dimensional value\n");
  Hostonly_ fprintf(fp, "\n\tNon-dimensional value\n");
  REAL_TYPE f_min, f_max, min_tmp, max_tmp, vec_min[4], vec_max[4];
  
  // Velocity
  fb_minmax_v_ (vec_min, vec_max, size, &guide, v00, d_v, &flop); // allreduceすること
  
  if ( numProc > 1 )
  {
    REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
    if( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
    if( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_ fprintf(stdout, "\t\tV : min=%13.6e / max=%13.6e\n", vec_min[0], vec_max[0]);
  Hostonly_ fprintf(fp, "\t\tV : min=%13.6e / max=%13.6e\n", vec_min[0], vec_max[0]);
  
  
  // Pressure
  fb_minmax_s_ (&f_min, &f_max, size, &guide, d_p, &flop);
  
  if ( numProc > 1 )
  {
    min_tmp = f_min;
    if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    max_tmp = f_max;
    if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_ fprintf(stdout, "\t\tP : min=%13.6e / max=%13.6e\n", f_min, f_max);
  Hostonly_ fprintf(fp, "\t\tP : min=%13.6e / max=%13.6e\n", f_min, f_max);
  
  
  // temperature
  if ( C.isHeatProblem() )
  {
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ie, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    Hostonly_ fprintf(stdout, "\t\tT : min=%13.6e / max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tT : min=%13.6e / max=%13.6e\n", f_min, f_max);
  }

}



// #################################################################
/**
 * @brief リスタート時の瞬時値ファイル読み込み
 * @param [in]  fp   ファイルポインタ
 * @param [out] flop 浮動小数点演算数
 */
void FFV::RestartInstantaneous(FILE* fp, double& flop)
{
  double time, r_time;
  const unsigned step = C.Restart_step;
  std::string fname;
  std::string fmt(C.file_fmt_ext);

  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  REAL_TYPE refD = C.RefDensity;
  REAL_TYPE refV = C.RefVelocity;
  
  // ガイド出力
  int gs = C.GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;

  
  const int* m_div = paraMngr->GetDivNum();
  
  // 自身の領域終点インデックス
  int tail[3];
  for (int i=0;i<3;i++) tail[i]=head[i]+size[i]-1;
  
  
  // Pressure
  if ( DFI_IN_PRS->ReadData(d_p,
                            step,
                            guide,
                            G_size,
                            (int *)m_div,
                            head,
                            tail,
                            r_time,
                            true,
                            i_dummy,
                            f_dummy) != CIO::E_CIO_SUCCESS ) Exit(0);
  
  if ( d_p == NULL ) Exit(0);
  time = r_time;
  
  // 有次元の場合，無次元に変換する
  if ( C.Unit.File == DIMENSIONAL )
  {
    REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
    U.convArrayPrsD2ND(d_p, size, guide, bp, C.RefDensity, C.RefVelocity, flop);
  }
  
  Hostonly_ printf     ("\tPressure has read :\tstep=%d  time=%e [%s]\n",
                        step, time, (C.Unit.File == DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\tPressure has read :\tstep=%d  time=%e [%s]\n",
                    step, time, (C.Unit.File == DIMENSIONAL)?"sec.":"-");

  
  // ここでタイムスタンプを得る
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  CurrentStep = step;
  CurrentTime = time;
  
  // v00[]に値をセット
  copyV00fromRF(time);
  
  
  
  if ( DFI_IN_VEL->ReadData(d_wo,
                            step,
                            guide,
                            G_size,
                            (int *)m_div,
                            head,
                            tail,
                            r_time,
                            true,
                            i_dummy,
                            f_dummy) != CIO::E_CIO_SUCCESS ) Exit(0);
  
  if( d_wo == NULL ) Exit(0);
  
  REAL_TYPE refv = (C.Unit.File == DIMENSIONAL) ? refV : 1.0;
  REAL_TYPE u0[4];
  u0[0] = v00[0];
  u0[1] = v00[1];
  u0[2] = v00[2];
  u0[3] = v00[3];
  
  Hostonly_ printf     ("\tVelocity has read :\tstep=%d  time=%e [%s]\n",
                        step, time, (C.Unit.File == DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\tVelocity has read :\tstep=%d  time=%e [%s]\n",
                    step, time, (C.Unit.File == DIMENSIONAL)?"sec.":"-");
  
  time = r_time;
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  if ( time != CurrentTime )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // indexの変換と無次元化
  fb_vin_nijk_(d_v, size, &guide, d_wo, u0, &refv, &flop);
  

  
  if ( !C.isHeatProblem() ) return;
  
  
  
  // Instantaneous Temperature fields
  if ( DFI_IN_TEMP->ReadData(d_ws,
                             C.Restart_step,
                             guide,
                             G_size,
                             (int *)m_div,
                             head,
                             tail,
                             r_time,
                             true,
                             i_dummy,
                             f_dummy) != CIO::E_CIO_SUCCESS ) Exit(0);
  
  if( d_ws == NULL ) Exit(0);
  
  time = r_time;
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  Hostonly_ printf     ("\tTemperature has read :\tstep=%d  time=%e [%s]\n",
                        step, time, (C.Unit.File == DIMENSIONAL)?"sec.":"-");
  Hostonly_ fprintf(fp, "\tTemperature has read :\tstep=%d  time=%e [%s]\n",
                    step, time, (C.Unit.File == DIMENSIONAL)?"sec.":"-");
  
  if ( time != CurrentTime )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  if (C.Unit.File == DIMENSIONAL)
  {
    U.convArrayTmp2IE(d_ie, size, guide, d_ws, d_bcd, mat_tbl, C.BaseTemp, C.DiffTemp, true, flop);
  }
  else
  {
    U.convArrayTmp2IE(d_ie, size, guide, d_ws, d_bcd, mat_tbl, C.BaseTemp, C.DiffTemp, false, flop);
  }
  
}


// #################################################################
/* @brief リスタートモードを判定
 */
void FFV::selectRestartMode()
{
  // エラーコード
  CIO::E_CIO_ERRORCODE cio_error;
  
  
  // 現在のセッションの領域分割数の取得
  int gdiv[3] = {1, 1, 1};
  
  if ( numProc > 1)
  {
    const int* p_div = paraMngr->GetDivNum();
    for (int i=0; i<3; i++ ) gdiv[i]=p_div[i];
  }
  
  
  // Instantaneous dataの初期化
  
  // Pressure
  DFI_IN_PRS = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_prs, G_size, gdiv, cio_error);
  if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
  
  
  // Velocity
  DFI_IN_VEL = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_vel, G_size, gdiv, cio_error);
  if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
  
  if ( DFI_IN_PRS == NULL || DFI_IN_VEL == NULL ) Exit(0);
  
  
  // Fvelocity
  DFI_IN_FVEL = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_fvel, G_size, gdiv, cio_error);
  if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
  if ( DFI_IN_FVEL == NULL ) Exit(0);
  
  // Temperature
  if ( C.isHeatProblem() )
  {
    DFI_IN_TEMP = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_temp, G_size, gdiv, cio_error);
    if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
    if ( DFI_IN_TEMP == NULL ) Exit(0);
  }
  
  
  ///CIO FileInfoの成分名の取得、成分名が登録されていないときは、空白が戻される
  /*
   std::string VCompVariable[3];
   VCompVariable[0]=DFI_IN_VEL->getComponentVariable(0);
   VCompVariable[1]=DFI_IN_VEL->getComponentVariable(1);
   VCompVariable[2]=DFI_IN_VEL->getComponentVariable(2);
   
   ///CIO TimeSliceのVectorMinMaxの取得、取得出来たときはCIO::E_CIO_SUCCESS
   double vec_minmax[2];
   cio_error = DFI_IN_VEL->getVectorMinMax(C.Restart_step,vec_minmax[0],vec_minmax[1]);
   
   ///CIO TimeSlice minmaxの取得、取得出来たときはCIO::E_CIO_SUCCESS
   double minmax[6];
   cio_error =  DFI_IN_VEL->getMinMax(C.Restart_step,0,minmax[0],minmax[1]);
   cio_error =  DFI_IN_VEL->getMinMax(C.Restart_step,1,minmax[2],minmax[3]);
   cio_error =  DFI_IN_VEL->getMinMax(C.Restart_step,2,minmax[4],minmax[5]);
   */
  
  
  
  /* Averaged dataの初期化
  if ( C.Mode.Average == ON && C.Interval[Control::tg_average].isStarted(CurrentStep, CurrentTime) )
  {
    DFI_IN_PRSA = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_prsa, G_size, gdiv, cio_error);
    if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
    
    DFI_IN_VELA = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_vela, G_size, gdiv, cio_error);
    if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
    
    if ( DFI_IN_PRSA == NULL || DFI_IN_VELA == NULL ) Exit(0);
    
    
    if ( C.isHeatProblem() )
    {
      DFI_IN_TEMPA = cio_DFI::ReadInit(MPI_COMM_WORLD, C.f_dfi_in_tempa, G_size, gdiv, cio_error);
      if ( cio_error != CIO::E_CIO_SUCCESS ) Exit(0);
      if ( DFI_IN_TEMPA == NULL ) Exit(0);
    }
  }
  */

  
  bool isSameDiv = true; // 同一分割数
  bool isSameRes = true; // 同一解像度
  
  
  // 前のセッションの領域分割数の取得
  int* DFI_div=NULL;
  
  if ( C.KindOfSolver != SOLID_CONDUCTION )
  {
    DFI_div = DFI_IN_PRS->GetDFIGlobalDivision();
  }
  else
  {
    DFI_div = DFI_IN_TEMP->GetDFIGlobalDivision();
  }
  
  
  
  // 前セッションと領域分割数が異なる場合
  for (int i=0; i<3; i++ )
  {
    if ( gdiv[i] != DFI_div[i] )
    {
      isSameDiv = false;
    }
  }
  
  // 前のセッションの全要素数の取得
  int* DFI_G_size = DFI_IN_PRS->GetDFIGlobalVoxel();
  
  
  // 前セッションと全要素数が異なる場合
  for (int i=0; i<3; i++ )
  {
    if ( G_size[i] != DFI_G_size[i] )
    {
      isSameRes = false;
    }
  }

  //  ボクセル数が2倍のチェック
  if ( !isSameRes )
  {
    for(int i=0; i<3; i++)
    {
      if ( G_size[i] != DFI_G_size[i]*2 )
      {
        printf("\tDimension size error (%d %d %d)\n", G_size[0], G_size[1], G_size[2]);
        Exit(0);
      }
    }
  }
  
  
  // モード判定と登録
  if ( isSameDiv )
  {
    if ( isSameRes ) // 同一解像度、同一分割数
    {
      C.Start = restart_sameDiv_sameRes;
    }
    else // Refinement、同一分割数
    {
      C.Start = restart_sameDiv_refinement;
    }
  }
  else
  {
    if ( isSameRes ) // 同一解像度、異なる分割数
    {
      C.Start = restart_diffDiv_sameRes;
    }
    else // Refinement、異なる分割数
    {
      C.Start = restart_diffDiv_refinement;
    }
  }
  
}

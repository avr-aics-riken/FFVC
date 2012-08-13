// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   ffv.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"

// コンストラクタ
FFV::FFV()
{
  ffv_procGrp = 0;
  ModeTiming = 0;
  G_Acell = 0;
  G_Fcell = 0;
  G_Wcell = 0;
  L_Acell = 0;
  L_Fcell = 0;
  L_Wcell = 0;
  
  CurrentTime = 0.0;
  CurrentTime_Avr = 0.0;
  Session_StartTime = 0.0;
  Session_CurrentTime = 0.0;
  
  Session_LastStep = 0;
  Session_CurrentStep = 0;
  Session_StartStep = 0;
  CurrentStep = 0;
  CurrentStep_Avr = 0;
  
  REAL_TYPE deltaT; ///< 時間積分幅（無次元）
  id_of_solid = 2;
  
  for (int i=0; i<3; i++) 
  {
    G_size[i]= 0;
    G_origin[i] = 0.0;
    G_region[i] = 0.0;
  }
  
  // dfi管理
  for (int i=0; i<var_END; i++) dfi_mng[i]=0;
  
  fp_b = NULL;
  fp_w = NULL;
  fp_c = NULL;
  fp_d = NULL;
  fp_i = NULL;
  fp_f = NULL;
  
  Ex = NULL;
  mat = NULL;
  cmp = NULL;
  paraMngr = NULL;
  
}


// デストラクタ
FFV::~FFV()
{
  
}



// 時間平均値のファイル出力
void FFV::AverageOutput(double& flop)
{
  // 出力ファイルの指定が有次元の場合
  REAL_TYPE timeAvr;
  if (C.Unit.File == DIMENSIONAL) 
  {
    timeAvr = CurrentTime_Avr * C.Tscale;
  }
  else 
  {
    timeAvr = CurrentTime_Avr;
  }
  
  // 出力用のヘッダ
  REAL_TYPE m_org[3], m_pit[3];
  
  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++) 
  {
    m_org[i] = origin[i] - pitch[i]*(REAL_TYPE)C.GuideOut;
    m_pit[i] = pitch[i];
  }
  
  // 平均操作の母数
  int stepAvr = (int)CurrentStep_Avr;
  REAL_TYPE scale = 1.0 / (REAL_TYPE)stepAvr;
  
  // ガイドセル出力
  int gc_out = C.GuideOut;
  
  // 出力ファイル名
  std::string tmp;
  
  int d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  // Pressure
  if (C.Unit.File == DIMENSIONAL) 
  {
    REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
    fb_prs_nd2d_(d_ws, d_ap, &d_length, &bp, &C.RefDensity, &C.RefVelocity, &scale, &flop);
  }
  else 
  {
    fb_xcopy_(d_ws, d_ap, &d_length, &scale, &flop);
  }
  
  tmp = DFI.Generate_FileName(C.f_AvrPressure, stepAvr, myRank, (bool)C.FIO.IO_Output);
  F.writeScalar(tmp, size, guide, d_ws, stepAvr, timeAvr, m_org, m_pit, gc_out);
  Hostonly_ if ( !DFI.Write_DFI_File(C.f_AvrPressure, stepAvr, dfi_mng[var_Pressure_Avr], (bool)C.FIO.IO_Output) ) Exit(0);
  
  // Velocity
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  fb_shift_refv_out_(d_wo, d_av, size, &guide, v00, &scale, &unit_velocity, &flop);
  
  tmp = DFI.Generate_FileName(C.f_AvrVelocity, stepAvr, myRank, (bool)C.FIO.IO_Output);
  F.writeVector(tmp, size, guide, d_wo, stepAvr, timeAvr, m_org, m_pit, gc_out);
  Hostonly_ if ( !DFI.Write_DFI_File(C.f_AvrVelocity, stepAvr, dfi_mng[var_Velocity_Avr], (bool)C.FIO.IO_Output) ) Exit(0);
  
  // Temperature
  if( C.isHeatProblem() )
  {
    if (C.Unit.File == DIMENSIONAL) 
    {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      fb_tmp_nd2d_(d_ws, d_at, &d_length, &C.BaseTemp, &C.DiffTemp, &klv, &scale, &flop);
    }
    else 
    {
      fb_xcopy_(d_ws, d_at, &d_length, &scale, &flop);
    }
    
    tmp = DFI.Generate_FileName(C.f_AvrTemperature, stepAvr, myRank, (bool)C.FIO.IO_Output);
    F.writeScalar(tmp, size, guide, d_ws, stepAvr, timeAvr, m_org, m_pit, gc_out);
    Hostonly_ if( !DFI.Write_DFI_File(C.f_AvrTemperature, stepAvr, dfi_mng[var_Temperature_Avr], (bool)C.FIO.IO_Output) ) Exit(0);
  }
}



// 時間平均操作を行う
void FFV::Averaging_Time(double& flop)
{
  CurrentStep_Avr++;
  CurrentTime_Avr += DT.get_DT();
  
  fb_average_s_(d_ap, size, &guide, d_p, &flop);
  fb_average_v_(d_av, size, &guide, d_v, &flop);
  
  if ( C.isHeatProblem() ) 
  {
    fb_average_s_(d_at, size, &guide, d_t, &flop);
  }
}



// 全ノードについて，ローカルノード1面・一層あたりの通信量の和を返す
double FFV::count_comm_size(const int sz[3], const int guide)
{
  double c = 0.0;
  
  // 内部面のみをカウントする
  for (int n=0; n<6; n++) 
  {
    if ( nID[n] >= 0 ) {
      
      switch (n) 
      {
        case X_MINUS:
        case X_PLUS:
          c += (double)(sz[1]*sz[2]);
          break;
          
        case Y_MINUS:
        case Y_PLUS:
          c += (double)(sz[0]*sz[2]);
          break;
          
        case Z_MINUS:
        case Z_PLUS:
          c += (double)(sz[0]*sz[1]);
          break;
      }
    }
  }
  
  if ( numProc > 1 )
  {
    double tmp = c;
    if ( paraMngr->Allreduce(&tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return c * sizeof(double); // Byte
}



// 外部計算領域の各面における総流量と対流流出速度を計算する
// 流出境界のみ和をとる，その他は既知
void FFV::DomainMonitor(BoundaryOuter* ptr, Control* R, double& flop)
{
  if ( !ptr ) Exit(0);
  BoundaryOuter* obc=NULL;
  
  obc = ptr;
  
  REAL_TYPE ddh;
  REAL_TYPE *vv, ec;
  REAL_TYPE u_sum, u_min, u_max, u_avr;
  REAL_TYPE tmp_sum, tmp_min, tmp_max;
  int ofv;
  
	ddh = deltaX * deltaX;
  
  for (int face=0; face<NOFACE; face++) 
  {
    
    // ofv (1-MINMAX, 2-AVERAGE) ゼロでなければ，流出境界
    ofv = obc[face].get_ofv();
    
    // OUTFLOW, SPEC_VELのときの有効セル数
    ec = (REAL_TYPE)obc[face].get_ValidCell();

    
    // 各プロセスの外部領域面の速度をvv[]にコピー
    vv = obc[face].getDomainV();
    
    
    // 流出境界のモード
    if (ofv == V_AVERAGE) // average
    {
      
      // 非境界面にはゼロなので，単に足し込むだけ
      u_sum = vv[0];
      
      if ( numProc > 1 ) 
      {
        tmp_sum = u_sum;
        if ( paraMngr->Allreduce(&tmp_sum, &u_sum, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      }

      u_avr = (ec != 0.0) ? u_sum / ec : 0.0;
      flop = flop + 1.0;
      
    }
    else if (ofv == V_MINMAX) // minmax
    {
      u_sum = vv[0];
      u_min = vv[1];
      u_max = vv[2];

      printf("rank=%d : u_min=%e umax=%e u_avr=%e\n", myRank, u_min, u_max, u_avr);
      
      if ( numProc > 1 )
      {
        tmp_sum = u_sum;
        tmp_min = u_min;
        tmp_max = u_max;
        if ( paraMngr->Allreduce(&tmp_sum, &u_sum, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        if ( paraMngr->Allreduce(&tmp_min, &u_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        if ( paraMngr->Allreduce(&tmp_max, &u_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }

      u_avr = 0.5*(u_min+u_max);
      flop = flop + 2.0;
      
    }
    else // 非OUTFLOW BCは速度がストアされている
    {
      u_sum = vv[0] * ec;
      u_avr = vv[0];
      flop = flop + 1.0;
    }
    
    // コントロールクラスにコピー
    R->V_Dface[face] = u_avr;
    R->Q_Dface[face] = u_sum * ddh;
    flop = flop + 1.0;
  }
  
}



// ファイル出力
void FFV::FileOutput(double& flop, const bool restart)
{
  REAL_TYPE scale = 1.0;
  int d_length;
  
  // d_p0をワークとして使用
  
  // 出力用のヘッダ
  REAL_TYPE m_org[3], m_pit[3];
  
  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++) 
  {
    m_org[i] = origin[i] - pitch[i]*(REAL_TYPE)C.GuideOut;
    m_pit[i] = pitch[i];
  }
  
  // 出力ファイルの指定が有次元の場合
  if ( C.Unit.File == DIMENSIONAL ) 
  {
    for (int i=0; i<3; i++) 
    {
      m_org[i] *= C.RefLength;
      m_pit[i] *= C.RefLength;
    }
  }
  
  // ステップ数
  int m_step = (int)CurrentStep;
  
  // 時間の次元変換
  REAL_TYPE m_time;
  if (C.Unit.File == DIMENSIONAL) 
  {
    m_time = (REAL_TYPE)CurrentTime * C.Tscale;
  }
  else 
  {
    m_time = (REAL_TYPE)CurrentTime;
  }
  
  // ガイドセル出力
  int gc_out = C.GuideOut;
  
  // 出力ファイル名
  std::string tmp;
  
  // 並列出力モード
  bool pout = ( C.FIO.IO_Output == IO_GATHER ) ? false : true;
  
  // Divergence デバッグ用なので無次元のみ
  if ( C.FIO.Div_Debug == ON ) 
  {
    
    REAL_TYPE coef = (REAL_TYPE)DT.get_DT()/(deltaX*deltaX); /// 発散値を計算するための係数　dt/h^2
    F.cnv_Div(d_ws, d_sq, size, guide, coef, flop);
    
    tmp = DFI.Generate_FileName(C.f_DivDebug, m_step, myRank, pout);
    F.writeScalar(tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.Write_DFI_File(C.f_DivDebug, m_step, dfi_mng[var_Divergence], pout) ) Exit(0);
    
  }
  
  
  // リスタート用
  std::string prs_restart("prs_restart_");
  std::string vel_restart("vel_restart_");
  std::string temp_restart("temp_restart_");
  
  // Pressure
  d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  if (C.Unit.File == DIMENSIONAL) 
  {
    REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
    fb_prs_nd2d_(d_ws, d_p, &d_length, &bp, &C.RefDensity, &C.RefVelocity, &scale, &flop);
  }
  else 
  {
    fb_xcopy_(d_ws, d_p, &d_length, &scale, &flop);
  }
  
  if ( !restart ) 
  {
    tmp = DFI.Generate_FileName(C.f_Pressure, m_step, myRank, pout);
  }
  else 
  {
    tmp = DFI.Generate_FileName(prs_restart, m_step, myRank, pout);
  }
  
  F.writeScalar(tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
  
  Hostonly_ if ( !DFI.Write_DFI_File(C.f_Pressure, m_step, dfi_mng[var_Pressure], pout) ) Exit(0);
  
  
  // Velocity
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  fb_shift_refv_out_(d_wo, d_v, size, &guide, v00, &scale, &unit_velocity, &flop);
  
  if ( !restart ) 
  {
    tmp = DFI.Generate_FileName(C.f_Velocity, m_step, myRank, pout);
  }
  else 
  {
    tmp = DFI.Generate_FileName(vel_restart, m_step, myRank, pout);
  }
  
  F.writeVector(tmp, size, guide, d_wo, m_step, m_time, m_org, m_pit, gc_out);
  
  Hostonly_ if ( !DFI.Write_DFI_File(C.f_Velocity, m_step, dfi_mng[var_Velocity], pout) ) Exit(0);
  
  
  // Tempearture
  if( C.isHeatProblem() )
  {
    
    d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
    
    if (C.Unit.File == DIMENSIONAL) 
    {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      fb_tmp_nd2d_(d_ws, d_t, &d_length, &C.BaseTemp, &C.DiffTemp, &klv, &scale, &flop);
    }
    else 
    {
      fb_xcopy_(d_ws, d_t, &d_length, &scale, &flop);
    }
    
    if ( !restart ) 
    {
      tmp = DFI.Generate_FileName(C.f_Temperature, m_step, myRank, pout);
    }
    else 
    {
      tmp = DFI.Generate_FileName(temp_restart, m_step, myRank, pout);
    }
    
    F.writeScalar(tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.Write_DFI_File(C.f_Temperature, m_step, dfi_mng[var_Temperature], pout) ) Exit(0);
  }

  
  
  // Total Pressure
  if (C.Mode.TP == ON ) {
    
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop);
    
    // convert non-dimensional to dimensional, iff file is dimensional
    if (C.Unit.File == DIMENSIONAL) 
    {
      F.cnv_TP_ND2D(d_ws, d_p0, size, guide, C.RefDensity, C.RefVelocity, flop);
    }
    else 
    {
      REAL_TYPE* tp;
      tp = d_ws; d_ws = d_p0; d_p0 = tp;
    }
    
    tmp = DFI.Generate_FileName(C.f_TotalP, m_step, myRank, pout);
    F.writeScalar(tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.Write_DFI_File(C.f_TotalP, m_step, dfi_mng[var_TotalP], pout) ) Exit(0);

  }
  
  
  // Vorticity
  if (C.Mode.VRT == ON ) 
  {
    
    rot_v_(d_wv, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    REAL_TYPE  vz[3];
    vz[0] = vz[1] = vz[2] = 0.0;
    unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity/C.RefLength : 1.0;
    fb_shift_refv_out_(d_wo, d_wv, size, &guide, vz, &scale, &unit_velocity, &flop);
    
    tmp = DFI.Generate_FileName(C.f_Vorticity, m_step, myRank, pout);
    F.writeVector(tmp, size, guide, d_wo, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.Write_DFI_File(C.f_Vorticity, m_step, dfi_mng[var_Vorticity], pout) ) Exit(0);

  }
  
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C.Mode.I2VGT == ON ) {
    
    i2vgt_ (d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    // 無次元で出力
    d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
    fb_xcopy_(d_ws, d_p0, &d_length, &scale, &flop);
    
    tmp = DFI.Generate_FileName(C.f_I2VGT, m_step, myRank, pout);
    F.writeScalar(tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.Write_DFI_File(C.f_I2VGT, m_step, dfi_mng[var_I2vgt], pout) ) Exit(0);

  }
  
  // Helicity
  if (C.Mode.Helicity == ON ) 
  {
    
    helicity_(d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    // 無次元で出力
    d_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
    fb_xcopy_(d_ws, d_p0, &d_length, &scale, &flop);
    
    tmp = DFI.Generate_FileName(C.f_Helicity, m_step, myRank, pout);
    F.writeScalar(tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.Write_DFI_File(C.f_Helicity, m_step, dfi_mng[var_Helicity], pout) ) Exit(0);

  }

}


// タイムステップループ
int FFV::MainLoop()
{
  int ret = 1;
  
  for (int i=1; i<=Session_LastStep; i++)
  {
    Session_CurrentStep = i;
    
    int loop_ret = Loop(i);
    
    switch (loop_ret) 
    {
      case -1: // error
        return -1;
        break;
        
      case 0: // forced terminated
        ret = 1;
        break;
        
      case 1: // normal
        ret = 1;
        break;
        
      default:
        return -1;
    }
    
    if( loop_ret == 0 ) break;
  }

  if ( !stepPost() ) return -1;

  return ret;
}



/**
 * @brief Poissonのノルムを計算する
 * @retval 収束値
 * @param [in] IC ItrCtlクラス
 */
REAL_TYPE FFV::Norm_Poisson(ItrCtl* IC)
{
  REAL_TYPE nrm;
  REAL_TYPE rms;
  REAL_TYPE convergence; /// 収束判定値
  REAL_TYPE coef = deltaT/(deltaX*deltaX); /// 発散値を計算するための係数　dt/h^2
  REAL_TYPE tmp; /// テンポラリバッファ
  double flop_count;
  
  //d_sqに発散がストアされている
  
  
  switch (IC->get_normType()) 
  {
    case ItrCtl::v_div_max:
      TIMING_start(tm_poi_itr_sct_5); // >>> Poisson Iteration subsection 5
      
      TIMING_start(tm_norm_div_max);
      flop_count=0.0;
      norm_v_div_max_(&nrm, size, &guide, d_sq, &coef, d_bcp, &flop_count);
      TIMING_stop(tm_norm_div_max, flop_count);
      
      if ( numProc > 1 ) 
      {
        TIMING_start(tm_norm_comm);
        tmp = nrm;
        if ( paraMngr->Allreduce(&tmp, &nrm, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0); // 最大値
        TIMING_stop(tm_norm_comm, 2.0*numProc*sizeof(REAL_TYPE) ); // 双方向 x ノード数
      }
      convergence = nrm;
      IC->set_normValue( convergence );
      
      TIMING_stop(tm_poi_itr_sct_5, 0.0); // <<< Poisson Iteration subsection 5
      break;
      
    case ItrCtl::v_div_l2:
      TIMING_start(tm_poi_itr_sct_5); // >>> Poisson Iteration subsection 5
      
      TIMING_start(tm_norm_div_l2);
      flop_count=0.0;
      norm_v_div_l2_(&rms, size, &guide, d_sq, &coef, d_bcp, &flop_count);
      TIMING_stop(tm_norm_div_l2, flop_count);
      
      if ( numProc > 1 ) 
      {
        TIMING_start(tm_norm_comm);
        tmp = rms; 
        if ( paraMngr->Allreduce(&tmp, &rms, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0); // 和
        TIMING_stop(tm_norm_comm, 2.0*numProc*sizeof(REAL_TYPE) ); // 双方向 x ノード数
      }
      convergence = sqrt(rms/(REAL_TYPE)numProc); //RMS
      IC->set_normValue( convergence );
      
      TIMING_stop(tm_poi_itr_sct_5, 0.0); // <<< Poisson Iteration subsection 5
      break;
      
    case ItrCtl::p_res_l2_a:
    case ItrCtl::p_res_l2_r:
      convergence = IC->get_normValue();
      break;
      
    case ItrCtl::v_div_max_dbg:
      TIMING_start(tm_poi_itr_sct_5); // >>> Poisson Iteration subsection 5
      
      TIMING_start(tm_norm_div_max_dbg);
      flop_count=0.0;
      int index[3];
      index[0] = 0;
      index[1] = 0;
      index[2] = 0;
      norm_v_div_dbg_(&nrm, &rms, index, size, &guide, d_sq, &coef, d_bcp, &flop_count);
      TIMING_stop(tm_norm_div_max_dbg, flop_count);
      
      //@todo ここで，最大値のグローバルなindexの位置を計算する
      
      if ( numProc > 1 ) 
      {
        TIMING_start(tm_norm_comm);
        tmp = nrm;
        if ( paraMngr->Allreduce(&tmp, &nrm, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0); // 最大値
        tmp = rms;
        if ( paraMngr->Allreduce(&tmp, &rms, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0); // 和
        TIMING_stop(tm_norm_comm, 2.0*numProc*sizeof(REAL_TYPE)*2.0 ); // 双方向 x ノード数
      }
      rms = sqrt(rms/(REAL_TYPE)numProc); // RMS
      convergence = nrm; // ノルムは最大値を返す
      IC->set_normValue( convergence );
      
      Hostonly_ {
        H->printHistoryItr(fp_i, IC->LoopCount, nrm, index);
        fflush(fp_i);
      }
      
      TIMING_stop(tm_poi_itr_sct_5, 0.0); // <<< Poisson Iteration subsection 5
      break;
      
    default:
      stamped_printf("\tInvalid convergence type\n");
      Exit(0);
  }
  
  return convergence;
}



// 圧力の引き戻し操作を行う
void FFV::Pressure_Shift()
{
  REAL_TYPE base=0.0;
  
  int n = C.Mode.Pshift;
  
  face_avr_sampling_(d_p, size, &guide, &n, &base);
  shift_pressure_   (d_p, size, &guide, &base);
}


// タイミング測定区間にラベルを与えるラッパー
void FFV::set_label(const int key, char* label, PerfMonitor::Type type, bool exclusive)
{
  // 文字数がTM_LABEL_MAX-1を超えるものはカット
  int len = strlen(label);
  char label_tmp[TM_LABEL_MAX];
  
  if ( len>TM_LABEL_MAX-1 ) 
  {
    strncpy(label_tmp, label, TM_LABEL_MAX-1);
    printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
  }
  else 
  {
    strcpy(label_tmp, label);
  }
  
  // Performance Monitorへの登録
  string tmp(label_tmp);
  PM.setProperties(key, tmp, type, exclusive);
  
  // FX用プロファイラの文字列登録
  strcpy(tm_label_ptr[key], label_tmp);
}


// タイミング測定区間にラベルを与える
void FFV::set_timing_label()
{
  // ラベルの設定
  set_label(tm_init_sct,           "Initialization_Section",  PerfMonitor::CALC, false);
  set_label(tm_init_alloc,         "Allocate_Arrays",         PerfMonitor::CALC);
  
  set_label(tm_voxel_prep_sct,     "Voxel_Prep_Section",      PerfMonitor::CALC, false);
  set_label(tm_voxel_load,         "Loading_Voxel_File",      PerfMonitor::CALC);
  set_label(tm_polygon_load,       "Loading_Polygon_File",    PerfMonitor::CALC);
  set_label(tm_cutinfo,            "Cut_Information",         PerfMonitor::CALC);
  set_label(tm_cmp_vertex8,        "Compo_Vertex8",           PerfMonitor::CALC);
  set_label(tm_cmp_subdivision,    "Compo_Subdivision",       PerfMonitor::CALC);
  // end of Voxel Prep. Section
  
  set_label(tm_restart,            "Restart_Process",         PerfMonitor::CALC);
  
  // end of Initialization Section
  
  
  // Loop section
  set_label(tm_loop_sct,           "Time_Step_Loop_Section",  PerfMonitor::CALC, false); 
  
  set_label(tm_vmax,               "Search_Vmax",             PerfMonitor::CALC);
  set_label(tm_vmax_comm,          "A_R_Vmax",                PerfMonitor::COMM);
  
  set_label(tm_flow_sct,           "Flow_Section",            PerfMonitor::CALC, false); 
  
  set_label(tm_frctnl_stp_sct,     "NS__F_Step_Section",      PerfMonitor::CALC, false); 
  
  set_label(tm_frctnl_stp_sct_1,   "NS__F_Step_Sct_1",        PerfMonitor::CALC, false);
  set_label(tm_spec_vel,           "Assign_BC_Velocity",      PerfMonitor::CALC);
  set_label(tm_WallFunc,           "Friction_Velocity",       PerfMonitor::CALC);
  // end of NS: F-Step Sct:1
  
  set_label(tm_frctnl_stp_sct_2,   "NS__F_Step_Sct_2",        PerfMonitor::CALC, false);
  set_label(tm_pseudo_vec,         "Pseudo_Velocity",         PerfMonitor::CALC);
  set_label(tm_pvec_flux,          "Pseudo_Vel_Flux_BC",      PerfMonitor::CALC);
  set_label(tm_pvec_ee,            "Pvec_Euler_Explicit",     PerfMonitor::CALC);
  set_label(tm_pvec_ab,            "Pvec_Adams_Bashforth",    PerfMonitor::CALC);
  set_label(tm_pvec_abcn,          "Pvec_AB_CN",              PerfMonitor::CALC);
  set_label(tm_pvec_abcn_df_ee,    "Pvec_AB_CN_Diff_EE",      PerfMonitor::CALC);
  set_label(tm_pvec_abcn_df_ee_BC, "Pvec_AB_CN_Diff_EE_BC",   PerfMonitor::CALC);
  // end of NS: F-Step Sct:2
  
  set_label(tm_frctnl_stp_sct_3,   "NS__F_Step_Sct_3",        PerfMonitor::CALC, false);
  set_label(tm_forcing,            "Forcing_Pseudo_Velocity", PerfMonitor::CALC);
  set_label(tm_buoyancy,           "Buoyancy",                PerfMonitor::CALC);
  set_label(tm_pvec_BC,            "Pseudo_Velocity_BC",      PerfMonitor::CALC);
  set_label(tm_pvec_comm,          "Sync_Pseudo_Velocity",    PerfMonitor::COMM);
  // end of NS: F-Step Sct:3
  
  set_label(tm_frctnl_stp_sct_4,   "NS__F_Step_Sct_4",        PerfMonitor::CALC, false);
  // end of NS: F-Step Sct:4
  
  // end of NS: F-Step Section
  
  set_label(tm_poi_src_sct,        "Poisson__Source_Section", PerfMonitor::CALC, false);
  set_label(tm_div_pvec,           "Divergence_of_Pvec",      PerfMonitor::CALC);
  set_label(tm_poi_src_vbc,        "Poisson_Src_VBC",         PerfMonitor::CALC);
  set_label(tm_poi_src_nrm,        "Poisson_Src_Norm",        PerfMonitor::CALC);
  set_label(tm_poi_src_comm,       "A_R_Poisson_Src",         PerfMonitor::COMM);
  // end of Poisson: Source Section
  
  set_label(tm_poi_itr_sct,        "Poisson__Iteration_Sct",  PerfMonitor::CALC, false);
  set_label(tm_hstry_itr,          "History_Iteration",       PerfMonitor::CALC);
  
  set_label(tm_poi_itr_sct_1,      "Poisson__Itr_Sct_1",      PerfMonitor::CALC, false);
  set_label(tm_force_src,          "Forcing Source",          PerfMonitor::CALC);
  // end of Poisson: Itr. Sct:1
  
  set_label(tm_poi_itr_sct_2,      "Poisson__Itr_Sct_2",      PerfMonitor::CALC, false); // LS_Binary(), LS_Planar()
  set_label(tm_poi_PSOR,           "Poisson_PSOR",            PerfMonitor::CALC);
  set_label(tm_poi_setup,          "Poisson_Setup_for_Itr",   PerfMonitor::CALC);
  set_label(tm_poi_SOR2SMA,        "Poisson_SOR2_SMA",        PerfMonitor::CALC);
  set_label(tm_poi_SOR2CMA,        "Poisson_SOR2_CMA",        PerfMonitor::CALC);
  set_label(tm_poi_BC,             "Poisson_BC",              PerfMonitor::CALC);
  // end of Poisson: Itr. Sct:2
  
  set_label(tm_poi_itr_sct_3,      "Poisson__Itr_Sct_3",      PerfMonitor::CALC, false); // LS_Binary(), LS_Planar()
  set_label(tm_poi_comm,           "Sync_Pressure",           PerfMonitor::COMM);
  set_label(tm_poi_res_comm,       "A_R_Poisson_Residual",    PerfMonitor::COMM);
  // end of Poisson: Itr. Sct:3
  
  set_label(tm_poi_itr_sct_4,      "Poisson__Itr_Sct_4",      PerfMonitor::CALC, false);
  set_label(tm_prj_vec,            "Projection_Velocity",     PerfMonitor::CALC);
  set_label(tm_prj_vec_bc,         "Projection_Velocity_BC",  PerfMonitor::CALC);
  set_label(tm_prj_vec_bc_comm,    "A_R_Projection_VBC",      PerfMonitor::COMM);
  set_label(tm_prj_frc_mod,        "Projection_Forcing",      PerfMonitor::CALC);
  set_label(tm_prj_frc_mod_comm,   "A_R_Projection_Forcing",  PerfMonitor::COMM);
  set_label(tm_prj_frc_dir,        "Forcing_Direction",       PerfMonitor::CALC);
  set_label(tm_vec_BC,             "Velocity_BC",             PerfMonitor::CALC);
  // end of Poisson: Itr. Sct:4
  
  set_label(tm_poi_itr_sct_5,      "Poisson__Itr_Sct_5",      PerfMonitor::CALC, false);
  set_label(tm_norm_div_max,       "Poisson_Norm_Div_max",    PerfMonitor::CALC);
  set_label(tm_norm_div_l2,        "Poisson_Norm_Div_L2",     PerfMonitor::CALC);
  set_label(tm_norm_div_max_dbg,   "Poi_Norm_Div_Max_Dbg",    PerfMonitor::CALC);
  set_label(tm_norm_comm,          "A_R_Poisson_Norm",        PerfMonitor::COMM);
  // end of Poisson: Itr. Sct:5
  
  // end of Poisson: Iteration Sct.
  
  set_label(tm_NS_loop_post_sct,   "NS__Loop_Post_Section",   PerfMonitor::CALC, false);
  set_label(tm_vectors_comm,       "Sync_Velocity",           PerfMonitor::COMM);
  set_label(tm_domain_monitor,     "Domain_Monitor",          PerfMonitor::CALC);
  set_label(tm_VBC_update,         "Velocity_BC_Update",      PerfMonitor::CALC);
  set_label(tm_LES_eddy,           "Eddy_Viscosity",          PerfMonitor::CALC);
  set_label(tm_LES_eddy_comm,      "Sync_Eddy_Viscosity",     PerfMonitor::COMM);
  set_label(tm_pressure_shift,     "Shift_Pressure",          PerfMonitor::COMM);
  // end of NS: Loop Post Section
  
  // end of Flow section
  
  
  set_label(tm_heat_sct,           "Heat_Section",            PerfMonitor::CALC, false);
  
  set_label(tm_heat_convection_sct,"Thermal_Convection_Sct",  PerfMonitor::CALC, false);
  set_label(tm_heat_spec_temp,     "Thermal_Assign_BC_Temp",  PerfMonitor::CALC);
  set_label(tm_heat_cnv,           "Thermal_Convection",      PerfMonitor::CALC);
  set_label(tm_heat_cnv_BC,        "Thermal_Convection_BC",   PerfMonitor::CALC);
  set_label(tm_heat_cnv_EE,        "Thermal_Convection_EE",   PerfMonitor::CALC);
  // end of Thermal Convection Sct.
  
  set_label(tm_heat_diffusion_sct, "Thermal_Diffusion_Sct",   PerfMonitor::CALC, false);
  set_label(tm_heat_diff_OBC,      "Thermal_Diff_Outer_BC",   PerfMonitor::CALC);
  // end of Thermal Diffusion Sct.
  
  set_label(tm_heat_diff_sct_1,    "Thermal_Diffusion_Sct_1", PerfMonitor::CALC, false);
  set_label(tm_heat_diff_IBC_vol,  "Thermal_Diff_IBC_Vol",    PerfMonitor::CALC);
  set_label(tm_heat_diff_comm,     "Sync_Thermal",            PerfMonitor::COMM);
  // end of Thermal Diffusion Sct:1
  
  set_label(tm_heat_diff_sct_2,    "Thermal_Diffusion_Sct_2", PerfMonitor::CALC, false);
  set_label(tm_heat_diff_IBC_face, "Thermal_Diff_IBC_Face",   PerfMonitor::CALC);
  set_label(tm_heat_diff_OBC_face, "Thermal_Diff_OBC_Face",   PerfMonitor::CALC);
  set_label(tm_heat_diff_QBC_comm, "Sync_Thermal_QBC",        PerfMonitor::COMM);
  // end of Thermal Diffusion Sct:2
  
  set_label(tm_heat_diff_sct_3,    "Thermal_Diffusion_Sct_3", PerfMonitor::CALC, false);
  set_label(tm_heat_diff_EE,       "Thermal_Diff_EE",         PerfMonitor::CALC);
  set_label(tm_heat_diff_PSOR,     "Thermal_Diff_PSOR",       PerfMonitor::CALC);
  set_label(tm_heat_update_comm,   "Sync_Thermal_Update",     PerfMonitor::COMM);
  set_label(tm_heat_diff_res_comm, "A_R_Thermal_Diff_Res",    PerfMonitor::COMM);
  // end of Thermal Diffusion Sct:3
  
  set_label(tm_heat_loop_post_sct, "Thermal_Loop_Post_Sct",   PerfMonitor::CALC, false);
  set_label(tm_heat_range,         "Thermal_Range_Cut",       PerfMonitor::CALC);
  // end of Thermal Loop Post Sct.
  
  // end of Heat Section
  
  
  set_label(tm_vof_sct,            "VOF_Section",             PerfMonitor::CALC, false);
  set_label(tm_vof_cnv,            "VOF_Convection",          PerfMonitor::CALC);
  set_label(tm_vof_cnv_comm,       "Sync_VOF_Convection",     PerfMonitor::COMM);
  // end of VOF section
  
  
  set_label(tm_loop_uty_sct,       "Loop_Utility_Section",    PerfMonitor::CALC, false);
  
  set_label(tm_loop_uty_sct_1,     "Loop_Utility_Sct_1",      PerfMonitor::CALC, false);
  set_label(tm_average_time,       "Averaging_Time",          PerfMonitor::CALC);
  set_label(tm_stat_space,         "Variation_Space",         PerfMonitor::CALC);
  set_label(tm_stat_space_comm,    "Sync_Variation",          PerfMonitor::COMM);
  // end of Loop Utility Sct:1
  
  set_label(tm_loop_uty_sct_2,     "Loop_Utility_Sct_2",      PerfMonitor::CALC, false);
  set_label(tm_hstry_stdout,       "History_Stdout",          PerfMonitor::CALC);
  set_label(tm_file_out,           "File_Output",             PerfMonitor::CALC);
  set_label(tm_hstry_base,         "History_Base",            PerfMonitor::CALC);
  set_label(tm_hstry_wall,         "History_Wall_Info",       PerfMonitor::CALC);
  set_label(tm_hstry_dmfx,         "History_Domain_Flux",     PerfMonitor::CALC);
  set_label(tm_total_prs,          "Total_Pressure",          PerfMonitor::CALC);
  set_label(tm_compo_monitor,      "Component_Monitoring",    PerfMonitor::CALC);
  set_label(tm_hstry_compo,        "History_Component",       PerfMonitor::CALC);
  set_label(tm_sampling,           "Sampling",                PerfMonitor::CALC);
  set_label(tm_hstry_sampling,     "History_Sampling",        PerfMonitor::CALC);
  set_label(tm_cal_force,          "Force_Calculation",       PerfMonitor::CALC);
  set_label(tm_hstry_force,        "History_Force",           PerfMonitor::CALC);
  // end of Loop Utility Sct:2
  
  // end of Loop Utility Section
  
  
  // end of Loop Section
  
  
  // 共通にまとめて利用
  set_label(tm_copy_array,         "Copy_Array",              PerfMonitor::CALC);
  set_label(tm_assign_const,       "assign_Const_to_Array",   PerfMonitor::CALC);
  
  // 統計処理
  set_label(tm_statistic,          "Statistic",               PerfMonitor::CALC, false);
  
}


// タイムステップループの後の処理
bool FFV::stepPost() 
{
  return true;
}



// 利用例の表示
void FFV::Usage()
{
  FBUtility::printVersion(stdout, "Frontflow/violet", FFV_VERS);
  
  cout << " Usage : ";
  cout << "ffv"
  << " parameter_file" << endl;
  cout << endl;
  
  cout << " \tparameter_file includes all parameters for simulation." << endl;
  cout << endl;

}



// 空間平均操作と変動量の計算を行う
// スカラ値は算術平均，ベクトル値は自乗和
void FFV::Variation_Space(double* avr, double* rms, double& flop)
{
  double m_var[2];
  
  // 速度
  fb_delta_v_(m_var, size, &guide, d_v, d_v0, d_bcd, &flop); // 速度反復でV_res_L2_を計算している場合はスキップすること
  rms[var_Velocity] = m_var[0];
  avr[var_Velocity] = m_var[1];
  
  // 圧力
  fb_delta_s_(m_var, size, &guide, d_p, d_p0, d_bcd, &flop);
  rms[var_Pressure] = m_var[0];
  avr[var_Pressure] = m_var[1];
  
  // 温度
  if ( C.isHeatProblem() ) 
  {
    fb_delta_s_(m_var, size, &guide, d_t, d_t0, d_bcd, &flop);
    rms[var_Temperature] = m_var[0];
    avr[var_Temperature] = m_var[1];
  }
  
}

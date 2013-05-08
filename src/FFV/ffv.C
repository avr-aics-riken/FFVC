// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
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
  
  deltaT = 0.0;
  
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
  
  // Vector3D
  d_v   = NULL;
  d_vf  = NULL;
  d_vc  = NULL;
  d_v0  = NULL;
  d_wv  = NULL;
  d_abf = NULL;
  d_av  = NULL;
  d_wo  = NULL;
  d_qbc = NULL;
  
  // Scalar3D
  d_mid = NULL;
  d_bcd = NULL;
  d_bcp = NULL;
  d_bcv = NULL;
  d_bh1 = NULL;
  d_bh2 = NULL;
  
  d_p   = NULL;
  d_p0  = NULL;
  d_ws  = NULL;
  d_sq  = NULL;
  d_dv  = NULL;
  d_b   = NULL;
  d_t   = NULL;
  d_t0  = NULL;
  d_vt  = NULL;
  d_vof = NULL;
  d_ap  = NULL;
  d_at  = NULL;
  d_cvf = NULL;
  
  // Coarse initial
  d_r_v = NULL;
  d_r_p = NULL;
  d_r_t = NULL;
  
  // GMRES
  d_wg  = NULL;
  d_res = NULL;
  d_vm  = NULL;
  d_zm  = NULL;

  
  // PCG & PBiCGSTAB
	d_pcg_r = NULL;
  d_pcg_p = NULL;
  
	// PCG
	d_pcg_q = NULL;
	d_pcg_z = NULL;
  
	// PBiCGSTAB
	d_pcg_r0 = NULL;
	d_pcg_p_ = NULL;
	d_pcg_q_ = NULL;
	d_pcg_s  = NULL;
	d_pcg_s_ = NULL;
	d_pcg_t_ = NULL;
  
  compo_global_bbox = NULL;
  
  cutPos = NULL;
  cutBid = NULL;
  
  // カット情報
  d_cut = NULL;
  d_bid = NULL;
  
  // SOR2のバッファ
  cf_x = NULL;
  cf_y = NULL;
  cf_z = NULL;
}


// デストラクタ
FFV::~FFV()
{
  
}


// #################################################################
/**
 * @brief 時間平均値のファイル出力
 * @param [in,out] flop 浮動小数点演算数
 */
void FFV::AverageOutput(double& flop)
{  
  // 出力用のヘッダ
  REAL_TYPE m_org[3], m_pit[3];
  
  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++) 
  {
    m_org[i] = origin[i] - pitch[i]*(REAL_TYPE)C.GuideOut;
    m_pit[i] = pitch[i];
  }
  
  // セルセンター位置を基点とする
  for (int i=0; i<3; i++)
  {
    m_org[i] += 0.5*m_pit[i];
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
  
  // 出力ファイルの指定が有次元の場合
  double timeAvr;
  if (C.Unit.File == DIMENSIONAL)
  {
    timeAvr = CurrentTime_Avr * C.Tscale;
  }
  else
  {
    timeAvr = CurrentTime_Avr;
  }
  
  // 平均操作の母数
  unsigned stepAvr = (unsigned)CurrentStep_Avr;
  REAL_TYPE scale = 1.0;
  
  // ガイドセル出力
  int gc_out = C.GuideOut;
  
  // ファイル出力のタイムスタンプに使うステップ数
  unsigned m_step = (unsigned)CurrentStep;
  
  // ファイル出力のタイムスタンプの次元変換
  REAL_TYPE m_time;
  if (C.Unit.File == DIMENSIONAL)
  {
    m_time = (REAL_TYPE)(CurrentTime * C.Tscale);
  }
  else
  {
    m_time = (REAL_TYPE)CurrentTime;
  }
  
  // 出力ファイル名
  std::string tmp;
  std::string dtmp = DFI.GenerateDirName(C.FIO.OutDirPath, m_step, C.FIO.Slice);
  
  // 出力ディレクトリの作成
  if ( !FBUtility::mkdirs(dtmp) ) {
    Hostonly_ printf("Error : create directory \"%s\"\n", dtmp.c_str());
    Exit(-1);
  }
  
  // 出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  
  // Pressure
  if (C.Unit.File == DIMENSIONAL) 
  {
    REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
    U.prs_array_ND2D(d_ws, size, guide, d_ap, bp, C.RefDensity, C.RefVelocity, scale, flop);
  }
  else 
  {
    U.xcopy(d_ws, size, guide, d_ap, scale, kind_scalar, flop);
  }
  
  // 最大値と最小値
  REAL_TYPE f_min, f_max, min_tmp, max_tmp;
  
  fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
  
  if ( numProc > 1 )
  {
    min_tmp = f_min;
    if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    max_tmp = f_max;
    if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  
  REAL_TYPE minmax[2] = {f_min, f_max};
  
  tmp = DFI.GenerateFileName(C.f_AvrPressure, C.file_fmt_ext, m_step, myRank, mio); // e.g., prsa_0000000000_id000000.sph
  F.writeScalar(dtmp+tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out, false, stepAvr, timeAvr);
  Hostonly_ if ( !DFI.WriteDFIindex(C.f_AvrPressure,
                                    C.FIO.OutDirPath,
                                    C.file_fmt_ext,
                                    m_step,
                                    m_time,
                                    dfi_mng[var_Pressure_Avr],
                                    "ijkn",
                                    1,
                                    minmax,
                                    mio,   // 分割
                                    false, // 平均値
                                    stepAvr,
                                    timeAvr) ) Exit(0);
  
  // Velocity
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  fb_shift_refv_out_(d_wo, d_av, size, &guide, v00, &scale, &unit_velocity, &flop);
  
  fb_minmax_vex_ (&f_min, &f_max, size, &guide, v00, d_wo, &flop);
  
  if ( numProc > 1 )
  {
    min_tmp = f_min;
    if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    max_tmp = f_max;
    if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  
  minmax[0] = f_min;
  minmax[1] = f_max;
  
  tmp = DFI.GenerateFileName(C.f_AvrVelocity, C.file_fmt_ext, m_step, myRank, mio);
  F.writeVector(dtmp+tmp, size, guide, d_wo, m_step, m_time, m_org, m_pit, gc_out, false, stepAvr, timeAvr);
  Hostonly_ if ( !DFI.WriteDFIindex(C.f_AvrVelocity,
                                    C.FIO.OutDirPath,
                                    C.file_fmt_ext,
                                    m_step,
                                    m_time,
                                    dfi_mng[var_Velocity_Avr],
                                    "nijk",
                                    3,
                                    minmax,
                                    mio,
                                    false,
                                    stepAvr,
                                    timeAvr) ) Exit(0);
  
  // Temperature
  if( C.isHeatProblem() )
  {
    if (C.Unit.File == DIMENSIONAL)
    {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      U.tmp_array_ND2D(d_ws, size, guide, d_at, C.BaseTemp, C.DiffTemp, klv, scale, flop);
    }
    else 
    {
      U.xcopy(d_ws, size, guide, d_at, scale, kind_scalar, flop);
    }
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    tmp = DFI.GenerateFileName(C.f_AvrTemperature, C.file_fmt_ext, m_step, myRank, mio);
    F.writeScalar(dtmp+tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out, false, stepAvr, timeAvr);
    Hostonly_ if( !DFI.WriteDFIindex(C.f_AvrTemperature,
                                     C.FIO.OutDirPath,
                                     C.file_fmt_ext,
                                     m_step,
                                     m_time,
                                     dfi_mng[var_Temperature_Avr],
                                     "ijkn",
                                     1,
                                     minmax,
                                     mio,
                                     false,
                                     stepAvr,
                                     timeAvr) ) Exit(0);
  }
}


// #################################################################
/**
 * @brief 時間平均操作を行う
 * @param [in,out] flop 浮動小数点演算数
 */
void FFV::Averaging_Time(double& flop)
{
  CurrentStep_Avr++;
  CurrentTime_Avr += DT.get_DT();
  REAL_TYPE nadd = (REAL_TYPE)CurrentStep_Avr;
  
  fb_average_s_(d_ap, size, &guide, d_p, &nadd, &flop);
  fb_average_v_(d_av, size, &guide, d_v, &nadd, &flop);
  
  if ( C.isHeatProblem() ) 
  {
    fb_average_s_(d_at, size, &guide, d_t, &nadd, &flop);
  }
}


// #################################################################
/**
 * @brief 全ノードについて，ローカルノード1面・一層あたりの通信量の和を返す
 * @retval 通信量(Byte)
 * @param [in] sz    配列サイズ
 * @param [in] guide ガイドセル
 */
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


// #################################################################
/**
 * @brief 外部計算領域の各面における総流量と対流流出速度を計算する
 * @param [in] ptr  BoundaryOuterクラスのポインタ
 * @param [in] R    Controlクラスのポインタ
 */
void FFV::DomainMonitor(BoundaryOuter* ptr, Control* R)
{
  if ( !ptr ) Exit(0);
  BoundaryOuter* obc=NULL;
  
  obc = ptr;
  
  REAL_TYPE ddh = deltaX * deltaX;
  REAL_TYPE u_sum, u_avr;

  
  for (int face=0; face<NOFACE; face++) 
  {
    
    // 有効セル数 => 外部境界でガイドセルと内側のセルで挟まれる面がFluidの場合のセル数
    REAL_TYPE ec = (REAL_TYPE)obc[face].get_ValidCell();
    
    // 各プロセスの外部領域面の速度をvv[]にコピー
    REAL_TYPE* vv = obc[face].getDomainV();
    
    // 特殊条件
    if ( (R->Mode.Example == id_Jet) && (face==0) )
    {
      REAL_TYPE q[2] = {0.0, 0.0};
      
      // 外部境界以外はゼロ
      if ( nID[face] < 0 )
      {
        q[0] = vv[0]; // 無次元流量
        q[1] = vv[1]; // セル数
      }
      
      if ( numProc > 1 )
      {
        REAL_TYPE tmp[2] = {q[0], q[1]};
        if ( paraMngr->Allreduce(tmp, q, 2, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      }
      
      R->V_Dface[face] = q[0]/q[1];  // 無次元平均流速
      R->Q_Dface[face] = q[0] * ddh; // 無次元流量
    }
    else // 標準
    {
      // 外部境界以外はゼロにする
      u_sum = ( nID[face] < 0 ) ? vv[0] : 0.0;
      
      if ( numProc > 1 )
      {
        REAL_TYPE tmp_sum = u_sum;
        if ( paraMngr->Allreduce(&tmp_sum, &u_sum, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      }
      
      u_avr = (ec != 0.0) ? u_sum / ec : 0.0;
      
      R->V_Dface[face] = u_avr;       // 無次元平均流速
      R->Q_Dface[face] = u_sum * ddh; // 無次元流量
    }

  }
  
}


// #################################################################
/**
 * @brief ファイル出力
 * @param [in,out] flop       浮動小数点演算数
 * @param [in]     refinement 粗格子を用いたリスタート時の出力指定（trueの場合出力、default=false, ファイル名に_restart_が含まれる）
 */
void FFV::FileOutput(double& flop, const bool refinement)
{
  REAL_TYPE scale = 1.0;
  
  // d_p0をワークとして使用
  
  // 出力用のヘッダ
  REAL_TYPE m_org[3], m_pit[3];
  
  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++) 
  {
    m_org[i] = origin[i] - pitch[i]*(REAL_TYPE)C.GuideOut;
    m_pit[i] = pitch[i];
  }
  
  // セルセンター位置を基点とする
  for (int i=0; i<3; i++)
  {
    m_org[i] += 0.5*m_pit[i];
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
  unsigned m_step = (unsigned)CurrentStep;
  
  // 時間の次元変換
  REAL_TYPE m_time;
  if (C.Unit.File == DIMENSIONAL) 
  {
    m_time = (REAL_TYPE)(CurrentTime * C.Tscale);
  }
  else 
  {
    m_time = (REAL_TYPE)CurrentTime;
  }
  
  // ガイドセル出力
  int gc_out = C.GuideOut;
  
  // 出力ファイル名
  std::string tmp;
  std::string dtmp = DFI.GenerateDirName(C.FIO.OutDirPath, m_step, C.FIO.Slice);
  
  // 出力ディレクトリの作成
  if ( !FBUtility::mkdirs(dtmp) ) {
    Hostonly_ printf("Error : create directory \"%s\"\n", dtmp.c_str());
    Exit(-1);
  }
  
  // 出力モード
  bool mio = (bool)C.FIO.IOmode;

  // 最大値と最小値
  REAL_TYPE f_min, f_max, min_tmp, max_tmp;
  REAL_TYPE minmax[2];
  
  // Divergence デバッグ用なので無次元のみ
  if ( C.FIO.Div_Debug == ON ) 
  {
    
    REAL_TYPE coef = (REAL_TYPE)DT.get_DT()/(deltaX*deltaX); /// 発散値を計算するための係数　dt/h^2
    U.cnv_Div(d_ws, d_dv, size, guide, coef, flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    tmp = DFI.GenerateFileName(C.f_DivDebug, C.file_fmt_ext, m_step, myRank, mio);
    F.writeScalar(dtmp+tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.WriteDFIindex(C.f_DivDebug,
                                      C.FIO.OutDirPath,
                                      C.file_fmt_ext,
                                      m_step, m_time,
                                      dfi_mng[var_Divergence],
                                      "ijkn",
                                      1,
                                      minmax,
                                      mio) ) Exit(0);
    
  }
  
  
  // Pressure
  if (C.Unit.File == DIMENSIONAL) 
  {
    REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
    U.prs_array_ND2D(d_ws, size, guide, d_p, bp, C.RefDensity, C.RefVelocity, scale, flop);
  }
  else 
  {
    U.xcopy(d_ws, size, guide, d_p, scale, kind_scalar, flop);
  }
  
  fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
  
  if ( numProc > 1 )
  {
    min_tmp = f_min;
    if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    max_tmp = f_max;
    if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  minmax[0] = f_min;
  minmax[1] = f_max;
  
  std::string prs_restart;
  prs_restart = ( !refinement ) ? C.f_Pressure : "prs_restart_";
  
  tmp = DFI.GenerateFileName(prs_restart, C.file_fmt_ext, m_step, myRank, mio);
  DFI_OUT_PRS->WriteData(m_step, guide, m_time, d_ws, minmax);
  

  // Velocity
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  fb_shift_refv_out_(d_wo, d_v, size, &guide, v00, &scale, &unit_velocity, &flop);
  fb_minmax_vex_ (&f_min, &f_max, size, &guide, v00, d_wo, &flop);
  
  if ( numProc > 1 )
  {
    min_tmp = f_min;
    if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    max_tmp = f_max;
    if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
  }
  minmax[0] = f_min;
  minmax[1] = f_max;
  
  std::string vel_restart;
  vel_restart = ( !refinement ) ? C.f_Velocity : "vel_restart_";

  tmp = DFI.GenerateFileName(vel_restart, C.file_fmt_ext, m_step, myRank, mio);
  DFI_OUT_VEL->WriteData(m_step, guide, m_time, d_wo, minmax);
  
  
  // Tempearture
  if( C.isHeatProblem() )
  {
    
    if (C.Unit.File == DIMENSIONAL) 
    {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      U.tmp_array_ND2D(d_ws, size, guide, d_t, C.BaseTemp, C.DiffTemp, klv, scale, flop);
    }
    else 
    {
      U.xcopy(d_ws, size, guide, d_t, scale, kind_scalar, flop);
    }
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    std::string temp_restart;
    temp_restart = ( !refinement ) ? C.f_Temperature : "temp_restart_";
    
    tmp = DFI.GenerateFileName(temp_restart, C.file_fmt_ext, m_step, myRank, mio);
    F.writeScalar(dtmp+tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.WriteDFIindex(C.f_Temperature,
                                      C.FIO.OutDirPath,
                                      C.file_fmt_ext,
                                      m_step,
                                      m_time,
                                      dfi_mng[var_Temperature],
                                      "ijkn",
                                      1,
                                      minmax,
                                      mio) ) Exit(0);
  }

  
  
  // Total Pressure
  if (C.Mode.TP == ON ) {
    
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop);
    
    // convert non-dimensional to dimensional, iff file is dimensional
    if (C.Unit.File == DIMENSIONAL) 
    {
      U.tp_array_ND2D(d_ws, d_p0, size, guide, C.RefDensity, C.RefVelocity, flop);
    }
    else 
    {
      REAL_TYPE* tp;
      tp = d_ws; d_ws = d_p0; d_p0 = tp;
    }
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    tmp = DFI.GenerateFileName(C.f_TotalP, C.file_fmt_ext, m_step, myRank, mio);
    F.writeScalar(dtmp+tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.WriteDFIindex(C.f_TotalP,
                                      C.FIO.OutDirPath,
                                      C.file_fmt_ext,
                                      m_step,
                                      m_time,
                                      dfi_mng[var_TotalP],
                                      "ijkn",
                                      1,
                                      minmax,
                                      mio) ) Exit(0);
  }
  
  
  // Vorticity
  if (C.Mode.VRT == ON ) 
  {
    
    rot_v_(d_wv, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    REAL_TYPE  vz[3];
    vz[0] = vz[1] = vz[2] = 0.0;
    unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity/C.RefLength : 1.0;
    fb_shift_refv_out_(d_wo, d_wv, size, &guide, vz, &scale, &unit_velocity, &flop);
    fb_minmax_vex_ (&f_min, &f_max, size, &guide, v00, d_wo, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    tmp = DFI.GenerateFileName(C.f_Vorticity, C.file_fmt_ext, m_step, myRank, mio);
    F.writeVector(dtmp+tmp, size, guide, d_wo, m_step, m_time, m_org, m_pit, gc_out, 1);
    
    Hostonly_ if ( !DFI.WriteDFIindex(C.f_Vorticity,
                                      C.FIO.OutDirPath,
                                      C.file_fmt_ext,
                                      m_step,
                                      m_time,
                                      dfi_mng[var_Vorticity],
                                      "nijk",
                                      3,
                                      minmax,
                                      mio) ) Exit(0);
  }
  
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C.Mode.I2VGT == ON ) {
    
    i2vgt_ (d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    // 無次元で出力
    U.xcopy(d_ws, size, guide, d_p0, scale, kind_scalar, flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    tmp = DFI.GenerateFileName(C.f_I2VGT, C.file_fmt_ext, m_step, myRank, mio);
    F.writeScalar(dtmp+tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.WriteDFIindex(C.f_I2VGT,
                                      C.FIO.OutDirPath,
                                      C.file_fmt_ext,
                                      m_step,
                                      m_time,
                                      dfi_mng[var_I2vgt],
                                      "ijkn",
                                      1,
                                      minmax,
                                      mio) ) Exit(0);
  }
  
  
  // Helicity
  if (C.Mode.Helicity == ON ) 
  {
    helicity_(d_p0, size, &guide, &deltaX, d_v, d_bcv, v00, &flop);
    
    // 無次元で出力
    U.xcopy(d_ws, size, guide, d_p0, scale, kind_scalar, flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    tmp = DFI.GenerateFileName(C.f_Helicity, C.file_fmt_ext, m_step, myRank, mio);
    F.writeScalar(dtmp+tmp, size, guide, d_ws, m_step, m_time, m_org, m_pit, gc_out);
    
    Hostonly_ if ( !DFI.WriteDFIindex(C.f_Helicity,
                                      C.FIO.OutDirPath,
                                      C.file_fmt_ext,
                                      m_step,
                                      m_time,
                                      dfi_mng[var_Helicity],
                                      "ijkn",
                                      1,
                                      minmax,
                                      mio) ) Exit(0);
  }

  
  // Face Velocity
  if (C.Mode.FaceV == ON )
  {
    fb_shift_refv_out_(d_wo, d_vf, size, &guide, v00, &scale, &unit_velocity, &flop);
    fb_minmax_vex_ (&f_min, &f_max, size, &guide, v00, d_wo, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    std::string fvel_restart;
    fvel_restart = ( !refinement ) ? C.f_Fvelocity : "fvel_restart_";
    
    tmp = DFI.GenerateFileName(fvel_restart, C.file_fmt_ext, m_step, myRank, mio);
    DFI_OUT_FVEL->WriteData(m_step, guide, m_time, d_wo, minmax);
    
  }
  
  
}


// #################################################################
/**
 * @brief シミュレーションの1ステップの処理
 * @note Loop() + stepPost()
 */
int FFV::MainLoop()
{
  int ret = 1;
  
  for (int i=1; i<=Session_LastStep; i++)
  {
    if ( FFV_TerminateCtrl::getTerminateFlag() )
    {
      return 0; // forced terminate
      break;
    }
    
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



// #################################################################
/**
 * @brief VP反復の発散値を計算する
 * @param [in] IC ItrCtlクラス
 * @retval 発散値の最大の場所のインデクス
 */
FB::Vec3i FFV::Norm_Div(ItrCtl* IC)
{
  double nrm;
  double flop_count;
  REAL_TYPE coef = 1.0/deltaX; /// 発散値を計算するための係数

  int index[3];
  index[0] = 0;
  index[1] = 0;
  index[2] = 0;
  
  switch (IC->get_normType())
  {

    case ItrCtl::v_div_max:
      TIMING_start(tm_norm_div_max);
      flop_count=0.0;
      norm_v_div_max_(&nrm, size, &guide, d_dv, &coef, d_bcp, &flop_count);
      TIMING_stop(tm_norm_div_max, flop_count);
      
      if ( numProc > 1 )
      {
        TIMING_start(tm_norm_comm);
        double tmp;
        tmp = nrm;
        if ( paraMngr->Allreduce(&tmp, &nrm, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0); // 最大値
        TIMING_stop(tm_norm_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
      }
      IC->set_normValue(nrm);
      break;
      
      
    case ItrCtl::v_div_dbg:
      TIMING_start(tm_poi_itr_sct_5); // >>> Poisson Iteration subsection 5
      
      TIMING_start(tm_norm_div_dbg);
      flop_count=0.0;
      
      norm_v_div_dbg_(&nrm, index, size, &guide, d_dv, &coef, d_bcp, &flop_count);
      TIMING_stop(tm_norm_div_dbg, flop_count);
      
      //@todo ここで，最大値のグローバルなindexの位置を計算する
      
      if ( numProc > 1 )
      {
        TIMING_start(tm_norm_comm);
        double tmp;
        tmp = nrm;
        if ( paraMngr->Allreduce(&tmp, &nrm, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0); // 最大値
      }
      IC->set_normValue(nrm);
      
      
      TIMING_stop(tm_poi_itr_sct_5, 0.0); // <<< Poisson Iteration subsection 5
      break;
      
      
    default:
      stamped_printf("\tInvalid convergence type\n");
      Exit(0);
  }
  
  FB::Vec3i idx ( index[0], index[1], index[2] );
  
  return idx;
}



// #################################################################
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


// #################################################################
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
  set_label(tm_norm_div_dbg,       "Poisson_Norm_Div_dbg",    PerfMonitor::CALC);
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
  
  set_label(tm_gmres_sor_sct,      "Poisson__GMRES_Sct",      PerfMonitor::CALC, false);
  set_label(tm_gmres_mvprod,       "GMRES_MatVec_Product",    PerfMonitor::CALC);
  set_label(tm_gmres_res_sample,   "GMRES_Res_sample",        PerfMonitor::CALC);
  set_label(tm_gmres_others,       "GMRES_Others",            PerfMonitor::CALC);
  set_label(tm_gmres_comm,         "A_R_GMRES",               PerfMonitor::COMM);
  
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


// #################################################################
// タイムステップループの後の処理
bool FFV::stepPost() 
{
  return true;
}


// #################################################################
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


// #################################################################
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

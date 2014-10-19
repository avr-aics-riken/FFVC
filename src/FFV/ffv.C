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
 * @file   ffv.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"

// コンストラクタ
FFV::FFV()
{
  EXEC_MODE = -1;
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

  Session_CurrentStep = 0;
  Session_LastStep = 0;
  CurrentStep = 0;
  CurrentStep_Avr = 0;
  
  CM_F.previous = 0.0;
  CM_F.rate     = 0.0;
  CM_H.previous = 0.0;
  CM_H.rate     = 0.0;
  
  deltaT = 0.0;
  poly_factor = 0.0;
  
  
  for (int i=0; i<3; i++) 
  {
    G_size[i]= 0;
    G_origin[i] = 0.0;
    G_region[i] = 0.0;
    poly_org[i] = 0.0;
    poly_dx[i] = 0.0;
    poly_gc[i] = 0;
    ensPeriodic[i] = 0;
  }
  
  mat_tbl = NULL;
  
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
  d_cdf = NULL;
  
  d_p   = NULL;
  d_p0  = NULL;
  d_ws  = NULL;
  d_sq  = NULL;
  d_dv  = NULL;
  d_b   = NULL;
  d_ie  = NULL;
  d_ie0 = NULL;
  d_vt  = NULL;
  d_vof = NULL;
  d_ap  = NULL;
  d_ae  = NULL;
  d_cvf = NULL;
  d_pvf = NULL;
  d_vrt = NULL;
  
  // Coarse initial
  d_r_v = NULL;
  d_r_p = NULL;
  d_r_t = NULL;
  
  // GMRES
  d_wg  = NULL;
  d_res = NULL;
  d_vm  = NULL;
  d_zm  = NULL;

  
  // PCG & BiCGSTAB
	d_pcg_r = NULL;
  d_pcg_p = NULL;
  
	// PCG
	d_pcg_z = NULL;
  
	// BiCGSTAB
	d_pcg_r0 = NULL;
	d_pcg_q  = NULL;
	d_pcg_s  = NULL;
	d_pcg_t  = NULL;
  
  // BiCGSTAB with Preconditioning
  d_pcg_p_ = NULL;
  d_pcg_s_ = NULL;
  d_pcg_t_ = NULL;
  
  cutPos = NULL;
  cutBid = NULL;
  
  // カット情報
  d_cut = NULL;
  d_bid = NULL;
  
  // SOR2のバッファ
  cf_x = NULL;
  cf_y = NULL;
  cf_z = NULL;
  
  
  // Poisson係数（Naive実装の実験）
  d_pni = NULL;
  
  
  // OBSTACLEの力の積算
  local_force = NULL;
  global_force = NULL;
  global_obstacle = NULL;
  num_obstacle = 0;
  buffer_force = NULL;
  
  // 発散の収束判定
  DivC.MaxIteration = 0;
  DivC.Iteration = 0;
  DivC.divType = 0;
  DivC.divEPS = 0.0;
  DivC.divergence = 0.0;
  
  
  // ファイル入出力
  DFI_IN_PRS   = NULL;
  DFI_IN_VEL   = NULL;
  DFI_IN_FVEL  = NULL;
  DFI_IN_TEMP  = NULL;
  DFI_IN_PRSA  = NULL;
  DFI_IN_VELA  = NULL;
  DFI_IN_TEMPA = NULL;
  DFI_OUT_PRS  = NULL;
  DFI_OUT_VEL  = NULL;
  DFI_OUT_FVEL = NULL;
  DFI_OUT_TEMP = NULL;
  DFI_OUT_PRSA = NULL;
  DFI_OUT_VELA = NULL;
  DFI_OUT_TEMPA= NULL;
  DFI_OUT_TP   = NULL;
  DFI_OUT_VRT  = NULL;
  DFI_OUT_I2VGT= NULL;
  DFI_OUT_HLT  = NULL;
  DFI_OUT_DIV  = NULL;
}



// デストラクタ
FFV::~FFV()
{
  if( DFI_IN_PRS    != NULL ) delete DFI_IN_PRS;
  if( DFI_IN_VEL    != NULL ) delete DFI_IN_VEL;
  if( DFI_IN_FVEL   != NULL ) delete DFI_IN_FVEL;
  if( DFI_IN_TEMP   != NULL ) delete DFI_IN_TEMP;
  if( DFI_IN_PRSA   != NULL ) delete DFI_IN_PRSA;
  if( DFI_IN_VELA   != NULL ) delete DFI_IN_VELA;
  if( DFI_IN_TEMPA  != NULL ) delete DFI_IN_TEMPA;
  if( DFI_OUT_PRS   != NULL ) delete DFI_OUT_PRS;
  if( DFI_OUT_VEL   != NULL ) delete DFI_OUT_VEL;
  if( DFI_OUT_FVEL  != NULL ) delete DFI_OUT_FVEL;
  if( DFI_OUT_TEMP  != NULL ) delete DFI_OUT_TEMP;
  if( DFI_OUT_PRSA  != NULL ) delete DFI_OUT_PRSA;
  if( DFI_OUT_VELA  != NULL ) delete DFI_OUT_VELA;
  if( DFI_OUT_TEMPA != NULL ) delete DFI_OUT_TEMPA;
  if( DFI_OUT_TP    != NULL ) delete DFI_OUT_TP;
  if( DFI_OUT_VRT   != NULL ) delete DFI_OUT_VRT;
  if( DFI_OUT_I2VGT != NULL ) delete DFI_OUT_I2VGT;
  if( DFI_OUT_HLT   != NULL ) delete DFI_OUT_HLT;
  if( DFI_OUT_DIV   != NULL ) delete DFI_OUT_DIV;
}



// #################################################################
/**
 * @brief 時間平均操作を行う
 * @param [in,out] flop 浮動小数点演算数
 */
void FFV::Averaging(double& flop)
{
  CurrentStep_Avr++;
  CurrentTime_Avr += DT.get_DT();
  REAL_TYPE nadd = (REAL_TYPE)CurrentStep_Avr;
  
  fb_average_s_(d_ap, size, &guide, d_p, &nadd, &flop);
  fb_average_v_(d_av, size, &guide, d_v, &nadd, &flop);
  
  if ( C.isHeatProblem() ) 
  {
    fb_average_s_(d_ae, size, &guide, d_ie, &nadd, &flop);
  }
}



// #################################################################
// 物体に働く力の計算
void FFV::calcForce(double& flop)
{
  REAL_TYPE dh = deltaX;
  int gd = guide;
  int st[3], ed[3];
  REAL_TYPE vec[3];
  
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].getType()==OBSTACLE )
    {
      cmp[n].getBbox(st, ed);
      int key = cmp[n].getMatodr();

      force_compo_(vec, size, &gd, &key, d_p, d_bid, &dh, st, ed, &flop);
      
      local_force[3*n+0] = vec[0];
      local_force[3*n+1] = vec[1];
      local_force[3*n+2] = vec[2];
    }
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
        case X_minus:
        case X_plus:
          c += (double)(sz[1]*sz[2]);
          break;
          
        case Y_minus:
        case Y_plus:
          c += (double)(sz[0]*sz[2]);
          break;
          
        case Z_minus:
        case Z_plus:
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
  
  return c * sizeof(REAL_TYPE); // Byte
}


// #################################################################
/**
 * @brief 外部計算領域の各面における総流量と対流流出速度を計算する
 * @param [in] ptr  BoundaryOuterクラスのポインタ
 * @param [in] R    Controlクラスのポインタ
 * @note 系への流入を正の符号とする
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
    // 外部境界面でない場合にはゼロが戻る
    u_sum = 0.0;
    vobc_face_massflow_(&u_sum, size, &guide, &face, d_vf, d_cdf, nID);
    
    
    // 有効セル数 => 外部境界でガイドセルと内側のセルで挟まれる面がFluidの場合のセル数
    REAL_TYPE ec = (REAL_TYPE)obc[face].getValidCell();
    
    // 各プロセスの外部領域面の速度をvv[]にコピー
    REAL_TYPE* vv = obc[face].getDomainV();
    
    REAL_TYPE q[2] = {0.0, 0.0};
    
    // 外部境界のみ値をもつ
    if ( nID[face] < 0 )
    {
      q[0] = u_sum; // 無次元流量
      q[1] = vv[1]; // セル数
    }
    
    if ( numProc > 1 )
    {
      REAL_TYPE tmp[2] = {q[0], q[1]};
      if ( paraMngr->Allreduce(tmp, q, 2, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    
    // 特殊条件
    if ( (R->Mode.Example == id_Jet) && (face==X_minus) )
    {
      R->V_Dface[face] = q[0]/q[1];  // 無次元平均流速
      R->Q_Dface[face] = q[0] * ddh; // 無次元流量
    }
    else // 標準
    {
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
 * @brief 外部計算領域の各面における総流量と対流流出速度を計算する
 * @param [in] ptr  BoundaryOuterクラスのポインタ
 * @param [in] R    Controlクラスのポインタ
 * @note 系への流入を正の符号とする
 *
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
    REAL_TYPE ec = (REAL_TYPE)obc[face].getValidCell();
    
    // 各プロセスの外部領域面の速度をvv[]にコピー
    REAL_TYPE* vv = obc[face].getDomainV();
    
    REAL_TYPE q[2] = {0.0, 0.0};
    
    // 外部境界のみ値をもつ
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
    
    
    // 特殊条件
    if ( (R->Mode.Example == id_Jet) && (face==X_minus) )
    {
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
  
}*/



// #################################################################
// 物体に働く力を各ランクから集めて積算する
void FFV::gatherForce(REAL_TYPE* m_frc)
{
  
  if( !m_frc ) Exit(0);
  
  for (int n=1; n<=C.NoCompo; n++)
  {
    
    if ( cmp[n].getType()==OBSTACLE )
    {
      REAL_TYPE vec[3];
      
      vec[0] = local_force[3*n+0];
      vec[1] = local_force[3*n+1];
      vec[2] = local_force[3*n+2];
      
      if ( numProc > 1 )
      {
        if ( paraMngr->Gather(vec, 3, m_frc, 3, 0) != CPM_SUCCESS ) Exit(0);
      }
      else
      {
        memcpy(m_frc, vec, 3*sizeof(REAL_TYPE));
      }
    }

    REAL_TYPE fx = 0.0;
    REAL_TYPE fy = 0.0;
    REAL_TYPE fz = 0.0;
    
    for (int i=0; i<numProc; i++)
    {
      fx = fx + m_frc[3*i+0];
      fy = fy + m_frc[3*i+1];
      fz = fz + m_frc[3*i+2];
    }
    
    global_force[3*n+0] = fx;
    global_force[3*n+1] = fy;
    global_force[3*n+2] = fz;
  }
}


// #################################################################
/**
 * @brief 時間平均値のファイル出力
 * @param [in,out] flop 浮動小数点演算数
 */
void FFV::OutputAveragedVarables(double& flop)
{
  REAL_TYPE f_min, f_max, min_tmp, max_tmp, vec_min[4], vec_max[4];
  REAL_TYPE minmax[2];
  REAL_TYPE cdm_minmax[8];
  
  // エラーコード
  CDM::E_CDM_ERRORCODE ret;
  
  
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
  
  
  if ( C.KindOfSolver != SOLID_CONDUCTION )
  {
    // Pressure
    if (C.Unit.File == DIMENSIONAL)
    {
      REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
      U.convArrayPrsND2D(d_ws, size, guide, d_ap, bp, C.RefDensity, C.RefVelocity, flop);
    }
    else
    {
      U.copyS3D(d_ws, size, guide, d_ap, scale);
    }
    
    // 最大値と最小値
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[0] = f_min;
    minmax[1] = f_max;

    
    if ( !DFI_OUT_PRSA )
    {
      printf("[%d] DFI_OUT_PRSA Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }

    ret = DFI_OUT_PRSA->WriteData(m_step,   // 出力step番号
                                  m_time,   // 出力時刻
                                  size,     // d_wsの実ボクセル数
                                  1,        // 成分数
                                  guide,    // 仮想セル数
                                  d_ws,     // フィールドデータポインタ
                                  minmax,   // 最小値と最大値
                                  false,    // 平均出力指示 false:出力あり
                                  stepAvr,  // 平均をとったステップ数
                                  timeAvr); // 平均をとった時刻

    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
    
    
    
    // Velocity
    REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
    
    if ( DFI_OUT_VELA->GetArrayShape() == CDM::E_CDM_NIJK ) // Velocityの型は CDM::E_CDM_NIJK
    {
      fb_vout_nijk_(d_wo, d_av, size, &guide, v00, &unit_velocity, &flop); // 配列並びを変換
      fb_minmax_vex_ (vec_min, vec_max, size, &guide, v00, d_wo, &flop);
    }
    else 
    {
      fb_vout_ijkn_(d_wo, d_av, size, &guide, v00, &unit_velocity, &flop); // 並び変換なし
      fb_minmax_v_ (vec_min, vec_max, size, &guide, v00, d_wo, &flop);
    }
    
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
      if ( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
      if ( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    
    if ( !DFI_OUT_VELA )
    {
      printf("[%d] DFI_OUT_VELA Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    cdm_minmax[0] = vec_min[1]; ///<<< vec_u min
    cdm_minmax[1] = vec_max[1]; ///<<< vec_u max
    cdm_minmax[2] = vec_min[2]; ///<<< vec_v min
    cdm_minmax[3] = vec_max[2]; ///<<< vec_v max
    cdm_minmax[4] = vec_min[3]; ///<<< vec_w min
    cdm_minmax[5] = vec_max[3]; ///<<< vec_w max
    cdm_minmax[6] = vec_min[0]; ///<<< u,v,wの合成値のmin
    cdm_minmax[7] = vec_max[0]; ///<<< u,v,wの合成値のmax
    
    DFI_OUT_VELA->setComponentVariable(0, "u");
    DFI_OUT_VELA->setComponentVariable(1, "v");
    DFI_OUT_VELA->setComponentVariable(2, "w");

    ret = DFI_OUT_VELA->WriteData(m_step,
                                  m_time,
                                  size,
                                  3,
                                  guide,
                                  d_wo,
                                  cdm_minmax,
                                  false,
                                  stepAvr,
                                  timeAvr);

    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }

  
  // Temperature
  if ( C.isHeatProblem() )
  {
    U.convArrayIE2Tmp(d_ws, size, guide, d_ae, d_bcd, mat_tbl, C.BaseTemp, C.DiffTemp, C.Unit.File, flop);
    
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_ws, &flop);
    
    if ( numProc > 1 )
    {
      min_tmp = f_min;
      if ( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if ( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    minmax[0] = f_min;
    minmax[1] = f_max;
    
    if ( !DFI_OUT_TEMPA )
    {
      printf("[%d] DFI_OUT_TEMPA Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    ret = DFI_OUT_TEMPA->WriteData(m_step,
                                   m_time,
                                   size,
                                   1,
                                   guide,
                                   d_ws,
                                   minmax,
                                   false,
                                   stepAvr,
                                   timeAvr);

    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
}



// #################################################################
/**
 * @brief 基本変数のファイル出力
 * @param [in,out] flop       浮動小数点演算数
 * @note d_p0をワークとして使用
 */
void FFV::OutputBasicVariables(double& flop)
{
  REAL_TYPE scale = 1.0;
  
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
  

  // 最大値と最小値
  REAL_TYPE f_min, f_max, min_tmp, max_tmp, vec_min[4], vec_max[4];
  REAL_TYPE minmax[2];
  REAL_TYPE cdm_minmax[8];
  
  
  // エラーコード
  CDM::E_CDM_ERRORCODE ret;
  
  
  
  if ( C.KindOfSolver != SOLID_CONDUCTION )
  {
    // Pressure
    if (C.Unit.File == DIMENSIONAL)
    {
      REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
      U.convArrayPrsND2D(d_ws, size, guide, d_p, bp, C.RefDensity, C.RefVelocity, flop);
    }
    else
    {
      U.copyS3D(d_ws, size, guide, d_p, scale);
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
    
    if( !DFI_OUT_PRS )
    {
      printf("[%d] DFI_OUT_PRS Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    ret = DFI_OUT_PRS->WriteData(m_step,
                                 m_time,
                                 size,
                                 1,
                                 guide,
                                 d_ws,
                                 minmax,
                                 true,
                                 0,
                                 0.0);

    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
    
    
    // Velocity
    REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
    
    if ( DFI_OUT_VEL->GetArrayShape() == CDM::E_CDM_NIJK ) // Velocityの型は CDM::E_CDM_NIJK
    {
      fb_vout_nijk_(d_wo, d_v, size, &guide, v00, &unit_velocity, &flop);
      fb_minmax_vex_ (vec_min, vec_max, size, &guide, v00, d_wo, &flop);
    }
    else 
    {
      fb_vout_ijkn_(d_wo, d_v, size, &guide, v00, &unit_velocity, &flop);
      fb_minmax_v_ (vec_min, vec_max, size, &guide, v00, d_wo, &flop);
    }
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( !DFI_OUT_VEL )
    {
      printf("[%d] DFI_OUT_VEL Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    cdm_minmax[0] = vec_min[1]; ///<<< vec_u min
    cdm_minmax[1] = vec_max[1]; ///<<< vec_u max
    cdm_minmax[2] = vec_min[2]; ///<<< vec_v min
    cdm_minmax[3] = vec_max[2]; ///<<< vec_v max
    cdm_minmax[4] = vec_min[3]; ///<<< vec_w min
    cdm_minmax[5] = vec_max[3]; ///<<< vec_w max
    cdm_minmax[6] = vec_min[0]; ///<<< u,v,wの合成値のmin
    cdm_minmax[7] = vec_max[0]; ///<<< u,v,wの合成値のmax
    
    DFI_OUT_VEL->setComponentVariable(0, "u");
    DFI_OUT_VEL->setComponentVariable(1, "v");
    DFI_OUT_VEL->setComponentVariable(2, "w");
    
    ret = DFI_OUT_VEL->WriteData(m_step,
                                 m_time,
                                 size,
                                 3,
                                 guide,
                                 d_wo,
                                 cdm_minmax,
                                 true,
                                 0,
                                 0.0);

    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
    
    
    // Face Velocity
    if ( DFI_OUT_VEL->GetArrayShape() == CDM::E_CDM_NIJK ) // FVelocityの型は CDM::E_CDM_NIJK
    {
      fb_vout_nijk_(d_wo, d_vf, size, &guide, v00, &unit_velocity, &flop);
      fb_minmax_vex_ (vec_min, vec_max, size, &guide, v00, d_wo, &flop);
    }
    else 
    {
      fb_vout_ijkn_(d_wo, d_vf, size, &guide, v00, &unit_velocity, &flop);
      fb_minmax_v_ (vec_min, vec_max, size, &guide, v00, d_wo, &flop);
    }
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    cdm_minmax[0] = vec_min[1]; ///<<< vec_u min
    cdm_minmax[1] = vec_max[1]; ///<<< vec_u max
    cdm_minmax[2] = vec_min[2]; ///<<< vec_v min
    cdm_minmax[3] = vec_max[2]; ///<<< vec_v max
    cdm_minmax[4] = vec_min[3]; ///<<< vec_w min
    cdm_minmax[5] = vec_max[3]; ///<<< vec_w max
    cdm_minmax[6] = vec_min[0]; ///<<< u,v,wの合成値のmin
    cdm_minmax[7] = vec_max[0]; ///<<< u,v,wの合成値のmax
    
    DFI_OUT_FVEL->setComponentVariable(0, "u");
    DFI_OUT_FVEL->setComponentVariable(1, "v");
    DFI_OUT_FVEL->setComponentVariable(2, "w");
    
    ret = DFI_OUT_FVEL->WriteData(m_step,
                                  m_time,
                                  size,
                                  3,
                                  guide,
                                  d_wo,
                                  cdm_minmax,
                                  true,
                                  0,
                                  0.0);
    
    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  
  if ( !C.isHeatProblem() ) return;
  
  // Tempearture
  U.convArrayIE2Tmp(d_ws, size, guide, d_ie, d_bcd, mat_tbl, C.BaseTemp, C.DiffTemp, C.Unit.File, flop);
  
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
  
  if ( !DFI_OUT_TEMP )
  {
    printf("[%d] DFI_OUT_TEMP Pointer Error\n", paraMngr->GetMyRankID());
    Exit(-1);
  }
  
  ret = DFI_OUT_TEMP->WriteData(m_step,
                                m_time,
                                size,
                                1,
                                guide,
                                d_ws,
                                minmax,
                                true,
                                0,
                                0.0);
  
  if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  
}


// #################################################################
/**
 * @brief 派生変数のファイル出力
 * @param [in,out] flop       浮動小数点演算数
 * @note d_p0をワークとして使用
 */
void FFV::OutputDerivedVariables(double& flop)
{
  REAL_TYPE scale = 1.0;
  
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
  
  
  // 最大値と最小値
  REAL_TYPE f_min, f_max, min_tmp, max_tmp, vec_min[4], vec_max[4];
  REAL_TYPE minmax[2];
  REAL_TYPE cdm_minmax[8];
  
  
  // エラーコード
  CDM::E_CDM_ERRORCODE ret;
  
  
  REAL_TYPE unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
  
  
  // Total Pressure
  if (C.varState[var_TotalP] == ON )
  {
    fb_totalp_ (d_p0, size, &guide, d_v, d_p, v00, &flop);
    
    // convert non-dimensional to dimensional, iff file is dimensional
    if (C.Unit.File == DIMENSIONAL)
    {
      U.convArrayTpND2D(d_ws, d_p0, size, guide, C.RefDensity, C.RefVelocity);
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
    
    if ( !DFI_OUT_TP )
    {
      printf("[%d] DFI_OUT_TP Pointer Error\n",paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    ret = DFI_OUT_TP->WriteData(m_step,
                                m_time,
                                size,
                                1,
                                guide,
                                d_ws,
                                minmax,
                                true,
                                0,
                                0.0);

    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // Vorticity
  if (C.varState[var_Vorticity] == ON )
  {
    rot_v_(d_wv, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    REAL_TYPE  vz[3];
    vz[0] = vz[1] = vz[2] = 0.0;
    unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity/C.RefLength : 1.0;
    
    if ( DFI_OUT_VRT->GetArrayShape() == CDM::E_CDM_NIJK ) // Vorticityの型は CDM::E_CDM_NIJK
    {
      fb_vout_nijk_(d_wo, d_wv, size, &guide, vz, &unit_velocity, &flop);
      fb_minmax_vex_ (vec_min, vec_max, size, &guide, v00, d_wo, &flop);
    }
    else 
    {
      fb_vout_ijkn_(d_wo, d_wv, size, &guide, vz, &unit_velocity, &flop);
      fb_minmax_v_ (vec_min, vec_max, size, &guide, v00, d_wo, &flop);
    }
    
    if ( numProc > 1 )
    {
      REAL_TYPE vmin_tmp[4] = {vec_min[0], vec_min[1], vec_min[2], vec_min[3]};
      if( paraMngr->Allreduce(vmin_tmp, vec_min, 4, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE vmax_tmp[4] = {vec_max[0], vec_max[1], vec_max[2], vec_max[3]};
      if( paraMngr->Allreduce(vmax_tmp, vec_max, 4, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( !DFI_OUT_VRT )
    {
      printf("[%d] DFI_OUT_VRT Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    cdm_minmax[0] = vec_min[1]; ///<<< vec_u min
    cdm_minmax[1] = vec_max[1]; ///<<< vec_u max
    cdm_minmax[2] = vec_min[2]; ///<<< vec_v min
    cdm_minmax[3] = vec_max[2]; ///<<< vec_v max
    cdm_minmax[4] = vec_min[3]; ///<<< vec_w min
    cdm_minmax[5] = vec_max[3]; ///<<< vec_w max
    cdm_minmax[6] = vec_min[0]; ///<<< u,v,wの合成値のmin
    cdm_minmax[7] = vec_max[0]; ///<<< u,v,wの合成値のmax

    ret = DFI_OUT_VRT->WriteData(m_step,
                                 m_time,
                                 size,
                                 3,
                                 guide,
                                 d_wo,
                                 cdm_minmax,
                                 true,
                                 0,
                                 0.0);

    if ( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C.varState[var_Qcr] == ON )
  {
    i2vgt_ (d_p0, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    // 無次元で出力
    U.copyS3D(d_ws, size, guide, d_p0, scale);
    
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
    
    if ( !DFI_OUT_I2VGT )
    {
      printf("[%d] DFI_OUT_I2VGT Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    ret = DFI_OUT_I2VGT->WriteData(m_step,
                                   m_time,
                                   size,
                                   1,
                                   guide,
                                   d_ws,
                                   minmax,
                                   true,
                                   0,
                                   0.0);

    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  
  // Helicity
  if (C.varState[var_Helicity] == ON )
  {
    helicity_(d_p0, size, &guide, &deltaX, d_v, d_cdf, v00, &flop);
    
    // 無次元で出力
    U.copyS3D(d_ws, size, guide, d_p0, scale);
    
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
    
    if ( !DFI_OUT_HLT )
    {
      printf("[%d] DFI_OUT_HLT Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    ret = DFI_OUT_HLT->WriteData(m_step,
                                 m_time,
                                 size,
                                 1,
                                 guide,
                                 d_ws,
                                 minmax,
                                 true,
                                 0,
                                 0.0);

    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
  // Divergence for Debug
  if (C.varState[var_Div] == ON )
  {
    REAL_TYPE coef = 1.0/deltaX; /// 発散値を計算するための係数　1/h
    U.cnv_Div(d_ws, d_dv, size, guide, coef);
    
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
    
    if( !DFI_OUT_DIV )
    {
      printf("[%d] DFI_OUT_DIV Pointer Error\n", paraMngr->GetMyRankID());
      Exit(-1);
    }
    
    ret = DFI_OUT_DIV->WriteData(m_step,
                                 m_time,
                                 size,
                                 1,
                                 guide,
                                 d_ws,
                                 minmax,
                                 true,
                                 0,
                                 0.0);
    
    if( ret != CDM::E_CDM_SUCCESS ) Exit(0);
  }
  
}


// #################################################################
/**
 * @brief シミュレーションの1ステップの処理
 * @note Loop() + stepPost()
 */
int FFV::MainLoop()
{
  //>> Graph Ploter
  C.Interval[Control::tg_compute].printInfo("tg_compute");
  C.Interval[Control::tg_console].printInfo("tg_console");
  C.Interval[Control::tg_history].printInfo("tg_history");
  C.Interval[Control::tg_basic].printInfo("tg_basic");
  C.Interval[Control::tg_average].printInfo("tg_average");
  C.Interval[Control::tg_derived].printInfo("tg_derived");
  C.Interval[Control::tg_accelra].printInfo("tg_accelra");
  C.Interval[Control::tg_sampled].printInfo("tg_sampled");
  C.Interval[Control::tg_END].printInfo("tg_END");
  //<< Graph Ploter
  
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
        ret = -1;
        break;
        
      case 0: // forced terminated
        ret = 1;
        break;
        
      case 1: // normal
        ret = 1;
        break;
        
      default:
        ret = -1;
    }
    
    if ( loop_ret == 0 ) break;
  }
  
  
  // サンプリングファイルのクローズ
  MO.closeFile();

  if ( fp_b ) fclose(fp_b);  ///< 基本情報
  if ( fp_w ) fclose(fp_w);  ///< 壁面情報
  if ( fp_c ) fclose(fp_c);  ///< コンポーネント情報
  if ( fp_d ) fclose(fp_d);  ///< 流量収支情報
  if ( fp_i ) fclose(fp_i);  ///< 反復履歴情報
  if ( fp_f ) fclose(fp_f);  ///< 力の履歴情報

  
  if ( !stepPost() ) ret = -1;

  return ret;
}



// #################################################################
/**
 * @brief 発散値を計算する
 * @param [in] LSd  線形ソルバクラスのDiv反復
 */
void FFV::NormDiv()
{
  REAL_TYPE dv;
  double flop_count, tmp;
  REAL_TYPE coef = 1.0/deltaX; /// 発散値を計算するための係数

  
  TIMING_start(tm_poi_itr_sct_5); // >>> Poisson Iteration subsection 5
  
  if ( DivC.divType == nrm_div_max )
  {
    TIMING_start(tm_norm_div_max);
    flop_count=0.0;
    norm_v_div_max_(&dv, size, &guide, d_dv, &coef, d_bcp, &flop_count);
    TIMING_stop(tm_norm_div_max, flop_count);
    
    if ( numProc > 1 )
    {
      TIMING_start(tm_norm_comm);
      REAL_TYPE tmp = dv;
      if ( paraMngr->Allreduce(&tmp, &dv, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0); // 最大値
      TIMING_stop(tm_norm_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
    }
    DivC.divergence = (double)dv;
  }
  else // nrm_div_l2
  {
    TIMING_start(tm_norm_div_max);
    flop_count=0.0;
    norm_v_div_l2_(&dv, size, &guide, d_dv, &coef, d_bcp, &flop_count);
    TIMING_stop(tm_norm_div_max, flop_count);
    
    if ( numProc > 1 )
    {
      TIMING_start(tm_norm_comm);
      REAL_TYPE tmp = dv;
      if ( paraMngr->Allreduce(&tmp, &dv, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0); // 自乗和
      TIMING_stop(tm_norm_comm, 2.0*numProc*sizeof(double));
    }
    tmp = (double)dv;
    DivC.divergence = sqrt( tmp );
  }
  
  TIMING_stop(tm_poi_itr_sct_5, 0.0); // <<< Poisson Iteration subsection 5

}



// #################################################################
/**
 * @brief タイミング測定区間にラベルを与えるラッパー
 * @param [in] key       キー番号
 * @param [in] label     ラベル
 * @param [in] type      測定対象タイプ(COMM or CALC)
 * @param [in] exclusive 排他測定フラグ(ディフォルトtrue)
 */
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
/**
 * @brief タイミング測定区間にラベルを与える
 */
void FFV::set_timing_label()
{
  // ラベルの設定
  set_label(tm_init_sct,           "Initialization_Section",  PerfMonitor::CALC, false);
  set_label(tm_init_alloc,         "Allocate_Arrays",         PerfMonitor::CALC);
  
  set_label(tm_voxel_prep_sct,     "Voxel_Prep_Section",      PerfMonitor::CALC, false);
  set_label(tm_voxel_load,         "Loading_Voxel_File",      PerfMonitor::CALC);
  set_label(tm_polygon_load,       "Loading_Polygon_File",    PerfMonitor::CALC);
  set_label(tm_cutinfo,            "Cut_Information",         PerfMonitor::CALC);
  set_label(tm_cut_min,            "Cut_Minimum_search",      PerfMonitor::CALC);
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
  set_label(tm_norm_div_max,       "Norm_Div_max",            PerfMonitor::CALC);
  set_label(tm_norm_comm,          "A_R_Poisson_Norm",        PerfMonitor::COMM);
  // end of Poisson: Itr. Sct:5
  
  // end of Poisson: Iteration Sct.
  
  set_label(tm_NS_loop_post_sct,   "NS__Loop_Post_Section",   PerfMonitor::CALC, false);
  set_label(tm_vectors_comm,       "Sync_Velocity",           PerfMonitor::COMM);
  set_label(tm_domain_monitor,     "Domain_Monitor",          PerfMonitor::CALC);
  set_label(tm_VBC_update,         "Velocity_BC_Update",      PerfMonitor::CALC);
  set_label(tm_LES_eddy,           "Eddy_Viscosity",          PerfMonitor::CALC);
  set_label(tm_LES_eddy_comm,      "Sync_Eddy_Viscosity",     PerfMonitor::COMM);
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
  set_label(tm_average_time,       "Averaging",               PerfMonitor::CALC);
  set_label(tm_stat_space,         "Variation_Space",         PerfMonitor::CALC);
  set_label(tm_stat_space_comm,    "Sync_Variation",          PerfMonitor::COMM);
  // end of Loop Utility Sct:1
  
  set_label(tm_loop_uty_sct_2,     "Loop_Utility_Sct_2",      PerfMonitor::CALC, false);
  set_label(tm_hstry_stdout,       "History_Stdout",          PerfMonitor::CALC);
  set_label(tm_file_out,           "File_Output",             PerfMonitor::CALC);
  set_label(tm_hstry_base,         "History_Base",            PerfMonitor::CALC, false);
  set_label(tm_cal_force,          "Force_Calculation",       PerfMonitor::CALC);
  set_label(tm_hstry_wall,         "History_Wall_Info",       PerfMonitor::CALC);
  set_label(tm_total_prs,          "Total_Pressure",          PerfMonitor::CALC);
  set_label(tm_sampling,           "Sampling",                PerfMonitor::CALC);
  set_label(tm_hstry_sampling,     "History_Sampling",        PerfMonitor::CALC);
  
  // end of Loop Utility Sct:2
  
  // end of Loop Utility Section
  
  
  // end of Loop Section
  
  
  // 共通にまとめて利用
  set_label(tm_copy_array,         "Copy_Array",              PerfMonitor::CALC);
  set_label(tm_assign_const,       "assign_Const_to_Array",   PerfMonitor::CALC);
  
  // Blas
  set_label(tm_blas_dot1,          "Blas_Dot1",               PerfMonitor::CALC);
  set_label(tm_blas_dot2,          "Blas_Dot2",               PerfMonitor::CALC);
  set_label(tm_blas_copy,          "Blas_Copy",               PerfMonitor::CALC);
  set_label(tm_blas_calc_r,        "Blas_Residual",           PerfMonitor::CALC);
  set_label(tm_blas_ax,            "Blas_Ax",                 PerfMonitor::CALC);
  set_label(tm_blas_z_xpay,        "Blas_Z=X+aY" ,            PerfMonitor::CALC);
  set_label(tm_blas_axpbypz,       "Blas_Z=Z+aX+bY",          PerfMonitor::CALC);
  set_label(tm_blas_bicg_update_p, "Blas_BiCGupdateP",        PerfMonitor::CALC);
  set_label(tm_blas_bicg_update_x, "Blas_BiCGupdateX",        PerfMonitor::CALC);
  set_label(tm_blas_comm,          "Blas_Comm",               PerfMonitor::COMM);
  
  set_label(tm_bicgstab_sct,       "BiCG",                    PerfMonitor::CALC, false);
  
  
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
/**
 * @brief コマンドラインヘルプ
 */
void FFV::Usage()
{
  FBUtility::printVersion(stdout, "Frontflow/violet", FFVC_VERSION_NO);
  
  cout << " Usage : ";
  cout << "ffv"
  << " parameter_file" << endl;
  cout << endl;
  
  cout << " \tparameter_file includes all parameters for simulation." << endl;
  cout << endl;

}


// #################################################################
/**
 * @brief 空間平均操作と変動量の計算を行う
 * @param [out]    rms  変動値の自乗和
 * @param [out]    avr  平均値の自乗和
 * @param [in,out] flop 浮動小数演算数
 */
void FFV::VariationSpace(double* rms, double* avr, double& flop)
{
  double m_var[2];
  
  // 速度
  fb_delta_v_(m_var, size, &guide, d_v, d_v0, d_bcd, &flop); // 速度反復でV_res_L2_を計算している場合はスキップすること
  rms[var_Velocity] = m_var[0];
  //avr[var_Velocity] = m_var[1]; 意味を持たない
  
  // 圧力
  fb_delta_s_(m_var, size, &guide, d_p, d_p0, d_bcd, &flop);
  rms[var_Pressure] = m_var[0];
  avr[var_Pressure] = m_var[1];
  
  // 温度
  if ( C.isHeatProblem() ) 
  {
    fb_delta_s_(m_var, size, &guide, d_ie, d_ie0, d_bcd, &flop);
    rms[var_Temperature] = m_var[0];
    avr[var_Temperature] = m_var[1];
  }
  
}

// #################################################################
/**
 * @brief sphファイルの書き出し（内部領域のみ）
 * @param [in] vf               スカラデータ
 * @param [in] sz               配列サイズ
 * @param [in] gc               ガイドセル
 * @param [in] gc_out           出力するガイドセル数
 * @param [in] org              基点
 * @param [in] ddx              ピッチ
 * @param [in] m_ModePrecision  浮動小数点の精度
 * @note 標記上，long 対応になっているが，ファイルフォーマットとの対応を確認のこと
 */
void FFV::writeRawSPH(const REAL_TYPE *vf, const int* sz, const int gc, const int gc_out, const REAL_TYPE* org, const REAL_TYPE* ddx, const int m_ModePrecision)
{
  int pad, dType, stp, svType;
  REAL_TYPE ox, oy, oz, dx, dy, dz, tm;
  long long szl[3], stpl;
  
  
  char sph_fname[512];
  
  if ( paraMngr->IsParallel() )
  {
    sprintf( sph_fname, "field%010d.sph", paraMngr->GetMyRankID() );
  }
  else
  {
    sprintf( sph_fname, "field.sph" );
  }
  
  ofstream ofs(sph_fname, ios::out | ios::binary);
  if (!ofs)
  {
    printf("\tCan't open %s file\n", sph_fname);
    Exit(0);
  }
  
  int m_sz[3];
  m_sz[0] = sz[0]+2*gc_out;
  m_sz[1] = sz[1]+2*gc_out;
  m_sz[2] = sz[2]+2*gc_out;
  int gd = gc;
  
  size_t nx = m_sz[0] * m_sz[1] * m_sz[2];
  
  ox = org[0]-ddx[0]*(REAL_TYPE)gc_out;
  oy = org[1]-ddx[1]*(REAL_TYPE)gc_out;
  oz = org[2]-ddx[2]*(REAL_TYPE)gc_out;
  dx = ddx[0];
  dy = ddx[1];
  dz = ddx[2];
  //printf("org: %f %f %f\n", ox, oy, oz);
  //printf("dx : %f %f %f\n", dx, dy, dz);
  
  svType = kind_scalar;
  if ( sizeof(REAL_TYPE) == sizeof(double) )
  {
    for (int i=0; i<3; i++)   szl[i] = (long long)m_sz[i];
  }
  
  REAL_TYPE *f = new REAL_TYPE[nx];
  
  size_t m, l;
  
  for (int k=1; k<=m_sz[2]; k++) {
    for (int j=1; j<=m_sz[1]; j++) {
      for (int i=1; i<=m_sz[0]; i++) {
        l = _F_IDX_S3D(i, j, k, m_sz[0], m_sz[1], m_sz[2], gc_out);
        m = _F_IDX_S3D(i, j, k, m_sz[0], m_sz[1], m_sz[2], gd);
        f[l] = (REAL_TYPE)vf[m];
      }
    }
  }
  
  // data property
  ( m_ModePrecision == sizeof(float) ) ? dType=1 : dType=2;
  pad = sizeof(int)*2;
  ofs.write( (char*)&pad, sizeof(int) );
  ofs.write( (char*)&svType, sizeof(int) );
  ofs.write( (char*)&dType, sizeof(int) );
  ofs.write( (char*)&pad, sizeof(int) );
  
  // voxel size
  if (dType == 1) {
    pad = sizeof(int)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&m_sz[0], sizeof(int) );
    ofs.write( (char*)&m_sz[1], sizeof(int) );
    ofs.write( (char*)&m_sz[2], sizeof(int) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(long long)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&szl[0], sizeof(long long) );
    ofs.write( (char*)&szl[1], sizeof(long long) );
    ofs.write( (char*)&szl[2], sizeof(long long) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // original point of domain
  if (dType == 1) {
    pad = sizeof(float)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(float) );
    ofs.write( (char*)&oy, sizeof(float) );
    ofs.write( (char*)&oz, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(double)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(double) );
    ofs.write( (char*)&oy, sizeof(double) );
    ofs.write( (char*)&oz, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // pitch of voxel
  if (dType == 1) {
    pad = sizeof(float)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(float) );
    ofs.write( (char*)&dy, sizeof(float) );
    ofs.write( (char*)&dz, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(double)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(double) );
    ofs.write( (char*)&dy, sizeof(double) );
    ofs.write( (char*)&dz, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // time stamp
  if (dType == 1) {
    stp = 0;
    tm = 0.0;
    pad = sizeof(int)+sizeof(float);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&stp, sizeof(int) );
    ofs.write( (char*)&tm, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    stpl =0;
    tm = 0.0;
    pad = sizeof(long long)+sizeof(double);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&stpl, sizeof(long long) );
    ofs.write( (char*)&tm, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  if (svType == kind_scalar) {
    int pad = (m_ModePrecision == sizeof(float)) ? nx * sizeof(float) : nx * sizeof(double);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)f,   pad );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  ofs.close();
  
  if (f) { delete [] f; f=NULL; }
}

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
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
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
  order_of_PM_key = 0;
  
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
  CurrentTimeStat = 0.0;

  Session_CurrentStep = 0;
  Session_LastStep = 0;
  CurrentStep = 0;
  CurrentStepStat = 0;
  
  CM_F.previous = 0.0;
  CM_F.rate     = 0.0;
  CM_H.previous = 0.0;
  CM_H.rate     = 0.0;
  
  deltaT = 0.0;
  poly_factor = 0.0;
  
  
  for (int i=0; i<3; i++) 
  {
    ensPeriodic[i] = 0;
  }
  
  mat_tbl = NULL;
  vec_tbl = NULL;
  
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
  
  
  // OBSTACLEの力の積算
  cmp_force_local = NULL;
  cmp_force_global = NULL;
  cmp_force_avr = NULL;
  global_obstacle = NULL;
  num_obstacle = 0;
  buffer_force = NULL;
  
  // 発散の収束判定
  DivC.MaxIteration = 0;
  DivC.Iteration = 0;
  DivC.divType = 0;
  DivC.divEPS = 0.0;
  DivC.divergence = 0.0;
}



// デストラクタ
FFV::~FFV() {}



// #################################################################
/**
 * @brief 時間平均操作を行う
 * @param [in,out] flop 浮動小数点演算数
 */
void FFV::Averaging(double& flop)
{
  CurrentStepStat++;
  CurrentTimeStat += DT.get_DT();
  REAL_TYPE nadd = (REAL_TYPE)CurrentStepStat;
  
  fb_average_s_(d_ap, size, &guide, d_p, &nadd, &flop);
  fb_average_v_(d_av, size, &guide, d_v, &nadd, &flop);
  
  if ( C.isHeatProblem() ) 
  {
    fb_average_s_(d_ae, size, &guide, d_ie, &nadd, &flop);
  }
}



// #################################################################
// 物体に働く力を計算し、各ランクから集めて積算する
void FFV::calcForce(double& flop)
{
  REAL_TYPE dh = deltaX;
  int gd = guide;
  int st[3], ed[3];
  REAL_TYPE vec[3];
  
  
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].getType()==OBSTACLE || cmp[n].getType()==SOLIDREV )
    {
      cmp[n].getBbox(st, ed);
      int key = cmp[n].getMatodr();

      // 力の計算
      force_compo_(vec, size, &gd, &key, d_p, d_bid, &dh, st, ed, &flop);
      
      cmp_force_local[3*n+0] = vec[0];
      cmp_force_local[3*n+1] = vec[1];
      cmp_force_local[3*n+2] = vec[2];
      
      
      // 集約
      if ( numProc > 1 )
      {
        if ( paraMngr->Gather(vec, 3, buffer_force, 3, 0) != CPM_SUCCESS ) Exit(0);
      }
      else
      {
        buffer_force[0] = vec[0];
        buffer_force[1] = vec[1];
        buffer_force[2] = vec[2];
      }
      
      REAL_TYPE fx = 0.0;
      REAL_TYPE fy = 0.0;
      REAL_TYPE fz = 0.0;
      
      for (int i=0; i<numProc; i++)
      {
        fx = fx + buffer_force[3*i+0];
        fy = fy + buffer_force[3*i+1];
        fz = fz + buffer_force[3*i+2];
      }
      
      cmp_force_global[3*n+0] = fx;
      cmp_force_global[3*n+1] = fy;
      cmp_force_global[3*n+2] = fz;
      
      
      // average マスターノードのみ、有効な値
      if ( C.Mode.Statistic == ON && C.Interval[Control::tg_statistic].isStarted(CurrentStep, CurrentTime) )
      {
        REAL_TYPE c2 = 1.0 / (REAL_TYPE)CurrentStepStat;
        REAL_TYPE c1 = 1.0 - c2;
        cmp_force_avr[3*n+0] = c1 * cmp_force_avr[3*n+0] + c2 * fx;
        cmp_force_avr[3*n+1] = c1 * cmp_force_avr[3*n+1] + c2 * fy;
        cmp_force_avr[3*n+2] = c1 * cmp_force_avr[3*n+2] + c2 * fz;
      }
      
    }
  }
  
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
  C.Interval[Control::tg_statistic].printInfo("tg_statistic");
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
 * @param [in] div  \sum{u}
 */
void FFV::NormDiv(REAL_TYPE* div)
{
  REAL_TYPE dv;
  double flop_count, tmp;
  REAL_TYPE coef = 1.0/deltaX; /// 発散値を計算するための係数

  
  if ( DivC.divType == nrm_div_max )
  {
    TIMING_start("Norm_Div_max");
    flop_count=0.0;
    norm_v_div_max_(&dv, size, &guide, div, &coef, d_bcp, &flop_count);
    TIMING_stop("Norm_Div_max", flop_count);
    
    if ( numProc > 1 )
    {
      TIMING_start("All_Reduce");
      REAL_TYPE tmp = dv;
      if ( paraMngr->Allreduce(&tmp, &dv, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0); // 最大値
      TIMING_stop("All_Reduce", 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
    }
    DivC.divergence = (double)dv;
  }
  else // nrm_div_l2
  {
    TIMING_start("Norm_Div_L2");
    flop_count=0.0;
    norm_v_div_l2_(&dv, size, &guide, div, &coef, d_bcp, &flop_count);
    TIMING_stop("Norm_Div_L2", flop_count);
    
    if ( numProc > 1 )
    {
      TIMING_start("All_Reduce");
      REAL_TYPE tmp = dv;
      if ( paraMngr->Allreduce(&tmp, &dv, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0); // 自乗和
      TIMING_stop("All_Reduce", 2.0*numProc*sizeof(double));
    }
    tmp = (double)dv;
    DivC.divergence = sqrt( tmp );
  }
  


}



void FFV::printCriteria(FILE* fp)
{
  // Criteria ------------------
  if ( EXEC_MODE == ffvc_solver )
  {
    fprintf(fp,"\n\tParameter of Linear Equation\n");
    LinearSolver* ICp1= &LS[ic_prs1];  /// 圧力のPoisson反復
    LinearSolver* ICp2= &LS[ic_prs2];  /// 圧力のPoisson反復　2回目
    LinearSolver* ICv = &LS[ic_vel1];  /// 粘性項のCrank-Nicolson反復
    LinearSolver* ICt = &LS[ic_tmp1];  /// 温度の拡散項の反復
    
    if ( C.Hide.PM_Test == ON )
    {
      fprintf(fp,"\t ### Performance Test Mode >> The iteration number is fixed by Iteration max.\n\n");
    }
    
    if ( C.KindOfSolver != SOLID_CONDUCTION )
    {
      // 1st iteration
      fprintf(fp,"\t     1st Pressure Iteration \n");
      printIteratoinParameter(fp, ICp1);
      
      if ( C.AlgorithmF == Flow_FS_RK_CN )
      {
        fprintf(fp,"\t     2nd Pressure Iteration \n");
        printIteratoinParameter(fp, ICp2);
      }
      
      // CN iteration
      if ( (C.AlgorithmF == Flow_FS_AB_CN) || (C.AlgorithmF == Flow_FS_RK_CN) )
      {
        fprintf(fp,"\n");
        fprintf(fp,"\t     Velocity CN Iteration \n");
        printIteratoinParameter(fp, ICv);
      }
      
      // Divergence
      fprintf(fp,"\n");
      fprintf(fp,"\t     Div Iteration \n");
      fprintf(fp,"\t       Iteration max          :   %d\n"  ,  DivC.MaxIteration);
      if ( DivC.divType == nrm_div_max)
      {
        fprintf(fp,"\t       Error    Norm type     :   Max divergence\n");
      }
      else
      {
        fprintf(fp,"\t       Error    Norm type     :   L2 divergence\n");
      }
      
      fprintf(fp,"\t       Threshold for Div.     :   %9.3e\n", DivC.divEPS);
    }
    
    // for Temperature
    if ( C.isHeatProblem() )
    {
      if ( C.AlgorithmH == Heat_EE_EI )
      {
        fprintf(fp,"\n");
        fprintf(fp,"\t     Temperature Iteration  \n");
        printIteratoinParameter(fp, ICt);
      }
    }
    
  } // End of Criteria
}


// #################################################################
/**
 * @brief 線形ソルバー種別のパラメータ表示
 * @param [in] fp ファイルポインタ
 * @param [in] IC LinearSOlver
 */
void FFV::printIteratoinParameter(FILE* fp, LinearSolver* IC)
{
  switch (IC->getLS())
  {
    case JACOBI:
      fprintf(fp,"\t       Linear Solver          :   Jacobi method\n");
      break;
      
    case SOR:
      fprintf(fp,"\t       Linear Solver          :   Point SOR method\n");
      break;
      
    case SOR2SMA:
      if (IC->getNaive()==OFF)
      {
        fprintf(fp,"\t       Linear Solver          :   2-colored SOR SMA (Stride Memory Access, Bit compressed 1-decode)\n");
      }
      else
      {
        fprintf(fp,"\t       Linear Solver          :   2-colored SOR SMA (Stride Memory Access, Naive Implementation)\n");
      }
      break;
      
    case GMRES:
      fprintf(fp,"\t       Linear Solver          :   GMRES\n");
      break;
      
    case PCG:
      fprintf(fp,"\t       Linear Solver          :   PCG\n");
      break;
      
    case BiCGSTAB:
      if (IC->getNaive()==OFF)
      {
        fprintf(fp,"\t       Linear Solver          :   BiCGstab");
        if ( IC->isPreconditioned() ) fprintf(fp," with Preconditioner\n");
      }
      else
      {
        fprintf(fp,"\t       Linear Solver          :   BiCGstab (Naive)");
        if ( IC->isPreconditioned() ) fprintf(fp," with Preconditioner\n");
      }
      break;
      
    default:
      stamped_printf("Error: Linear Solver section\n");
  }
  
  fprintf(fp,"\t       Iteration max          :   %d\n"  ,  IC->getMaxIteration());
  fprintf(fp,"\t       Residual Norm type     :   %s\n",    IC->getResNormString().c_str());
  fprintf(fp,"\t       Threshold for residual :   %9.3e\n", IC->getResCriterion());
  fprintf(fp,"\t       Error    Norm type     :   %s\n",    IC->getErrNormString().c_str());
  fprintf(fp,"\t       Threshold for error    :   %9.3e\n", IC->getErrCriterion());
  
  switch (IC->getLS())
  {
    case JACOBI:
      fprintf(fp,"\t       Coef. of Relaxation    :   %9.3e\n", IC->getOmega());
      break;
      
    case SOR:
      fprintf(fp,"\t       Coef. of Acceleration  :   %9.3e\n", IC->getOmega());
      break;
      
    case SOR2SMA:
      fprintf(fp,"\t       Coef. of Acceleration  :   %9.3e\n", IC->getOmega());
      fprintf(fp,"\t       Communication Mode     :   %s\n",   (IC->getSyncMode()==comm_sync) ? "SYNC" : "ASYNC");
      break;
      
    case GMRES:
      fprintf(fp,"\t       Linear Solver          :   GMRES\n");
      break;
      
    case PCG:
      fprintf(fp,"\t       Linear Solver          :   PCG\n");
      break;
      
    case BiCGSTAB:
      if (IC->isPreconditioned() == true)
      {
        fprintf(fp,"\t       Inner Iteration        :   %d\n"  ,  IC->getInnerItr());
        fprintf(fp,"\t       Coef. of Acceleration  :   %9.3e\n", IC->getOmega());
        fprintf(fp,"\t       Communication Mode     :   %s\n",   (IC->getSyncMode()==comm_sync) ? "SYNC" : "ASYNC");
      }
      break;
      
    default:
      stamped_printf("Error: Linear Solver section\n");
  }
  
}



// #################################################################
/**
 * @brief タイミング測定区間にラベルを与えるラッパー
 * @param [in] label     ラベル
 * @param [in] type      測定対象タイプ(COMM or CALC)
 * @param [in] exclusive 排他測定フラグ(ディフォルトtrue)
 */
void FFV::set_label(const string label, PerfMonitor::Type type, bool exclusive)
{
  // 登録個数のチェック
  order_of_PM_key++;
  
  if ( order_of_PM_key > PM_NUM_MAX )
  {
    Hostonly_ printf("\tThe number of labels for Performance monitor goes over limit.\n");
    Exit(0);
  }
  
  // 文字数がTM_LABEL_MAX-1を超えるものはカット
  if ( strlen(label.c_str()) > TM_LABEL_MAX-1 )
  {
    printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
  }
  
  // Performance Monitorへの登録
  PM.setProperties(label, type, exclusive);
}


// #################################################################
/**
 * @brief タイミング測定区間にラベルを与える
 */
void FFV::set_timing_label()
{
  // common
  set_label("Allocate_Arrays",         PerfMonitor::CALC);
  set_label("Copy_Array",              PerfMonitor::CALC);
  set_label("assign_Const_to_Array",   PerfMonitor::CALC);
  set_label("All_Reduce",              PerfMonitor::COMM);
  
  
  
  set_label("Initialization_Section",  PerfMonitor::CALC, false);
  
  set_label("Voxel_Prep_Section",      PerfMonitor::CALC, false);
  
  set_label("Polylib_Section",         PerfMonitor::CALC, false);
  set_label("Loading_Polygon_File",    PerfMonitor::CALC);
  set_label("Write_Polygon_File",      PerfMonitor::CALC);
  // Polylib_Section
  
  set_label("Cut_Section",             PerfMonitor::CALC, false);
  set_label("Cut_Information",         PerfMonitor::CALC);
  set_label("Cut_Minimum_search",      PerfMonitor::CALC);
  // Cut_Section
  
  set_label("Geometry_Section",        PerfMonitor::CALC, false);
  set_label("SeedFilling",             PerfMonitor::CALC);
  set_label("SubSampling",             PerfMonitor::CALC);
  set_label("Fill",                    PerfMonitor::CALC);
  // Geometry_Section
  
  set_label("Compo_Vertex8",           PerfMonitor::CALC);
  set_label("Compo_Subdivision",       PerfMonitor::CALC);
  set_label("Compo_Fraction",          PerfMonitor::CALC);
  
  set_label("Encode_BCindex",          PerfMonitor::CALC);
  
  set_label("Gather_DomainInfo",       PerfMonitor::CALC);
  
  set_label("Generate_Glyph",          PerfMonitor::CALC);
  
  // Voxel_Prep_Section
  // Initialization_Section
  
  
  set_label("Restart_Process",         PerfMonitor::CALC);
  
  
  set_label("Time_Step_Loop_Section",  PerfMonitor::CALC, false);
  
  set_label("Search_Vmax",             PerfMonitor::CALC);

  set_label("Flow_Section",            PerfMonitor::CALC, false);
  
  set_label("NS__F_Step_Section",      PerfMonitor::CALC, false);
  set_label("Pvec_MUSCL",              PerfMonitor::CALC);
  set_label("Pvec_MUSCL_LES",          PerfMonitor::CALC);
  set_label("Pvec_Central",            PerfMonitor::CALC);
  set_label("Pvec_Central_LES",        PerfMonitor::CALC);
  set_label("Pvec_Flux_BC",            PerfMonitor::CALC);
  set_label("Pvec_Euler_Explicit",     PerfMonitor::CALC);
  set_label("Pvec_Adams_Bashforth",    PerfMonitor::CALC);
  set_label("Pvec_AB_CN",              PerfMonitor::CALC);
  set_label("Pvec_Forcing",            PerfMonitor::CALC);
  set_label("Pvec_Buoyancy",           PerfMonitor::CALC);
  set_label("Pvec_BC",                 PerfMonitor::CALC);
  set_label("Sync_Pvec",               PerfMonitor::COMM);
  // NS__F_Step_Section
  
  
  set_label("Poisson__Source_Section", PerfMonitor::CALC, false);
  set_label("Divergence_of_Pvec",      PerfMonitor::CALC);
  set_label("Poisson_Src_VBC",         PerfMonitor::CALC);
  set_label("Poisson_Src_Norm",        PerfMonitor::CALC);
  set_label("A_R_Poisson_Src_L2",      PerfMonitor::COMM);
  set_label("Poisson_Init_Res",        PerfMonitor::CALC);
  set_label("A_R_Poisson_Init_Res_L2", PerfMonitor::COMM);
  // Poisson__Source_Section
  
  
  set_label("VP-Iteration_Section",    PerfMonitor::CALC, false);
  set_label("Point_SOR",               PerfMonitor::CALC, false);
  set_label("2-colored_SOR_stride",    PerfMonitor::CALC, false);
  set_label("PBiCGstab",               PerfMonitor::CALC, false);
  set_label("Projection_Velocity",     PerfMonitor::CALC);
  set_label("Projection_Velocity_BC",  PerfMonitor::CALC);
  set_label("A_R_Projection_VBC",      PerfMonitor::COMM);
  set_label("Projection_Forcing",      PerfMonitor::CALC);
  set_label("A_R_Projection_Forcing",  PerfMonitor::COMM);
  set_label("Forcing Source",          PerfMonitor::CALC);
  set_label("Sync_Face_Velocity",      PerfMonitor::COMM);
  set_label("Velocity_BC",             PerfMonitor::CALC);
  set_label("Norm_Div_max",            PerfMonitor::CALC);
  set_label("Norm_Div_L2",             PerfMonitor::CALC);
  // VP-Iteration_Section

  
  set_label("NS__Loop_Post_Section",   PerfMonitor::CALC, false);
  set_label("Sync_Velocity",           PerfMonitor::COMM);
  set_label("Domain_Monitor",          PerfMonitor::CALC);
  // NS__Loop_Post_Section
  
  // Flow_Section
  
  
  set_label("Heat_Section",            PerfMonitor::CALC, false);
  
  set_label("Thermal_Convection",      PerfMonitor::CALC);
  set_label("Thermal_Convection_BC",   PerfMonitor::CALC);
  set_label("Thermal_Convection_EE",   PerfMonitor::CALC);
  set_label("Thermal_Diff_Outer_BC",   PerfMonitor::CALC);
  set_label("Thermal_Diff_IBC_Vol",    PerfMonitor::CALC);
  set_label("Sync_Thermal",            PerfMonitor::COMM);
  set_label("Thermal_Diff_IBC_Face",   PerfMonitor::CALC);
  set_label("Thermal_Diff_OBC_Face",   PerfMonitor::CALC);
  set_label("Sync_Thermal_QBC",        PerfMonitor::COMM);
  set_label("Thermal_Diff_EE",         PerfMonitor::CALC);
  set_label("A_R_Thermal_Diff_Res",    PerfMonitor::COMM);
  set_label("Thermal_Range_Cut",       PerfMonitor::CALC);
  
  
  set_label("Thermal_Diff_PSOR",       PerfMonitor::CALC); //
  
  
  
  
  
  set_label("VOF_Section",             PerfMonitor::CALC, false);
  set_label("VOF_Convection",          PerfMonitor::CALC);
  set_label("Sync_VOF_Convection",     PerfMonitor::COMM);
  
  
  set_label("Loop_Utility_Section",    PerfMonitor::CALC, false);
  set_label("Averaging",               PerfMonitor::CALC);
  set_label("Turbulence Statistic",    PerfMonitor::CALC);
  set_label("Variation_Space",         PerfMonitor::CALC);
  set_label("A_R_variation_space",     PerfMonitor::COMM);
  set_label("File_Output",             PerfMonitor::CALC);
  set_label("Total_Pressure",          PerfMonitor::CALC);
  set_label("Sampling",                PerfMonitor::CALC);
  set_label("History_out",             PerfMonitor::CALC);
  set_label("Force_Calculation",       PerfMonitor::CALC);
  // Loop_Utility_Section
  
  // Time_Step_Loop_Section
 
  set_label("Statistic",               PerfMonitor::CALC, false);
  
  
  // LinearSolver class
  set_label("A_R_Convergence",         PerfMonitor::COMM);
  set_label("Dot",                     PerfMonitor::CALC);
  set_label("Dot1",                    PerfMonitor::CALC);
  set_label("Dot2",                    PerfMonitor::CALC);
  set_label("A_R_Dot",                 PerfMonitor::COMM);
  set_label("Poisson_PSOR",            PerfMonitor::CALC);
  set_label("Poisson_BC",              PerfMonitor::CALC);
  set_label("Sync_Poisson",            PerfMonitor::COMM);
  set_label("Poisson_SOR2_SMA",        PerfMonitor::CALC);
  set_label("Blas_Clear",              PerfMonitor::CALC);
  set_label("Blas_Copy",               PerfMonitor::CALC);
  set_label("Blas_Residual",           PerfMonitor::CALC);
  set_label("Blas_BiCG_1" ,            PerfMonitor::CALC);
  set_label("Blas_BiCG_2",             PerfMonitor::CALC);
  set_label("Blas_AX",                 PerfMonitor::CALC);
  set_label("Blas_TRIAD",              PerfMonitor::CALC);

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


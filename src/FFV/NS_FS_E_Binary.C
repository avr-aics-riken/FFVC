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
 * @file   NS_FS_E_Binary.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"


// Fractional Step法でNavier-Stokes方程式を解く．バイナリ近似．
void FFV::NS_FS_E_Binary()
{
  // local variables
  double flop;                         /// 浮動小数演算数
  double b_l2 = 0.0;                   /// 反復解法での定数項ベクトルのL2ノルム
  double res0_l2 = 0.0;                /// 反復解法での初期残差ベクトルのL2ノルム
  
  REAL_TYPE dt = deltaT;               /// 時間積分幅
  REAL_TYPE Re = C.Reynolds;           /// レイノルズ数
  REAL_TYPE rei = C.getRcpReynolds();  /// レイノルズ数の逆数
  REAL_TYPE half = 0.5;                /// 定数
  REAL_TYPE one = 1.0;                 /// 定数
  REAL_TYPE zero = 0.0;                /// 定数
  int cnv_scheme = C.CnvScheme;        /// 対流項スキーム
  

  // 境界処理用
  Gemini_R* m_buf = new Gemini_R [C.NoCompo+1];
  REAL_TYPE* m_snd = new REAL_TYPE [(C.NoCompo+1)*2];
  REAL_TYPE* m_rcv = new REAL_TYPE [(C.NoCompo+1)*2];

  
  LinearSolver* LSp = &LS[ic_prs1];  /// 圧力のPoisson反復
  LinearSolver* LSv = &LS[ic_vel1];  /// 粘性項のCrank-Nicolson反復

  // #### タイムステップ間で保持されるデータ ####
  // d_v   セルセンタ速度 v^n -> v^{n+1}
  // d_vf  セルフェイス速度
  // d_p   圧力 p^n -> p^{n+1}
  // d_dv  発散値, div(u)の値を保持
  // d_bcd IDのビットフラグ
  // d_bcp 圧力のビットフラグ
  // d_cdf Component Directional BC Flag
  // d_ie0  内部エネルギー t^n
  
  // #### タイムステップ間で保持されるないテンポラリ ####
  // d_v0  セルセンタ速度 v^nの保持
  // d_vc  疑似速度ベクトル
  // d_wv  陰解法の時の疑似速度ベクトル，射影ステップの境界条件
  // d_p0  圧力 p^nの保持
  // d_ws  Poissonのソース項0　速度境界を考慮
  // d_sq  Poissonのソース項1　反復毎に変化するソース項，摩擦速度，
  // d_b   反復の右辺ベクトル
  // d_cvf コンポーネントの体積率
  // d_vt  LES計算の渦粘性係数
  // d_ab0 Adams-Bashforth用のワーク
  
  
  // >>> Fractional step section
  TIMING_start("NS__F_Step_Section");
  
  
  // n stepの値を保持 >> In use (d_v0, d_p0)
  TIMING_start("Copy_Array");
  U.copyS3D(d_p0, size, guide, d_p, one);
  U.copyV3D(d_v0, size, guide, d_v, one);
  TIMING_stop("Copy_Array", 0.0, 2);
  
  
  
  // 対流項と粘性項の評価 >> In use (d_vc, d_wv)
  switch (C.AlgorithmF)
  {
    case Flow_FS_EE_EE:
    case Flow_FS_AB2:
      switch ( cnv_scheme )
      {
        case Control::O1_upwind:
        case Control::O3_muscl:
          if ( C.LES.Calc == ON )
          {
            TIMING_start("Pvec_MUSCL_LES");
            flop = 0.0;
            pvec_muscl_les_ (d_vc, size, &guide, pitch, &cnv_scheme, v00, &rei, d_v0, d_vf, d_cdf, d_bid, &one, &C.LES.Cs, &C.LES.Model, &C.RefKviscosity, &C.RefDensity, &flop);
            TIMING_stop("Pvec_MUSCL_LES", flop);
          }
          else
          {
            TIMING_start("Pvec_MUSCL");
            flop = 0.0;
            pvec_muscl_(d_vc, size, &guide, pitch, &cnv_scheme, v00, &rei, d_v0, d_vf, d_cdf, d_bid, &one, &flop);
            TIMING_stop("Pvec_MUSCL", flop);
          }
          break;
          
        case Control::O2_central:
        case Control::O4_central:
          if ( C.LES.Calc == ON )
          {
            TIMING_start("Pvec_Central_LES");
            flop = 0.0;
            pvec_central_les_(d_vc, size, &guide, pitch, &cnv_scheme, v00, &rei, d_v0, d_vf, d_cdf, d_bid, &one, &C.LES.Cs, &C.LES.Model, &C.RefKviscosity, &C.RefDensity, &flop);
            TIMING_stop("Pvec_Central_LES", flop);
          }
          else
          {
            TIMING_start("Pvec_Central");
            flop = 0.0;
            pvec_central_(d_vc, size, &guide, pitch, &cnv_scheme, v00, &rei, d_v0, d_vf, d_cdf, d_bid, &one, &flop);
            TIMING_stop("Pvec_Central", flop);
          }
          break;
      }

      TIMING_start("Pvec_Flux_BC");
      flop = 0.0;
      BC.modPvecFlux(d_vc, d_v0, d_cdf, CurrentTime, &C, v00, flop);
      TIMING_stop("Pvec_Flux_BC", flop);
      break;
      
      
    case Flow_FS_AB_CN:
      switch ( cnv_scheme )
      {
        case Control::O1_upwind:
        case Control::O3_muscl:
          if ( C.LES.Calc == ON )
          {
            TIMING_start("Pvec_MUSCL_LES");
            flop = 0.0;
            //pvec_les_(wv, sz, &guide, dh, (int*)&C.CnvScheme, v00, &rei, v0, vf, (int*)bcv, vt, &flop);
            TIMING_stop("Pvec_MUSCL_LES", flop);
          }
            else
          {
            TIMING_start("Pvec_MUSCL");
            flop = 0.0;
            pvec_muscl_(d_wv, size, &guide, pitch, &cnv_scheme, v00, &rei, d_v0, d_vf, d_cdf, d_bid, &half, &flop);
            TIMING_stop("Pvec_MUSCL", flop);
          }
          break;
        
        case Control::O2_central:
        case Control::O4_central:
          if ( C.LES.Calc == ON )
          {
            TIMING_start("Pvec_Central_LES");
            flop = 0.0;
            //pvec_les_(wv, sz, &guide, dh, (int*)&C.CnvScheme, v00, &rei, v0, vf, (int*)bcv, vt, &flop);
            TIMING_stop("Pvec_Central_LES", flop);
          }
          else
          {
            TIMING_start("Pvec_Central");
            flop = 0.0;
            pvec_central_(d_wv, size, &guide, pitch, &cnv_scheme, v00, &rei, d_v0, d_vf, d_cdf, d_bid, &half, &flop);
            TIMING_stop("Pvec_Central", flop);
          }
          break;
      }
      
      TIMING_start("Pvec_Flux_BC");
      flop = 0.0;
      BC.modPvecFlux(d_wv, d_v0, d_cdf, CurrentTime, &C, v00, flop);
      TIMING_stop("Pvec_Flux_BC", flop);
      break;
      
    default:
      Exit(0);
  }

  
  // 時間積分
  switch (C.AlgorithmF) 
  {
    case Flow_FS_EE_EE:
      TIMING_start("Pvec_Euler_Explicit");
      flop = 0.0;
      euler_explicit_ (d_vc, size, &guide, &dt, d_v0, d_bcd, &flop);
      TIMING_stop("Pvec_Euler_Explicit", flop);
      break;
      
    case Flow_FS_AB2:
      TIMING_start("Pvec_Adams_Bashforth");
      flop = 0.0;
      if ( Session_CurrentStep == 1 ) // 初期とリスタート後，1ステップめ
      {
        euler_explicit_ (d_vc, size, &guide, &dt, d_v0, d_bcd, &flop);
      }
      else 
      {
        ab2_(d_vc, size, &guide, &dt, d_v0, d_abf, d_bcd, v00, &flop);
      }
      TIMING_stop("Pvec_Adams_Bashforth", flop);
      break;
      
    case Flow_FS_AB_CN:
      TIMING_start("Pvec_AB_CN");
      flop = 0.0;
      if ( Session_CurrentStep == 1 ) 
      {
        euler_explicit_ (d_wv, size, &guide, &dt, d_v0, d_bcd, &flop);
      }
      else 
      {
        ab2_(d_wv, size, &guide, &dt, d_v0, d_abf, d_bcd, v00, &flop);
      }
      TIMING_stop("Pvec_AB_CN", flop);
      
// 陰解法部分
      break;
      
    default:
      Exit(0);
  }


  // FORCINGコンポーネントの疑似速度ベクトルの方向修正と力の加算
  if ( C.EnsCompo.forcing == ON ) 
  {
    TIMING_start("Pvec_Forcing");
    flop = 0.0;
    BC.mod_Pvec_Forcing(d_vc, d_v, d_bcd, d_cvf, v00, dt, flop);
    TIMING_stop("Pvec_Forcing", flop);
  }

  
  // 浮力項
  if ( C.isHeatProblem() && (C.Mode.Buoyancy == BOUSSINESQ) ) 
  {
    TIMING_start("Pvec_Buoyancy");
    REAL_TYPE dgr = dt*C.Grashof*rei*rei * v00[0];
    flop = 0.0;
    //Buoyancy(d_vc, dgr, d_ie0, d_bcd, flop);
    ps_buoyancy_(d_vc, size, &guide, &dgr, d_ie0, d_bcd, &C.NoCompo, mat_tbl, &flop);
    TIMING_stop("Pvec_Buoyancy", flop);
  }

  
  // 疑似ベクトルの境界条件
  TIMING_start("Pvec_BC");
  BC.OuterVBCfacePrep (d_vc, d_v0, d_cdf, dt, &C, ensPeriodic, Session_CurrentStep);
  BC.InnerVBCperiodic(d_vc, d_bcd);
  TIMING_stop("Pvec_BC");

  
  // 疑似ベクトルの同期
  if ( numProc > 1 )
  {
    TIMING_start("Sync_Pvec");
    if ( paraMngr->BndCommV3D(d_vc, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("Sync_Pvec", face_comm_size*3.0*guide*sizeof(REAL_TYPE)); // ガイドセル数 x ベクトル
  }
  

  
  /* Crank-Nicolson Iteration
  if ( C.AlgorithmF == Flow_FS_AB_CN ) 
  {
    TIMING_start(tm_copy_array);;
    U.copyV3D(d_wv, size, guide, d_vc, one);
    TIMING_stop(tm_copy_array, 0.0);
    
    for (LSv->setLoopCount(0); LSv->getLoopCount() < LSv->getMaxIteration(); LSv->incLoopCount())
    {
      //CN_Itr(LSv);
      if ( LSv->isErrConverged() || LSv->isResConverged() ) break;
    }
  }
   */
  
  
  TIMING_stop("NS__F_Step_Section");
  // <<< Fractional step section
  

  
  // Poissonのソース部分
  // >>> Poisson Source section
  TIMING_start("Poisson__Source_Section");

  
  // 非VBC面に対してのみ，セルセンターの値から div{u^*} を計算
  TIMING_start("Divergence_of_Pvec");
  flop = 0.0;
  divergence_cc_(d_ws, size, &guide, pitch, d_vc, d_cdf, d_bid, &flop);
  TIMING_stop("Divergence_of_Pvec", flop);
  
  
  // Poissonソース項の速度境界条件（VBC）面による修正
  TIMING_start("Poisson_Src_VBC");
  flop = 0.0;
  BC.modPsrcVBC(d_ws, d_cdf, CurrentTime, &C, v00, d_vf, d_vc, d_v0, dt, flop);
  TIMING_stop("Poisson_Src_VBC", flop);
  
  
  
  // (Neumann_BCType_of_Pressure_on_solid_wall == grad_NS)　のとき，\gamma^{N2}の処理
  //hogehoge
  

  // ソース項のコピー
  //TIMING_start(tm_copy_array);
  //flop = 0.0;
  //U.copyS3D(d_sq, size, guide, d_ws, one);
  //TIMING_stop(tm_copy_array, flop);

  // 反復ソースの初期化
  //TIMING_start(tm_assign_const);
  U.initS3D(d_sq, size, guide, zero);
  //TIMING_stop(tm_assign_const, 0.0);
  
  // Forcingコンポーネントによるソース項の寄与分
  if ( C.EnsCompo.forcing == ON )
  {
    //TIMING_start(tm_force_src);
    flop=0.0;
    BC.mod_Psrc_Forcing(d_sq, d_v0, d_bcd, d_cvf, v00, component_array, flop);
    //TIMING_stop(tm_force_src, flop);
  }
  
  
  // 内部周期境界部分のディリクレソース項
  //TIMING_start(tm_prdc_src);
  //BC.InnerPrdc_Src(d_sq, d_p, d_bcd);
  //TIMING_stop(tm_prdc_src, flop);
  
  
  // 定数項bの自乗和　b_l2
  TIMING_start("Poisson_Src_Norm");
  b_l2 = 0.0;
  flop = 0.0;
  blas_calc_b_(&b_l2, d_b, d_ws, d_bcp, size, &guide, pitch, &dt, &flop);
  TIMING_stop("Poisson_Src_Norm", flop);
  
  
  if ( numProc > 1 )
  {
    TIMING_start("A_R_Poisson_Src_L2");
    double m_tmp = b_l2;
    if ( paraMngr->Allreduce(&m_tmp, &b_l2, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("A_R_Poisson_Src_L2", 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
  }
  
  // L2 norm of b vector
  b_l2 = sqrt(b_l2);
  
  
  
  // Initial residual
  if ( LSp->getResType() == nrm_r_r0 )
  {
    TIMING_start("Poisson_Init_Res");
    res0_l2 = 0.0;
    flop = 0.0;
    blas_calc_r2_(&res0_l2, d_p, d_b, d_bcp, size, &guide, pitch, &flop);
    TIMING_stop("Poisson_Init_Res", flop);
    
    if ( numProc > 1 )
    {
      TIMING_start("A_R_Poisson_Init_Res_L2");
      double m_tmp = res0_l2;
      if ( paraMngr->Allreduce(&m_tmp, &res0_l2, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      TIMING_stop("A_R_Poisson_Init_Res_L2", 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
    }
    
    res0_l2 = sqrt(res0_l2);
  }

  
  TIMING_stop("Poisson__Source_Section");
  // <<< Poisson Source section
  

  
  // VP-Iteration
  // >>> Poisson Iteration section
  TIMING_start("VP-Iteration_Section");


  // 反復回数の積算
  int loop_p  = 0;
  int loop_vp;
  
  
  for (loop_vp=1; loop_vp<DivC.MaxIteration; loop_vp++)
  {
    // 線形ソルバー
    switch (LSp->getLS())
    {
      case SOR:
        TIMING_start("Point_SOR");
        if ( (loop_p += LSp->PointSOR(d_p, d_b, b_l2, res0_l2)) < 0 ) Exit(0);
        TIMING_stop("Point_SOR");
        //if ( (loop_p += LSp->PointSOR_4th(d_p, d_b, d_ws, d_p0, d_sq, dt, dh, b_l2, res0_l2)) < 0 ) Exit(0);
        break;
        
      case SOR2SMA:
        TIMING_start("2-colored_SOR_stride");
        if ( (loop_p += LSp->SOR2_SMA(d_p, d_b, LSp->getMaxIteration(), b_l2, res0_l2)) < 0 ) Exit(0);
        TIMING_stop("2-colored_SOR_stride");
        break;
        
        //case GMRES:
        //  Fgmres(LSp, b_l2, res0_l2);
        //  break;
        
        
      case BiCGSTAB:
        TIMING_start("PBiCGstab");
        if ( (loop_p += LSp->PBiCGstab(d_p, d_b, b_l2, res0_l2)) < 0 ) Exit(0);
        TIMING_stop("PBiCGstab");
        break;
        
      default:
        printf("\tInvalid Linear Solver for Pressure\n");
        Exit(0);
        break;
    }
    
    
    // スカラポテンシャルによる射影と速度の発散の計算 d_dvはdiv(u)のテンポラリ保持に利用
    TIMING_start("Projection_Velocity");
    flop = 0.0;
    update_vec_(d_v, d_vf, d_dv, size, &guide, &dt, pitch, d_vc, d_p, d_bcp, d_cdf, &flop);
    //update_vec4_(d_v, d_vf, d_dv, size, &guide, &dt, pitch, d_vc, d_p, d_bcp, d_cdf, d_bid, &flop, &cnv_scheme);
    TIMING_stop("Projection_Velocity", flop);
    
    
    // 速度の流束形式の境界条件による発散値の修正
    TIMING_start("Projection_Velocity_BC");
    flop=0.0;
    BC.modDivergence(d_dv, d_cdf, CurrentTime, &C, v00, d_vf, d_v, m_buf, flop);
    TIMING_stop("Projection_Velocity_BC", flop);
    

    
    // セルフェイス速度の境界条件の通信部分
    if ( C.EnsCompo.outflow )
    {
      if ( numProc > 1 )
      {
        for (int n=1; n<=C.NoCompo; n++)
        {
          m_snd[2*n]   = m_rcv[2*n]   = m_buf[n].p0; // 積算速度
          m_snd[2*n+1] = m_rcv[2*n+1] = m_buf[n].p1; // 積算回数
        }
        
        TIMING_start("A_R_Projection_VBC");
        if ( paraMngr->Allreduce(m_snd, m_rcv, 2*(C.NoCompo+1), MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop("A_R_Projection_VBC", 2.0*C.NoCompo*numProc*sizeof(REAL_TYPE)*2.0 ); // 双方向 x ノード数 x 変数
        
        for (int n=1; n<=C.NoCompo; n++)
        {
          m_buf[n].p0 = m_rcv[2*n];
          m_buf[n].p1 = m_rcv[2*n+1];
        }
      }
      
      for (int n=1; n<=C.NoCompo; n++)
      {
        if ( cmp[n].getType() == OUTFLOW )
        {
          cmp[n].val[var_Velocity] = m_buf[n].p0 / m_buf[n].p1; // 無次元平均流速
          //printf("%e  = %e  / %e\n", cmp[n].val[var_Velocity], m_buf[n].p0, m_buf[n].p1);
        }
      }
    }
    
    // Forcingコンポーネントによる速度と発散値の修正
    if ( C.EnsCompo.forcing == ON )
    {
      TIMING_start("Projection_Forcing");
      flop=0.0;
      BC.mod_Vdiv_Forcing(d_v, d_bcd, d_cvf, d_dv, dt, v00, m_buf, component_array, flop);
      TIMING_stop("Projection_Forcing", flop);
      
      // 通信部分
      if ( numProc > 1 )
      {
        for (int n=1; n<=C.NoCompo; n++)
        {
          m_snd[2*n]   = m_rcv[2*n]   = m_buf[n].p0; // 積算速度
          m_snd[2*n+1] = m_rcv[2*n+1] = m_buf[n].p1; // 積算圧力損失
        }
        
        TIMING_start("A_R_Projection_Forcing");
        if ( paraMngr->Allreduce(m_snd, m_rcv, 2*(C.NoCompo+1), MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop("A_R_Projection_Forcing", 2.0*(C.NoCompo+1)*numProc*sizeof(REAL_TYPE)*2.0);
        
        for (int n=1; n<=C.NoCompo; n++)
        {
          m_buf[n].p0 = m_rcv[2*n];
          m_buf[n].p1 = m_rcv[2*n+1];
        }
      }
      
      for (int n=1; n<=C.NoCompo; n++)
      {
        if ( cmp[n].isFORCING() )
        {
          REAL_TYPE aa = (REAL_TYPE)cmp[n].getElement();
          cmp[n].val[var_Velocity] = m_buf[n].p0 / aa; // 平均速度
          cmp[n].val[var_Pressure] = m_buf[n].p1 / aa; // 平均圧力損失量
        }
      }
    }
    
    // 反復ソース項
    if ( C.EnsCompo.forcing == ON )
    {
      TIMING_start("Forcing Source");
      flop=0.0;
      BC.mod_Psrc_Forcing(d_sq, d_v, d_bcd, d_cvf, v00, component_array, flop);
      TIMING_stop("Forcing Source", flop);
      
      TIMING_start("Poisson_Src_Norm");
      b_l2 = 0.0;
      flop = 0.0;
      blas_calc_b_(&b_l2, d_b, d_ws, d_bcp, size, &guide, pitch, &dt, &flop);
      TIMING_stop("Poisson_Src_Norm", flop);
      
      if ( numProc > 1 )
      {
        TIMING_start("A_R_Poisson_Src_L2");
        double m_tmp = b_l2;
        if ( paraMngr->Allreduce(&m_tmp, &b_l2, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop("A_R_Poisson_Src_L2", 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
      }
      
      b_l2 = sqrt(b_l2);
    }
    
    
    /* トラクションフリーの場合 >> not neccesarry
    if ( C.EnsCompo.tfree )
    {
      if ( numProc > 1 )
      {
        TIMING_start("Sync_Face_Velocity");
        if ( paraMngr->BndCommV3D(d_vf, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
        TIMING_stop("Sync_Face_Velocity", face_comm_size*guide*3.0*sizeof(REAL_TYPE));
      }
    }
     */

    
    // 速度境界条件　値を代入する境界条件
    TIMING_start("Velocity_BC");
    BC.OuterVBC(d_v, d_vf, d_cdf, CurrentTime, &C, v00, ensPeriodic);
    BC.InnerVBCperiodic(d_v, d_bcd);
    TIMING_stop("Velocity_BC");
    
    
    
    // \nabla {}^f u^{n+1})の計算
    NormDiv(d_dv);
    
    
    /* Forcingコンポーネントによる速度の方向修正(収束判定から除外)  >> TEST
     TIMING_start(tm_prj_frc_dir);
     flop=0.0;
     BC.mod_Dir_Forcing(d_v, d_bcd, d_cvf, v00, flop);
     TIMING_stop(tm_prj_frc_dir, flop);
     */
    
    LSp->setLoopCount(loop_p);
    
    // 収束判定
    if ( DivC.divergence <= DivC.divEPS ) break;
  }

  
  // 総反復回数を代入
  DivC.Iteration = loop_vp;


  TIMING_stop("VP-Iteration_Section", 0.0);
  // <<< Poisson Iteration section
  
  
  
  
  /// >>> NS Loop post section
  TIMING_start("NS__Loop_Post_Section");
  
  
  // 同期
  if ( numProc > 1 )
  {
    TIMING_start("Sync_Velocity");
    if ( paraMngr->BndCommV3D(d_v, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("Sync_Velocity", face_comm_size*guide*3.0*sizeof(REAL_TYPE));
  }
  
  
  // 外部領域境界面での速度や流量を計算 > 外部流出境界条件の移流速度に利用
  TIMING_start("Domain_Monitor");
  DomainMonitor(BC.exportOBC(), &C);
  TIMING_stop("Domain_Monitor");
  
  
  
  /* 非同期にして隠す
  if (C.LES.Calc==ON) 
  {
    TIMING_start(tm_LES_eddy);
    flop = 0.0;
    eddy_viscosity_(d_vt, size, &guide, &dh, &C.Reynolds, &C.LES.Cs, d_v, d_cdf, range_Ut, range_Yp, v00);
    TIMING_stop(tm_LES_eddy, flop);
    
    if ( numProc > 1 ) 
    {
      TIMING_start(tm_LES_eddy_comm);
      if ( paraMngr->BndCommS3D(d_vt, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      TIMING_stop(tm_LES_eddy_comm, face_comm_size*guide*sizeof(REAL_TYPE));
    }
  }
   */


  TIMING_stop("NS__Loop_Post_Section", 0.0);
  // >>> NS loop post section
  
  // 後始末
  if ( m_buf ) delete [] m_buf;
  if ( m_snd ) delete [] m_snd;
  if ( m_rcv ) delete [] m_rcv;

}

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
 * @file   NS_FS_E_Binary.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"


// Fractional Step法でNavier-Stokes方程式を解く．バイナリ近似．
void FFV::NS_FS_E_Binary()
{
  // local variables
  double flop;                         /// 浮動小数演算数
  double rhs_nrm = 0.0;                /// 反復解法での定数項ベクトルのL2ノルム
  double res_init = 0.0;               /// 反復解法での初期残差ベクトルのL2ノルム
  double comm_size;                    /// 通信面1面あたりの通信量
  double convergence=0.0;              /// 定常収束モニター量
  
  REAL_TYPE dt = deltaT;               /// 時間積分幅
  REAL_TYPE dh = (REAL_TYPE)deltaX;    /// 空間格子幅
  REAL_TYPE coef = deltaX/dt;          /// Poissonソース項の係数
  REAL_TYPE Re = C.Reynolds;           /// レイノルズ数
  REAL_TYPE rei = C.getRcpReynolds();  /// レイノルズ数の逆数
  REAL_TYPE half = 0.5;                /// 定数
  REAL_TYPE one = 1.0;                 /// 定数
  REAL_TYPE zero = 0.0;                /// 定数
  int wall_prof = C.Mode.Wall_profile; /// 壁面条件（slip/noslip）
  int cnv_scheme = C.CnvScheme;        /// 対流項スキーム
  
  
  // 境界処理用
  Gemini_R* m_buf = new Gemini_R [C.NoBC];
  REAL_TYPE* m_snd = new REAL_TYPE [C.NoBC*2];
  REAL_TYPE* m_rcv = new REAL_TYPE [C.NoBC*2];
  
  comm_size = count_comm_size(size, 1);
  
  int v_mode=0;
  
  ItrCtl* ICp = &IC[ItrCtl::ic_prs_pr];  /// 圧力のPoisson反復
  ItrCtl* ICv = &IC[ItrCtl::ic_vis_cn];  /// 粘性項のCrank-Nicolson反復
  ItrCtl* ICd = &IC[ItrCtl::ic_div];     /// 圧力-速度反復

  // point Data
  // d_v   セルセンタ速度 v^n -> v^{n+1}
  // d_vf  セルフェイス速度
  // d_v0  セルセンタ速度 v^nの保持
  // d_vc  疑似速度ベクトル
  // d_wv  ワーク　陰解法の時の疑似速度ベクトル，射影ステップの境界条件
  // d_p   圧力 p^n -> p^{n+1}
  // d_p0  圧力 p^nの保持
  // d_ws  Poissonのソース項0　速度境界を考慮
  // d_sq  Poissonのソース項1　反復毎に変化するソース項，摩擦速度，
  // d_dv  発散値, div(u)の値を保持
  // d_b   反復の右辺ベクトル
  // d_bcd IDのビットフラグ
  // d_bcp 圧力のビットフラグ
  // d_bcv 速度のビットフラグ
  // d_wo  ワーク　壁関数利用時のWSS，ベクトル出力時のテンポラリ
  // d_cvf コンポーネントの体積率
  // d_t0  温度 t^n 
  // d_vt  LES計算の渦粘性係数
  // d_ab0 Adams-Bashforth用のワーク

  
  // >>> Fractional step section
  TIMING_start(tm_frctnl_stp_sct); 
  
  // >>> Fractional step sub-section 1
  TIMING_start(tm_frctnl_stp_sct_1); 
  
  
  // n stepの値を保持 >> In use (d_v0, d_p0)
  TIMING_start(tm_copy_array);
  flop = 0.0;
  U.xcopy(d_p0, size, guide, d_p, one, kind_scalar, flop);
  U.xcopy(d_v0, size, guide, d_v, one, kind_vector, flop);
  TIMING_stop(tm_copy_array, flop, 2);
  

  // 壁関数指定時の摩擦速度の計算 src0をテンポラリのワークとして利用
  if ( C.Mode.Wall_profile == Control::Log_Law ) 
  {
    TIMING_start(tm_WallFunc);
    flop = 0.0;
    friction_velocity_(d_ws, size, &guide, &dh, &Re, d_v0, d_bcp, range_Yp, range_Ut, v00, &flop);
    TIMING_stop(tm_WallFunc, flop);
  }

  TIMING_stop(tm_frctnl_stp_sct_1, 0.0);
  // <<< Fractional step subsection 1


  
  // >>> Fractional step sub-section 2
  TIMING_start(tm_frctnl_stp_sct_2);
  
  // 対流項と粘性項の評価 >> In use (d_vc, d_wv)
  switch (C.AlgorithmF) {
    case Control::Flow_FS_EE_EE:
    case Control::Flow_FS_AB2:
      TIMING_start(tm_pseudo_vec);     
      
      flop = 0.0;
      v_mode = (C.Mode.Wall_profile == Control::Log_Law) ? 2 : 1;

      if ( C.LES.Calc == ON )
      {
        //pvec_les_(vc, sz, &guide, dh, (int*)&C.CnvScheme, v00, &rei, v0, vf, (int*)bcv, vt, &flop);
        Exit(0);
      }
      else
      {
        pvec_muscl_(d_vc, size, &guide, &dh, &cnv_scheme, v00, &rei, d_v0, d_vf, d_bcv, d_bcp, &v_mode, d_ws, &wall_prof, d_bcd, d_cvf, &one, &flop);
      }
      TIMING_stop(tm_pseudo_vec, flop);

      TIMING_start(tm_pvec_flux);      
      flop = 0.0;
      BC.mod_Pvec_Flux(d_vc, d_v0, d_vf, d_bcv, CurrentTime, &C, v_mode, v00, flop);
      TIMING_stop(tm_pvec_flux, flop);
      break;
      
    case Control::Flow_FS_AB_CN:
      TIMING_start(tm_pseudo_vec);
      flop = 0.0;
      v_mode = 0;
      if ( C.LES.Calc == ON )
      {
        //pvec_les_(wv, sz, &guide, &dh, (int*)&C.CnvScheme, v00, &rei, v0, vf, (int*)bcv, vt, &flop);
      }
      else
      {
        pvec_muscl_(d_wv, size, &guide, &dh, &cnv_scheme, v00, &rei, d_v0, d_vf, d_bcv, d_bcp, &v_mode, d_ws, &wall_prof, d_bcd, d_cvf, &half, &flop);
      }
      TIMING_stop(tm_pseudo_vec, flop);
      
      TIMING_start(tm_pvec_flux);
      flop = 0.0;
      BC.mod_Pvec_Flux(d_wv, d_v0, d_vf, d_bcv, CurrentTime, &C, v_mode, v00, flop);
      TIMING_stop(tm_pvec_flux, flop);
      break;
      
    default:
      Exit(0);
  }

  // 時間積分
  switch (C.AlgorithmF) 
  {
    case Control::Flow_FS_EE_EE:
      TIMING_start(tm_pvec_ee);
      flop = 0.0;
      euler_explicit_ (d_vc, size, &guide, &dt, d_v0, d_bcd, &flop);
      TIMING_stop(tm_pvec_ee, flop);
      break;
      
    case Control::Flow_FS_AB2:
      TIMING_start(tm_pvec_ab);
      flop = 0.0;
      if ( Session_CurrentStep == 1 ) // 初期とリスタート後，1ステップめ
      {
        euler_explicit_ (d_vc, size, &guide, &dt, d_v0, d_bcd, &flop);
      }
      else 
      {
        ab2_(d_vc, size, &guide, &dt, d_v0, d_abf, d_bcd, v00, &flop);
      }
      TIMING_stop(tm_pvec_ab, flop);
      break;
      
    case Control::Flow_FS_AB_CN:
      TIMING_start(tm_pvec_abcn);
      flop = 0.0;
      if ( Session_CurrentStep == 1 ) 
      {
        euler_explicit_ (d_wv, size, &guide, &dt, d_v0, d_bcd, &flop);
      }
      else 
      {
        ab2_(d_wv, size, &guide, &dt, d_v0, d_abf, d_bcd, v00, &flop);
      }
      TIMING_stop(tm_pvec_abcn, flop);
      
      TIMING_start(tm_pvec_abcn_df_ee);
      flop = 0.0;
      vis_ee_(d_vc, size, &guide, &dh, &dt, v00, &rei, d_wv, d_v0, d_bcv, &half, &flop);
      TIMING_stop(tm_pvec_abcn_df_ee, flop);
      
      TIMING_start(tm_pvec_abcn_df_ee_BC);
      flop = 0.0;
      BC.mod_Vis_EE(d_vc, d_v0, half, d_bcv, CurrentTime, dt, v00, flop);
      TIMING_stop(tm_pvec_abcn_df_ee_BC, flop);
      break;
      
    default:
      Exit(0);
  }

  TIMING_stop(tm_frctnl_stp_sct_2, 0.0);
  // <<< Fractional step subsection 2
  

  
  // >>> Fractional step sub-section 3
  TIMING_start(tm_frctnl_stp_sct_3);
  
  // FORCINGコンポーネントの疑似速度ベクトルの方向修正と力の加算
  if ( C.isForcing() == ON ) 
  {
    TIMING_start(tm_forcing);
    flop = 0.0;
    BC.mod_Pvec_Forcing(d_vc, d_v, d_bcd, d_cvf, v00, dt, flop);
    TIMING_stop(tm_forcing, flop);
  }

  // 浮力項
  if ( C.isHeatProblem() && (C.Mode.Buoyancy == BOUSSINESQ) ) 
  {
    TIMING_start(tm_buoyancy);
    REAL_TYPE dgr = dt*C.Grashof*rei*rei * v00[0];
    flop = 3.0;
    //Buoyancy(d_vc, dgr, d_t0, d_bcd, flop);
    ps_buoyancy_(d_vc, size, &guide, &dgr, d_t0, d_bcd, &flop);
    TIMING_stop(tm_buoyancy, flop);
  }

  // 疑似ベクトルの境界条件
  TIMING_start(tm_pvec_BC);
  flop = 0.0;
  BC.OuterVBC_Pseudo(d_vc, d_v0, CurrentTime, dt, &C, d_bcv, flop);
  BC.OuterVBC_Periodic(d_vc);
  BC.InnerVBC_Periodic(d_vc, d_bcd);
  TIMING_stop(tm_pvec_BC, flop);

  // 疑似ベクトルの同期
  if ( numProc > 1 ) 
  {
    TIMING_start(tm_pvec_comm);
    if ( paraMngr->BndCommV3D(d_vc, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    TIMING_stop(tm_pvec_comm, comm_size*1.0*3.0); // ガイドセル数 x ベクトル
  }

  TIMING_stop(tm_frctnl_stp_sct_3, 0.0);
  // <<< Fractional step subsection 3
  

  
  // >>> Fractional step sub-section 4
  TIMING_start(tm_frctnl_stp_sct_4);
  
  // Crank-Nicolson Iteration
  if ( C.AlgorithmF == Control::Flow_FS_AB_CN ) 
  {
    TIMING_start(tm_copy_array);
    flop = 0.0;
    U.xcopy(d_wv, size, guide, d_vc, one, kind_vector, flop);
    TIMING_stop(tm_copy_array, flop);
    
    for (ICv->LoopCount=0; ICv->LoopCount< ICv->get_ItrMax(); ICv->LoopCount++) {
      //CN_Itr(ICv);
      if (  ICv->get_normValue() < ICv->get_eps() ) break;
    }
  }
  
  TIMING_stop(tm_frctnl_stp_sct_4, 0.0);
  // <<< Fractional step subsection 4
  
  TIMING_stop(tm_frctnl_stp_sct, 0.0);
  // <<< Fractional step section
  

  
  // Poissonのソース部分
  // >>> Poisson Source section
  TIMING_start(tm_poi_src_sct);
  
  
  // vの初期値をvcにしておく
  TIMING_start(tm_copy_array);
  flop = 0.0;
  U.xcopy(d_v, size, guide, d_vc, one, kind_vector, flop);
  TIMING_stop(tm_copy_array, flop);
  
  
  // 非VBC面に対してのみ，セルセンターの値から \sum{u^*} を計算
  TIMING_start(tm_div_pvec);
  flop = 0.0;
  divergence_(d_ws, size, &guide, d_vc, d_bcv, &flop);
  TIMING_stop(tm_div_pvec, flop);
  
  
  // Poissonソース項の速度境界条件（VBC）面による修正
  TIMING_start(tm_poi_src_vbc);
  flop = 0.0;
  BC.mod_Psrc_VBC(d_ws, d_vc, d_v0, d_vf, d_bcv, CurrentTime, dt, &C, v00, flop);
  TIMING_stop(tm_poi_src_vbc, flop);

  
  // (Neumann_BCType_of_Pressure_on_solid_wall == grad_NS)　のとき，\gamma^{N2}の処理
  //hogehoge
  

  // ソース項のコピー
  //TIMING_start(tm_copy_array);
  //flop = 0.0;
  //U.xcopy(d_sq, size, guide, d_ws, one, kind_scalar, flop);
  //TIMING_stop(tm_copy_array, flop);

  // 反復ソースの初期化
  TIMING_start(tm_assign_const);
  U.xset(d_sq, size, guide, zero, kind_scalar);
  TIMING_stop(tm_assign_const, 0.0);
  
  // Forcingコンポーネントによるソース項の寄与分
  if ( C.isForcing() == ON )
  {
    TIMING_start(tm_force_src);
    flop=0.0;
    BC.mod_Psrc_Forcing(d_sq, d_v0, d_bcd, d_cvf, v00, component_array, flop);
    TIMING_stop(tm_force_src, flop);
  }
  
  
  // 内部周期境界部分のディリクレソース項
  //TIMING_start(tm_prdc_src);
  //BC.InnerPrdc_Src(d_sq, d_p, d_bcd);
  //TIMING_stop(tm_prdc_src, flop);
  
  
  // 定数項のL2ノルム　rhs_nrm
  TIMING_start(tm_poi_src_nrm);
  rhs_nrm = 0.0;
  flop = 0.0;
  poi_rhs_(&rhs_nrm, d_b, size, &guide, d_ws, d_sq, d_bcp, &dh, &dt, &flop);
  TIMING_stop(tm_poi_src_nrm, flop);
  
  if ( ICp->get_LS() == RBGS ||
       ICp->get_LS() == PCG  ||
       ICp->get_LS() == PBiCGSTAB )
  {
		blas_calcb_(d_b, d_ws, d_sq, d_bcp, &dh, &dt, size, &guide);
		REAL_TYPE bb = 0.0;
		blas_dot_(&bb, d_b, d_b, size, &guide);
		rhs_nrm = bb;
	}
  
  if ( numProc > 1 )
  {
    TIMING_start(tm_poi_src_comm);
    double m_tmp = rhs_nrm;
    if ( paraMngr->Allreduce(&m_tmp, &rhs_nrm, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    TIMING_stop(tm_poi_src_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
  }
  
  rhs_nrm = sqrt(rhs_nrm);
  
  
  // Initial residual
  if ( ICp->get_normType() == ItrCtl::r_r0 )
  {
    TIMING_start(tm_poi_src_nrm);
    res_init = 0.0;
    flop = 0.0;
    poi_residual_(&res_init, size, &guide, d_p, d_b, d_bcp, &flop);
    TIMING_stop(tm_poi_src_nrm, flop);
    
    if ( numProc > 1 )
    {
      TIMING_start(tm_poi_src_comm);
      double m_tmp = res_init;
      if ( paraMngr->Allreduce(&m_tmp, &res_init, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      TIMING_stop(tm_poi_src_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
    }
  }
  

  TIMING_stop(tm_poi_src_sct, 0.0);
  // <<< Poisson Source section
  

  
  
  // VP-Iteration
  // >>> Poisson Iteration section
  TIMING_start(tm_poi_itr_sct);
  
  
  // 流出境界条件のアップデート（u^{n+1}を先に）
  for (int face=0; face<NOFACE; face++)
  {
    if( nID[face] < 0 )
    {
      BoundaryOuter* T = BC.export_OBC(face);
      
      if ( T->get_Class() == OBC_OUTFLOW )
      {
        vobc_update_(d_v, size, &guide, d_vc, &face);
      }
      
    }
  }
  
  
  if ( C.Mode.Log_Itr == ON )
  {
    TIMING_start(tm_hstry_itr);
    Hostonly_ H->printHistoryItrTitle(fp_i);
    TIMING_stop(tm_hstry_itr, 0.0);
  }
  
  // 反復回数の積算
  int loop_p = 0;
  ICp->LoopCount = 0;
  ICd->LoopCount = 0;
  
  for (int loop_d=0; loop_d<ICd->get_ItrMax(); loop_d++)
  {

    // 線形ソルバー
    switch (ICp->get_LS())
    {
      case SOR:
        loop_p += Point_SOR(ICp, d_p, d_b, rhs_nrm, res_init);
        break;
        
      case SOR2SMA:
        loop_p += SOR_2_SMA(ICp, d_p, d_b, rhs_nrm, res_init);
        break;
        
      case GMRES:
        Fgmres(ICp, rhs_nrm, res_init);
        break;
        
      case RBGS:
        loop_p += Frbgs(ICp, d_p, d_b, rhs_nrm, res_init);
        break;
        
      case PCG:
        loop_p += Fpcg(ICp, d_p, d_b, rhs_nrm, res_init);
        break;
        
      case PBiCGSTAB:
        loop_p += Fpbicgstab(ICp, d_p, d_b, rhs_nrm, res_init);
        break;
        
      default:
        printf("\tInvalid Linear Solver for Pressure\n");
        Exit(0);
        break;
    }

    
    // >>> Poisson Iteration subsection 4
    TIMING_start(tm_poi_itr_sct_4);

    // 速度のスカラポテンシャルによる射影と速度の発散の計算 d_dvはdiv(u)のテンポラリ保持に利用
    TIMING_start(tm_prj_vec);
    flop = 0.0;
    update_vec_(d_v, d_dv, size, &guide, &dt, &dh, d_vc, d_vf, d_p, d_bcp, d_bcv, &flop);
    TIMING_stop(tm_prj_vec, flop);


    // 速度の境界条件による修正
    TIMING_start(tm_prj_vec_bc);
    flop=0.0;
    BC.mod_div(d_dv, d_bcv, CurrentTime, v00, m_buf, d_vf, d_v, &C, flop);
    TIMING_stop(tm_prj_vec_bc, flop);

    
    // セルフェイス速度の境界条件の通信部分
    if ( C.isOutflow() == ON ) 
    {
      if ( !C.isCDS() ) // Binary
      {
        if ( numProc > 1 ) 
        {
          for (int n=0; n<C.NoBC; n++) {
            m_snd[2*n]   = m_rcv[2*n]   = m_buf[n].p0; // 積算速度
            m_snd[2*n+1] = m_rcv[2*n+1] = m_buf[n].p1; // 積算回数
          }
          
          TIMING_start(tm_prj_vec_bc_comm);
          if ( paraMngr->Allreduce(m_snd, m_rcv, 2*C.NoBC, MPI_SUM) != CPM_SUCCESS ) Exit(0);
          TIMING_stop(tm_prj_vec_bc_comm, 2.0*C.NoBC*numProc*sizeof(REAL_TYPE)*2.0 ); // 双方向 x ノード数 x 変数
          
          for (int n=0; n<C.NoBC; n++) {
            m_buf[n].p0 = m_rcv[2*n];
            m_buf[n].p1 = m_rcv[2*n+1];
          }
        }
        
        for (int n=1; n<=C.NoBC; n++) {
          if ( cmp[n].getType() == OUTFLOW ) 
          {
            cmp[n].val[var_Velocity] = m_buf[n-1].p0 / m_buf[n-1].p1; // 無次元平均流速
          }
        }
      }
      else // Cut-Distance
      {
        ;
      }
    }
    
    // Forcingコンポーネントによる速度と発散値の修正
    if ( C.isForcing() == ON ) 
    {
      TIMING_start(tm_prj_frc_mod);
      flop=0.0;
      BC.mod_Vdiv_Forcing(d_v, d_bcd, d_cvf, d_dv, dt, v00, m_buf, component_array, flop);
      TIMING_stop(tm_prj_frc_mod, flop);

      // 通信部分
      if ( numProc > 1 ) 
      {
        for (int n=0; n<C.NoBC; n++) {
          m_snd[2*n]   = m_rcv[2*n]   = m_buf[n].p0; // 積算速度
          m_snd[2*n+1] = m_rcv[2*n+1] = m_buf[n].p1; // 積算圧力損失
        }
        
        TIMING_start(tm_prj_frc_mod_comm);
        if ( paraMngr->Allreduce(m_snd, m_rcv, 2*C.NoBC, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop(tm_prj_frc_mod_comm, 2.0*C.NoBC*numProc*sizeof(REAL_TYPE)*2.0);
        
        for (int n=0; n<C.NoBC; n++) {
          m_buf[n].p0 = m_rcv[2*n];
          m_buf[n].p1 = m_rcv[2*n+1];
        }
      }
      
      for (int n=1; n<=C.NoBC; n++) {
        if ( cmp[n].isFORCING() ) 
        {
          REAL_TYPE aa = (REAL_TYPE)cmp[n].getElement();
          cmp[n].val[var_Velocity] = m_buf[n-1].p0 / aa; // 平均速度
          cmp[n].val[var_Pressure] = m_buf[n-1].p1 / aa; // 平均圧力損失量
        }
      }
    }

    // 反復ソース項
    if ( C.isForcing() == ON )
    {
      TIMING_start(tm_force_src);
      flop=0.0;
      BC.mod_Psrc_Forcing(d_sq, d_v, d_bcd, d_cvf, v00, component_array, flop);
      TIMING_stop(tm_force_src, flop);
      
      TIMING_start(tm_poi_src_nrm);
      rhs_nrm = 0.0;
      flop = 0.0;
      poi_rhs_(&rhs_nrm, d_b, size, &guide, d_ws, d_sq, d_bcp, &dh, &dt, &flop);
      TIMING_stop(tm_poi_src_nrm, flop);
      
      if ( numProc > 1 )
      {
        TIMING_start(tm_poi_src_comm);
        double m_tmp = rhs_nrm;
        if ( paraMngr->Allreduce(&m_tmp, &rhs_nrm, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop(tm_poi_src_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
      }
      
      rhs_nrm = sqrt(rhs_nrm);
    }


    // トラクションフリーの場合
    if ( C.isTfree() ) {
      if ( numProc > 1 )
      {
        TIMING_start(tm_vectors_comm);
        if ( paraMngr->BndCommV3D(d_vf, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
        TIMING_stop(tm_vectors_comm, 2*comm_size*guide*3.0);
      }
    }
    
    // 速度境界条件　値を代入する境界条件
    TIMING_start(tm_vec_BC);
    flop=0.0;
    
    for (int face=0; face<NOFACE; face++) {
      REAL_TYPE vsum=0.0;
      
      if( nID[face] < 0 )
      {
        BoundaryOuter* T = BC.export_OBC(face);
        
        switch (T->get_Class())
        {
          case OBC_TRC_FREE:
            vobc_tfree_(d_v, size, &guide, &face, d_vf, d_bcv, &vsum, &flop);
            T->setDomainV(vsum);
            break;
            
          case OBC_FAR_FIELD:
            vobc_neumann_(d_v, size, &guide, &face, &vsum);
            T->setDomainV(vsum);
            break;
            
          case OBC_OUTFLOW:
            vobc_get_massflow_(size, &guide, &face, &vsum, d_v, d_bcv, &flop);
            T->setDomainV(vsum);
            break;
            
          default:
            break;
        }
        
      }
    }
    
    BC.OuterVBC_Periodic(d_v);
    BC.InnerVBC_Periodic(d_v, d_bcd);
    TIMING_stop(tm_vec_BC, flop);
    
    TIMING_stop(tm_poi_itr_sct_4, 0.0);
    // <<< Poisson Iteration subsection 4

    
    // div(u^{n+1})の計算
    FB::Vec3i divmaxidx = Norm_Div(ICd);
    
    // 総反復回数を代入
    ICp->LoopCount = loop_p;
    ICd->LoopCount++;
    
    if ( C.Mode.Log_Itr == ON )
    {
      TIMING_start(tm_hstry_itr);
      Hostonly_ {
        H->printHistoryItr(fp_i, IC, divmaxidx);
        fflush(fp_i);
      }
      TIMING_stop(tm_hstry_itr, 0.0);
    }
    

    /* Forcingコンポーネントによる速度の方向修正(収束判定から除外)  >> TEST
    TIMING_start(tm_prj_frc_dir);
    flop=0.0;
    BC.mod_Dir_Forcing(d_v, d_bcd, d_cvf, v00, flop);
    TIMING_stop(tm_prj_frc_dir, flop);
    */
    
    
    // 収束判定　性能測定モードのときは収束判定を行わない
    if ( (C.Hide.PM_Test == OFF) && (ICd->get_normValue() < ICd->get_eps()) ) break;
  } // end of iteration

  
  
  TIMING_stop(tm_poi_itr_sct, 0.0);
  // <<< Poisson Iteration section

  
  
  /// >>> NS Loop post section
  TIMING_start(tm_NS_loop_post_sct);
  
  // 同期
  if ( numProc > 1 )
  {
    TIMING_start(tm_vectors_comm);
    if ( paraMngr->BndCommV3D(d_v, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    TIMING_stop(tm_vectors_comm, 2*comm_size*guide*3.0);
  }
  
  
  // 外部領域境界面での速度や流量を計算 > 外部流出境界条件の移流速度に利用
  TIMING_start(tm_domain_monitor);
  flop=0.0;
  DomainMonitor(BC.export_OBC(), &C);
  TIMING_stop(tm_domain_monitor, flop);
  
  
  // 流出境界のガイドセル値の更新とガイドセル代入
  TIMING_start(tm_VBC_update);
  flop=0.0;
  BC.InnerVBC(d_v, d_bcv, CurrentTime, v00, flop);
  BC.OuterVBC(d_v, d_vc, d_bcv, CurrentTime, dt, &C, v00, flop);
  TIMING_stop(tm_VBC_update, flop);
  
  
  // 非同期にして隠す
  if (C.LES.Calc==ON) 
  {
    TIMING_start(tm_LES_eddy);
    flop = 0.0;
    eddy_viscosity_(d_vt, size, &guide, &dh, &C.Reynolds, &C.LES.Cs, d_v, d_bcv, range_Ut, range_Yp, v00);
    TIMING_stop(tm_LES_eddy, flop);
    
    if ( numProc > 1 ) 
    {
      TIMING_start(tm_LES_eddy_comm);
      if ( paraMngr->BndCommS3D(d_vt, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      TIMING_stop(tm_LES_eddy_comm, comm_size*guide);
    }
  }
  
  
  // ノルムの増加率が規定値をこえたら，終了
  if (convergence_prev != 0.0 ) 
  {
    convergence_rate = convergence / convergence_prev;
  }
  else 
  {
    convergence_rate = 1.0;
  }
  convergence_prev = convergence;
  
  
  // 圧力値の引き戻しオプション
  if ( C.Mode.Pshift != -1 ) 
  {
    TIMING_start(tm_pressure_shift);
    flop = 0.0;
    Pressure_Shift();
    TIMING_stop(tm_pressure_shift, flop);
  }

  TIMING_stop(tm_NS_loop_post_sct, 0.0);
  // >>> NS loop post section
  
  // 後始末
  if ( m_buf ) delete [] m_buf;
  if ( m_snd ) delete [] m_snd;
  if ( m_rcv ) delete [] m_rcv;

}

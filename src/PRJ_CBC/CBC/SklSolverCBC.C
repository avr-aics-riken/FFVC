/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file SklSolverCBC.h
//@brief SklSolverCBC class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <string.h>
#include "SklGloval.h"
#include "SklSolverCBC.h"

//@fn SklSolverCBC::SklSolverCBC()
//@brief コンストラクタ
SklSolverCBC::SklSolverCBC() {
  para_key = -1; /// parallel managerの識別ID，-1は計算領域全体のID
  cm_mode = 0;   /// 同期通信をデフォルト
  
  id_of_solid = 2;
  
  m_condition = m_log = NULL;
  cmp = NULL;
  mat = NULL;
  Ex  = NULL;
  compo_global_bbox= NULL;
  
  // Polylib
  PL = NULL;

  // Cutlib
  cutPos = NULL;
  cutBid = NULL;
  
  // default constructor for inheritance
  ix = jx = kx = gc = NULL;
  dh = dh0 = x0 = y0 = z0 = NULL;
  for (int i=0; i<4; i++) v00[i] = 0.0;
  
  ixc = jxc = kxc = 0;
  guide = 0;
  sz[0] = sz[1] = sz[2] = 0;
  size[0] = size[1] = size[2] = 0;
  
  G_Fcell = G_Acell = G_Wcell = 0;
  G_size[0] = G_size[1] = G_size[2] = 0;
  G_Lbx[0] = G_Lbx[1] = G_Lbx[2] = 0.0;
  G_org[0] = G_org[1] = G_org[2] = 0.0;
  
  range_Ut[0] = range_Ut[1] = 0.0;
  range_Yp[0] = range_Yp[1] = 0.0;
  
  // 履歴情報のファイルポインタ
  mp = stdout;
  fp_b = fp_w = fp_c = fp_d = fp_i = fp_f = NULL;
  
  checkTime = 0.0;
  convergence_prev = 1.0;
  convergence_rate = 0.0;
  
  // (3, ix+guide*2, jx+guide*2, kx+guide*2)
  dc_v   = NULL;
  dc_vc  = NULL;
  dc_v0  = NULL;
  dc_wv  = NULL;
  dc_abf = NULL;
  dc_vf0 = NULL;
  
  dc_av  = NULL;
  
  // (3, ix+guide*2, jx+guide*2, kx+guide*2)
  dc_wvex = NULL;
  dc_qbc  = NULL;
  
  // (ix+guide*2, jx+guide*2, kx+guide*2)
  dc_mid  = NULL;
  dc_bcd  = NULL;
  dc_bcp  = NULL;
  dc_bcv  = NULL;
  dc_bh1  = NULL;
  dc_bh2  = NULL;
  dc_ws   = NULL;
  dc_p    = NULL;
  dc_wk2  = NULL;
  dc_dp   = NULL;
  dc_p0   = NULL;
  dc_t    = NULL;
  dc_t0   = NULL;
  dc_vt   = NULL;
  dc_vof  = NULL;
  
  dc_cvf  = NULL;
  
  dc_ap  = NULL;
  dc_at  = NULL;
  
  dc_r_p = NULL;
  dc_r_v = NULL;
  dc_r_t = NULL;
  
  dc_cut = NULL;
  dc_bid = NULL;
  
  // Cut
  cut = NULL;
  cut_id = NULL;
  
  // コンポーネント配列
  component_array=NULL;
  
  m_outPrs = m_outUVW = m_outTmp = m_outVrt = m_outTP = NULL;
  m_outAvrPrs = m_outAvrUVW = m_outAvrTmp = m_outVOF = NULL;
  m_outI2VGT = m_outHlcty = m_outDiv = NULL;
  
  // for tuning
  cf_sz[0] = cf_sz[1] = cf_sz[2] = 0;
  cf_x = cf_y = cf_z = NULL;
  
  pn.procGrp = 0;
  pn.ID = 0; 
  for (int i=0; i<NOFACE; i++) pn.nID[i] = -1;
  ModeTiming = OFF;

}

//@fn SklSolverCBC::SklSolverCBC()
//@brief コンストラクタ
SklSolverCBC::SklSolverCBC(int sType) {
  para_key = -1; /// parallel managerの識別ID，-1は計算領域全体のID
  cm_mode = 0;   /// 同期通信をデフォルト
  
  id_of_solid = 2;
  
  m_condition = m_log = NULL;
  cmp = NULL;
  mat = NULL;
  Ex  = NULL;
  compo_global_bbox= NULL;
  
  // Polylib
  PL = NULL;
  
  // Cutlib
  cutPos = NULL;
  cutBid = NULL;
  
  // default constructor for inheritance
  ix = jx = kx = gc = NULL;
  dh = dh0 = x0 = y0 = z0 = NULL;
  for (int i=0; i<4; i++) v00[i] = 0.0;
  
  ixc = jxc = kxc = 0;
  guide = 0;
  sz[0] = sz[1] = sz[2] = 0;
  size[0] = size[1] = size[2] = 0;
  
  G_Fcell = G_Acell = G_Wcell = 0;
  G_size[0] = G_size[1] = G_size[2] = 0;
  G_Lbx[0] = G_Lbx[1] = G_Lbx[2] = 0.0;
  G_org[0] = G_org[1] = G_org[2] = 0.0;
  
  range_Ut[0] = range_Ut[1] = 0.0;
  range_Yp[0] = range_Yp[1] = 0.0;
  
  // 履歴情報のファイルポインタ
  mp = stdout;
  fp_b = fp_w = fp_c = fp_d = fp_i = fp_f = NULL;
  
  checkTime = 0.0;
  convergence_prev = 1.0;
  convergence_rate = 0.0;
  
  // (3, ix+guide*2, jx+guide*2, kx+guide*2)
  dc_v   = NULL;
  dc_vc  = NULL;
  dc_v0  = NULL;
  dc_wv  = NULL;
  dc_abf = NULL;
  dc_vf0 = NULL;
  
  dc_av  = NULL;
  
  // (3, ix+guide*2, jx+guide*2, kx+guide*2)
  dc_wvex = NULL;
  dc_qbc  = NULL;
  
  // (ix+guide*2, jx+guide*2, kx+guide*2)
  dc_mid  = NULL;
  dc_bcd  = NULL;
  dc_bcp  = NULL;
  dc_bcv  = NULL;
  dc_bh1  = NULL;
  dc_bh2  = NULL;
  dc_ws   = NULL;
  dc_p    = NULL;
  dc_wk2  = NULL;
  dc_dp   = NULL;
  dc_p0   = NULL;
  dc_t    = NULL;
  dc_t0   = NULL;
  dc_vt   = NULL;
  dc_vof  = NULL;
  
  dc_cvf  = NULL;
  
  dc_ap  = NULL;
  dc_at  = NULL;
  
  dc_r_p = NULL;
  dc_r_v = NULL;
  dc_r_t = NULL;
  
  dc_cut = NULL;
  dc_bid = NULL;
  
  // Cut
  cut = NULL;
  cut_id = NULL;
  
  // コンポーネント配列
  component_array=NULL;
  
  m_outPrs = m_outUVW = m_outTmp = m_outVrt = m_outTP = NULL;
  m_outAvrPrs = m_outAvrUVW = m_outAvrTmp = m_outVOF = NULL;
  m_outI2VGT = m_outHlcty = m_outDiv = NULL;
  
  // for tuning
  cf_sz[0] = cf_sz[1] = cf_sz[2] = 0;
  cf_x = cf_y = cf_z = NULL;
  
  pn.procGrp = 0;
  pn.ID = 0; 
  for (int i=0; i<6; i++) pn.nID[i] = -1;
  ModeTiming = OFF;

}

//@fn SklSolverCBC::~SklSolverCBC()
//@brief デストラクタ
SklSolverCBC::~SklSolverCBC() {
  if( m_condition )  delete [] m_condition;
  if( m_log )        delete [] m_log;
  if( compo_global_bbox )        delete [] compo_global_bbox;
  if( cmp )          delete [] cmp;
  if( mat )          delete [] mat;
  
  if ( mp ) fclose(mp);
  
  if (C.Mode.Log_Base == ON) {
    if (fp_b) fclose(fp_b);
    if (fp_c) fclose(fp_c);
    if (fp_d) fclose(fp_d);
    if (fp_f) fclose(fp_f);
  }
  
  if (C.Mode.Log_Itr == ON) {
    if (fp_i) fclose(fp_i);
  }
    
  if (C.Mode.Log_Wall == ON) {
    if (fp_w) fclose(fp_w);
  }

}

/**
 @fn void SklSolverCBC::DomainMonitor(BoundaryOuter* ptr, Control* R, REAL_TYPE& flop)
 @brief 外部計算領域の各面における総流量と対流流出速度を計算する
 @param ptr BoundaryOuterクラスのポインタ
 @param R Controlクラスのポインタ
 @param flop
 @note 流出境界のみ和をとる，その他は既知
 */
void SklSolverCBC::DomainMonitor(BoundaryOuter* ptr, Control* R, REAL_TYPE& flop)
{
  if ( !ptr ) Exit(0);
  BoundaryOuter* obc=NULL;
  
  obc = ptr;
  
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  REAL_TYPE dh, ddh;
  REAL_TYPE vv[3], ec;
  REAL_TYPE u_sum, u_min, u_max, u_avr;
  REAL_TYPE tmp_sum, tmp_min, tmp_max;
  int ofv;
  
	dh = R->dh;
	ddh= dh*dh;
  
  for (int face=0; face<NOFACE; face++) {
    
    // ofv (1-MINMAX, 2-AVERAGE) ゼロでなければ，流出境界
    ofv = obc[face].get_ofv();
    
    // OUTFLOW, SPEC_VELのときの有効セル数
    ec = (REAL_TYPE)obc[face].get_ValidCell();
    
    // 各プロセスの外部領域面の速度をvv[]にコピー
    obc[face].getDomainV(vv);
    
    // 流出境界のモード
    if (ofv == V_AVERAGE) { // average
      
      // 非境界面はゼロなので，単に足し込むだけ
      if ( para_mng->IsParallel() ) {
        u_sum = tmp_sum = vv[0];
        para_mng->Allreduce(&tmp_sum, &u_sum, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
      } 
      u_avr = (ec != 0.0) ? u_sum / ec : 0.0;
      flop = flop + 1.0;
    }
    else if (ofv == V_MINMAX) { // minmax
      
      if ( para_mng->IsParallel() ) {
        u_sum = tmp_sum = vv[0];
        u_min = tmp_min = vv[1];
        u_max = tmp_max = vv[2];
        para_mng->Allreduce(&tmp_sum, &u_sum, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
        para_mng->Allreduce(&tmp_min, &u_min, 1, SKL_ARRAY_DTYPE_REAL, SKL_MIN, pn.procGrp);
        para_mng->Allreduce(&tmp_max, &u_max, 1, SKL_ARRAY_DTYPE_REAL, SKL_MAX, pn.procGrp);
      }
      u_avr = 0.5*(u_min+u_max);
      flop = flop + 2.0;
    }
    else { // 非OUTFLOW BCは速度がストアされている
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

//@fn void SklSolverCBC::set_label(void)
//@brief タイミング測定区間にラベルを与えるラッパー
//@param[in] key キー番号
//@param[in] label ラベル
//@param[in] type  測定対象タイプ(COMM or CALC)
//@param[in] exclusive 排他測定フラグ(ディフォルトtrue)
void SklSolverCBC::set_label(unsigned key, char* label, PerfMonitor::Type type, bool exclusive)
{
  // 文字数がTM_LABEL_MAX-1を超えるものはカット
  int len = strlen(label);
  char label_tmp[TM_LABEL_MAX];
  
  if ( len>TM_LABEL_MAX-1 ) {
    strncpy(label_tmp, label, TM_LABEL_MAX-1);
    printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
  }
  else {
    strcpy(label_tmp, label);
  }
  
  // Performance Monitorへの登録
  string tmp(label_tmp);
  PM.setProperties(key, tmp, type, exclusive);

  // FX用プロファイラの文字列登録
  strcpy(tm_label_ptr[key], label_tmp);
}

//@fn void SklSolverCBC::set_timing_label(void)
//@brief タイミング測定区間にラベルを与える
void SklSolverCBC::set_timing_label(void)
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

//@fn void SklSolverCBC::Averaging_Time(REAL_TYPE& flop)
//@brief 時間平均操作を行う
void SklSolverCBC::Averaging_Time(REAL_TYPE& flop)
{
  SklIncrementAverageStep();
  SklIncrementAverageTime();
  
  REAL_TYPE *p, *ap;
  REAL_TYPE *v, *av;
  REAL_TYPE *t, *at;
  p = ap = NULL;
  v = av = NULL;
  t = at = NULL;
  
  if( !(v  = dc_v->GetData()) )   Exit(0);
  if( !(av = dc_av->GetData()) )  Exit(0);
  if( !(p  = dc_p->GetData()) )   Exit(0);
  if( !(ap = dc_ap->GetData()) )  Exit(0);
  
  int d_length;

  d_length = (int)dc_ap->GetArrayLength();
  fb_average_(ap, p, &d_length, &flop);
  
  d_length = (int)dc_av->GetArrayLength();
  fb_average_(av, v, &d_length, &flop);

  if ( C.isHeatProblem() ) {
    if( !(t  = dc_t->GetData()) )  Exit(0);
    if( !(at = dc_at->GetData()) ) Exit(0);
    
    d_length = (int)dc_at->GetArrayLength();
    fb_average_(at, t, &d_length, &flop);
  }
}

//@fn void SklSolverCBC::Variation_Space(REAL_TYPE* avr, REAL_TYPE& flop)
//@brief 空間平均操作と変動量の計算を行う
//@param avr[out] 平均値と変動値
//@param flop[out] 浮動小数演算数
//@note スカラ値は算術平均，ベクトル値は自乗和
void SklSolverCBC::Variation_Space(REAL_TYPE* avr, REAL_TYPE& flop)
{
  REAL_TYPE *v =NULL;
  REAL_TYPE *v0=NULL;
  REAL_TYPE *p =NULL;
  REAL_TYPE *p0=NULL;
  REAL_TYPE *t =NULL;
  REAL_TYPE *t0=NULL;
  unsigned *bd=NULL;
  
  if( !(v  = dc_v->GetData()) )   Exit(0);
  if( !(v0 = dc_v0->GetData()) )  Exit(0);
  if( !(p  = dc_p->GetData()) )   Exit(0);
  if( !(p0 = dc_p0->GetData()) )  Exit(0);
  if( !(bd = dc_bcd->GetData()) ) Exit(0);
  
  if ( C.isHeatProblem() ) {
    if( !(t  = dc_t->GetData()) )  Exit(0);
    if( !(t0 = dc_t0->GetData()) ) Exit(0);
  }
  
  REAL_TYPE m_var[2];
  
  // 速度
  fb_delta_v_(m_var, sz, gc, v, v0, (int*)bd, &flop); // 速度反復でV_res_L2_を計算している場合はスキップすること
  avr[var_Velocity]   = m_var[0]; // 0
  avr[var_Velocity+3] = m_var[1]; // 3

  // 圧力
  fb_delta_s_(m_var, sz, gc, p, p0, (int*)bd, &flop);
  avr[var_Pressure]   = m_var[0]; // 1
  avr[var_Pressure+3] = m_var[1]; // 4

  // 温度
  if ( C.isHeatProblem() ) {
    fb_delta_s_(m_var, sz, gc, t, t0, (int*)bd, &flop);
    avr[var_Temperature]   = m_var[0]; // 2
    avr[var_Temperature+3] = m_var[1]; // 5
  }
  
}

//@fn void SklSolverCBC::AverageOutput (unsigned mode, REAL_TYPE& flop)
//@brief 時間平均値のファイル出力
//@param mode 強制出力
//@param flop 浮動小数演算数
void SklSolverCBC::AverageOutput (unsigned mode, REAL_TYPE& flop)
{
  bool forceFlag=false;
  bool correct_flag=false; // no correction
  unsigned guide_zero=0;  // for averaged fields
  int stepAvr = (int)SklGetAverageTotalStep();
  REAL_TYPE timeAvr, unit_velocity, scale;
  int d_length;
  
  REAL_TYPE *vo = NULL;
  REAL_TYPE *av = NULL;
  REAL_TYPE *ws = NULL;
  REAL_TYPE *ap = NULL;
  if( !(vo = dc_wvex->GetData()) ) Exit(0);
  if( !(av = dc_av->GetData()) )   Exit(0);
  if( !(ws = dc_ws->GetData()) )   Exit(0);
  if( !(ap = dc_ap->GetData()) )   Exit(0);
  
  
  if (C.Unit.File == DIMENSIONAL) {
    timeAvr = SklGetAverageTotalTime() * C.Tscale;
  }
  else {
    timeAvr = SklGetAverageTotalTime();
  }
  
  if ( mode == Control::IO_forced ) forceFlag = true;
  
  scale = 1.0 / (REAL_TYPE)stepAvr;
  
  // Pressure
  if( m_outAvrPrs ){
    d_length = (int)dc_ws->GetArrayLength();
    
    if (C.Unit.File == DIMENSIONAL) {
      REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
      fb_prs_nd2d_(ws, ap, &d_length, &bp, &C.RefDensity, &C.RefVelocity, &scale, &flop);
    }
    else {
      fb_xcopy_(ws, ap, &d_length, &scale, &flop);
    }
    
    if( !WriteFile(m_outAvrPrs, stepAvr, timeAvr, pn.procGrp, forceFlag, guide_zero, correct_flag) ) {
      stamped_printf("\tFile output was failed.\n");
      Exit(0);
    }
  }
  
  // Velocity
  if( m_outAvrUVW ){
    unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
    fb_shift_refv_out_(vo, sz, gc, av, v00, &scale, &unit_velocity, &flop);
    
    if( !WriteFile(m_outAvrUVW, stepAvr, timeAvr, pn.procGrp, forceFlag, guide_zero, correct_flag) ) {
      stamped_printf("\tFile output was failed.\n");
      Exit(0);
    }
  }
  
  // Temperature
  if( C.isHeatProblem() && m_outAvrTmp ){
    REAL_TYPE *at = NULL;
    if( !(at = dc_ap->GetData()) )   Exit(0);
    d_length = (int)dc_ws->GetArrayLength();
    
    if (C.Unit.File == DIMENSIONAL) {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      fb_tmp_nd2d_(ws, at, &d_length, &C.BaseTemp, &C.DiffTemp, &klv, &scale, &flop);
    }
    else {
      fb_xcopy_(ws, at, &d_length, &scale, &flop);
    }
    
    if( !WriteFile(m_outAvrTmp, stepAvr, timeAvr, pn.procGrp, forceFlag, guide_zero, correct_flag) ) {
      stamped_printf("\tFile output was failed.\n");
      Exit(0);
    }
  }
}

/**
 @fn void SklSolverCBC::FileOutput (unsigned mode, REAL_TYPE& flop)
 @brief ファイル出力
 @param mode 強制出力
 @param flop 浮動小数点演算数
 @note dc_p0をワークとして使用
 */
void SklSolverCBC::FileOutput (unsigned mode, REAL_TYPE& flop)
{
  SklParaComponent* para_cmp = ParaCmpo;
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  bool forceFlag=false;
  bool correct_flag=false; // no correction
  unsigned step = SklGetTotalStep();
  REAL_TYPE time, unit_velocity;
  REAL_TYPE scale = 1.0;
  int d_length;
  
  REAL_TYPE *vo = NULL;
  REAL_TYPE *v  = NULL;
  REAL_TYPE *ws = NULL;
  REAL_TYPE *p  = NULL;
  
  if( !(vo = dc_wvex->GetData()) ) Exit(0);
  if( !(v  = dc_v->GetData()) )    Exit(0);
  if( !(ws = dc_ws->GetData()) )   Exit(0);
  if( !(p  = dc_p->GetData()) )    Exit(0);
  

  if (C.Unit.File == DIMENSIONAL) {
    time = SklGetTotalTime() * C.Tscale;
  }
  else {
    time = SklGetTotalTime();
  }

  if ( mode == Control::IO_forced ) forceFlag = true;
  
  // Divergence デバッグ用なので無次元のみ
  if( m_outDiv ){
    REAL_TYPE coef = SklGetDeltaT()/(C.dh*C.dh); /// 発散値を計算するための係数　dt/h^2
    F.cnv_Div(dc_ws, dc_wk2, coef, flop);
    
    if( !WriteFile(m_outDiv, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
      stamped_printf("\tFile output was failed.\n");
      Exit(0);
    }
  }
  
  // Pressure
  if( m_outPrs ){
    d_length = (int)dc_ws->GetArrayLength();
    
    if (C.Unit.File == DIMENSIONAL) {
      REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
      fb_prs_nd2d_(ws, p, &d_length, &bp, &C.RefDensity, &C.RefVelocity, &scale, &flop);
    }
    else {
      fb_xcopy_(ws, p, &d_length, &scale, &flop);
    }
    
    if( !WriteFile(m_outPrs, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
      stamped_printf("\tFile output was failed.\n");
      Exit(0);
    }
  }

  // Velocity
  if( m_outUVW ){
    unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity : 1.0;
    fb_shift_refv_out_(vo, sz, gc, v, v00, &scale, &unit_velocity, &flop);
    
    if( !WriteFile(m_outUVW, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
      stamped_printf("\tFile output was failed.\n");
      Exit(0);
    }
  }

  // Tempearture
  if( C.isHeatProblem() && m_outTmp ){
    REAL_TYPE *t  = NULL;
    if( !(t  = dc_t->GetData()) )    Exit(0);
    d_length = (int)dc_ws->GetArrayLength();
    
    if (C.Unit.File == DIMENSIONAL) {
      REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
      fb_tmp_nd2d_(ws, t, &d_length, &C.BaseTemp, &C.DiffTemp, &klv, &scale, &flop);
    }
    else {
      fb_xcopy_(ws, t, &d_length, &scale, &flop);
    }
    if( !WriteFile(m_outTmp, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
      stamped_printf("\tFile output was failed.\n");
      Exit(0);
    }
  }
  
  // Total Pressure
  if (C.Mode.TP == ON ) {
    if( m_outTP ){
      REAL_TYPE  *p, *tp;
      if( !(p  = dc_p->GetData()) )   Exit(0);
      if( !(tp = dc_p0->GetData()) )  Exit(0);
      
      fb_totalp_ (tp, sz, gc, v, p, v00, &flop);
      
      // convert non-dimensional to dimensional, iff file is dimensional
      if (C.Unit.File == DIMENSIONAL) {
        F.cnv_TP_ND2D(dc_ws, dc_p0, C.RefDensity, C.RefVelocity, flop);
      }
      else {
        if( !SklUtil::cpyS3D(dc_ws, dc_p0) ) Exit(0);
      }
      if( !WriteFile(m_outTP, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
        stamped_printf("\tFile output was failed.\n");
        Exit(0);
      }
    }
  }
  
  // Vorticity
  if (C.Mode.VRT == ON ) {
    if( m_outVrt ) {
      REAL_TYPE *vrt= NULL;
      unsigned *bcv=NULL;
      
      if( !(bcv = dc_bcv->GetData()) ) Exit(0);
      if( !(vrt = dc_wv->GetData()) )  Exit(0);
      
      cbc_rot_v_(vrt, sz, gc, dh, v, (int*)bcv, v00, &flop);
      
      REAL_TYPE  vz[3];
      vz[0] = vz[1] = vz[2] = 0.0;
      unit_velocity = (C.Unit.File == DIMENSIONAL) ? C.RefVelocity/C.RefLength : 1.0;
      fb_shift_refv_out_(vo, sz, gc, vrt, vz, &scale, &unit_velocity, &flop);
      
      if( !WriteFile(m_outVrt, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
        stamped_printf("\tFile output was failed.\n");
        Exit(0);
      }
    }
  }
  
  // 2nd Invariant of Velocity Gradient Tensor
  if (C.Mode.I2VGT == ON ) {
    if( m_outI2VGT ){
      REAL_TYPE  *v=NULL, *q=NULL;
      unsigned *bcv=NULL;
      
      if( !(v  = dc_v->GetData()) )   Exit(0);
      if( !(q  = dc_p0->GetData()) )  Exit(0);
      if( !(bcv= dc_bcv->GetData()) ) Exit(0);
      
      cbc_i2vgt_ (q, sz, gc, dh, v, (int*)bcv, v00, &flop);
      
      // 無次元で出力
      d_length = (int)dc_ws->GetArrayLength();
      fb_xcopy_(ws, q, &d_length, &scale, &flop);

      if( !WriteFile(m_outI2VGT, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
        stamped_printf("\tFile output was failed.\n");
        Exit(0);
      }
    }
  }
  
  // Helicity
  if (C.Mode.Helicity == ON ) {
    if( m_outHlcty ){
      REAL_TYPE  *v=NULL, *q=NULL;
      unsigned *bcv=NULL;
      
      if( !(v  = dc_v->GetData()) )   Exit(0);
      if( !(q  = dc_p0->GetData()) )  Exit(0);
      if( !(bcv= dc_bcv->GetData()) ) Exit(0);
      
      cbc_helicity_(q, sz, gc, dh, v, (int*)bcv, v00, &flop);
      
      // 無次元で出力
      d_length = (int)dc_ws->GetArrayLength();
      fb_xcopy_(ws, q, &d_length, &scale, &flop);
      
      if( !WriteFile(m_outHlcty, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
        stamped_printf("\tFile output was failed.\n");
        Exit(0);
      }
    }
  }
  
  // Interface Function
  if ( C.BasicEqs == INCMP_2PHASE ) {
    if( m_outVOF ){
      if( !SklUtil::cpyS3D(dc_ws, dc_vof) ) Exit(0);
      if( !WriteFile(m_outVOF, (int)step, time, pn.procGrp, forceFlag, C.GuideOut, correct_flag) ) {
        stamped_printf("\tFile output was failed.\n");
        Exit(0);
      }
    }
  }
}

/**
 @fn void SklSolverCBC::LS_Binary(ItrCtl* IC, REAL_TYPE b2)
 @brief 線形ソルバーの選択実行
 @param IC ItrCtlクラス
 @param b2 ソースベクトルの自乗和
 */
void SklSolverCBC::LS_Binary(ItrCtl* IC, REAL_TYPE b2)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
	
	REAL_TYPE omg, r;
	REAL_TYPE *p, *src0, *p0, *src1;
	unsigned *bcp;
  REAL_TYPE flop_count=0.0;
  REAL_TYPE np_f = (REAL_TYPE)para_mng->GetNodeNum(pn.procGrp); /// 全ノード数
  REAL_TYPE comm_size;              /// 通信面1面あたりの通信量
  unsigned int wait_num=0;
  int req[12];
  REAL_TYPE clear_value=0.0;
	
	p = src0 = src1 = p0 = NULL;
	bcp = NULL;
  omg = IC->get_omg();
	r = 0.0;
  comm_size = count_comm_size(size, guide);
  
	if( !(p   = dc_p->GetData()) )   Exit(0); // 圧力 p^{n+1}
  if( !(p0  = dc_p0->GetData()) )  Exit(0); // 圧力 p^n
	if( !(src0= dc_ws->GetData()) )  Exit(0); // 非反復のソース項
  if( !(src1= dc_wk2->GetData()) ) Exit(0); // 反復毎に変化するソース項
	if( !(bcp = dc_bcp->GetData()) ) Exit(0); // ビットフラグ
  
  // 反復処理
  switch (IC->get_LS()) {
      
    case SOR:
      // >>> Poisson Iteration section 2
      TIMING_start(tm_poi_itr_sct_2);
      
      // 反復処理
      TIMING_start(tm_poi_PSOR);
      flop_count = 0.0;
      cbc_psor_(p, sz, gc, &omg, &r, src0, src1, (int*)bcp, &flop_count);
      //r = PSOR(p, src0, src1, bcp, IC, flop_count); //実装速度比較
      TIMING_stop(tm_poi_PSOR, flop_count);
      
      // 境界条件
      TIMING_start(tm_poi_BC);
      BC.OuterPBC(dc_p);
      if ( C.isPeriodic() == ON ) BC.InnerPBC_Periodic(dc_p, dc_bcd);
      TIMING_stop(tm_poi_BC, 0.0);
      
      TIMING_stop(tm_poi_itr_sct_2, 0.0);
      // <<< Poisson Iteration subsection 2
      
      
      
      // >>> Poisson Iteration section 3
      TIMING_start(tm_poi_itr_sct_3);

      // 同期処理
      if ( para_mng->IsParallel() ) {
        TIMING_start(tm_poi_comm);
        if (cm_mode == 0 ) {
          if( !dc_p->CommBndCell(1) ) Exit(0); // 1 layer communication
        }
        else {
          if( !dc_p->CommBndCell2(1, wait_num, req) ) Exit(0); // 1 layer communication
          para_mng->WaitAll(wait_num, req);
        }
        TIMING_stop(tm_poi_comm, comm_size);
      }
      
      // 残差の集約
      if ( para_mng->IsParallel() ) {
        TIMING_start(tm_poi_res_comm);
        REAL_TYPE tmp = r;
        para_mng->Allreduce(&tmp, &r, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
        TIMING_stop(tm_poi_res_comm, 2.0*np_f*(REAL_TYPE)sizeof(REAL_TYPE) ); // 双方向 x ノード数
      }

      TIMING_stop(tm_poi_itr_sct_3, 0.0);
      // <<< Poisson Iteration subsection 3
      break;
      
    case SOR2SMA:
      // 2色のマルチカラーのセットアップ
      TIMING_start(tm_poi_setup);

      int nID[6], sidx[3], ip;
      
      if ( para_mng->IsParallel() ){
        sidx[0] = para_mng->GetVoxelHeadIndex(pn.ID, 0, pn.procGrp) + 1;
        sidx[1] = para_mng->GetVoxelHeadIndex(pn.ID, 1, pn.procGrp) + 1;
        sidx[2] = para_mng->GetVoxelHeadIndex(pn.ID, 2, pn.procGrp) + 1;
        ip = (sidx[0]+sidx[1]+sidx[2]) % 2;
      } 
      else {
        sidx[0] = sidx[1] = sidx[2] = 1;
        ip = 0;
      }
      TIMING_stop(tm_poi_setup, 0.0);
      
      // 各カラー毎の間に同期
      r = 0.0;          // 色間で積算する
      for (int color=0; color<2; color++) {
        // >>> Poisson Iteration section 2
        TIMING_start(tm_poi_itr_sct_2);
        
        TIMING_start(tm_poi_SOR2SMA);
        flop_count = 0.0; // 色間で積算しない
        cbc_psor2sma_core_(p, sz, gc, &ip, &color, &omg, &r, src0, src1, (int*)bcp, &flop_count);
        //PSOR2sma_core(p, ip, color, src0, src1, bcp, IC, flop_count);
        TIMING_stop(tm_poi_SOR2SMA, flop_count);
        
        // 境界条件
        TIMING_start(tm_poi_BC);
        BC.OuterPBC(dc_p);
        if ( C.isPeriodic() == ON ) BC.InnerPBC_Periodic(dc_p, dc_bcd);
        TIMING_stop(tm_poi_BC, 0.0);

        TIMING_stop(tm_poi_itr_sct_2, 0.0);
        // <<< Poisson Iteration subsection 2
        
        
        
        // >>> Poisson Iteration section 3
        TIMING_start(tm_poi_itr_sct_3);
        
        // 残差の集約と同期処理
        if ( para_mng->IsParallel() ) {
          TIMING_start(tm_poi_res_comm);
          REAL_TYPE tmp = r;
          para_mng->Allreduce(&tmp, &r, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp);
          TIMING_stop(tm_poi_res_comm, 2.0*np_f*(REAL_TYPE)sizeof(REAL_TYPE)*0.5 ); // 双方向 x ノード数 check
          
          TIMING_start(tm_poi_comm);
          if (cm_mode == 0 ) {
            if( !dc_p->CommBndCell(1) ) Exit(0); // 1 layer communication
          }
          else {
            cbc_sma_comm_     (p, sz, gc, &color, &ip, cf_sz, cf_x, cf_y, cf_z, req, &para_key);
            cbc_sma_comm_wait_(p, sz, gc, &color, &ip, cf_sz, cf_x, cf_y, cf_z, req);
          }
          TIMING_stop(tm_poi_comm, comm_size*0.5);
        }
        
        TIMING_stop(tm_poi_itr_sct_3, 0.0);
        // <<< Poisson Iteration subsection 3
      }
      break;
      
    default:
      printf("\tInvalid Linear Solver for Pressure\n");
      Exit(0);
      break;
  }

  // 残差の保存 
  /// @note Control::p_res_l2_r/aのときのみ，それ以外はNorm_Poissonで計算
  if ( (IC->get_normType() == ItrCtl::p_res_l2_r) || (IC->get_normType() == ItrCtl::p_res_l2_a) ) {
    REAL_TYPE res_a = sqrt( r /(REAL_TYPE)G_Acell ); // 残差のRMS
    REAL_TYPE bb  = sqrt( b2/(REAL_TYPE)G_Acell ); // ソースベクトルのRMS
    REAL_TYPE res_r = ( bb == 0.0 ) ? res_a : res_a/bb;
    REAL_TYPE nrm = ( IC->get_normType() == ItrCtl::p_res_l2_r ) ? res_r : res_a;
    IC->set_normValue( nrm );
  }
}

/**
 @fn void SklSolverCBC::LS_Planar(ItrCtl* IC, REAL_TYPE b2)
 @brief 線形ソルバーの選択実行
 @param IC ItrCtlクラス
 @param b2 ソースベクトルのL2ノルム
 */
void SklSolverCBC::LS_Planar(ItrCtl* IC, REAL_TYPE b2)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
	REAL_TYPE omg, r;
	REAL_TYPE *p, *ws, *p0;
  unsigned *bcv;
	
	p = ws = cut = p0 = NULL;
  bcv = NULL;
  omg = IC->get_omg();
	r = 0.0;

  unsigned int wait_num=0;
  int req[12];
	
	if( !(p   = dc_p->GetData()) )   Exit(0);
	if( !(ws  = dc_ws->GetData()) )  Exit(0);
  if( !(p0  = dc_p0->GetData()) )  Exit(0);
  if( !(bcv = dc_bcv->GetData()) ) Exit(0);
  
  // 反復処理
  switch (IC->get_LS()) {
    case SOR:
      //cpc3d_psor_(p, sz, gc, &omg, &r, ws, bnd, cut, &para_key);
      break;
      
    default:
      stamped_printf("\tInvalid Linear Solver for Pressure\n");
      Exit(0);
      break;
  }
  
  // 外部境界条件
  BC.OuterPBC(dc_p);
  if ( C.isPeriodic() == ON ) BC.InnerPBC_Periodic(dc_p, dc_bcd);
  
  // 同期処理
  switch (IC->get_LS()) {
    case SOR:
      if (cm_mode == 0 ) {
        if( !dc_p->CommBndCell(guide) ) Exit(0);
      }
      else {
        if( !dc_p->CommBndCell2(guide, wait_num, req) ) Exit(0);
        para_mng->WaitAll(wait_num, req);
      }
      break;
  }
  

}

/**
 @fn REAL_TYPE SklSolverCBC::count_comm_size(unsigned sz[3], unsigned guide) const
 @brief 全ノードについて，ローカルノード1面・一層あたりの通信量の和を返す
 @retval 通信量(Byte)
 @param sz 配列サイズ
 @param guide ガイドセル
 */
REAL_TYPE SklSolverCBC::count_comm_size(unsigned sz[3], unsigned guide) const
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  REAL_TYPE c = 0.0;
  
  // 内部面のみをカウントする
  for (unsigned n=0; n<6; n++) {
    if ( pn.nID[n] >= 0 ) {
      
      switch (n) {
        case X_MINUS:
        case X_PLUS:
          c += (REAL_TYPE)(sz[1]*sz[2]);
          break;
          
        case Y_MINUS:
        case Y_PLUS:
          c += (REAL_TYPE)(sz[0]*sz[2]);
          break;
          
        case Z_MINUS:
        case Z_PLUS:
          c += (REAL_TYPE)(sz[0]*sz[1]);
          break;
      }
    }
  }
  
  if( para_mng->IsParallel() ){
    REAL_TYPE tmp = c;
    if ( !para_mng->Allreduce(&tmp, &c, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp) ) Exit(0);
  }
  
  return c*(REAL_TYPE)sizeof(REAL_TYPE); // Byte
}

/**
 @fn void SklSolverCBC::CN_Itr(ItrCtl* IC)
 @brief 線形ソルバーの選択実行
 @param IC ItrCtlクラス
 
void SklSolverCBC::CN_Itr(ItrCtl* IC)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
	
	REAL_TYPE *vc, *wv, *wk, *v, *p0, *v0, *vf0;
	unsigned *bcv, *bcd;
  unsigned int wait_num=0;
  int req[12];
	
  vc = wv = wk = v = p0 = v0 = vf0 = NULL;
	bcv = bcd = NULL;
  REAL_TYPE tm = SklGetTotalTime();   /// 計算開始からの積算時刻
  REAL_TYPE dt = SklGetDeltaT();      /// 時間積分幅
  REAL_TYPE omg = IC->get_omg();      /// 加速係数
  REAL_TYPE rei = C.getRcpReynolds(); /// レイノルズ数の逆数
	REAL_TYPE r = 0.0;                  /// 残差
  REAL_TYPE flop_count=0.0;           /// 浮動小数演算数
  REAL_TYPE half=0.5;                 /// 定数
  REAL_TYPE machine_epsilon = (C.Mode.Precision == SPH_SINGLE) ? SINGLE_EPSILON : DOUBLE_EPSILON;
  REAL_TYPE np_f = (REAL_TYPE)para_mng->GetNodeNum(pn.procGrp); /// 全ノード数
  
	if( !(vc  = dc_vc->GetData()) )  Exit(0);
	if( !(wv  = dc_wv->GetData()) )  Exit(0);
  if( !(wk  = dc_vf0->GetData()) ) Exit(0); // assign vf0 as working array
  if( !(v   = dc_v->GetData()) )   Exit(0);
  if( !(bcv = dc_bcv->GetData()) ) Exit(0);
  if( !(bcd = dc_bcd->GetData()) ) Exit(0);
  if( !(p0  = dc_p0->GetData()) )  Exit(0);
  if( !(v0  = dc_v0->GetData()) )  Exit(0);
  if( !(vf0 = dc_vf0->GetData()) ) Exit(0);
  
  
  
  // >>> Fractional step sub-section 7
  TIMING_start(tm_frctnl_stp_sct_7);
  
  // 反復処理
  switch (IC->get_LS()) {
      
    case JACOBI:
      TIMING_start(tm_pvec_cn_Jacobi);
      flop_count = 0.0;
      //SklInitializeREAL_TYPE(dc_vf0->GetData(), 0.0, dc_vf0->GetArrayLength());
      fb_set_value_real_(dc_vf0->GetData(), dc_vf0->GetArrayLength(), 0.0);
      cbc_vis_cn_jcb_(vc, sz, gc, dh, &dt, v00, &rei, &omg, wv, (int*)bcv, wk, &half, &r, &flop_count); 
      TIMING_stop(tm_pvec_cn_Jacobi, flop_count);
      
      TIMING_start(tm_pvec_cn_mod);
      BC.mod_Vis_CN(vc, wv, half, bcv, wk, tm, dt, omg, &C, &r, IC->get_LS(), v00, flop_count);
      TIMING_stop(tm_pvec_cn_mod, flop_count);
      
      TIMING_start(tm_copy_array);
      CU.copy_REAL_TYPE(vc, vf0, dc_vc->GetArrayLength());
      TIMING_stop(tm_copy_array, 0.0);
      break;
      
    case SOR:
      TIMING_start(tm_pvec_cn_PSOR);
      flop_count = 0.0;
      cbc_vis_cn_sor_(vc, sz, gc, dh, &dt, v00, &rei, &omg, wv, (int*)bcv, &half, &r, &flop_count);
      TIMING_stop(tm_pvec_cn_PSOR, flop_count);
      
      TIMING_start(tm_pvec_cn_mod);
      BC.mod_Vis_CN(vc, wv, half, bcv, wk, tm, dt, omg, &C, &r, IC->get_LS(), v00, flop_count);
      TIMING_stop(tm_pvec_cn_mod, flop_count);
      break;
      
    case SOR2SMA: // SOR2は同期処理も行う
      TIMING_start(tm_pvec_cn_SOR2SMA);
      
      TIMING_stop(tm_pvec_cn_SOR2SMA, flop_count);
      break;
      
    default:
      printf("\tInvalid Linear Solver for Velocity_CN\n");
      Exit(0);
  }
  
  // 境界条件
  TIMING_start(tm_pvec_BC);
  flop_count=0.0;
  BC.OuterVBC(dc_vc, p0, bcv, &C, v00, flop_count);
  BC.assign_Velocity(v, bcv, tm, &C, v00);
  TIMING_stop(tm_pvec_BC, flop_count);
  
  TIMING_stop(tm_frctnl_stp_sct_7, 0.0);
  // <<< Fractional step subsection 7
  
  
  
  // >>> Fractional step sub-section 8
  TIMING_start(tm_frctnl_stp_sct_8);
  
  // 同期処理
  if ( para_mng->IsParallel() ) {
    switch (IC->get_LS()) {
      case JACOBI:
      case SOR:
        TIMING_start(tm_pvec_cn_comm);
        if (cm_mode == 0 ) {
          if( !dc_vc->CommBndCell(guide) ) Exit(0);
        }
        else {
          if( !dc_vc->CommBndCell2(guide, wait_num, req) ) Exit(0);
          para_mng->WaitAll(wait_num, req);
        }
        TIMING_stop(tm_pvec_cn_comm, (REAL_TYPE)sizeof(REAL_TYPE)*2.0, 1);
        break;
    }
  }
  
  // Residual reduction
  if ( para_mng->IsParallel() ) {
    TIMING_start(tm_pvec_cn_res_comm);
    REAL_TYPE tmp = r;
    para_mng->Allreduce(&tmp, &r, 1, SKL_ARRAY_DTYPE_REAL, SKL_MAX, pn.procGrp); // In fact, CN is MAX norm
    TIMING_stop(tm_pvec_cn_res_comm, 2.0*np_f*(REAL_TYPE)sizeof(REAL_TYPE) );
  }
  
  // 残差の保存
  IC->set_normValue( sqrt(r/(REAL_TYPE)G_Acell) );
  
  TIMING_stop(tm_frctnl_stp_sct_8, 0.0);
  // <<< Fractional step subsection 8
}*/

/**
 @fn REAL_TYPE SklSolverCBC::Norm_Poisson(ItrCtl* IC)
 @brief Poissonのノルムを計算する
 @retval 収束値
 @param IC ItrCtlクラス
 @note src1に発散がストアされている
 */
REAL_TYPE SklSolverCBC::Norm_Poisson(ItrCtl* IC)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  REAL_TYPE nrm, rms, convergence, flop_count;
  REAL_TYPE coef = SklGetDeltaT()/(C.dh*C.dh); /// 発散値を計算するための係数　dt/h^2
  REAL_TYPE np_f = (REAL_TYPE)para_mng->GetNodeNum(pn.procGrp); /// 全ノード数
  REAL_TYPE tmp;
  REAL_TYPE *src1=NULL;  /// 発散値
  unsigned *bcp=NULL;
  
  if( !(bcp = dc_bcp->GetData()) )  Exit(0);
  if( !(src1 = dc_wk2->GetData()) ) Exit(0);
  
  switch (IC->get_normType()) {
    case ItrCtl::v_div_max:
      // >>> Poisson Iteration subsection 5
      TIMING_start(tm_poi_itr_sct_5);
      
      TIMING_start(tm_norm_div_max);
      flop_count=0.0;
      cbc_norm_v_div_max_(&nrm, sz, gc, src1, &coef, (int*)bcp, &flop_count);
      TIMING_stop(tm_norm_div_max, flop_count);
      
      if ( para_mng->IsParallel() ) {
        TIMING_start(tm_norm_comm);
        tmp = nrm;
        para_mng->Allreduce(&tmp, &nrm, 1, SKL_ARRAY_DTYPE_REAL, SKL_MAX, pn.procGrp); // 最大値
        TIMING_stop(tm_norm_comm, 2.0*np_f*(REAL_TYPE)sizeof(REAL_TYPE) ); // 双方向 x ノード数
      }
      convergence = nrm;
      IC->set_normValue( convergence );
      
      TIMING_stop(tm_poi_itr_sct_5, 0.0);
      // <<< Poisson Iteration subsection 5
      break;
      
    case ItrCtl::v_div_l2:
      // >>> Poisson Iteration subsection 5
      TIMING_start(tm_poi_itr_sct_5);
      
      TIMING_start(tm_norm_div_l2);
      flop_count=0.0;
      cbc_norm_v_div_l2_(&rms, sz, gc, src1, &coef, (int*)bcp, &flop_count);
      TIMING_stop(tm_norm_div_l2, flop_count);
      
      if ( para_mng->IsParallel() ) {
        TIMING_start(tm_norm_comm);
        tmp = rms;
        para_mng->Allreduce(&tmp, &rms, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp); // 和
        TIMING_stop(tm_norm_comm, 2.0*np_f*(REAL_TYPE)sizeof(REAL_TYPE) ); // 双方向 x ノード数
      }
      convergence = sqrt(rms/np_f); //RMS
      IC->set_normValue( convergence );
      
      TIMING_stop(tm_poi_itr_sct_5, 0.0);
      // <<< Poisson Iteration subsection 5
      break;
      
    case ItrCtl::p_res_l2_a:
    case ItrCtl::p_res_l2_r:
      convergence = IC->get_normValue();
      break;
      
    case ItrCtl::v_div_max_dbg:
      // >>> Poisson Iteration subsection 5
      TIMING_start(tm_poi_itr_sct_5);
      
      TIMING_start(tm_norm_div_max_dbg);
      flop_count=0.0;
      int index[3];
      index[0] = 0;
      index[1] = 0;
      index[2] = 0;
      cbc_norm_v_div_dbg_(&nrm, &rms, index, sz, gc, src1, &coef, (int*)bcp, &flop_count);
      TIMING_stop(tm_norm_div_max_dbg, flop_count);
      
      //@todo ここで，最大値のグローバルなindexの位置を計算する
      
      if ( para_mng->IsParallel() ) {
        TIMING_start(tm_norm_comm);
        tmp = nrm;
        para_mng->Allreduce(&tmp, &nrm, 1, SKL_ARRAY_DTYPE_REAL, SKL_MAX, pn.procGrp); // 最大値
        tmp = rms;
        para_mng->Allreduce(&tmp, &rms, 1, SKL_ARRAY_DTYPE_REAL, SKL_SUM, pn.procGrp); // 和
        TIMING_stop(tm_norm_comm, 2.0*np_f*(REAL_TYPE)sizeof(REAL_TYPE)*2.0 ); // 双方向 x ノード数
      }
      rms = sqrt(rms/np_f); // RMS
      convergence = nrm; // ノルムは最大値を返す
      IC->set_normValue( convergence );
      
      Hostonly_ {
        H->printHistoryItr(fp_i, IC->LoopCount, nrm, index);
        fflush(fp_i);
      }
      
      TIMING_stop(tm_poi_itr_sct_5, 0.0);
      // <<< Poisson Iteration subsection 5
      break;
      
    default:
      stamped_printf("\tInvalid convergence type\n");
      Exit(0);
  }
  
  return convergence;
}

//@fn void SklSolverCBC::Pressure_Shift(REAL_TYPE* p)
//@brief 圧力の引き戻し操作を行う
void SklSolverCBC::Pressure_Shift(REAL_TYPE* p)
{
  REAL_TYPE base=0.0;
  
  int n = C.Mode.Pshift;
  
  cbc_face_avr_sampling_(p, sz, gc, &n, &base);
  cbc_shift_pressure_(p, sz, gc, &base);
}

/* test用のコア fortranのcbc_psor()と等価
 インラインメソッドとマクロの実行性能比較テスト
*/
REAL_TYPE SklSolverCBC::PSOR(REAL_TYPE* p, REAL_TYPE* src0, REAL_TYPE* src1, unsigned* bp, ItrCtl* IC, REAL_TYPE& flop)
{
  int i,j,k;
  unsigned m_p, m_w, m_e, m_s, m_n, m_b, m_t;
  REAL_TYPE a_p, g_w, g_e, g_s, g_n, g_b, g_t;
  REAL_TYPE t_p, t_w, t_e, t_s, t_n, t_b, t_t;
  REAL_TYPE dd, dp, ss;
  REAL_TYPE omg;
  REAL_TYPE res; // 残差の自乗和
  unsigned register s;

  omg = IC->get_omg();
  res = 0.0;
  flop += (REAL_TYPE)(ixc*jxc*kxc)* 22.0;
  
  for (k=1; k<=kxc; k++) {
    for (j=1; j<=jxc; j++) {
      for (i=1; i<=ixc; i++) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
        m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
        m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
        m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
        m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
        m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
        
        //m_p = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j  , k  );
        //m_w = F_INDEX_S3D(ixc, jxc, kxc, guide, i-1, j  , k  );
        //m_e = F_INDEX_S3D(ixc, jxc, kxc, guide, i+1, j  , k  );
        //m_s = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j-1, k  );
        //m_n = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j+1, k  );
        //m_b = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j  , k-1);
        //m_t = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j  , k+1);
        
        t_p = p[m_p];
        t_w = p[m_w];
        t_e = p[m_e];
        t_s = p[m_s];
        t_n = p[m_n];
        t_b = p[m_b];
        t_t = p[m_t];
        
        s   = bp[m_p];
        a_p = GET_SHIFT_F(s, ACTIVE_BIT);
        g_w = GET_SHIFT_F(s, BC_NDAG_W);
        g_e = GET_SHIFT_F(s, BC_NDAG_E);
        g_s = GET_SHIFT_F(s, BC_NDAG_S);
        g_n = GET_SHIFT_F(s, BC_NDAG_N);
        g_b = GET_SHIFT_F(s, BC_NDAG_B);
        g_t = GET_SHIFT_F(s, BC_NDAG_T);
        
        dd = 1.0 / (REAL_TYPE)( (s>>BC_DIAG) & 0x7 ); // 3bitを取り出す
        ss = g_w * t_w  // west  
           + g_e * t_e  // east  
           + g_s * t_s  // south 
           + g_n * t_n  // north 
           + g_b * t_b  // bottom
           + g_t * t_t  // top
           + src0[m_p] + src1[m_p];

        dp = (dd * ss - t_p) * a_p;
        p[m_p] += omg * dp;
        res += dp * dp;
      }
    }
  }
  
	return res;
}


REAL_TYPE SklSolverCBC::PSOR2sma_core(REAL_TYPE* p, int ip, int color, REAL_TYPE* src0, REAL_TYPE* src1, unsigned* bp, ItrCtl* IC, REAL_TYPE& flop)
{
  int i,j,k;
  unsigned m_p, m_w, m_e, m_s, m_n, m_b, m_t;
  REAL_TYPE a_p, g_w, g_e, g_s, g_n, g_b, g_t;
  REAL_TYPE t_p, t_w, t_e, t_s, t_n, t_b, t_t;
  REAL_TYPE dd, dp, ss;
  REAL_TYPE omg;
  REAL_TYPE res; // 残差の自乗和
  unsigned register s;
  
  omg = IC->get_omg();
  res = 0.0;
  flop += (REAL_TYPE)(ixc*jxc*kxc)* 22.0;
  
  for (k=1; k<=kxc; k++) {
    for (j=1; j<=jxc; j++) {
      for (i=1+(k+j+color+ip)%2; i<=ixc; i+=2) {
        m_p = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        m_w = FBUtility::getFindexS3D(size, guide, i-1, j  , k  );
        m_e = FBUtility::getFindexS3D(size, guide, i+1, j  , k  );
        m_s = FBUtility::getFindexS3D(size, guide, i  , j-1, k  );
        m_n = FBUtility::getFindexS3D(size, guide, i  , j+1, k  );
        m_b = FBUtility::getFindexS3D(size, guide, i  , j  , k-1);
        m_t = FBUtility::getFindexS3D(size, guide, i  , j  , k+1);
        
        //m_p = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j  , k  );
        //m_w = F_INDEX_S3D(ixc, jxc, kxc, guide, i-1, j  , k  );
        //m_e = F_INDEX_S3D(ixc, jxc, kxc, guide, i+1, j  , k  );
        //m_s = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j-1, k  );
        //m_n = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j+1, k  );
        //m_b = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j  , k-1);
        //m_t = F_INDEX_S3D(ixc, jxc, kxc, guide, i  , j  , k+1);
        
        t_p = p[m_p];
        t_w = p[m_w];
        t_e = p[m_e];
        t_s = p[m_s];
        t_n = p[m_n];
        t_b = p[m_b];
        t_t = p[m_t];
        
        s   = bp[m_p];
        a_p = GET_SHIFT_F(s, ACTIVE_BIT);
        g_w = GET_SHIFT_F(s, BC_NDAG_W);
        g_e = GET_SHIFT_F(s, BC_NDAG_E);
        g_s = GET_SHIFT_F(s, BC_NDAG_S);
        g_n = GET_SHIFT_F(s, BC_NDAG_N);
        g_b = GET_SHIFT_F(s, BC_NDAG_B);
        g_t = GET_SHIFT_F(s, BC_NDAG_T);
        
        dd = 1.0 / (REAL_TYPE)( (s>>BC_DIAG) & 0x7 ); // 3bitを取り出す
        ss = g_w * t_w  // west  
           + g_e * t_e  // east  
           + g_s * t_s  // south 
           + g_n * t_n  // north 
           + g_b * t_b  // bottom
           + g_t * t_t  // top
           - src0[m_p] + src1[m_p];
        
        dp = dd * ss - t_p;
        p[m_p] = t_p + omg * dp;
        res += dp * dp * a_p;
      }
    }
  }
  
	return res;
}

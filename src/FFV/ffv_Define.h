#ifndef _FFV_DEFINE_H_
#define _FFV_DEFINE_H_

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
 * @file   ffv_Define.h
 * @brief  FFV Class Definition Header
 * @author aics
 */

#define TM_LABEL_MAX 24
  
/** 計算性能測定キー番号 */
enum timing_key {
  tm_init_sct,
    tm_init_alloc,
  
    tm_voxel_prep_sct,
      tm_voxel_load,
      tm_polygon_load,
      tm_cutinfo,
      tm_cut_min,
  
    tm_cmp_vertex8,
    tm_cmp_subdivision,
  
    tm_restart,
  
  tm_loop_sct,
  
  tm_vmax,          
  tm_vmax_comm,
  
  // common use
  tm_copy_array,
  tm_assign_const,
  tm_blas_dot1,
  tm_blas_dot2,
  tm_blas_copy,
  tm_blas_calc_r,
  tm_blas_ax,
  tm_blas_z_xpay,
  tm_blas_axpbypz,
  tm_blas_bicg_update_x,
  tm_blas_bicg_update_p,
  tm_blas_comm,
  
  tm_bicgstab_sct,
  
  tm_flow_sct,
    tm_frctnl_stp_sct,
  
      tm_frctnl_stp_sct_1,       
        tm_spec_vel,
        tm_WallFunc,
  
      tm_frctnl_stp_sct_2,
        tm_pseudo_vec,
        tm_pvec_flux,
        tm_pvec_ee,
        tm_pvec_ab,
        tm_pvec_abcn,
        tm_pvec_abcn_df_ee,
        tm_pvec_abcn_df_ee_BC,
  
      tm_frctnl_stp_sct_3,
        tm_forcing,
        tm_buoyancy,
        tm_pvec_BC,
        tm_pvec_comm,
  
      tm_frctnl_stp_sct_4,
  
    tm_poi_src_sct,
      tm_div_pvec,
      tm_poi_src_vbc,
      tm_poi_src_nrm,
      tm_poi_src_comm,
  
    tm_poi_itr_sct,
      tm_hstry_itr,
  
      tm_poi_itr_sct_1,
        tm_force_src,
  
      tm_poi_itr_sct_2,     
        tm_poi_PSOR,
        tm_poi_setup,
        tm_poi_SOR2SMA,
        tm_poi_BC,
  
      tm_poi_itr_sct_3,
        tm_poi_comm,
        tm_poi_res_comm,
  
      tm_poi_itr_sct_4,
        tm_prj_vec,
        tm_prj_vec_bc,
        tm_prj_vec_bc_comm,
        tm_prj_frc_mod,
        tm_prj_frc_mod_comm,
        tm_prj_frc_dir,
        tm_vec_BC,
  
      tm_poi_itr_sct_5,
        tm_norm_div_max,
        tm_norm_comm,
  
      tm_gmres_sor_sct,
        tm_gmres_mvprod,
        tm_gmres_res_sample,
        tm_gmres_others,
        tm_gmres_comm,
  
    tm_NS_loop_post_sct,
      tm_vectors_comm,
      tm_VBC_update,
      tm_LES_eddy,
      tm_LES_eddy_comm,
      tm_domain_monitor,


  tm_heat_sct,
  
    tm_heat_convection_sct,
      tm_heat_spec_temp,
      tm_heat_cnv,
      tm_heat_cnv_BC,
      tm_heat_cnv_EE,
  
    tm_heat_diffusion_sct,
      tm_heat_diff_OBC,
  
      tm_heat_diff_sct_1,
        tm_heat_diff_IBC_vol,
        tm_heat_diff_comm,
  
      tm_heat_diff_sct_2,
        tm_heat_diff_IBC_face,
        tm_heat_diff_OBC_face,
        tm_heat_diff_QBC_comm,
  
      tm_heat_diff_sct_3,
        tm_heat_diff_EE,
        tm_heat_diff_PSOR,
        tm_heat_update_comm,
        tm_heat_diff_res_comm,
  
      tm_heat_loop_post_sct,
        tm_heat_range,
  
  
  tm_vof_sct,
    tm_vof_cnv,
    tm_vof_cnv_comm,
  
  tm_loop_uty_sct,

    tm_loop_uty_sct_1,
      tm_average_time,
      tm_stat_space,
      tm_stat_space_comm,
  
    tm_loop_uty_sct_2,
      tm_hstry_stdout,
      tm_file_out,
      tm_hstry_base,
        tm_cal_force,
      tm_hstry_wall,
      tm_total_prs,
      tm_sampling,
      tm_hstry_sampling,
  
  tm_statistic,
  
  tm_END
};  

#endif // _FFV_DEFINE_H_

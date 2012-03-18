/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file FortranFuncCBC.h
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

#include "FBDefine.h"

#ifndef _SKL_FORTRAN_FUNC_CBC_H_
#define _SKL_FORTRAN_FUNC_CBC_H_

#ifdef _WIN32
// cbc_poisson.f90
#define cbc_div_cnst_           CBC_DIV_CNST
#define cbc_psor_               CBC_PSOR
#define cbc_psor2sma_core_      CBC_PSOR2SMA_CORE
#define cbc_sma_comm_           CBC_SMA_COMM
#define cbc_sma_comm_wait_      CBC_SMA_COMM_WAIT

// cbc_utility.f90
#define cbc_norm_v_div_dbg_     CBC_NORM_V_DIV_DBG
#define cbc_norm_v_div_l2_      CBC_NORM_V_DIV_L2
#define cbc_norm_v_div_max_     CBC_NORM_V_DIV_MAX
#define cbc_helicity_           CBC_HELICITY
#define cbc_i2vgt_              CBC_I2VGT
#define cbc_rot_v_              CBC_ROT_V
#define cbc_vmax_               CBC_VMAX

// BCvec_cc.f90
#define cbc_vibc_drchlt_        CBC_VIBC_DRCHLT
#define cbc_vibc_outflow_       CBC_VIBC_OUTFLOW
#define cbc_pvec_vibc_specv_    CBC_PVEC_VIBC_SPECV
#define cbc_pvec_vibc_oflow_    CBC_PVEC_VIBC_OFLOW
#define cbc_pvec_vobc_specv_    CBC_PVEC_VOBC_SPECV
#define cbc_pvec_vobc_oflow_    CBC_PVEC_VOBC_OFLOW
#define cbc_pvec_vobc_wall_     CBC_PVEC_VOBC_WALL
#define cbc_pvec_vobc_symtrc_   CBC_PVEC_VOBC_SYMTRC
#define cbc_vobc_drchlt_        CBC_VOBC_DRCHLT
#define cbc_vobc_outflow_       CBC_VOBC_OUTFLOW
#define cbc_vobc_tfree_         CBC_VOBC_TFREE
#define cbc_vobc_update_        CBC_VOBC_UPDATE
#define cbc_div_ibc_drchlt_     CBC_DIV_IBC_DRCHLT
#define cbc_div_ibc_oflow_pvec_ CBC_DIV_IBC_OFLOW_PVEC
#define cbc_div_ibc_oflow_vec_  CBC_DIV_IBC_OFLOW_VEC
#define cbc_div_obc_drchlt_     CBC_DIV_OBC_DRCHLT
#define cbc_div_obc_oflow_pvec_ CBC_DIV_OBC_OFLOW_PVEC
#define cbc_div_obc_oflow_vec_  CBC_DIV_OBC_OFLOW_VEC

// BCprs_cc.f90
#define cbc_pobc_drchlt_        CBC_POBC_DRCHLT
#define cbc_pobc_neumann_       CBC_POBC_NEUMANN

// cbc_3d.f90
#define cbc_ab2_                CBC_AB2
#define cbc_div_                CBC_DIV
#define cbc_eddy_viscosity_     CBC_EDDY_VISCOSITY
#define cbc_ee_                 CBC_EE
#define cbc_friction_velocity_  CBC_FRICTION_VELOCITY
#define cbc_pvec_muscl_         CBC_PVEC_MUSCL
#define cbc_pvec_les_           CBC_PVEC_LES
#define cbc_update_vec_         CBC_UPDATE_VEC
#define cbc_vis_cn_jcb_         CBC_VIS_CN_JCB
#define cbc_vis_cn_mod_jcb_     CBC_VIS_CN_MOD_JCB
#define cbc_vis_cn_mod_sor_     CBC_VIS_CN_MOD_SOR
#define cbc_vis_cn_sor_         CBC_VIS_CN_SOR
#define cbc_vis_ee_             CBC_VIS_EE
#define cbc_vis_ee_vbc_         CBC_VIS_EE_VBC

// cbc_forcing.f90
#define cbc_force_keep_vec_     CBC_FORCE_KEEP_VEC
#define cbc_hex_psrc_           CBC_HEX_PSRC
#define cbc_hex_force_pvec_     CBC_HEX_FORCE_PVEC
#define cbc_hex_force_vec_      CBC_HEX_FORCE_VEC

// cbc_pscalar.f90
#define cbc_ps_muscl_           CBC_PS_MUSCL
#define cbc_ps_buoyancy_        CBC_PS_BUOYANCY
#define cbc_ps_diff_ee_         CBC_PS_DIFF_EE
#define cbc_hbc_drchlt_         CBC_HBC_DRCHLT

// cds_poisson.f90
#define cds_div_cc_             CDS_DIV_CC
#define cds_div_cf_             CDS_DIV_CF
#define cds_psor_               CDS_PSOR

// c3d_vof.f90
#define vof_uwd_                VOF_UWD
#define vof_muscl_              VOF_MUSCL

// cds_vector.f90
#define cds_div_                CDS_DIV
#define cds_pvec_muscl          CDS_PVEC_MUSCL
#define cds_update_vec_         CDS_UPDATE_VEC

#endif // _WIN32


extern "C" {
  // cbc_poisson.f90
  void cbc_div_cnst_        (REAL_TYPE* dv, int* sz, int* g, REAL_TYPE* b2, int* bp, REAL_TYPE* flop);
  void cbc_psor_            (REAL_TYPE* p,  int* sz, int* g, REAL_TYPE* omg, REAL_TYPE* res, REAL_TYPE* s0, REAL_TYPE* s1, 
                             int* bp, REAL_TYPE* flop);
  void cbc_psor2sma_core_   (REAL_TYPE* p,  int* sz, int* g, int* ip, int* color, REAL_TYPE* omg, REAL_TYPE* res, REAL_TYPE* s0, 
                             REAL_TYPE* s1, int* bp, REAL_TYPE* flop);
  void cbc_sma_comm_        (REAL_TYPE* p, int* sz, int* g, int* col, int* ip, int* cf_sz, REAL_TYPE* cf_x, REAL_TYPE* cf_y, 
                             REAL_TYPE* cf_z, int* key, int* para_key);
  void cbc_sma_comm_wait_   (REAL_TYPE* p, int* sz, int* g, int* col, int* ip, int* cf_sz, REAL_TYPE* cf_x, REAL_TYPE* cf_y, 
                             REAL_TYPE* cf_z, int* key);
  
  // cbc_utility.f90
  void cbc_norm_v_div_dbg_ (REAL_TYPE* ds, REAL_TYPE* rm, int* index, int* sz, int* g, REAL_TYPE* div, REAL_TYPE* coef, int* bp, REAL_TYPE* flop);
  void cbc_norm_v_div_l2_  (REAL_TYPE* rms,   int* sz, int* g, REAL_TYPE* div, REAL_TYPE* coef, int* bp, REAL_TYPE* flop);
  void cbc_norm_v_div_max_ (REAL_TYPE* ds,    int* sz, int* g, REAL_TYPE* div, REAL_TYPE* coef, int* bp, REAL_TYPE* flop);
  void cbc_helicity_       (REAL_TYPE* ht,    int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v, int* bv, REAL_TYPE* v00, REAL_TYPE* flop);
  void cbc_i2vgt_          (REAL_TYPE* q,     int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v, int* bv, REAL_TYPE* v00, REAL_TYPE* flop);
  void cbc_rot_v_          (REAL_TYPE* rot,   int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v, int* bv, REAL_TYPE* v00, REAL_TYPE* flop);
  void cbc_vmax_           (REAL_TYPE* v_max, int* sz, int* g, REAL_TYPE* v00, REAL_TYPE* v, REAL_TYPE* flop);
  
  // BCvec_cc.f90
  void cbc_vibc_drchlt_       (REAL_TYPE* v, int* sz, int* g, int* st, int* ed, REAL_TYPE* v00, int* bv, int* odr, REAL_TYPE* vec);
  void cbc_vibc_outflow_      (REAL_TYPE* vc,int* sz, int* g, int* st, int* ed, int* bv, int* odr);
  void cbc_pvec_vibc_specv_   (REAL_TYPE* wv, int* sz, int* g, int* st, int* ed, REAL_TYPE* dh, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v, 
                               int* bv, int* odr, REAL_TYPE* vec, int* v_mode, REAL_TYPE* flop);
  void cbc_pvec_vibc_oflow_   (REAL_TYPE* wv, int* sz, int* g, int* st, int* ed, REAL_TYPE* dh, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v, 
                               int* bv, int* odr, REAL_TYPE* vec, int* v_mode, REAL_TYPE* flop);
  void cbc_pvec_vobc_specv_   (REAL_TYPE* wv, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v,
                               int* bv, REAL_TYPE* vec, int* v_mode, int* face, REAL_TYPE* flop);
  void cbc_pvec_vobc_oflow_   (REAL_TYPE* wv, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v,
                               int* bv, REAL_TYPE* vec, int* v_mode, int* face, REAL_TYPE* flop);
  void cbc_pvec_vobc_wall_    (REAL_TYPE* wv, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v,
                               int* bv, REAL_TYPE* vec, int* v_mode, int* face, REAL_TYPE* flop);
  void cbc_pvec_vobc_symtrc_  (REAL_TYPE* wv, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v,
                               int* bv, REAL_TYPE* vec, int* v_mode, int* face, REAL_TYPE* flop);
  void cbc_vobc_drchlt_       (REAL_TYPE* v, int* sz, int* g, REAL_TYPE* v00, int* bv, int* face, REAL_TYPE* vec);
  void cbc_vobc_outflow_      (REAL_TYPE* v, int* sz, int* g, REAL_TYPE* c, int* bv, int* face, REAL_TYPE* v0, REAL_TYPE* flop);
  void cbc_vobc_tfree_        (REAL_TYPE* v, int* sz, int* g, int* face, int* bv, REAL_TYPE* v00, REAL_TYPE* flop);
  void cbc_vobc_update_       (REAL_TYPE* v, int* sz, int* g, REAL_TYPE*vc, int* face);
  void cbc_div_ibc_drchlt_    (REAL_TYPE* div, int* sz, int* g, int* st, int* ed, REAL_TYPE* v00, REAL_TYPE* coef, int* bv, int* odr, 
                               REAL_TYPE* vec, REAL_TYPE* flop);
  void cbc_div_ibc_oflow_pvec_(REAL_TYPE* div, int* sz, int* g, int* st, int* ed, REAL_TYPE* v00, REAL_TYPE* v_cnv, REAL_TYPE* coef, 
                               int* bv, int* odr, REAL_TYPE* v0, REAL_TYPE* flop);
  void cbc_div_ibc_oflow_vec_ (REAL_TYPE* div, int* sz, int* g, int* st, int* ed, REAL_TYPE* v00, REAL_TYPE* coef, int* bv, int* odr, 
                               REAL_TYPE* av, REAL_TYPE* flop);
  void cbc_div_obc_drchlt_    (REAL_TYPE* div, int* sz, int* g, int* face, REAL_TYPE* v00, REAL_TYPE* coef, int* bv, REAL_TYPE* vec, REAL_TYPE* flop);
  void cbc_div_obc_oflow_pvec_(REAL_TYPE* div, int* sz, int* g, int* face, REAL_TYPE* v00, REAL_TYPE* v_out, REAL_TYPE* coef, int* bv, 
                               REAL_TYPE* v0, REAL_TYPE* flop);
  //void cbc_div_obc_oflow_pvec_(REAL_TYPE* div, int* sz, int* g, int* face, REAL_TYPE* v00, REAL_TYPE* coef, int* bv, REAL_TYPE* v0, REAL_TYPE* flop);
  void cbc_div_obc_oflow_vec_ (REAL_TYPE* div, int* sz, int* g, int* face, REAL_TYPE* v00, REAL_TYPE* coef, int* bv, REAL_TYPE* aa, REAL_TYPE* flop);
  
  // BCprs_cc.f90
  void cbc_pobc_drchlt_  (REAL_TYPE* p, int* sz, int* g, int* face, REAL_TYPE* pv);
  void cbc_pobc_neumann_ (REAL_TYPE* p, int* sz, int* g, int* face);
  
  // cbc_3d.f90
  void cbc_ab2_               (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* v, REAL_TYPE* ab, int* bd, REAL_TYPE* v00, REAL_TYPE* flop);
  void cbc_div_               (REAL_TYPE* dv, int* sz, int* g, REAL_TYPE* coef, REAL_TYPE* vc, int* bv, REAL_TYPE* v00, REAL_TYPE* flop);
  void cbc_eddy_viscosity_    (REAL_TYPE* vt, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* re, REAL_TYPE* cs, REAL_TYPE* v, int* bx, 
                               REAL_TYPE* vt_range, REAL_TYPE* yp_range, REAL_TYPE* v00);
  void cbc_ee_                (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* v, int* bd, REAL_TYPE* flop);
  void cbc_friction_velocity_ (REAL_TYPE* ut, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* re, REAL_TYPE* v, int* bp, 
                               REAL_TYPE* range_Yp, REAL_TYPE* range_Ut, REAL_TYPE* v00, REAL_TYPE* flop);
  void cbc_pvec_muscl_        (REAL_TYPE* wv, int* sz, int* g, REAL_TYPE* dh, int* c_scheme, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v, 
                               int* bv, int* bp, int* v_mode, REAL_TYPE* ut, int* wall_type, int* bd, float* cvf, REAL_TYPE* flop);
  void cbc_update_vec_        (REAL_TYPE* v, REAL_TYPE* div, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* vc, REAL_TYPE* p, int* bp, int* bv, 
                               REAL_TYPE* v00, REAL_TYPE* flop);
  void cbc_vis_cn_jcb_        (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                               REAL_TYPE* omg, REAL_TYPE* wv, int* bx, REAL_TYPE* wk, REAL_TYPE* coef, REAL_TYPE* res, REAL_TYPE* flop);
  void cbc_vis_cn_mod_jcb_    (REAL_TYPE* vc, int* sz, int* g, int* st, int* ed, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                               REAL_TYPE* omg, REAL_TYPE* wv, int* bx, REAL_TYPE* wk, REAL_TYPE* coef, REAL_TYPE* res, REAL_TYPE* vec, REAL_TYPE* flop);
  void cbc_vis_cn_mod_sor_    (REAL_TYPE* vc, int* sz, int* g, int* st, int* ed, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                               REAL_TYPE* omg, REAL_TYPE* wv, int* bx, REAL_TYPE* coef, REAL_TYPE* res, REAL_TYPE* vec, REAL_TYPE* flop);
  void cbc_vis_cn_sor_        (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                               REAL_TYPE* omg, REAL_TYPE* wv, int* bx, REAL_TYPE* coef, REAL_TYPE* res, REAL_TYPE* flop);
  void cbc_vis_ee_            (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* wv, 
                               REAL_TYPE* v, int* bx, REAL_TYPE* coef, REAL_TYPE* flop);
  void cbc_vis_ee_vbc_        (REAL_TYPE* vc, int* sz, int* g, int* st, int* ed, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                               REAL_TYPE* v0, int* bx, int* odr, REAL_TYPE* coef, REAL_TYPE* vec, REAL_TYPE* flop);

  // cbc_forcing.f90
  void cbc_force_keep_vec_    (REAL_TYPE* wk, int* c_sz, int* st, int* ed, REAL_TYPE* v, int* sz, int* g);
  void cbc_hex_psrc_          (REAL_TYPE* src, int* sz, int* g, int* st, int* ed, int* bd, float* vf, REAL_TYPE* wk, int* c_sz, int* odr, 
                               REAL_TYPE* v00, REAL_TYPE* dh, REAL_TYPE* nv, REAL_TYPE* c, REAL_TYPE* flop);
  void cbc_hex_force_pvec_    (REAL_TYPE* vc,  int* sz, int* g, int* st, int* ed, int* bd, REAL_TYPE* vf, REAL_TYPE* v, int* odr, 
                               REAL_TYPE* v00, REAL_TYPE* dt, REAL_TYPE* nv, REAL_TYPE* c, REAL_TYPE* flop);
  void cbc_hex_force_vec_     (REAL_TYPE* v, REAL_TYPE* div, int* sz, int* g, int* st, int* ed, int* bd, REAL_TYPE* vf, REAL_TYPE* wk, int* c_sz, int* odr, 
                               REAL_TYPE* v00, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* nv, REAL_TYPE* c, REAL_TYPE* am, REAL_TYPE* flop);
  
  // cbc_pscalar.f90
  void cbc_ps_muscl_      (REAL_TYPE* ws, int* sz, int* g, REAL_TYPE* dh, int* c_scheme, REAL_TYPE* v00, REAL_TYPE* v, REAL_TYPE* t, int* bv, 
                           int* bh1, int* bh2, int* swt, REAL_TYPE* flop);
  void cbc_ps_buoyancy_   (REAL_TYPE* v, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* gr, REAL_TYPE* rei, REAL_TYPE* t, int* bd, REAL_TYPE* flop);
  void cbc_ps_diff_ee_    (REAL_TYPE* t, int* sz, int* g, REAL_TYPE* res, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* pei, REAL_TYPE* qbc, 
                           int* bh, REAL_TYPE* ws, REAL_TYPE* flop);
  void cbc_hbc_drchlt_    (REAL_TYPE* t,  int* sz, int* g, int* st, int* ed, int* bh, int* odr, REAL_TYPE* tc);
  
  // cds_poisson.f90
  void cds_div_cc_        (REAL_TYPE* div, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* b2, REAL_TYPE* vc, REAL_TYPE* bnd, 
                          REAL_TYPE* cut, REAL_TYPE* v00, int* mode);
  void cds_div_cf_        (REAL_TYPE* div, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* b2, REAL_TYPE* bnd, 
                          REAL_TYPE* cut, REAL_TYPE* v00, int* mode);
  void cds_psor_          (REAL_TYPE* p,   int* sz, int* g, REAL_TYPE* omg, REAL_TYPE* res, REAL_TYPE* div, REAL_TYPE* bnd, REAL_TYPE* cut, 
                          REAL_TYPE* epsilon, int* para_key);
  
  // c3d_vof.f90
  void vof_uwd_    (REAL_TYPE* f, int* sz, int* g, REAL_TYPE* v00, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* v, REAL_TYPE* q, int* bx, REAL_TYPE* flop);
  void vof_muscl_  (REAL_TYPE* f, int* sz, int* g, REAL_TYPE* v00, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* v, REAL_TYPE* q, int* bx, REAL_TYPE* flop);
  
  // cds_vector.f90
  void cds_div_           (REAL_TYPE* div, int* sz, int* g, REAL_TYPE* coef, REAL_TYPE* v, int* bv, float* cut, REAL_TYPE* v00, REAL_TYPE* flop);
  void cds_pvec_muscl_    (REAL_TYPE* wv, int* sz, int* g, REAL_TYPE* dh, int* c_scheme, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* v, 
                           int* bv, int* bp, int* v_mode, float* cut, REAL_TYPE* flop);
  void cds_update_vec_    (REAL_TYPE* v, REAL_TYPE* div, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* vc, REAL_TYPE* p, 
                           int* bp, int* bv, float* cut, REAL_TYPE* v00, REAL_TYPE* coef, REAL_TYPE* flop);
}

#endif // _SKL_FORTRAN_FUNC_CBC_H_

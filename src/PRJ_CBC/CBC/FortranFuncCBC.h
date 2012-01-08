/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file FortranFuncCBC.h
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"

#ifndef _SKL_FORTRAN_FUNC_CBC_H_
#define _SKL_FORTRAN_FUNC_CBC_H_

#ifdef _WIN32
// cbc_poisson.f90
#define cbc_div_cnst_           CBC_DIV_CNST
#define cbc_jacobi_             CBC_JACOBI
#define cbc_jacobi_if_          CBC_JACOBI_IF
#define cbc_psor_               CBC_PSOR
#define cbc_psor_if_            CBC_PSOR_IF
#define cbc_psor_index3_        CBC_PSOR_INDEX3
#define cbc_psor_index_         CBC_PSOR_INDEX
#define cbc_psor2sma_core_if_   CBC_PSOR2SMA_CORE_IF
#define cbc_psor2sma_core_      CBC_PSOR2SMA_CORE
#define cbc_sma_comm_           CBC_SMA_COMM
#define cbc_sma_comm_wait_      CBC_SMA_COMM_WAIT
#define cbc_psrc_pbc_           CBC_PSRC_PBC

// cbc_utility.f90
#define cbc_norm_v_div_max_     CBC_NORM_V_DIV_MAX
#define cbc_helicity_           CBC_HELICITY
#define cbc_i2vgt_              CBC_I2VGT
#define cbc_rot_v_              CBC_ROT_V
#define cbc_vmax_               CBC_VMAX

// BCvec_cc.f90
#define cbc_vibc_drchlt_        CBC_VIBC_DRCHLT
#define cbc_vibc_outflow_       CBC_VIBC_OUTFLOW
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
#define cbc_pvec_vibc_          CBC_PVEC_VIBC
#define cbc_pvec_vobc_          CBC_PVEC_VOBC
#define cbc_update_vec_         CBC_UPDATE_VEC
#define cbc_vis_cn_jcb_         CBC_VIS_CN_JCB
#define cbc_vis_cn_mod_jcb_     CBC_VIS_CN_MOD_JCB
#define cbc_vis_cn_mod_sor_     CBC_VIS_CN_MOD_SOR
#define cbc_vis_cn_sor_         CBC_VIS_CN_SOR
#define cbc_vis_ee_             CBC_VIS_EE
#define cbc_vis_ee_vbc_         CBC_VIS_EE_VBC

// cbc_forcing.f90
#define cbc_psrc_hex_           CBC_PSRC_HEX
#define cbc_pvec_hex_           CBC_PVEC_HEX
#define cbc_update_vcc_hex_     CBC_UPDATE_VCC_HEX
#define cbc_update_vcf_hex_     CBC_UPDATE_VCF_HEX

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
  void cbc_div_cnst_        (SKL_REAL* dv, int* sz, int* g, SKL_REAL* b2, int* bp, SKL_REAL* flop);
  void cbc_jacobi_          (SKL_REAL* p,  int* sz, int* g, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* s0, SKL_REAL* s1, 
                             int* bp, SKL_REAL* wk, SKL_REAL* flop);
  void cbc_jacobi_if_       (SKL_REAL* p,  int* sz, int* g, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* s0, SKL_REAL* s1, 
                             int* bp, SKL_REAL* wk, SKL_REAL* flop);
  void cbc_psor_            (SKL_REAL* p,  int* sz, int* g, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* s0, SKL_REAL* s1, 
                             int* bp, SKL_REAL* flop);
  void cbc_psor_if_         (SKL_REAL* p,  int* sz, int* g, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* s0, SKL_REAL* s1, 
                             int* bp, SKL_REAL* flop);
  void cbc_psor_index3_     (SKL_REAL* p,  int* sz, int* g, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* s0, SKL_REAL* s1, 
                             int* bp, int* index, int* idx_sz, SKL_REAL* flop);
  void cbc_psor_index_      (SKL_REAL* p,  int* sz, int* g, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* s0, SKL_REAL* s1, 
                             int* bp, int* index, int* idx_sz, SKL_REAL* flop);
  void cbc_psor2sma_core_   (SKL_REAL* p,  int* sz, int* g, int* ip, int* color, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* s0, 
                             SKL_REAL* s1, int* bp, SKL_REAL* flop);
	void cbc_psor2sma_core_if_(SKL_REAL* p,  int* sz, int* g, int* ip, int* color, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* s0, 
                             SKL_REAL* s1, int* bp, SKL_REAL* flop);
  void cbc_sma_comm_        (SKL_REAL* p, int* sz, int* g, int* col, int* ip, int* cf_sz, SKL_REAL* cf_x, SKL_REAL* cf_y, 
                             SKL_REAL* cf_z, int* key, int* para_key);
  void cbc_sma_comm_wait_   (SKL_REAL* p, int* sz, int* g, int* col, int* ip, int* cf_sz, SKL_REAL* cf_x, SKL_REAL* cf_y, 
                             SKL_REAL* cf_z, int* key);
  
  // cbc_utility.f90
  void cbc_norm_v_div_max_ (SKL_REAL* ds,    int* sz, int* g, SKL_REAL* div, SKL_REAL* coef, int* bp, SKL_REAL* flop);
  void cbc_helicity_       (SKL_REAL* ht,    int* sz, int* g, SKL_REAL* dh, SKL_REAL* v, int* bv, SKL_REAL* v00, SKL_REAL* flop);
  void cbc_i2vgt_          (SKL_REAL* q,     int* sz, int* g, SKL_REAL* dh, SKL_REAL* v, int* bv, SKL_REAL* v00, SKL_REAL* flop);
  void cbc_rot_v_          (SKL_REAL* rot,   int* sz, int* g, SKL_REAL* dh, SKL_REAL* v, int* bv, SKL_REAL* v00, SKL_REAL* flop);
  void cbc_vmax_           (SKL_REAL* v_max, int* sz, int* g, SKL_REAL* v00, SKL_REAL* v, SKL_REAL* flop);
  
  // BCvec_cc.f90
  void cbc_vibc_drchlt_       (SKL_REAL* v, int* sz, int* g, int* st, int* ed, SKL_REAL* v00, int* bv, int* odr, SKL_REAL* vec);
  void cbc_vibc_outflow_      (SKL_REAL* vc,int* sz, int* g, int* st, int* ed, int* bv, int* odr);
  void cbc_vobc_drchlt_       (SKL_REAL* v, int* sz, int* g, SKL_REAL* v00, int* bv, int* face, SKL_REAL* vec);
  void cbc_vobc_outflow_      (SKL_REAL* v, int* sz, int* g, SKL_REAL* c, int* bv, int* face, SKL_REAL* v0, SKL_REAL* flop);
  void cbc_vobc_tfree_        (SKL_REAL* v, int* sz, int* g, int* face, int* bv, SKL_REAL* v00, SKL_REAL* flop);
  void cbc_vobc_update_       (SKL_REAL* v, int* sz, int* g, SKL_REAL*vc, int* face);
  void cbc_div_ibc_drchlt_    (SKL_REAL* div, int* sz, int* g, int* st, int* ed, SKL_REAL* v00, SKL_REAL* coef, int* bv, int* odr, 
                               SKL_REAL* vec, SKL_REAL* flop);
  void cbc_div_ibc_oflow_pvec_(SKL_REAL* div, int* sz, int* g, int* st, int* ed, SKL_REAL* v00, SKL_REAL* v_cnv, SKL_REAL* coef, 
                               int* bv, int* odr, SKL_REAL* v0, SKL_REAL* flop);
  void cbc_div_ibc_oflow_vec_ (SKL_REAL* div, int* sz, int* g, int* st, int* ed, SKL_REAL* v00, SKL_REAL* coef, int* bv, int* odr, 
                               SKL_REAL* avr, SKL_REAL* flop);
  void cbc_div_obc_drchlt_    (SKL_REAL* div, int* sz, int* g, int* face, SKL_REAL* v00, SKL_REAL* coef, int* bv, SKL_REAL* vec, SKL_REAL* flop);
  void cbc_div_obc_oflow_pvec_(SKL_REAL* div, int* sz, int* g, int* face, SKL_REAL* v00, SKL_REAL* v_out, SKL_REAL* coef, int* bv, 
                               SKL_REAL* v0, SKL_REAL* flop);
  //void cbc_div_obc_oflow_pvec_(SKL_REAL* div, int* sz, int* g, int* face, SKL_REAL* v00, SKL_REAL* coef, int* bv, SKL_REAL* v0, SKL_REAL* flop);
  void cbc_div_obc_oflow_vec_ (SKL_REAL* div, int* sz, int* g, int* face, SKL_REAL* v00, SKL_REAL* coef, int* bv, SKL_REAL* aa, SKL_REAL* flop);
  
  // BCprs_cc.f90
  void cbc_pobc_drchlt_  (SKL_REAL* p, int* sz, int* g, int* face, SKL_REAL* pv);
  void cbc_pobc_neumann_ (SKL_REAL* p, int* sz, int* g, int* face);
  
  // cbc_3d.f90
  void cbc_ab2_               (SKL_REAL* vc, int* sz, int* g, SKL_REAL* dt, SKL_REAL* v, SKL_REAL* ab, int* bd, SKL_REAL* v00, SKL_REAL* flop);
  void cbc_div_               (SKL_REAL* dv, int* sz, int* g, SKL_REAL* coef, SKL_REAL* vc, int* bv, SKL_REAL* v00, SKL_REAL* flop);
  void cbc_eddy_viscosity_    (SKL_REAL* vt, int* sz, int* g, SKL_REAL* dh, SKL_REAL* re, SKL_REAL* cs, SKL_REAL* v, int* bx, 
                               SKL_REAL* vt_range, SKL_REAL* yp_range, SKL_REAL* v00);
  void cbc_ee_                (SKL_REAL* vc, int* sz, int* g, SKL_REAL* dt, SKL_REAL* v, int* bd, SKL_REAL* flop);
  void cbc_friction_velocity_ (SKL_REAL* ut, int* sz, int* g, SKL_REAL* dh, SKL_REAL* re, SKL_REAL* v, int* bp, 
                               SKL_REAL* range_Yp, SKL_REAL* range_Ut, SKL_REAL* v00, SKL_REAL* flop);
  void cbc_pvec_muscl_        (SKL_REAL* wv, int* sz, int* g, SKL_REAL* dh, int* c_scheme, SKL_REAL* v00, SKL_REAL* rei, SKL_REAL* v, 
                               int* bv, int* bp, int* v_mode, SKL_REAL* ut, int* wall_type, SKL_REAL* flop);
  void cbc_pvec_vibc_         (SKL_REAL* wv, int* sz, int* g, int* st, int* ed, SKL_REAL* dh, SKL_REAL* v00, SKL_REAL* rei, SKL_REAL* v, 
                               int* bv, int* odr, SKL_REAL* vec, int* v_mode, int* ofi, SKL_REAL* flop);
  void cbc_pvec_vobc_         (SKL_REAL* wv, int* sz, int* g, SKL_REAL* dh, SKL_REAL* v00, SKL_REAL* rei, SKL_REAL* v,
                               int* bv, SKL_REAL* vec, int* v_mode, int* ofi, int* face, SKL_REAL* flop);
  void cbc_update_vec_        (SKL_REAL* v, SKL_REAL* div, int* sz, int* g, SKL_REAL* dt, SKL_REAL* dh, SKL_REAL* vc, SKL_REAL* p, int* bp, int* bv, 
                               SKL_REAL* v00, SKL_REAL* coef, SKL_REAL* flop);
  void cbc_vis_cn_jcb_        (SKL_REAL* vc, int* sz, int* g, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* v00, SKL_REAL* rei, 
                               SKL_REAL* omg, SKL_REAL* wv, int* bx, SKL_REAL* wk, SKL_REAL* coef, SKL_REAL* res, SKL_REAL* flop);
  void cbc_vis_cn_mod_jcb_    (SKL_REAL* vc, int* sz, int* g, int* st, int* ed, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* v00, SKL_REAL* rei, 
                               SKL_REAL* omg, SKL_REAL* wv, int* bx, SKL_REAL* wk, SKL_REAL* coef, SKL_REAL* res, SKL_REAL* vec, SKL_REAL* flop);
  void cbc_vis_cn_mod_sor_    (SKL_REAL* vc, int* sz, int* g, int* st, int* ed, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* v00, SKL_REAL* rei, 
                               SKL_REAL* omg, SKL_REAL* wv, int* bx, SKL_REAL* coef, SKL_REAL* res, SKL_REAL* vec, SKL_REAL* flop);
  void cbc_vis_cn_sor_        (SKL_REAL* vc, int* sz, int* g, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* v00, SKL_REAL* rei, 
                               SKL_REAL* omg, SKL_REAL* wv, int* bx, SKL_REAL* coef, SKL_REAL* res, SKL_REAL* flop);
  void cbc_vis_ee_            (SKL_REAL* vc, int* sz, int* g, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* v00, SKL_REAL* rei, SKL_REAL* wv, 
                               SKL_REAL* v, int* bx, SKL_REAL* coef, SKL_REAL* flop);
  void cbc_vis_ee_vbc_        (SKL_REAL* vc, int* sz, int* g, int* st, int* ed, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* v00, SKL_REAL* rei, 
                               SKL_REAL* v0, int* bx, int* odr, SKL_REAL* coef, SKL_REAL* vec, SKL_REAL* flop);

  // cbc_forcing.f90
  void cbc_psrc_hex_          (SKL_REAL* src,int* sz, int* g, int* st, int* ed, SKL_REAL* dh, int* bd, int* odr, SKL_REAL* v00, 
                               SKL_REAL* nv, SKL_REAL* c, SKL_REAL* v, SKL_REAL* flop);
  void cbc_pvec_hex_          (SKL_REAL* v,  int* sz, int* g, int* st, int* ed, int* bd, int* odr, SKL_REAL* v00, SKL_REAL* vec, SKL_REAL* flop);
  void cbc_update_vcc_hex_    (SKL_REAL* v,  int* sz, int* g, int* st, int* ed, SKL_REAL* dt, int* bd, int* odr, SKL_REAL* v00, 
                               SKL_REAL* nv, SKL_REAL* c, SKL_REAL* vm, SKL_REAL* flop);
  void cbc_update_vcf_hex_    (SKL_REAL* vf, int* sz, int* g, int* st, int* ed, SKL_REAL* dt, int* bd, int* odr, SKL_REAL* v00, 
                               SKL_REAL* nv, SKL_REAL* c, SKL_REAL* v, SKL_REAL* flop);
  
  // cbc_pscalar.f90
  void cbc_ps_muscl_      (SKL_REAL* ws, int* sz, int* g, SKL_REAL* dh, int* c_scheme, SKL_REAL* v00, SKL_REAL* v, SKL_REAL* t, int* bv, 
                           int* bh1, int* bh2, int* swt, SKL_REAL* flop);
  void cbc_ps_buoyancy_   (SKL_REAL* v, int* sz, int* g, SKL_REAL* dt, SKL_REAL* gr, SKL_REAL* rei, SKL_REAL* t, int* bd, SKL_REAL* flop);
  void cbc_ps_diff_ee_    (SKL_REAL* t, int* sz, int* g, SKL_REAL* res, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* pei, SKL_REAL* qbc, 
                           int* bh, SKL_REAL* ws, SKL_REAL* flop);
  void cbc_hbc_drchlt_    (SKL_REAL* t,  int* sz, int* g, int* st, int* ed, int* bh, int* odr, SKL_REAL* tc);
  
  // cds_poisson.f90
  void cds_div_cc_        (SKL_REAL* div, int* sz, int* g, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* b2, SKL_REAL* vc, SKL_REAL* bnd, 
                          SKL_REAL* cut, SKL_REAL* v00, int* mode);
  void cds_div_cf_        (SKL_REAL* div, int* sz, int* g, SKL_REAL* dh, SKL_REAL* dt, SKL_REAL* b2, SKL_REAL* bnd, 
                          SKL_REAL* cut, SKL_REAL* v00, int* mode);
  void cds_psor_          (SKL_REAL* p,   int* sz, int* g, SKL_REAL* omg, SKL_REAL* res, SKL_REAL* div, SKL_REAL* bnd, SKL_REAL* cut, 
                          SKL_REAL* epsilon, int* para_key);
  
  // c3d_vof.f90
  void vof_uwd_    (SKL_REAL* f, int* sz, int* g, SKL_REAL* v00, SKL_REAL* dt, SKL_REAL* dh, SKL_REAL* v, SKL_REAL* q, int* bx, SKL_REAL* flop);
  void vof_muscl_  (SKL_REAL* f, int* sz, int* g, SKL_REAL* v00, SKL_REAL* dt, SKL_REAL* dh, SKL_REAL* v, SKL_REAL* q, int* bx, SKL_REAL* flop);
  
  // cds_vector.f90
  void cds_div_           (SKL_REAL* div, int* sz, int* g, SKL_REAL* coef, SKL_REAL* v, int* bv, float* cut, SKL_REAL* v00, SKL_REAL* flop);
  void cds_pvec_muscl_    (SKL_REAL* wv, int* sz, int* g, SKL_REAL* dh, int* c_scheme, SKL_REAL* v00, SKL_REAL* rei, SKL_REAL* v, 
                           int* bv, int* bp, int* v_mode, float* cut, SKL_REAL* flop);
  void cds_update_vec_    (SKL_REAL* v, SKL_REAL* div, int* sz, int* g, SKL_REAL* dt, SKL_REAL* dh, SKL_REAL* vc, SKL_REAL* p, 
                           int* bp, int* bv, float* cut, SKL_REAL* v00, SKL_REAL* coef, SKL_REAL* flop);
}

#endif // _SKL_FORTRAN_FUNC_CBC_H_
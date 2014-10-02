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
 * @file   ffv_Ffunc.h
 * @brief  FFV Fortran function Header
 * @author aics
 */

#include "cpm_Define.h"
#include "../FB/FB_Define.h"

#ifndef _FFV_F_FUNC_H_
#define _FFV_F_FUNC_H_

#ifdef _WIN32

// ffv_forcing.f90
#define hex_dir_            HEX_DIR
#define force_keep_vec_     FORCE_KEEP_VEC
#define hex_psrc_           HEX_PSRC
#define hex_force_pvec_     HEX_FORCE_PVEC
#define hex_force_vec_      HEX_FORCE_VEC


// ffv_pbc.f90
#define pobc_drchlt_        POBC_DRCHLT
#define pobc_neumann_       POBC_NEUMANN


// ffv_pscalar.f90
#define ps_muscl_           PS_MUSCL
#define ps_buoyancy_        PS_BUOYANCY
#define ps_diff_ee_         PS_DIFF_EE
#define hbc_drchlt_         HBC_DRCHLT

// ffv_vbc_inner.f90
#define pvec_vibc_specv_    PVEC_VIBC_SPECV
#define pvec_vibc_oflow_    PVEC_VIBC_OFLOW
#define div_ibc_drchlt_     DIV_IBC_DRCHLT
#define div_ibc_oflow_pvec_ DIV_IBC_OFLOW_PVEC
#define div_ibc_oflow_vec_  DIV_IBC_OFLOW_VEC

// ffv_vbc_outer.f90
#define vobc_pv_specv_      VOBC_PV_SPECV
#define vobc_pv_wall_       VOBC_PV_WALL
#define vobc_drchlt_        VOBC_DRCHLT
#define vobc_face_drchlt_   VOBC_FACE_DRCHLT
#define vobc_outflow_       VOBC_OUTFLOW
#define vobc_tfree1_        VOBC_TFREE1
#define vobc_tfree2_        VOBC_TFREE2
#define vobc_update_        VOBC_UPDATE
#define vobc_div_drchlt_    VOBC_DIV_DRCHLT
#define vobc_get_massflow_  VOBC_GET_MASSFLOW
#define vobc_neumann_       VOBC_NEUMANN

// ffv_velocity_binary.f90
#define ab2_                AB2
#define divergence_         DIVERGENCE
#define eddy_viscosity_     EDDY_VISCOSITY
#define euler_explicit_     EULER_EXPLICIT
#define friction_velocity_  FRICTION_VELOCITY
#define pvec_muscl_         PVEC_MUSCL
#define pvec_les_           PVEC_LES
#define update_vec_         UPDATE_VEC
#define vis_cn_jcb_         VIS_CN_JCB
#define vis_cn_mod_jcb_     VIS_CN_MOD_JCB
#define vis_cn_mod_sor_     VIS_CN_MOD_SOR
#define vis_cn_sor_         VIS_CN_SOR
#define vis_ee_             VIS_EE
#define vis_ee_vbc_         VIS_EE_VBC

// ffv_utility.f90
#define norm_v_div_dbg_     NORM_V_DIV_DBG
#define norm_v_div_max_     NORM_V_DIV_MAX
#define helicity_           HELICITY
#define i2vgt_              I2VGT
#define rot_v_              ROT_V
#define find_vmax_          FIND_VMAX
#define force_compo_        FORCE_COMPO


#define cds_pvec_vibc_specv_    CDS_PVEC_VIBC_SPECV


// c3d_vof.f90
#define vof_uwd_                VOF_UWD
#define vof_muscl_              VOF_MUSCL

// cds_vector.f90
#define pvec_muscl_cds_     PVEC_MUSCL_CDS
#define update_vec_cds_     UPDATE_VEC_CDS
#define divergence_cds_     DIVERGENCE_CDS
#define force_cds_          FORCE_CDS


// FB_util.f90
#define fb_average_s_            FB_AVERAGE_S
#define fb_average_v_            FB_AVERAGE_V
#define fb_delta_s_              FB_DELTA_S
#define fb_delta_v_              FB_DELTA_V
#define fb_interp_coarse_s_      FB_INTERP_COARSE_S
#define fb_interp_coarse_v_      FB_INTERP_COARSE_V
#define fb_limit_scalar_         FB_LIMIT_SCALAR
#define fb_minmax_s_             FB_MINMAX_S
#define fb_minmax_v_             FB_MINMAX_V
#define fb_minmax_vex_           FB_MINMAX_VEX
#define fb_set_vector_           FB_SET_VECTOR
#define fb_set_fvector_          FB_SET_FVECTOR
#define fb_vin_nijk_             FB_VIN_NIJK
#define fb_vout_nijk_            FB_VOUT_NIJK
#define fb_vout_ijkn_            FB_VOUT_IJKN
#define fb_totalp_               FB_TOTALP

#endif // _WIN32


extern "C" {
  //***********************************************************************************************
  // ffv_forcing.f90
  void hex_dir_(REAL_TYPE* v, int* sz, int* g, int* st, int* ed, int* bd, REAL_TYPE* vf, int* odr, REAL_TYPE* v00, REAL_TYPE* nv, double* flop);
  void force_keep_vec_(REAL_TYPE* wk, int* c_sz, int* st, int* ed, REAL_TYPE* v, int* sz, int* g);
  
  void hex_psrc_ (REAL_TYPE* src,
                  int* sz,
                  int* g,
                  int* st,
                  int* ed,
                  int* bd,
                  REAL_TYPE* vf,
                  REAL_TYPE* wk,
                  int* c_sz,
                  int* odr,
                  REAL_TYPE* v00,
                  REAL_TYPE* nv,
                  REAL_TYPE* c,
                  double* flop);
  
  void hex_force_pvec_ (REAL_TYPE* vc,
                        int* sz,
                        int* g,
                        int* st,
                        int* ed,
                        int* bd,
                        REAL_TYPE* vf,
                        REAL_TYPE* v,
                        int* odr,
                        REAL_TYPE* v00,
                        REAL_TYPE* dt,
                        REAL_TYPE* nv,
                        REAL_TYPE* c,
                        double* flop);
  
  void hex_force_vec_ (REAL_TYPE* v,
                       REAL_TYPE* dv,
                       int* sz,
                       int* g,
                       int* st,
                       int* ed,
                       int* bd,
                       REAL_TYPE* vf,
                       REAL_TYPE* wk,
                       int* c_sz,
                       int* odr,
                       REAL_TYPE* v00,
                       REAL_TYPE* dt,
                       REAL_TYPE* dh,
                       REAL_TYPE* nv,
                       REAL_TYPE* c,
                       REAL_TYPE* am,
                       double* flop);
  
  
  
  
  
  //***********************************************************************************************
  // ffv_pbc.f90
  void pobc_drchlt_(REAL_TYPE* p, int* sz, int* g, int* m_face, REAL_TYPE* pv, int* nID);
  
  void pobc_neumann_(REAL_TYPE* p, int* sz, int* g, int* m_face, int* nID);
  
  
  
  
  //***********************************************************************************************
  // ffv_pscalar.f90
  void ps_muscl_ (REAL_TYPE* ws,
                  int* sz,
                  int* g,
                  REAL_TYPE* dh,
                  int* scheme,
                  REAL_TYPE* v00,
                  REAL_TYPE* v,
                  REAL_TYPE* ie,
                  int* bp,
                  int* cdf,
                  int* bh,
                  int* swt,
                  double* flop);
  
  void ps_buoyancy_ (REAL_TYPE* v,
                     int* sz,
                     int* g,
                     REAL_TYPE* dgr,
                     REAL_TYPE* ie,
                     int* bd,
                     int* ncompo,
                     double* mtbl,
                     double* flop);
  
  void ps_diff_ee_ (REAL_TYPE* ie,
                    int* sz,
                    int* g,
                    double* res,
                    REAL_TYPE* dh,
                    REAL_TYPE* dt,
                    REAL_TYPE* qbc,
                    int* bh,
                    REAL_TYPE* ws,
                    int* ncompo,
                    double* mtbl,
                    int* h_mode,
                    double* flop);
  
  
  
  //***********************************************************************************************
  // ffv_vbc_inner.f90
  void pvec_vibc_oflow_   (REAL_TYPE* wv,
                           int* sz,
                           int* g,
                           int* st,
                           int* ed,
                           REAL_TYPE* dh,
                           REAL_TYPE* rei,
                           REAL_TYPE* v,
                           int* bv,
                           int* odr,
                           REAL_TYPE* vec,
                           double* flop);
  
  void pvec_vibc_specv_   (REAL_TYPE* wv,
                           int* sz,
                           int* g,
                           int* st,
                           int* ed,
                           REAL_TYPE* dh,
                           REAL_TYPE* v00,
                           REAL_TYPE* rei,
                           REAL_TYPE* v,
                           int* bv,
                           int* odr,
                           REAL_TYPE* vec,
                           double* flop);
  
  void div_ibc_drchlt_    (REAL_TYPE* div,
                           int* sz,
                           int* g,
                           int* st,
                           int* ed,
                           REAL_TYPE* v00,
                           int* bv,
                           int* odr,
                           REAL_TYPE* vec,
                           double* flop);
  
  void div_ibc_oflow_pvec_(REAL_TYPE* div,
                           int* sz,
                           int* g,
                           int* st,
                           int* ed,
                           REAL_TYPE* v00,
                           REAL_TYPE* v_cnv,
                           int* bv,
                           int* odr,
                           REAL_TYPE* v0,
                           REAL_TYPE* vf,
                           double* flop);
  
  void div_ibc_oflow_vec_ (REAL_TYPE* div,
                           int* sz,
                           int* g,
                           int* st,
                           int* ed,
                           int* bv,
                           int* odr,
                           REAL_TYPE* av,
                           double* flop);
  
  
  //***********************************************************************************************
  // ffv_vbc_outer.f90
  
  void vobc_pv_specv_ (REAL_TYPE* wv,
                       int* sz,
                       int* g,
                       int* m_face,
                       REAL_TYPE* dh,
                       REAL_TYPE* rei,
                       REAL_TYPE* v,
                       int* bv,
                       REAL_TYPE* vec,
                       int* nID,
                       double* flop);
  
  void vobc_pv_wall_ (REAL_TYPE* wv,
                      int* sz,
                      int* g,
                      int* m_face,
                      REAL_TYPE* dh,
                      REAL_TYPE* rei,
                      REAL_TYPE* v,
                      REAL_TYPE* vec,
                      int* nID,
                      double* flop);
  
  void vobc_drchlt_ (REAL_TYPE* v,
                     int* sz,
                     int* g,
                     int* m_face,
                     int* bv,
                     REAL_TYPE* vec,
                     int* nID);
  
  void vobc_neumann_ (REAL_TYPE* v,
                      int* sz,
                      int* g,
                      int* m_face,
                      REAL_TYPE* sum,
                      int* nID);

  void vobc_tfree2_ (REAL_TYPE* v,
                     int* sz,
                     int* g,
                     int* m_face,
                     REAL_TYPE* vf,
                     REAL_TYPE* sum,
                     int* nID);
  
  void vobc_tfree1_ (REAL_TYPE* vf,
                     int* sz,
                     int* g,
                     int* m_face,
                     int* nID);
  
  void vobc_update_ (REAL_TYPE* v,
                     int* sz,
                     int* g,
                     int* m_face,
                     REAL_TYPE* vc,
                     int* nID);
  
  void vobc_div_drchlt_ (REAL_TYPE* div,
                         int* sz,
                         int* g,
                         int* m_face,
                         int* bv,
                         REAL_TYPE* vec,
                         int* nID);
  
  void vobc_get_massflow_ (REAL_TYPE* sum,
                           int* sz,
                           int* g,
                           int* m_face,
                           REAL_TYPE* v,
                           int* bv,
                           int* nID);
  
  void vobc_face_drchlt_  (REAL_TYPE* vf,
                           int* sz,
                           int* g,
                           int* m_face,
                           int* bv,
                           REAL_TYPE* vec,
                           REAL_TYPE* vsum,
                           int* nID);

  
  //***********************************************************************************************
  // ffv_velocity_binary.f90
  void ab2_               (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* v, REAL_TYPE* ab, int* bd, REAL_TYPE* v00, double* flop);
  
  void divergence_        (REAL_TYPE* dv,
                           int* sz,
                           int* g,
                           REAL_TYPE* vc,
                           int* bv,
                           int* bp,
                           double* flop);
  
  void eddy_viscosity_    (REAL_TYPE* vt,
                           int* sz,
                           int* g,
                           REAL_TYPE* dh,
                           REAL_TYPE* re,
                           REAL_TYPE* cs,
                           REAL_TYPE* v,
                           int* bx,
                           REAL_TYPE* vt_range,
                           REAL_TYPE* yp_range,
                           REAL_TYPE* v00);
  
  void euler_explicit_    (REAL_TYPE* vc,
                           int* sz,
                           int* g,
                           REAL_TYPE* dt,
                           REAL_TYPE* v,
                           int* bd,
                           double* flop);
  
  void pvec_muscl_        (REAL_TYPE* wv,
                           int* sz,
                           int* g,
                           REAL_TYPE* dh,
                           int* c_scheme,
                           REAL_TYPE* v00,
                           REAL_TYPE* rei,
                           REAL_TYPE* v,
                           REAL_TYPE* vf,
                           int* bv,
                           int* bp,
                           REAL_TYPE* vcs_coef,
                           double* flop);
  
  void pvec_central_      (REAL_TYPE* wv,
                           int* sz,
                           int* g,
                           REAL_TYPE* dh,
                           int* c_scheme,
                           REAL_TYPE* v00,
                           REAL_TYPE* rei,
                           REAL_TYPE* v,
                           REAL_TYPE* vf,
                           int* bv,
                           int* bp,
                           REAL_TYPE* vcs_coef,
                           double* flop);
  
  void update_vec_        (REAL_TYPE* v,
                           REAL_TYPE* div,
                           int* sz,
                           int* g,
                           REAL_TYPE* dt,
                           REAL_TYPE* dh,
                           REAL_TYPE* vc,
                           REAL_TYPE* vf,
                           REAL_TYPE* p,
                           int* bp,
                           int* bv,
                           double* flop);
  
  void vis_cn_jcb_        (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                           REAL_TYPE* omg, REAL_TYPE* wv, int* bx, REAL_TYPE* wk, REAL_TYPE* coef, REAL_TYPE* res, double* flop);
  void vis_cn_mod_jcb_    (REAL_TYPE* vc, int* sz, int* g, int* st, int* ed, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                           REAL_TYPE* omg, REAL_TYPE* wv, int* bx, REAL_TYPE* wk, REAL_TYPE* coef, REAL_TYPE* res, REAL_TYPE* vec, double* flop);
  void vis_cn_mod_sor_    (REAL_TYPE* vc, int* sz, int* g, int* st, int* ed, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                           REAL_TYPE* omg, REAL_TYPE* wv, int* bx, REAL_TYPE* coef, REAL_TYPE* res, REAL_TYPE* vec, double* flop);
  void vis_cn_sor_        (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                           REAL_TYPE* omg, REAL_TYPE* wv, int* bx, REAL_TYPE* coef, REAL_TYPE* res, double* flop);
  void vis_ee_            (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, REAL_TYPE* wv, 
                           REAL_TYPE* v, int* bx, REAL_TYPE* coef, double* flop);
  void vis_ee_vbc_        (REAL_TYPE* vc, int* sz, int* g, int* st, int* ed, REAL_TYPE* dh, REAL_TYPE* dt, REAL_TYPE* v00, REAL_TYPE* rei, 
                           REAL_TYPE* v0, int* bx, int* odr, REAL_TYPE* coef, REAL_TYPE* vec, double* flop);
  
  


  
  //***********************************************************************************************
  // ffv_utility.f90
  void norm_v_div_dbg_ (double* ds,
                        int* index,
                        int* sz,
                        int* g,
                        REAL_TYPE* div,
                        REAL_TYPE* coef,
                        int* bp,
                        double* flop);
  
  void norm_v_div_max_ (double* ds,
                        int* sz,
                        int* g,
                        REAL_TYPE* div,
                        REAL_TYPE* coef,
                        int* bp,
                        double* flop);
  
  void helicity_          (REAL_TYPE* ht,    int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v, int* bv, REAL_TYPE* v00, double* flop);
  void i2vgt_             (REAL_TYPE* q,     int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v, int* bv, REAL_TYPE* v00, double* flop);
  void rot_v_             (REAL_TYPE* rot,   int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* v, int* bv, REAL_TYPE* v00, double* flop);
  void find_vmax_         (REAL_TYPE* v_max, int* sz, int* g, REAL_TYPE* v00, REAL_TYPE* v, double* flop);
  void face_avr_sampling_ (REAL_TYPE* p, int* sz, int* g, int* face, REAL_TYPE* avr);
  void shift_pressure_    (REAL_TYPE* p, int* sz, int* g, REAL_TYPE* avr);
  void force_compo_       (REAL_TYPE* frc, int* sz, int* g, int* tgt, REAL_TYPE* p, int* bid, REAL_TYPE* dh, int* st, int* ed, double* flop);

  
  
  // c3d_vof.f90
  void vof_uwd_    (REAL_TYPE* f, int* sz, int* g, REAL_TYPE* v00, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* v, REAL_TYPE* q, int* bx, double* flop);
  void vof_muscl_  (REAL_TYPE* f, int* sz, int* g, REAL_TYPE* v00, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* v, REAL_TYPE* q, int* bx, double* flop);
  
  
  
  //***********************************************************************************************
  // ffv_velocity_cds.f90
  void pvec_muscl_cds_ (REAL_TYPE* wv,
                        int* sz,
                        int* g,
                        REAL_TYPE* dh,
                        int* c_scheme,
                        REAL_TYPE* v00,
                        REAL_TYPE* rei,
                        REAL_TYPE* v,
                        REAL_TYPE* vf,
                        int* bv,
                        int* bp,
                        int* v_mode,
                        float* cut,
                        double* flop);
  
  void update_vec_cds_    (REAL_TYPE* v, REAL_TYPE* div, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* dh, REAL_TYPE* vc, REAL_TYPE* p,
                           int* bp, int* bv, float* cut, REAL_TYPE* v00, double* flop);
  void divergence_cds_    (REAL_TYPE* div, int* sz, int* g, REAL_TYPE* coef, REAL_TYPE* v, int* bv, float* cut, REAL_TYPE* v00, double* flop);
  void force_cds_         (REAL_TYPE* force, int* sz, int* g, REAL_TYPE* p, int* bp, int* bid, int* id, REAL_TYPE* dh, double* flop);
  void eddy_viscosity_cds_(REAL_TYPE* vt, int* sz, int* g, REAL_TYPE* dh, REAL_TYPE* re, REAL_TYPE* cs, REAL_TYPE* v, float* cut,
                           REAL_TYPE* vt_range, REAL_TYPE* yp_range, REAL_TYPE* v00, double* flop);
  
  
  
  //***********************************************************************************************
  // FB_util.f90
  
  void fb_delta_s_        (double* d,
                           int* sz,
                           int* g,
                           REAL_TYPE* sn,
                           REAL_TYPE* so,
                           int* bx,
                           double* flop);
  
  void fb_delta_v_        (double* d,
                           int* sz,
                           int* g,
                           REAL_TYPE* vn,
                           REAL_TYPE* vo,
                           int* bx,
                           double* flop);
  
  void fb_limit_scalar_   (REAL_TYPE* t,
                           int* sz,
                           int* g);
  
  void fb_minmax_s_       (REAL_TYPE* f_min,
                           REAL_TYPE* f_max,
                           int* sz,
                           int* g,
                           REAL_TYPE* s,
                           double* flop);
  
  void fb_minmax_v_       (REAL_TYPE* f_min,
                           REAL_TYPE* f_max,
                           int* sz,
                           int* g,
                           REAL_TYPE* v00,
                           REAL_TYPE* v,
                           double* flop);
  
  void fb_minmax_vex_     (REAL_TYPE* f_min,
                           REAL_TYPE* f_max,
                           int* sz,
                           int* g,
                           REAL_TYPE* v00,
                           REAL_TYPE* v,
                           double* flop);
  
  void fb_average_s_      (REAL_TYPE* avr,
                           int* sz,
                           int* g,
                           REAL_TYPE* s,
                           REAL_TYPE* nadd,
                           double* flop);
  
  void fb_average_v_      (REAL_TYPE* avr,
                           int* sz,
                           int* g,
                           REAL_TYPE* v,
                           REAL_TYPE* nadd,
                           double* flop);
  
  void fb_interp_coarse0_s_(REAL_TYPE* dst,
                            int* sz,
                            int* g,
                            REAL_TYPE* src,
                            int* st,
                            int* bk);
  
  void fb_interp_coarse0_v_(REAL_TYPE* dst,
                            int* sz,
                            int* g,
                            REAL_TYPE* src,
                            int* st,
                            int* bk);
  
  void fb_interp_coarse_s_(REAL_TYPE* dst,
                           int* sz,
                           int* g,
                           REAL_TYPE* src,
                           int* st,
                           int* bk);
  
  void fb_interp_coarse_v_(REAL_TYPE* dst,
                           int* sz,
                           int* g,
                           REAL_TYPE* src,
                           int* st,
                           int* bk);
  
  void fb_set_vector_     (REAL_TYPE* var,
                           int* sz,
                           int* g,
                           REAL_TYPE* val,
                           int* bv);
  
  void fb_set_fvector_    (REAL_TYPE* var,
                           int* sz,
                           int* g,
                           REAL_TYPE* val,
                           int* bv);
  
  void fb_vin_nijk_       (REAL_TYPE* vo,
                           int* sz,
                           int* g,
                           REAL_TYPE* vi,
                           REAL_TYPE* v00,
                           REAL_TYPE* refv,
                           double* flop);
  
  void fb_vout_nijk_      (REAL_TYPE* vout,
                           REAL_TYPE* vin,
                           int* sz,
                           int* g,
                           REAL_TYPE* v00,
                           REAL_TYPE* unit_v,
                           double* flop);
  
  
  void fb_vout_ijkn_      (REAL_TYPE* vout,
                           REAL_TYPE* vin,
                           int* sz,
                           int* g,
                           REAL_TYPE* v00,
                           REAL_TYPE* unit_v,
                           double* flop);
  
  void fb_totalp_         (REAL_TYPE* tp,
                           int* sz,
                           int* g,
                           REAL_TYPE* v,
                           REAL_TYPE* p,
                           REAL_TYPE* v00,
                           double* flop);
  
}

#endif // _FFV_F_FUNC_H_

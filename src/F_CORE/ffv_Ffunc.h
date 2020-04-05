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
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   ffv_Ffunc.h
 * @brief  FFV Fortran function Header
 * @author aics
 */

#include "cpm_Define.h"
#include "FB_Define.h"

#ifndef _FFV_F_FUNC_H_
#define _FFV_F_FUNC_H_


extern "C" {
  //***********************************************************************************************
  // io_sph.f90
 void sph_write_tensor_(
                  int* sz,
                  int* gc,
                  REAL_TYPE* org,
                  REAL_TYPE* pch,
                  int* stp,
                  REAL_TYPE* tm,
                  REAL_TYPE* R,
                  int* compo,
                  char* fname,
                  int* dtype);


  //***********************************************************************************************
  // ffv_forcing.f90
  void hex_dir_(REAL_TYPE* v, int* sz, int* g, int* st, int* ed, int* bd, REAL_TYPE* vf, int* odr, REAL_TYPE* nv, double* flop);
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
                  REAL_TYPE* vf,
                  REAL_TYPE* ie,
                  int* bid,
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

  void ps_convection_ee_ (REAL_TYPE* ie,
                          int* sz,
                          int* g,
                          REAL_TYPE* dt,
                          int* bd,
                          REAL_TYPE* ie0,
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

  void perturb_u_y_ (REAL_TYPE* v,
                     int* sz,
                     int* g,
                     REAL_TYPE* dh,
                     REAL_TYPE* origin,
                     REAL_TYPE* width,
                     REAL_TYPE* Re_tau,
                     REAL_TYPE* Ubar,
                     REAL_TYPE* nu,
                     int* mode);

  void perturb_u_z_ (REAL_TYPE* v,
                     int* sz,
                     int* g,
                     REAL_TYPE* dh,
                     REAL_TYPE* origin,
                     REAL_TYPE* width,
                     REAL_TYPE* Re_tau,
                     REAL_TYPE* Ubar,
                     REAL_TYPE* nu,
                     int* mode);

  void pvec_ibc_oflow_ (REAL_TYPE* wv,
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

  void pvec_ibc_specv_fvm_ (REAL_TYPE* wv,
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

  void pvec_ibc_specv_fdm_ (REAL_TYPE* wv,
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

  void div_ibc_drchlt_    (REAL_TYPE* div,
                           int* sz,
                           int* g,
                           REAL_TYPE* dh,
                           int* st,
                           int* ed,
                           int* bv,
                           int* odr,
                           REAL_TYPE* vec,
                           double* flop);

  void div_ibc_oflow_pvec_(REAL_TYPE* div,
                           int* sz,
                           int* g,
                           int* st,
                           int* ed,
                           REAL_TYPE* v_cnv,
                           REAL_TYPE* dt,
                           REAL_TYPE* dh,
                           int* bv,
                           int* odr,
                           REAL_TYPE* v0,
                           REAL_TYPE* vf,
                           double* flop);

  void div_ibc_oflow_vec_ (REAL_TYPE* div,
                           int* sz,
                           int* g,
                           REAL_TYPE* dh,
                           int* st,
                           int* ed,
                           int* bv,
                           int* odr,
                           REAL_TYPE* av,
                           double* flop);

  void pvec_ibc_sldrev_fvm_ (REAL_TYPE* wv,
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
                             REAL_TYPE* ctr,
                             REAL_TYPE* org,
                             double* flop);

  void div_ibc_sldrev_ (REAL_TYPE* div,
                        int* sz,
                        int* g,
                        int* st,
                        int* ed,
                        REAL_TYPE* pch,
                        int* bv,
                        int* odr,
                        REAL_TYPE* vec,
                        REAL_TYPE* ctr,
                        REAL_TYPE* org,
                        double* flop);


  //***********************************************************************************************
  // ffv_vbc_outer.f90

  void vobc_cc_drchlt_ (REAL_TYPE* v,
                        int* sz,
                        int* g,
                        int* m_face,
                        REAL_TYPE* vec,
                        int* nID);

  void vobc_cc_neumann_ (REAL_TYPE* v,
                         int* sz,
                         int* g,
                         int* m_face,
                         int* nID);

  void vobc_cc_outflow_ (REAL_TYPE* vc,
                         REAL_TYPE* v0,
                         int* sz,
                         int* g,
                         REAL_TYPE* dh,
                         REAL_TYPE* dt,
                         int* bv,
                         REAL_TYPE* v_cnv,
                         int* m_face,
                         int* nID);

  void vobc_cc_copy_ (REAL_TYPE* v,
                      int* sz,
                      int* g,
                      int* m_face,
                      REAL_TYPE* vc,
                      int* nID);

  void vobc_cc_tfree_ (REAL_TYPE* v,
                       int* sz,
                       int* g,
                       int* m_face,
                       REAL_TYPE* vf,
                       int* nID);

  void vobc_cf_tfree_ (REAL_TYPE* vf,
                       int* sz,
                       int* g,
                       int* m_face,
                       int* nID);



  void vobc_div_drchlt_ (REAL_TYPE* div,
                         int* sz,
                         int* g,
                         REAL_TYPE* dh,
                         int* m_face,
                         int* bv,
                         REAL_TYPE* vec,
                         int* nID);



  //***********************************************************************************************
  // ffv_vbc_outer_face.f90

  void vobc_face_drchlt_  (REAL_TYPE* vf,
                           int* sz,
                           int* g,
                           int* m_face,
                           int* bv,
                           REAL_TYPE* vec,
                           int* nID);

  void vobc_face_massflow_ (REAL_TYPE* sum,
                            int* sz,
                            int* g,
                            int* m_face,
                            REAL_TYPE* vf,
                            int* bv,
                            int* nID);


  //***********************************************************************************************
  // ffv_vbc_outer_flux.f90

  void vobc_pv_specv_fdm_ (REAL_TYPE* wv,
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

  void vobc_pv_specv_fvm_ (REAL_TYPE* wv,
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


  //***********************************************************************************************
  // ffv_velocity_binary.f90
  void ab2_ (REAL_TYPE* vc, int* sz, int* g, REAL_TYPE* dt, REAL_TYPE* v, REAL_TYPE* ab, int* bd, double* flop);

  void divergence_cc_ (REAL_TYPE* dv,
                       int* sz,
                       int* g,
                       REAL_TYPE* dh,
                       REAL_TYPE* vc,
                       int* bid,
                       int* bcd,
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
                           REAL_TYPE* yp_range);

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
                           REAL_TYPE* rei,
                           REAL_TYPE* v,
                           REAL_TYPE* vf,
                           int* bid,
                           int* bcd,
                           REAL_TYPE* vcs_coef,
                           double* flop);

  void pvec_central_      (REAL_TYPE* wv,
                           int* sz,
                           int* g,
                           REAL_TYPE* dh,
                           int* c_scheme,
                           REAL_TYPE* rei,
                           REAL_TYPE* v,
                           REAL_TYPE* vf,
                           int* bid,
                           int* bcd,
                           REAL_TYPE* vcs_coef,
                           double* flop);

  void update_vec_ (REAL_TYPE* v,
                    REAL_TYPE* vf,
                    REAL_TYPE* div,
                    int* sz,
                    int* g,
                    REAL_TYPE* dt,
                    REAL_TYPE* dh,
                    REAL_TYPE* vc,
                    REAL_TYPE* p,
                    int* bp,
                    int* bcd,
                    double* flop);

  void update_vec4_ (REAL_TYPE* v,
                     REAL_TYPE* vf,
                     REAL_TYPE* div,
                     int* sz,
                     int* g,
                     REAL_TYPE* dt,
                     REAL_TYPE* dh,
                     REAL_TYPE* vc,
                     REAL_TYPE* p,
                     int* bp,
                     int* bid,
                     int* bcd,
                     double* flop,
                     int* c_scheme);

  void pvec_muscl_les_ (REAL_TYPE* wv,
                        int* sz,
                        int* g,
                        REAL_TYPE* dh,
                        int* c_scheme,
                        REAL_TYPE* rei,
                        REAL_TYPE* v,
                        REAL_TYPE* vf,
                        int* bid,
                        int* bcd,
                        REAL_TYPE* vcs_coef,
                        REAL_TYPE* Cs,
                        int* imodel,
                        REAL_TYPE* nu,
                        REAL_TYPE* rho,
                        double* flop);

  void pvec_central_les_ (REAL_TYPE* wv,
                          int* sz,
                          int* g,
                          REAL_TYPE* dh,
                          int* c_scheme,
                          REAL_TYPE* rei,
                          REAL_TYPE* v,
                          REAL_TYPE* vf,
                          int* bid,
                          int* bcd,
                          REAL_TYPE* vcs_coef,
                          REAL_TYPE* Cs,
                          int* imodel,
                          REAL_TYPE* nu,
                          REAL_TYPE* rho,
                          double* flop);

  void stabilize_ (REAL_TYPE* vc,
                   int* sz,
                   int* g,
                   REAL_TYPE* dt,
                   REAL_TYPE* v,
                   int* bcd,
                   REAL_TYPE* vref,
                   REAL_TYPE* st,
                   REAL_TYPE* ed,
                   REAL_TYPE* penalty,
                   int* count,
                   double* flop);


  //***********************************************************************************************
  // ffv_utility.f90
  void norm_v_div_l2_ (REAL_TYPE* ds,
                        int* sz,
                        int* g,
                        REAL_TYPE* div,
                        int* bp,
                        double* flop);

  void norm_v_div_l2_lc_ (REAL_TYPE* ds,
                          int* sz,
                          int* g,
                          REAL_TYPE* div,
                          int* bp,
                          REAL_TYPE* p,
                          REAL_TYPE* p0,
                          REAL_TYPE* cm,
                          double* flop);

  void norm_v_div_max_ (REAL_TYPE* ds,
                        int* sz,
                        int* g,
                        REAL_TYPE* div,
                        int* bp,
                        double* flop);

  void norm_v_div_max_lc_ (REAL_TYPE* ds,
                           int* sz,
                           int* g,
                           REAL_TYPE* div,
                           int* bp,
                           REAL_TYPE* p,
                           REAL_TYPE* p0,
                           REAL_TYPE* cm,
                           double* flop);

  void helicity_ (REAL_TYPE* ht,
                  int* sz, int* g,
                  REAL_TYPE* dh,
                  REAL_TYPE* v,
                  int* bv,
                  double* flop);

  void i2vgt_ (REAL_TYPE* q,
               int* sz,
               int* g,
               REAL_TYPE* dh,
               REAL_TYPE* v,
               int* bv,
               double* flop);

  void rot_v_ (REAL_TYPE* rot,
               int* sz,
               int* g,
               REAL_TYPE* dh,
               REAL_TYPE* v,
               int* bv,
               double* flop);

  void find_vmax_         (REAL_TYPE* v_max, int* sz, int* g, REAL_TYPE* v, double* flop);
  void face_avr_sampling_ (REAL_TYPE* p, int* sz, int* g, int* face, REAL_TYPE* avr);
  void shift_pressure_    (REAL_TYPE* p, int* sz, int* g, REAL_TYPE* avr);
  void force_compo_       (REAL_TYPE* frc, int* sz, int* g, int* tgt, REAL_TYPE* p, int* bid, REAL_TYPE* dh, int* st, int* ed, double* flop);

  void calc_rms_v_ (REAL_TYPE* rms,
                    REAL_TYPE* rmsmean,
                    int* sz,
                    int* g,
                    REAL_TYPE* v,
                    REAL_TYPE* av,
                    REAL_TYPE* accum,
                    double* flop);

  void calc_rms_s_ (REAL_TYPE* rms,
                    REAL_TYPE* rmsmean,
                    int* sz,
                    int* g,
                    REAL_TYPE* s,
                    REAL_TYPE* as,
                    REAL_TYPE* accum,
                    double* flop);

  void output_vtk_(int* step,
                   REAL_TYPE* G_origin,
                   int* G_division,
                   int* G_size,
                   int* myRank,
                   int* sz,
                   REAL_TYPE* dx,
                   REAL_TYPE* dy,
                   REAL_TYPE* dz,
                   int* g,
                   REAL_TYPE* v,
                   REAL_TYPE* p
                   );



  void calc_reynolds_stress_ (REAL_TYPE* R,
                              REAL_TYPE* R_ave,
                              int* sz,
                              int* g,
                              REAL_TYPE* v,
                              REAL_TYPE* v_ave,
                              REAL_TYPE* nadd,
                              double* flop);

  void calc_production_rate_ (REAL_TYPE* P_ave,
                              int* sz,
                              REAL_TYPE* dh,
                              int* g,
                              REAL_TYPE* v_ave,
                              REAL_TYPE* R,
                              int* bid,
                              REAL_TYPE* nadd,
                              double* flop);

  void calc_dissipation_rate_ (REAL_TYPE* E_ave,
                               int* sz,
                               REAL_TYPE* dh,
                               int* g,
                               REAL_TYPE* nu,
                               REAL_TYPE* v,
                               REAL_TYPE* v_ave,
                               int* bid,
                               REAL_TYPE* nadd,
                               double* flop);

  void calc_turb_transport_rate_ (REAL_TYPE* T_ave,
                                  int* sz,
                                  REAL_TYPE* dh,
                                  int* g,
                                  REAL_TYPE* v,
                                  REAL_TYPE* v_ave,
                                  REAL_TYPE* R,
                                  int* bid,
                                  REAL_TYPE* nadd,
                                  double* flop);

  void calc_vel_pregrad_term_ (REAL_TYPE* PI_ave,
                               int* sz,
                               REAL_TYPE* dh,
                               int* g,
                               REAL_TYPE* v,
                               REAL_TYPE* v_ave,
                               REAL_TYPE* p,
                               REAL_TYPE* p_ave,
                               int* bp,
                               REAL_TYPE* nadd,
                               double* flop);

  void averaging_xz_plane_ (REAL_TYPE* vmean,
                            REAL_TYPE* rmsmean,
                            REAL_TYPE* Rmean,
                            REAL_TYPE* Pmean,
                            REAL_TYPE* Emean,
                            REAL_TYPE* Tmean,
                            REAL_TYPE* PImean,
                            int*       sz,
                            int*       g,
                            REAL_TYPE* v_ave,
                            REAL_TYPE* rms_ave,
                            REAL_TYPE* R_ave,
                            REAL_TYPE* P_ave,
                            REAL_TYPE* E_ave,
                            REAL_TYPE* T_ave,
                            REAL_TYPE* PI_ave,
                            double* flop
                            );

  void src_trnc_ (REAL_TYPE* rhs,
                  int* sz,
                  int* g,
                  REAL_TYPE* dt,
                  REAL_TYPE* dh,
                  REAL_TYPE* u_sum,
                  double* flop);

  void src_1st_ (REAL_TYPE* rhs,
                 int* sz,
                 int* g,
                 REAL_TYPE* dt,
                 REAL_TYPE* dh,
                 REAL_TYPE* u_sum,
                 REAL_TYPE* b,
                 double* flop);

  void src_2nd_ (REAL_TYPE* rhs,
                 int* sz,
                 int* g,
                 REAL_TYPE* dh,
                 REAL_TYPE* b,
                 REAL_TYPE* psi,
                 double* flop);

  void generate_iblank_ (int* iblk, int* sz, int* g, int* bcd, int* g_out);

  void pack_scalar_   (REAL_TYPE* dst, int* sz, int* g, REAL_TYPE* src, int* g_out);
  void pack_vector_   (REAL_TYPE* dst, int* sz, int* g, REAL_TYPE* src, int* g_out);
  void unpack_scalar_ (REAL_TYPE* dst, int* sz, int* g, REAL_TYPE* src, int* g_out);
  void unpack_vector_ (REAL_TYPE* dst, int* sz, int* g, REAL_TYPE* src, int* g_out);

  void write_plot3d_ (REAL_TYPE* buf, int* sz, int* g, int* nx, char* fname);
  void read_plot3d_  (REAL_TYPE* buf, int* sz, int* g, int* nx, char* fname);


  //***********************************************************************************************
  // c3d_vof.f90
  void vof_uwd_    (REAL_TYPE* f,
                    int* sz,
                    int* g,
                    REAL_TYPE* dt,
                    REAL_TYPE* dh,
                    REAL_TYPE* v,
                    REAL_TYPE* q,
                    int* bx,
                    double* flop);

  void vof_muscl_  (REAL_TYPE* f,
                    int* sz,
                    int* g,
                    REAL_TYPE* dt,
                    REAL_TYPE* dh,
                    REAL_TYPE* v,
                    REAL_TYPE* q,
                    int* bx,
                    double* flop);


  //***********************************************************************************************
  // ffv_velocity_cds.f90
  void pvec_muscl_cds_ (REAL_TYPE* wv,
                        int* sz,
                        int* g,
                        REAL_TYPE* dh,
                        int* c_scheme,
                        REAL_TYPE* rei,
                        REAL_TYPE* v,
                        REAL_TYPE* vf,
                        int* bv,
                        int* bp,
                        int* bcd,
                        int* v_mode,
                        int* cut_l,
                        int* cut_u,
                        double* flop);

  void update_vec_cds_ (REAL_TYPE* v,
                        REAL_TYPE* div,
                        int* sz,
                        int* g,
                        REAL_TYPE* dt,
                        REAL_TYPE* dh,
                        REAL_TYPE* vc,
                        REAL_TYPE* p,
                        int* bp,
                        int* bv,
                        int* bcd,
                        int* cut_l,
                        int* cut_u,
                        double* flop);

  void divergence_cds_ (REAL_TYPE* div,
                        int* sz,
                        int* g,
                        REAL_TYPE* coef,
                        REAL_TYPE* v,
                        int* bv,
                        int* bcd,
                        int* cut_l,
                        int* cut_u,
                        double* flop);

  void force_cds_ (REAL_TYPE* force,
                   int* sz,
                   int* g,
                   REAL_TYPE* p,
                   int* bp,
                   int* bid,
                   int* id,
                   REAL_TYPE* dh,
                   double* flop);

  void eddy_viscosity_cds_(REAL_TYPE* vt,
                           int* sz,
                           int* g,
                           REAL_TYPE* dh,
                           REAL_TYPE* re,
                           REAL_TYPE* cs,
                           REAL_TYPE* v,
                           int* cut_l,
                           int* cut_u,
                           REAL_TYPE* vt_range,
                           REAL_TYPE* yp_range,
                           double* flop);

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
                           REAL_TYPE* v,
                           double* flop);

  void fb_minmax_vex_     (REAL_TYPE* f_min,
                           REAL_TYPE* f_max,
                           int* sz,
                           int* g,
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

  void fb_vin_ijkn_       (REAL_TYPE* vo,
                           int* sz,
                           int* g,
                           REAL_TYPE* refv,
                           double* flop);

  void fb_vin_nijk_       (REAL_TYPE* vo,
                           int* sz,
                           int* g,
                           REAL_TYPE* vi,
                           REAL_TYPE* refv,
                           double* flop);

  void fb_vout_nijk_      (REAL_TYPE* vout,
                           REAL_TYPE* vin,
                           int* sz,
                           int* g,
                           REAL_TYPE* unit_v,
                           double* flop);


  void fb_vout_ijkn_      (REAL_TYPE* vout,
                           REAL_TYPE* vin,
                           int* sz,
                           int* g,
                           REAL_TYPE* unit_v,
                           double* flop);

  void fb_totalp_         (REAL_TYPE* tp,
                           int* sz,
                           int* g,
                           REAL_TYPE* v,
                           REAL_TYPE* p,
                           double* flop);

}

#endif // _FFV_F_FUNC_H_

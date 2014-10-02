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
 * @file   ffv_LSfunc.h
 * @brief  FFV Linear Solver function Header
 * @author aics
 */

#include "cpm_Define.h"
#include "../FB/FB_Define.h"

#ifndef _FFV_LS_FUNC_H_
#define _FFV_LS_FUNC_H_

#ifdef _WIN32


// ffv_poisson.f90
#define div_cnst_           DIV_CNST
#define psor_               PSOR
#define psor2sma_core_      PSOR2SMA_CORE
#define psor2sma_naive_     PSOR2SMA_NAIVE
#define sma_comm_           SMA_COMM
#define sma_comm_wait_      SMA_COMM_WAIT
#define cds_psor_           CDS_PSOR
#define poi_residual_       POI_RESIDUAL
#define poi_rhs_            POI_RHS

// ffv_blas.f90
#define blas_clear_          BLAS_CLEAR
#define blas_copy_           BLAS_COPY
#define blas_xpay_           BLAS_XPAY
#define blas_axpy_           BLAS_AXPY
#define blas_triad_          BLAS_TRIAD
#define blas_axpbypz_        BLAS_AXPBYPZ
#define blas_dot1_           BLAS_DOT1
#define blas_dot2_           BLAS_DOT2
#define blas_calcr_          BLAS_CALCR
#define blas_calcr2_         BLAS_CALCR2
#define blas_calcax_         BLAS_CALCAX
#define blas_smoother_core_  BLAS_SMOOTHER_CORE
#define blas_preconditioner_ BLAS_PRECONDITIONER
#define blas_calcr_naive_    BLAS_CALCR_NAIVE
#define blas_calcax_naive_   BLAS_CALCAX_NAIVE

// ffv_cg.f90
#define bicg_update_p_  BICG_UPDATE_P
#define bicg_update_x_  BICG_UPDATE_X

#endif // _WIN32


extern "C" {
  
  //***********************************************************************************************
  // ffv_poisson.f90
  void poi_residual_ (double* res,
                      int* sz,
                      int* g,
                      REAL_TYPE* p,
                      REAL_TYPE* b,
                      int* bp,
                      double* flop);
  
  void poi_rhs_ (double* rhs,
                 REAL_TYPE* b,
                 int* sz,
                 int* g,
                 REAL_TYPE* s_0,
                 REAL_TYPE* s_1,
                 int* bp,
                 REAL_TYPE* dh,
                 REAL_TYPE* dt,
                 double* flop);
  
  void psor_ (REAL_TYPE* p,
              int* sz,
              int* g,
              REAL_TYPE* omg,
              double* cnv,
              REAL_TYPE* b,
              int* bp,
              double* flop);
  
  void psor2sma_core_ (REAL_TYPE* p,
                       int* sz,
                       int* g,
                       int* ip,
                       int* color,
                       REAL_TYPE* omg,
                       double* cnv,
                       REAL_TYPE* b,
                       int* bp,
                       double* flop);
  
  void psor2sma_naive_ (REAL_TYPE* p,
                        int* sz,
                        int* g,
                        int* ip,
                        int* color,
                        REAL_TYPE* omg,
                        double* cnv,
                        REAL_TYPE* b,
                        int* bp,
                        REAL_TYPE* pn,
                        double* flop);
  
  void sma_comm_      (REAL_TYPE* p,
                       int* sz,
                       int* g,
                       int* col,
                       int* ip,
                       int* cf_sz,
                       REAL_TYPE* cf_x,
                       REAL_TYPE* cf_y,
                       REAL_TYPE* cf_z,
                       int* key,
                       int* nID);
  
  void sma_comm_wait_ (REAL_TYPE* p,
                       int* sz,
                       int* g,
                       int* col,
                       int* ip,
                       int* cf_sz,
                       REAL_TYPE* cf_x,
                       REAL_TYPE* cf_y,
                       REAL_TYPE* cf_z,
                       int* key);
  
  //***********************************************************************************************
  // ffv_blas.f90
  void blas_clear_    (REAL_TYPE* x,
                       int* sz,
                       int* g);
  
  void blas_copy_     (REAL_TYPE* y,
                       REAL_TYPE* x,
                       int* sz,
                       int* g);
  
  void blas_xpay_     (REAL_TYPE* y,
                       REAL_TYPE* x,
                       REAL_TYPE* a,
                       int* sz,
                       int* g);
  
  void blas_axpy_     (REAL_TYPE* y,
                       REAL_TYPE* x,
                       REAL_TYPE* a,
                       int* sz,
                       int* g);
  
  void blas_triad_    (REAL_TYPE* z,
                       REAL_TYPE* x,
                       REAL_TYPE* y,
                       double* a,
                       int* sz,
                       int* g,
                       double* flop);
  
  void blas_axpbypz_  (REAL_TYPE* z,
                       REAL_TYPE* x,
                       REAL_TYPE* y,
                       REAL_TYPE* a,
                       REAL_TYPE* b,
                       int* sz,
                       int* g,
                       double* flop);
  
  void blas_dot1_     (double* r,
                       REAL_TYPE* p,
                       int* sz,
                       int* g,
                       double* flop);
  
  void blas_dot2_     (double* r,
                       REAL_TYPE* p,
                       REAL_TYPE* q,
                       int* sz,
                       int* g,
                       double* flop);
  
  void blas_calcr_    (REAL_TYPE* r,
                       REAL_TYPE* p,
                       REAL_TYPE* b,
                       int* bp,
                       int* sz,
                       int* g,
                       double* flop);
  
  void blas_calcr_naive_ (REAL_TYPE* r,
                          REAL_TYPE* p,
                          REAL_TYPE* b,
                          int* bp,
                          int* sz,
                          int*g,
                          REAL_TYPE* pn,
                          double* flop);
  
  void blas_calcr2_   (REAL_TYPE* rr,
                       REAL_TYPE* p,
                       REAL_TYPE* b,
                       int* bp,
                       int* sz,
                       int* g);
  
  void blas_calcax_   (REAL_TYPE* ap,
                       REAL_TYPE* p,
                       int* bp,
                       int* sz,
                       int* g,
                       double* flop);
  
  void blas_calcax_naive_ (REAL_TYPE* ap,
                           REAL_TYPE* p,
                           int* bp,
                           int* sz,
                           int* g,
                           REAL_TYPE* pn,
                           double* flop);
  
  void blas_smoother_core_  (REAL_TYPE* x,
                             REAL_TYPE* b,
                             int* bp,
                             int* ip,
                             int* color,
                             REAL_TYPE* omg,
                             int* sz,
                             int* g);
  
  void blas_preconditioner_  (REAL_TYPE* y,
                              REAL_TYPE* x,
                              int* bp,
                              REAL_TYPE* omg,
                              int* sz,
                              int* g);

  
  //***********************************************************************************************
  // ffv_cg.f90
  
  void bicg_update_p_ (REAL_TYPE* p,
                       REAL_TYPE* r,
                       REAL_TYPE* q,
                       double* beta,
                       double* omega,
                       int* sz,
                       int* g,
                       double* flop);
  
  void bicg_update_x_ (REAL_TYPE* x,
                       REAL_TYPE* p,
                       REAL_TYPE* s,
                       int* bp,
                       double* alpha,
                       double* omega,
                       double* x_l2,
                       double* err,
                       int* sz,
                       int* g,
                       double* flop);

  //***********************************************************************************************
  // ffv_poisson2.f90
  void matvec_p_   (REAL_TYPE* ax,  int* sz, int* g, REAL_TYPE* p, int* bp, double* flop);
  void residual_   (REAL_TYPE* rs,  int* sz, int* g, double* res_a, REAL_TYPE* div, REAL_TYPE* src, int* bp, double* flop);
  void orth_basis_ (REAL_TYPE* dst, int* sz, int* g, int* nc, int* l, double* s, REAL_TYPE* src, double* flop);
  void cp_orth_basis_ (REAL_TYPE* z, int* sz, int* g, REAL_TYPE* v);
  void copy_1_   (REAL_TYPE* dst, int* sz, int* g, int* nc, REAL_TYPE* src, int* im);
  void copy_2_   (REAL_TYPE* dst, int* sz, int* g, int* nc, REAL_TYPE* src, int* im);
  void ml_add_1_ (double* ac,     int* sz, int* g, int* nc, REAL_TYPE* s4, REAL_TYPE* s3, int* lm, double* flop);
  void ml_add_2_ (double* ac,     int* sz, int* g, REAL_TYPE* s3, double* flop);
  void ml_add_3_ (REAL_TYPE* s3,  int* sz, int* g, int* nc, double* s, REAL_TYPE* s4, int* lm, double* flop);
  void ml_add_4_ (REAL_TYPE* s3,  int* sz, int* g, int* nc, double* s, REAL_TYPE* s4, int* lm, double* flop);
  
   // cds_poisson.f90
  void cds_psor_          (REAL_TYPE* p,   int* sz, int* g, REAL_TYPE* omg, REAL_TYPE* res, REAL_TYPE* div, REAL_TYPE* bnd, REAL_TYPE* cut,
                           REAL_TYPE* epsilon, int* para_key);
  
}

#endif // _FFV_LS_FUNC_H_

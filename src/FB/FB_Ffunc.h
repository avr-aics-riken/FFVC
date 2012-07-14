#ifndef _FB_F_FUNC_H_
#define _FB_F_FUNC_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/**
 @file   FB_Ffunc.h
 @brief  FlowBase Fortran function Header
 @author kero
 */


#include "cpm_Define.h"
#include "FB_Define.h"

#ifdef _WIN32
// FB_util.f90
#define fb_average_s_         FB_AVERAGE_S
#define fb_average_v_         FB_AVERAGE_V
#define fb_copy_real_s_       FB_COPY_REAL_S
#define fb_copy_real_v_       FB_COPY_REAL_V
#define fb_delta_s_           FB_DELTA_S
#define fb_delta_v_           FB_DELTA_V
#define fb_interp_coarse_s_   FB_INTERP_COARSE_S
#define fb_interp_coarse_v_   FB_INTERP_COARSE_V
#define fb_limit_scalar_      FB_LIMIT_SCALAR
#define fb_minmax_s_          FB_MINMAX_S
#define fb_minmax_v_          FB_MINMAX_V
#define fb_set_int_s_         FB_SET_INT_S
#define fb_set_real_s_        FB_SET_REAL_S
#define fb_set_real_v_        FB_SET_REAL_V
#define fb_set_vector_        FB_SET_VECTOR
#define fb_shift_refv_in_     FB_SHIFT_REFV_IN
#define fb_shift_refv_out_    FB_SHIFT_REFV_OUT
#define fb_totalp_            FB_TOTALP

#define fb_copy_real_         FB_COPY_REAL
#define fb_copy_int_          FB_COPY_INT
#define fb_set_int_           FB_SET_INT
#define fb_set_real_          FB_SET_REAL
#define fb_prs_d2nd_          FB_PRS_D2ND
#define fb_prs_nd2d_          FB_PRS_ND2D
#define fb_tmp_d2nd_          FB_TMP_D2ND
#define fb_tmp_nd2d_          FB_TMP_ND2D
#define fb_xcopy_             FB_XCOPY

#endif // _WIN32

extern "C" {
	// FB_util.f90
  void fb_copy_real_s_    (REAL_TYPE* dst,
                           REAL_TYPE* src,
                           int* sz,
                           int* g);
  
  void fb_copy_real_v_    (REAL_TYPE* dst,
                           REAL_TYPE* src,
                           int* sz,
                           int* g);
  
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
  
  void fb_set_int_s_      (int* var,
                           int* sz,
                           int* g,
                           int* init);
  
  void fb_set_real_s_     (REAL_TYPE* var,
                           int* sz,
                           int* g,
                           REAL_TYPE* init);
  
  void fb_set_real_v_     (REAL_TYPE* var,
                           int* sz,
                           int* g,
                           REAL_TYPE* vec);
  
  void fb_copy_real_      (REAL_TYPE* dst,
                           REAL_TYPE* src,
                           int* sz);
  
  void fb_copy_int_       (int* dst,
                           int* src,
                           int* sz);
  
  void fb_set_int_        (int* var,
                           int* sz,
                           int* init);
  
  void fb_set_real_       (REAL_TYPE* var,
                           int* sz,
                           REAL_TYPE* init);
  
  void fb_xcopy_          (REAL_TYPE* dst,
                           REAL_TYPE* src,
                           int* sz,
                           REAL_TYPE* scale,
                           double* flop);
 
  void fb_average_s_      (REAL_TYPE* avr, 
                           int* sz, 
                           int* g, 
                           REAL_TYPE* s, 
                           double* flop);
  
  void fb_average_v_      (REAL_TYPE* avr, 
                           int* sz, 
                           int* g, 
                           REAL_TYPE* v, 
                           double* flop);
  
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
  
  void fb_prs_d2nd_       (REAL_TYPE* s, 
                           int* sz, 
                           REAL_TYPE* Base_prs, 
                           REAL_TYPE* Ref_rho, 
                           REAL_TYPE* Ref_v, 
                           REAL_TYPE* scale, 
                           double* flop);
  
  void fb_prs_nd2d_       (REAL_TYPE* dst, 
                           REAL_TYPE* src, 
                           int* sz, 
                           REAL_TYPE* Base_prs, 
                           REAL_TYPE* Ref_rho, 
                           REAL_TYPE* Ref_v, 
                           REAL_TYPE* scale, 
                           double* flop);
  
  void fb_read_sph_s_     (REAL_TYPE* s, 
                           int* sz, 
                           int* g, 
                           char* fname, 
                           int* step, 
                           REAL_TYPE* time, 
                           int* gs, 
                           int* avs,
                           int* step_avr, 
                           REAL_TYPE* time_avr);
  
  void fb_read_sph_v_     (REAL_TYPE* v, 
                           int* sz, 
                           int* g, 
                           char* fname, 
                           int* step, 
                           REAL_TYPE* time, 
                           int* gs,
                           int* avs,
                           int* step_avr, 
                           REAL_TYPE* time_avr);

  void fb_set_vector_     (REAL_TYPE* var, 
                           int* sz, 
                           int* g, 
                           REAL_TYPE* val, 
                           int* bv);
  
  void fb_shift_refv_in_  (REAL_TYPE* v, 
                           int* sz, 
                           int* g, 
                           REAL_TYPE* v00, 
                           REAL_TYPE* scale, 
                           REAL_TYPE* refv, 
                           double* flop);
  
  void fb_shift_refv_out_ (REAL_TYPE* vout, 
                           REAL_TYPE* vin, 
                           int* sz, 
                           int* g, 
                           REAL_TYPE* v00, 
                           REAL_TYPE* scale, 
                           REAL_TYPE* unit_v, 
                           double* flop);
  
  void fb_tmp_d2nd_       (REAL_TYPE* t, 
                           int* sz, 
                           REAL_TYPE* Base_tmp, 
                           REAL_TYPE* Diff_tm, 
                           REAL_TYPE* klv, 
                           REAL_TYPE* scale, 
                           double* flop);
  
  void fb_tmp_nd2d_       (REAL_TYPE* dst, 
                           REAL_TYPE* src, 
                           int* sz, 
                           REAL_TYPE* Base_tmp, 
                           REAL_TYPE* Diff_tm, 
                           REAL_TYPE* klv, 
                           REAL_TYPE* scale, 
                           double* flop);
  
  void fb_totalp_         (REAL_TYPE* tp, 
                           int* sz, 
                           int* g, 
                           REAL_TYPE* v, 
                           REAL_TYPE* p, 
                           REAL_TYPE* v00, 
                           double* flop);
  
  void fb_write_sph_s_    (REAL_TYPE* s, 
                           int* sz, 
                           int* g, 
                           char* fname, 
                           int* step, 
                           REAL_TYPE* time, 
                           REAL_TYPE* org, 
                           REAL_TYPE* pit, 
                           int* d_type, 
                           int* gs, 
                           int* avs,
                           int* step_avr, 
                           REAL_TYPE* time_avr);
  
  void fb_write_sph_v_    (REAL_TYPE* v, 
                           int* sz, 
                           int* g, 
                           char* fname, 
                           int* step, 
                           REAL_TYPE* time, 
                           REAL_TYPE* org, 
                           REAL_TYPE* pit, 
                           int* d_type, 
                           int* gs, 
                           int* avs,
                           int* step_avr, 
                           REAL_TYPE* time_avr);
  
  void fb_mulcpy_         (REAL_TYPE* dst, 
                           REAL_TYPE* src, 
                           int* sz, 
                           int* g, 
                           REAL_TYPE* cf, 
                           double* flop);
}

#endif // _FB_F_FUNC_H_

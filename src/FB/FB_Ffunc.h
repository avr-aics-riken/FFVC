#ifndef _FB_F_FUNC_H_
#define _FB_F_FUNC_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   FB_Ffunc.h
 * @brief  FlowBase Fortran function Header
 * @author kero
 */


#include "cpm_Define.h"
#include "FB_Define.h"

#ifdef _WIN32
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
#define fb_set_vector_           FB_SET_VECTOR
#define fb_vin_nijk_             FB_VIN_NIJK
#define fb_vout_nijk_            FB_VOUT_NIJK
#define fb_vout_ijkn_            FB_VOUT_IJKN
#define fb_totalp_               FB_TOTALP

#endif // _WIN32

extern "C" {
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
  
  // for CIOlib debugging
  void s3dwrite_          (int* ID,
                           REAL_TYPE* v,
                           int* sz,
                           int* g);
  
  void v3dwrite_          (int* ID,
                           REAL_TYPE* v,
                           int* sz,
                           int* g);
  
  void nv3dwrite_         (int* ID,
                           REAL_TYPE* v,
                           int* sz,
                           int* g);
}

#endif // _FB_F_FUNC_H_

#ifndef _SKL_FORTRAN_FUNC_FB_H_
#define _SKL_FORTRAN_FUNC_FB_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file FB_Ffunc.h
//@brief FlowBase Fortran function Header
//@author keno, FSI Team, VCAD, RIKEN

#include "FBDefine.h"

#ifdef _WIN32
// FB_util.f90
#define fb_average_s_         FB_AVERAGE_S
#define fb_average_v_         FB_AVERAGE_V
#define fb_copy_real_s_       FB_COPY_REAL_S
#define fb_copy_real_v_       FB_COPY_REAL_V
#define fb_delta_s_           FB_DELTA_S
#define fb_delta_v_           FB_DELTA_V
#define fb_limit_scalar_      FB_LIMIT_SCALAR
#define fb_minmax_s_          FB_MINMAX_S
#define fb_minmax_v_          FB_MINMAX_V
#define fb_set_int_s_         FB_SET_INT_S
#define fb_set_real_s_        FB_SET_REAL_S
#define fb_set_real_v_        FB_SET_REAL_V
#define fb_set_vector_        FB_SET_VECTOR
#define fb_totalp_            FB_TOTALP

#endif // _WIN32

extern "C" {
	// FB_util.f90
  void fb_average_s_     (REAL_TYPE* avr, int* sz, int* g, REAL_TYPE* s, REAL_TYPE* flop);
  void fb_average_v_     (REAL_TYPE* avr, int* sz, int* g, REAL_TYPE* v, REAL_TYPE* flop);
  void fb_copy_real_s_   (REAL_TYPE* dst, REAL_TYPE* src, int* sz, int* g);
  void fb_copy_real_v_   (REAL_TYPE* dst, REAL_TYPE* src, int* sz, int* g);
  void fb_delta_s_       (REAL_TYPE* d, int* sz, int* g, REAL_TYPE* sn, REAL_TYPE* so, int* bx, REAL_TYPE* flop);
  void fb_delta_v_       (REAL_TYPE* d, int* sz, int* g, REAL_TYPE* vn, REAL_TYPE* vo, int* bx, REAL_TYPE* flop);
  void fb_limit_scalar_  (REAL_TYPE* t, int* sz, int* g);
  void fb_minmax_s_      (REAL_TYPE* f_min, REAL_TYPE* f_max, int* sz, int* g, REAL_TYPE* s, REAL_TYPE* flop);
  void fb_minmax_v_      (REAL_TYPE* f_min, REAL_TYPE* f_max, int* sz, int* g, REAL_TYPE* v00, REAL_TYPE* v, REAL_TYPE* flop);
  void fb_set_int_s_     (int* var, int* sz, int* g, int* init);
  void fb_set_real_s_    (REAL_TYPE* var, int* sz, int* g, REAL_TYPE* init);
  void fb_set_real_v_    (REAL_TYPE* var, int* sz, int* g, REAL_TYPE* init);
  void fb_set_vector_    (REAL_TYPE* var, int* sz, int* g, REAL_TYPE* val);
  void fb_totalp_        (REAL_TYPE* tp,  int* sz, int* g, REAL_TYPE* v, REAL_TYPE* p, REAL_TYPE* v00, REAL_TYPE* flop);
}

#endif // _SKL_FORTRAN_FUNC_FB_H_

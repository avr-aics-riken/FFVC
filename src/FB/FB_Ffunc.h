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

#include "Skl.h"

#ifdef _WIN32
// FB_util.f90
#define fb_average_s_         FB_AVERAGE_S
#define fb_average_v_         FB_AVERAGE_V
#define fb_copy_real_         FB_COPY_REAL
#define fb_limit_scalar_      FB_LIMIT_SCALAR
#define fb_minmax_v_          FB_MINMAX_V
#define fb_minmax_s_          FB_MINMAX_S
#define fb_set_value_int_     FB_SET_VALUE_INT
#define fb_set_value_real_    FB_SET_VALUE_REAL
#define fb_set_vector_        FB_SET_VECTOR
#define fb_totalp_            FB_TOTALP

#endif // _WIN32

extern "C" {
	// FB_util.f90
  void fb_average_s_     (SKL_REAL* avr, int* sz, int* g, SKL_REAL* s, SKL_REAL* flop);
  void fb_average_v_     (SKL_REAL* avr, int* sz, int* g, SKL_REAL* v, SKL_REAL* flop);
  void fb_copy_real_     (SKL_REAL* dst, SKL_REAL* src, int* sz);
  void fb_limit_scalar_  (SKL_REAL* t, int* sz, int* g);
  void fb_minmax_s_      (SKL_REAL* f_min, SKL_REAL* f_max, int* sz, int* g, SKL_REAL* s, SKL_REAL* flop);
  void fb_minmax_v_      (SKL_REAL* f_min, SKL_REAL* f_max, int* sz, int* g, SKL_REAL* v00, SKL_REAL* v, SKL_REAL* flop);
  void fb_set_value_int_ (int* var, int* sz, int* init);
  void fb_set_value_real_(SKL_REAL* var, int* sz, SKL_REAL* init);
  void fb_set_vector_    (SKL_REAL* var, int* sz, int* g, SKL_REAL* val);
  void fb_totalp_        (SKL_REAL* tp,  int* sz, int* g, SKL_REAL* v, SKL_REAL* p, SKL_REAL* v00, int* gs, SKL_REAL* flop);
}

#endif // _SKL_FORTRAN_FUNC_FB_H_

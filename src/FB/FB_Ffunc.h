#ifndef _SKL_FORTRAN_FUNC_FB_H_
#define _SKL_FORTRAN_FUNC_FB_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file FB_Ffunc.h
//@brief FlowBase Fortran function Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"

#ifdef _WIN32
// FB_util.f90
#define fb_limit_scalar_      FB_LIMIT_SCALAR
#define fb_minmax_v_          FB_MINMAX_V
#define fb_minmax_s_          FB_MINMAX_S
#define fb_setdt_cp_          FB_SETDT_CP
#define fb_set_vector_        FB_SET_VECTOR
#define fb_totalp_            FB_TOTALP

#endif // _WIN32

extern "C" {
	// FB_util.f90
  void fb_limit_scalar_  (SKL_REAL* t, int* sz, int* g);
  void fb_minmax_s_      (SKL_REAL* f_min, SKL_REAL* f_max, int* sz, int* g, SKL_REAL* s);
  void fb_minmax_v_      (SKL_REAL* f_min, SKL_REAL* f_max, int* sz, int* g, SKL_REAL* v00, SKL_REAL* v);
  void fb_setdt_cp_      (SKL_REAL* q, SKL_REAL* gd, SKL_REAL* dx, SKL_REAL* dy, SKL_REAL* dz, SKL_REAL* cfl, SKL_REAL* dt, SKL_REAL* dxmin, int* sz, int* jm, int* km, int* lm);
  void fb_set_vector_    (SKL_REAL* var, int* sz, int* g, SKL_REAL* val);
  void fb_totalp_        (SKL_REAL* tp,  int* sz, int* g, SKL_REAL* v, SKL_REAL* p, SKL_REAL* v00, int* gs, SKL_REAL* flop);
}

#endif // _SKL_FORTRAN_FUNC_FB_H_

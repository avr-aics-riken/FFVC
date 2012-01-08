#ifndef _SKL_FB_CORE_UTIL_H_
#define _SKL_FB_CORE_UTIL_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file Core_Util.h
//@brief FlowBase Core_Utility class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <math.h>

#include "Skl.h"
#include "SklSolverBase.h"
#include "FBDefine.h"
#include "FB_Ffunc.h"
#include "Control.h"
#include "FBUtility.h"
#include "Parallel_node.h"
//#include "BndOuter.h"

class Core_Utility : public Parallel_Node {
  
public:
  Core_Utility() {}
  ~Core_Utility() {}
 
public:
  static bool shiftVin3D  (SklVector3D<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr=1);
  static bool shiftVout3D (SklVector3DEx<SKL_REAL>* dst, const SklVector3D<SKL_REAL>* src, SKL_REAL v00[3], unsigned stepAvr=1);
  
  SKL_REAL count_comm_size (unsigned sz[3], unsigned guide) const;
  SKL_REAL norm_v_div_l2   (unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const;
  SKL_REAL norm_v_div_max  (unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const;
  
  void average_Vector  (SKL_REAL* av, unsigned sz[3], unsigned guide, SKL_REAL* v, unsigned* bx, SKL_REAL* avr, SKL_REAL& flop) const;
	void delta_Scalar    (unsigned sz[3], unsigned guide, SKL_REAL* sn, SKL_REAL* so, unsigned* bd, SKL_REAL* var, SKL_REAL& flop) const;
  void delta_Vector    (unsigned sz[3], unsigned guide, SKL_REAL* vn, SKL_REAL* vo, unsigned* bd, SKL_REAL* var, SKL_REAL& flop) const;
  void CutOffRange     (SKL_REAL* t, CompoList* compo, unsigned* bx, Control* C) const;
  void norm_v_div_dbg  (SKL_REAL& nrm, SKL_REAL& rm, int* index, unsigned sz[3], unsigned guide, SKL_REAL coef, SKL_REAL* src, unsigned* bp, SKL_REAL& flop) const;
  void TotalPressure   (SKL_REAL* tp, SKL_REAL* v, SKL_REAL* p, Control* C, const char* arrangement, SKL_REAL* v00, SKL_REAL& flop) const;
  
  /**
   @fn void copy_SKL_REAL(SKL_REAL* dst, SKL_REAL* src, unsigned nx) const
   @brief SKL_REAL型のデータをコピーする
   @param dst ベクトル コピー先
   @param src ベクトル コピー元
   @param nx データサイズ
   */
  void copy_SKL_REAL(SKL_REAL* dst, SKL_REAL* src, unsigned nx) const
  {
    for (unsigned i=0; i<nx; i++)
      *dst++ = *src++;
  }
  
  /**
   @fn static inline void CalcIndex(cosnt unsigned dst_ilen, const unsigned dst_jlen, const unsigned dst_klen, cosnt unsigned dst_gc, const unsigned src_gc, 
   unsigned& dst_ix, unsigned& dst_jx, unsigned& dst_kx, int& diff, unsigned& sta)
   @brief calculation of index, this method is taken from SklUtil class
   @param dst_ilen dimension size /w guide cell for dst
   @param dst_jlen 
   @param dst_klen 
   @param dst_gc size of guide cell
   @param src_gc 
   @param dst_ix return value
   @param dst_jx 
   @param dst_kx 
   @param diff 
   @param sta 
   */
  static inline void CalcIndex(const unsigned dst_ilen,
                        const unsigned dst_jlen,
                        const unsigned dst_klen,
                        const unsigned dst_gc,
                        const unsigned src_gc,
                        unsigned& dst_ix,
                        unsigned& dst_jx,
                        unsigned& dst_kx,
                        int& diff,
                        unsigned& sta) {
    diff = src_gc - dst_gc;
    if( src_gc >= dst_gc ){
      sta = 0;
      dst_ix = dst_ilen; dst_jx = dst_jlen; dst_kx = dst_klen;
    }
    else {
      sta = abs(diff);
      dst_ix = dst_ilen+diff; 
      dst_jx = dst_jlen+diff; 
      dst_kx = dst_klen+diff;
    }
  }
};

#endif // _SKL_FB_CORE_UTIL_H_
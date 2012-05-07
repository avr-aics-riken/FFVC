#ifndef _SKL_FB_FILE_IO_H_
#define _SKL_FB_FILE_IO_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file FileIO.h
//@brief FlowBase FileIO class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"
#include "FBDefine.h"
#include "SklUtil.h"
#include <fstream>
#include "Parallel_node.h"
#include "fileio/SklSbxDataSet.h"
#include "parallel/SklParaComponent.h"
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include "FBUtility.h"
#include "FB_Ffunc.h"

#ifndef _WIN32
#include <unistd.h>
#include <strings.h>
#else
#include "sph_win32_util.h"
#endif
#include <sys/types.h>

#if defined(IA32_LINUX) || defined(IA64_LINUX) || defined(SGI_ALTIX)
#include <sys/stat.h>
#endif

#ifdef MacOSX
#include <sys/uio.h>
#endif

class FileIO : public Parallel_Node {
  
public:
  FileIO() {}
  ~FileIO() {}
  
protected:
  /// sphファイルのスカラ/ベクトルの種別
  enum sv_type {
    kind_scalar=1,
    kind_vector
  };
  
  void loadSBXfile (SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, SklScalar3D<unsigned char>* data);
  void loadSBXfile (SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, SklScalar3D<int>* mid_data, 
                    SklScalar3D<REAL_TYPE>* vol_data = NULL);
  
public:
  void cnv_Div           (SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE coef, REAL_TYPE& flop);
  void cnv_TP_ND2D       (SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, REAL_TYPE& flop);
  void loadPressure      (SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, SklScalar3D<REAL_TYPE>* dc_p, 
                          int& step, REAL_TYPE& time, const unsigned Dmode, const REAL_TYPE BasePrs, const REAL_TYPE RefDensity, 
                          const REAL_TYPE RefVelocity, REAL_TYPE& flop, const bool mode=true);
  void loadTemperature   (SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, SklScalar3D<REAL_TYPE>* dc_t, 
                          int& step, REAL_TYPE& time, const unsigned Dmode, const REAL_TYPE Base_tmp, const REAL_TYPE Diff_tmp, 
                          const REAL_TYPE klv, REAL_TYPE& flop, const bool mode=true);
  void loadVelocity      (SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, SklVector3DEx<REAL_TYPE>* dc_v, 
                          int& step, REAL_TYPE& time, const REAL_TYPE *v00, const unsigned Dmode, const REAL_TYPE RefVelocity, 
                          REAL_TYPE& flop, const bool mode=true);
  void readSBX           (SklSolverBase* obj, FILE* fp, const char* file_attr, unsigned* size, unsigned guide, 
                          SklScalar3D<unsigned char>* dc_mid);
  void readSBX           (SklSolverBase* obj, FILE* fp, const char* mid_str, unsigned* size, unsigned guide, 
                          SklScalar3D<int>* dc_mid, SklScalar3D<REAL_TYPE>* dc_vol = NULL);
  void readSVX           (SklSolverBase* obj, FILE* fp, const char* fname, unsigned* size, unsigned guide, 
                          SklScalar3D<int>* dc_mid, bool vf_mode=false, SklScalar3D<REAL_TYPE>* dc_ws=NULL);
  void writeRawSPH       (const REAL_TYPE *vf, const unsigned* size, const unsigned gc, const REAL_TYPE* org, const REAL_TYPE* ddx, 
                          const unsigned m_ModePrecision);
  
  void loadVector(FILE* fp, char* fname, const unsigned* size, const unsigned gc, 
                  REAL_TYPE* v, int& step, REAL_TYPE& time, const REAL_TYPE *v00, const unsigned Dmode, 
                  const REAL_TYPE RefVelocity, REAL_TYPE& flop, const bool mode=true);
  void loadScalar(FILE* fp, char* fname, const unsigned* size, const unsigned gc,
                  REAL_TYPE* s, int& step, REAL_TYPE& time, 
                  const unsigned Dmode, 
                  const REAL_TYPE BasePrs, 
                  const REAL_TYPE RefDensity, 
                  const REAL_TYPE RefVelocity, REAL_TYPE& flop, const bool mode=true);
  void loadTemp(FILE* fp, char* fname, const unsigned* size, const unsigned gc, 
                REAL_TYPE* t, int& step, REAL_TYPE& time, 
                const unsigned Dmode, 
                const REAL_TYPE Base_tmp, 
                const REAL_TYPE Diff_tmp, 
                const REAL_TYPE Kelvin, REAL_TYPE& flop, const bool mode=true);
  
  /**
   @fn inline void CalcIndex(cosnt unsigned dst_ilen, const unsigned dst_jlen, const unsigned dst_klen, cosnt unsigned dst_gc, const unsigned src_gc, 
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
  inline void CalcIndex(const unsigned dst_ilen,
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
#endif // _SKL_FB_FILE_IO_H_

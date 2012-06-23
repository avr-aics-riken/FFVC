#ifndef _FB_FILE_IO_H_
#define _FB_FILE_IO_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file FileIO.h
 * @brief FlowBase FileIO class Header
 * @author kero
 */

#include <fstream>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <string>

#include "cpm_Define.h"
#include "cpm_ParaManager.h"

#include "FB_Define.h"
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

class FileIO {
  
public:
  /** コンストラクタ */
  FileIO() 
  {
    paraMngr = NULL;
  }
  
  /**　デストラクタ */
  ~FileIO() {}
  
protected:
  cpm_ParaManager *paraMngr; ///< Cartesian Partition Maneger
  
  /** sphファイルのスカラ/ベクトルの種別 */
  enum sv_type {
    kind_scalar=1,
    kind_vector
  };
  
  
public:
  void cnv_Div           (SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE coef, REAL_TYPE& flop);
  void cnv_TP_ND2D       (SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, REAL_TYPE& flop);

  void writeRawSPH       (const REAL_TYPE *vf, const unsigned* size, const unsigned gc, const REAL_TYPE* org, const REAL_TYPE* ddx, 
                          const unsigned m_ModePrecision);
  
  
  void readPressure(FILE* fp, 
                    const std::string fname, 
                    const unsigned* size, 
                    const unsigned gc,
                    REAL_TYPE* s, 
                    int& step, 
                    REAL_TYPE& time, 
                    const unsigned Dmode, 
                    const REAL_TYPE BasePrs, 
                    const REAL_TYPE RefDensity, 
                    const REAL_TYPE RefVelocity, 
                    REAL_TYPE& flop, 
                    const int guide_out,
                    const bool mode,
                    int& step_avr,
                    REAL_TYPE& time_avr);
  
  void readVelocity(FILE* fp, 
                    const std::string fname, 
                    const unsigned* size, 
                    const unsigned gc, 
                    REAL_TYPE* v, 
                    int& step, 
                    REAL_TYPE& time, 
                    const REAL_TYPE *v00, 
                    const unsigned Dmode, 
                    const REAL_TYPE RefVelocity, 
                    REAL_TYPE& flop, 
                    const int guide_out,
                    const bool mode,
                    int& step_avr,
                    REAL_TYPE& time_avr);
  
  void readTemperature(FILE* fp, 
                       const std::string fname, 
                       const unsigned* size, 
                       const unsigned gc, 
                       REAL_TYPE* t, 
                       int& step, 
                       REAL_TYPE& time, 
                       const unsigned Dmode, 
                       const REAL_TYPE Base_tmp, 
                       const REAL_TYPE Diff_tmp, 
                       const REAL_TYPE Kelvin, 
                       REAL_TYPE& flop, 
                       const int guide_out,
                       const bool mode,
                       int& step_avr,
                       REAL_TYPE& time_avr);
  
  void writeScalar(const std::string fname, 
                   const unsigned* size, 
                   const unsigned gc,
                   REAL_TYPE* s, 
                   const int step, 
                   const REAL_TYPE time, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* pit, 
                   const int guide_out,
                   const bool mode=true,
                   const int step_avr=0,
                   const REAL_TYPE time_avr=0.0);
  
  void writeVector(const std::string fname, 
                   const unsigned* size, 
                   const unsigned gc, 
                   REAL_TYPE* v, 
                   const int step, 
                   const REAL_TYPE time, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* pit, 
                   const int guide_out,
                   const bool mode=true,
                   const int step_avr=0,
                   const REAL_TYPE time_avr=0.0);
  
  
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
  
  /** CPMlibのポインタをセット 
   * @param[in] m_paraMngr  
   */
  void setPartitionManager(cpm_ParaManager* m_paraMngr)
  
};
#endif // _FB_FILE_IO_H_

#ifndef _FB_FILE_IO_H_
#define _FB_FILE_IO_H_

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
 * @file   FileIO.h
 * @brief  FlowBase FileIO class Header
 * @author kero
 */

#include <fstream>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <string>

#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "FB_Ffunc.h"
#include "FBUtility.h"


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

using namespace std;

class FileIO : public DomainInfo {
  
public:
  /** コンストラクタ */
  FileIO() {}
  
  /**　デストラクタ */
  ~FileIO() {}
  
  
public:
  
  FBUtility U;

  
  // sphファイルの書き出し（内部領域のみ）
  void writeRawSPH(const REAL_TYPE *vf, 
                   const int* sz, 
                   const int gc, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* ddx, 
                   const int m_ModePrecision);
  
  // 圧力のファイルをロードする
  void readPressure(FILE* fp, 
                    const string fname, 
                    int* sz, 
                    int gc,
                    REAL_TYPE* s, 
                    unsigned& step, 
                    double& time, 
                    const int Dmode, 
                    const REAL_TYPE BasePrs, 
                    const REAL_TYPE RefDensity, 
                    const REAL_TYPE RefVelocity, 
                    double& flop, 
                    const int guide_out,
                    const bool mode,
                    unsigned& step_avr,
                    double& time_avr);
  
  
  // 速度のファイルをロードする
  void readVelocity(FILE* fp, 
                    const string fname, 
                    int* sz, 
                    int gc, 
                    REAL_TYPE* v,
                    REAL_TYPE* v_buf,
                    unsigned& step, 
                    double& time, 
                    const REAL_TYPE *v00, 
                    const int Dmode, 
                    const REAL_TYPE RefVelocity, 
                    double& flop, 
                    const int guide_out,
                    const bool mode,
                    unsigned& step_avr,
                    double& time_avr);
  
  
  // 温度のファイルをロードする
  void readTemperature(FILE* fp, 
                       const string fname, 
                       int* sz, 
                       int gc, 
                       REAL_TYPE* t, 
                       unsigned& step, 
                       double& time, 
                       const int Dmode, 
                       const REAL_TYPE Base_tmp, 
                       const REAL_TYPE Diff_tmp, 
                       const REAL_TYPE Kelvin, 
                       double& flop, 
                       const int guide_out,
                       const bool mode,
                       unsigned& step_avr,
                       double& time_avr);
  
  
  // スカラーファイルを出力する
  void writeScalar(const string fname, 
                   int* sz, 
                   int gc,
                   REAL_TYPE* s, 
                   const unsigned step, 
                   const double time, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* pit, 
                   const int guide_out,
                   const bool mode=true,
                   const unsigned step_avr=0,
                   const double time_avr=0.0);
  
  
  // ベクトルファイルを出力する
  void writeVector(const string fname, 
                   int* sz, 
                   int gc, 
                   REAL_TYPE* v, 
                   const unsigned step, 
                   const double time, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* pit, 
                   const int guide_out,
                   const bool mode=true,
                   const unsigned step_avr=0,
                   const double time_avr=0.0);
  
};
#endif // _FB_FILE_IO_H_

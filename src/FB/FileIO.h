#ifndef _FB_FILE_IO_H_
#define _FB_FILE_IO_H_

//##################################################################################
//
// Flow Base class
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
 * @file   FileIO.h
 * @brief  FlowBase FileIO class Header
 * @author kero
 */

#include "cpm_ParaManager.h"
#include <fstream>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <string>
#include "DomainInfo.h"
#include "FB_Define.h"
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

using namespace std;

class FileIO : public DomainInfo {
  
public:
  /** コンストラクタ */
  FileIO() {}
  
  /**　デストラクタ */
  ~FileIO() {}
  
  
public:

  /**
   * @brief sphファイルの書き出し（内部領域のみ）
   * @param [in] vf               スカラデータ
   * @param [in] sz               配列サイズ
   * @param [in] gc               ガイドセル
   * @param [in] org              基点
   * @param [in] ddx              ピッチ
   * @param [in] m_ModePrecision  浮動小数点の精度
   * @note 標記上，long 対応になっているが，ファイルフォーマットとの対応を確認のこと
   */
  void writeRawSPH(const REAL_TYPE *vf, 
                   const int* sz, 
                   const int gc, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* ddx, 
                   const int m_ModePrecision);
  
  
  /**
   * @brief スカラーファイルを出力する
   * @param [in] fname     ファイル名
   * @param [in] sz        分割数
   * @param [in] gc        ガイドセル数
   * @param [in] s         スカラー場
   * @param [in] step      ステップ
   * @param [in] time      時刻
   * @param [in] org       領域の基点
   * @param [in] pit       セル幅
   * @param [in] guide_out ガイドセル数
   * @param [in] mode      平均値出力指示（瞬時値のときtrue，平均値のときfalse）
   * @param [in] step_avr  平均操作したステップ数
   * @param [in] time_avr  平均操作した時間
   */
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
  
  
  /**
   * @brief ベクトルファイルを出力する
   * @param [in] fname     ファイル名
   * @param [in] sz        分割数
   * @param [in] gc        ガイドセル数
   * @param [in] v         ベクトル場
   * @param [in] step      ステップ
   * @param [in] time      時刻
   * @param [in] org       領域の基点
   * @param [in] pit       セル幅
   * @param [in] guide_out ガイドセル数
   * @param [in] mode      平均値出力指示（瞬時値のときtrue，平均値のときfalse）
   * @param [in] step_avr  平均操作したステップ数
   * @param [in] time_avr  平均操作した時間
   */
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

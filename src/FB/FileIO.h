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

#include "cpm_ParaManager.h"

#include "FB_Define.h"
#include "FBUtility.h"
#include "FB_Ffunc.h"
#include "mydebug.h"


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
  
  /**
   * @brief ファイル出力時，発散値を計算する
   * @param [out]    dst  単位変換後のデータ
   * @param [in]     src  単位変換前のデータ
   * @param [in]     sz   分割数
   * @param [in]     gc   ガイドセル数
   * @param [in]     coef 係数
   * @param [in/out] flop 浮動小数演算数
   */
  void cnv_Div(REAL_TYPE* dst, 
               REAL_TYPE* src, 
               const int* size, 
               const int guide, 
               const REAL_TYPE coef, 
               REAL_TYPE& flop);
  
  
  /**
   * @brief 全圧データについて，無次元から有次元単位に変換する
   * @param [out]    dst     単位変換後のデータ
   * @param [in]     src     単位変換前のデータ
   * @param [in]     sz      分割数
   * @param [in]     gc      ガイドセル数
   * @param [in]     Ref_rho 代表密度(kg/m^3)
   * @param [in]     Ref_v   代表速度(m/s)
   * @param [in/out] flop    浮動小数演算数
   */
  void cnv_TP_ND2D(REAL_TYPE* dst, 
                   REAL_TYPE* src, 
                   const int* size, 
                   const int guide, 
                   const REAL_TYPE Ref_rho, 
                   const REAL_TYPE Ref_v, 
                   REAL_TYPE& flop);

  
  /**
   * @brief sphファイルの書き出し（内部領域のみ）
   * @param [in] vf               スカラデータ
   * @param [in] size             配列サイズ
   * @param [in] gc               ガイドセル
   * @param [in] org              基点
   * @param [in] ddx              ピッチ
   * @param [in] m_ModePrecision  浮動小数点の精度
   * @note 標記上，long 対応になっているが，ファイルフォーマットとの対応を確認のこと
   */
  void writeRawSPH(const REAL_TYPE *vf, 
                   const int* size, 
                   const int gc, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* ddx, 
                   const int m_ModePrecision);
  
  /**
   * @brief 圧力のファイルをロードする
   * @param [in]     fp          ファイルポインタ（ファイル出力）
   * @param [in]     fname       ファイル名
   * @param [in]     size        サイズ
   * @param [in]     guide       ガイドセルサイズ
   * @param [out]    p           圧力データポインタ
   * @param [out]    step        ステップ
   * @param [out]    time        時刻
   * @param [in]     Dmode       次元（無次元-0 / 有次元-1）
   * @param [in]     BasePrs     基準圧力
   * @param [in]     RefDensity　代表密度
   * @param [in]     RefVelocity 代表速度
   * @param [in/out] flop        浮動小数点演算数
   * @param [in]     guide_out   出力ガイドセル数
   * @param [in]     mode        平均値出力指示（瞬時値のときtrue，平均値のときfalse）
   * @param [out]    step_avr    平均操作したステップ数
   * @param [out]    time_avr    平均操作した時間
   */
  void readPressure(FILE* fp, 
                    const std::string fname, 
                    int* size, 
                    int gc,
                    REAL_TYPE* s, 
                    int& step, 
                    REAL_TYPE& time, 
                    const int Dmode, 
                    const REAL_TYPE BasePrs, 
                    const REAL_TYPE RefDensity, 
                    const REAL_TYPE RefVelocity, 
                    REAL_TYPE& flop, 
                    const int guide_out,
                    const bool mode,
                    int& step_avr,
                    REAL_TYPE& time_avr);
  
  
  /**
   * @brief 速度のファイルをロードする
   * @param [in]     fp          ファイルポインタ（ファイル出力）
   * @param [in]     fname       ファイル名
   * @param [in]     size        サイズ
   * @param [in]     guide       ガイドセルサイズ
   * @param [out]    v           速度データポインタ
   * @param [out]    step        ステップ
   * @param [out]    time        時刻
   * @param [in]     Dmode       次元（無次元-0 / 有次元-1）
   * @param [in]     BasePrs     基準圧力
   * @param [in]     RefVelocity 代表速度
   * @param [in/out] flop        浮動小数点演算数
   * @param [in]     guide_out   出力ガイドセル数
   * @param [in]     mode        平均値出力指示（瞬時値のときtrue，平均値のときfalse）
   * @param [out]    step_avr    平均操作したステップ数
   * @param [out]    time_avr    平均操作した時間
   */
  void readVelocity(FILE* fp, 
                    const std::string fname, 
                    int* size, 
                    int gc, 
                    REAL_TYPE* v, 
                    int& step, 
                    REAL_TYPE& time, 
                    const REAL_TYPE *v00, 
                    const int Dmode, 
                    const REAL_TYPE RefVelocity, 
                    REAL_TYPE& flop, 
                    const int guide_out,
                    const bool mode,
                    int& step_avr,
                    REAL_TYPE& time_avr);
  
  
  /**
   * @brief 温度のファイルをロードする
   * @param [in]     fp          ファイルポインタ（ファイル出力）
   * @param [in]     fname       ファイル名
   * @param [in]     size        サイズ
   * @param [in]     guide       ガイドセルサイズ
   * @param [out]    t           温度データポインタ
   * @param [out]    step        ステップ
   * @param [out]    time        時刻
   * @param [in]     Dmode       次元（無次元-0 / 有次元-1）
   * @param [in]     Base_tmp    基準温度
   * @param [in]     Diff_tmp  　代表温度差
   * @param [in]     Kelvin      定数
   * @param [in/out] flop        浮動小数点演算数
   * @param [in]     guide_out   出力ガイドセル数
   * @param [in]     mode        平均値出力指示（瞬時値のときtrue，平均値のときfalse）
   * @param [out]    step_avr    平均操作したステップ数
   * @param [out]    time_avr    平均操作した時間
   */
  void readTemperature(FILE* fp, 
                       const std::string fname, 
                       int* size, 
                       int gc, 
                       REAL_TYPE* t, 
                       int& step, 
                       REAL_TYPE& time, 
                       const int Dmode, 
                       const REAL_TYPE Base_tmp, 
                       const REAL_TYPE Diff_tmp, 
                       const REAL_TYPE Kelvin, 
                       REAL_TYPE& flop, 
                       const int guide_out,
                       const bool mode,
                       int& step_avr,
                       REAL_TYPE& time_avr);
  
  
  /** 
   * @brief スカラーファイルを出力する
   * @param [in] fname     ファイル名
   * @param [in] size      分割数
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
  void writeScalar(const std::string fname, 
                   int* size, 
                   int gc,
                   REAL_TYPE* s, 
                   const int step, 
                   const REAL_TYPE time, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* pit, 
                   const int guide_out,
                   const bool mode=true,
                   const int step_avr=0,
                   const REAL_TYPE time_avr=0.0);
  
  
  /** 
   * @brief ベクトルファイルを出力する
   * @param [in] fname     ファイル名
   * @param [in] size      分割数
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
  void writeVector(const std::string fname, 
                   int* size, 
                   int gc, 
                   REAL_TYPE* v, 
                   const int step, 
                   const REAL_TYPE time, 
                   const REAL_TYPE* org, 
                   const REAL_TYPE* pit, 
                   const int guide_out,
                   const bool mode=true,
                   const int step_avr=0,
                   const REAL_TYPE time_avr=0.0);
  
  
  /** CPMlibのポインタをセット 
   * @param [in] m_paraMngr  CPMクラスのポインタ
   */
  void importCPM(cpm_ParaManager* m_paraMngr);
  
};
#endif // _FB_FILE_IO_H_

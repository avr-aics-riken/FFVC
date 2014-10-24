#ifndef _FFV_IO_BASE_H_
#define _FFV_IO_BASE_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
//

/**
 * @file   ffv_io_base.h
 * @brief  File IO base Class Header
 * @author aics
 */

#include "cpm_ParaManager.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>

#include "omp.h"
#include "../FB/DomainInfo.h"
#include "../FB/FB_Define.h"
#include "../FB/mydebug.h"
#include "../FB/Control.h"
#include "../FB/FBUtility.h"
#include "../FB/IntervalManager.h"
#include "../FFV/ffv_Define.h"

#include "../F_LS/ffv_LSfunc.h"
#include "../F_CORE/ffv_Ffunc.h"

#include "TextParser.h"

// CDMlib
#include "cdm_DFI.h"

using namespace std;

class IO_BASE : public DomainInfo {

protected:
  
  int IOmode;          ///< 逐次 or 並列
  int IO_Voxel;        ///< デバッグ用にボクセルを出力
  int Format;          ///< ファイル入出力モード（sph, bov, plot3d）
  int Slice;           ///< タイムスライス毎にまとめる
  int Iblank;          ///< PLOT3DのIBLANKオプション
  string OutDirPath;   ///< 出力ディレクトリパス
  string InDirPath;    ///< 入力ディレクトリパス
  string file_fmt_ext; ///< フォーマット識別子
  
  
  // variables
  REAL_TYPE* d_p;    ///< pressure
  REAL_TYPE* d_v;    ///< velocity
  REAL_TYPE* d_vf;   ///< face velocity
  REAL_TYPE* d_ie;   ///< internal energy
  REAL_TYPE* d_iobuf;///< IO buffer
  REAL_TYPE* d_ws;   ///< work for scalar
  REAL_TYPE* d_p0;   ///< work for scalar
  REAL_TYPE* d_wv;   ///< work for vector
  REAL_TYPE* d_vc;   ///< work for vector
  REAL_TYPE* d_ap;   ///< averaged pressure
  REAL_TYPE* d_av;   ///< averaged velocity
  REAL_TYPE* d_ae;   ///< averaged internal energy
  REAL_TYPE* d_dv;   ///< Divergence
  int* d_bcd;        ///< BCindex D
  int* d_cdf;        ///< BCindex C
  double* mat_tbl;   ///< material table
  int* d_iblk;       ///< IBLANK
  
  
  // class pointer
  Control* C;
  ReferenceFrame* RF;
  TextParser* tpCntl;
  
  // Utility class
  FBUtility U;
  
  
public:
  
  // default constructor
  IO_BASE() {
    
    IOmode   = 0;
    IO_Voxel = 0;
    Format   = 0;
    Slice    = 0;
    Iblank   = 0;
    
    // 変数
    d_p = NULL;
    d_v = NULL;
    d_vf = NULL;
    d_ie = NULL;
    d_ws = NULL;
    d_p0 = NULL;
    d_iobuf = NULL;
    d_wv = NULL;
    d_vc = NULL;
    d_ap = NULL;
    d_av = NULL;
    d_ae = NULL;
    d_dv = NULL;
    d_bcd = NULL;
    mat_tbl = NULL;
    d_bcd = NULL;
    d_cdf = NULL;
    d_iblk = NULL;
  }
  
  
  ~IO_BASE() {}
  
  
protected:
  
  /*
   * @brief フォーマットのオプションを指定
   * @param [in] form    format
   */
  void getFormatOption(const string form);

  
  
  
public:
  
  // @brief 出力ファイルをチェック
  //int checkOutFile();
  
  
  /*
   * @brief ファイル入出力に関するパラメータを取得し，出力の並列モードを指定する．
   * @note インターバルパラメータは，setParameters()で無次元して保持
   * @pre getTimeControl()
   */
  void getFIOparams();

  
  
  // formatを返す
  int getFormat() const
  {
    return Format;
  }
  
  
  // Iblankの利用の有無を返す
  int getIblank() const
  {
    return Iblank;
  }
  
  
  // リスタートに必要なDFIファイルを取得
  virtual void getRestartDFI() {}
  
  
  // ステージングオプション
  void getStagingOption();
  
  
  /*
   * @brief 初期値
   * @see Control::getTimeControl()
   */
  void getStartCondition();
  
  
  /**
   * @brief 外部クラスのポインタを受け取る
   * @param [in]  tp     TextParser class
   * @param [in]  m_RF   ReeferenceFrame class
   * @param [in]  m_C    Control class
   */
  void importExtClass(TextParser* tp, ReferenceFrame* m_RF, Control* m_C)
  {
    if ( !tp || !m_RF || !m_C ) Exit(0);
    
    tpCntl = tp;
    RF = m_RF;
    C = m_C;
  }
  
  
  // @brief ファイル出力の初期化
  virtual void initFileOut() {}
  
  
  /**
   * @brief 時間平均値のファイル出力
   * @param [in]     m_CurrentStep     CurrentStep
   * @param [in]     m_CurrentTime     CurrentTime
   * @param [in]     m_CurrentStepAvr  CurrentStepAvr
   * @param [in]     m_CurrentTimeAvr  CurrentTimeAvr
   * @param [in,out] flop              浮動小数点演算数
   */
  virtual void OutputAveragedVarables(const unsigned m_CurrentStep,
                                      const double m_CurrentTime,
                                      const unsigned m_CurrentStepAvr,
                                      const double m_CurrentTimeAvr,
                                      double& flop) {
  }
  
  
  /**
   * @brief 基本変数のファイル出力
   * @param [in]     m_CurrentStep     CurrentStep
   * @param [in]     m_CurrentTime     CurrentTime
   * @param [in,out] flop              浮動小数点演算数
   */
  virtual void OutputBasicVariables(const unsigned m_CurrentStep,
                                    const double m_CurrentTime,
                                    double& flop) {
  }
  
  
  /**
   * @brief 派生変数のファイル出力
   * @param [in]     m_CurrentStep  CurrentStep
   * @param [in]     m_CurrentTime  CurrentTime
   * @param [in,out] flop           浮動小数点演算数
   * @note d_p0をワークとして使用
   */
  virtual void OutputDerivedVariables(const unsigned m_CurrentStep,
                                      const double m_CurrentTime,
                                      double& flop) {
  }
  
  
  // 制御パラメータSTEERの表示
  void printSteerConditions(FILE* fp);
  
  
  /**
   * @brief リスタートプロセス
   * @param [in]     fp                ファイルポインタ
   * @param [out]    m_CurrentStep     CurrentStep
   * @param [out]    m_CurrentTime     CurrentTime
   */
  virtual void Restart(FILE* fp,
                       unsigned& m_CurrentStep,
                       double& m_CurrentTime) {
  }
  
  
  /**
   * @brief リスタート時の平均値ファイル読み込み
   * @param [in]  fp                ファイルポインタ
   * @param [in]  m_CurrentStep     CurrentStep
   * @param [in]  m_CurrentTime     CurrentTime
   * @param [out] m_CurrentStepAvr  CurrentStepAvr
   * @param [out] m_CurrentTimeAvr  CurrentTimeAvr
   * @param [out] flop              浮動小数点演算数
   */
  virtual void RestartAvrerage(FILE* fp,
                               const unsigned m_CurrentStep,
                               const double m_CurrentTime,
                               unsigned& m_CurrentStepAvr,
                               double& m_CurrentTimeAvr,
                               double& flop) {
  }
  

  /**
   * @brief リスタートの最大値と最小値の表示
   * @param [in]  fp   ファイルポインタ
   * @param [out] flop 浮動小数点演算数
   */
  void RestartDisplayMinmax(FILE* fp, double& flop);
  
  
  
  /**
   * @brief リスタート時の瞬時値ファイル読み込み
   * @param [in]  fp             ファイルポインタ
   * @param [out] m_CurrentStep  CurrentStep
   * @param [out] m_CurrentTime  CurrentTime
   * @param [out] flop           浮動小数点演算数
   */
  virtual void RestartInstantaneous(FILE* fp,
                                    unsigned& m_CurrentStep,
                                    double& m_CurrentTime,
                                    double& flop) {
  }
  
  
  // formatを登録
  void setFormat(const int key)
  {
    Format = key;
  }
  

  // 必要なポインタをセット
  void setVarPointers(REAL_TYPE* m_d_p,
                      REAL_TYPE* m_d_v,
                      REAL_TYPE* m_d_vf,
                      REAL_TYPE* m_d_ie,
                      REAL_TYPE* m_d_ws,
                      REAL_TYPE* m_d_p0,
                      REAL_TYPE* m_d_iob,
                      REAL_TYPE* m_d_wv,
                      REAL_TYPE* m_d_vc,
                      REAL_TYPE* m_d_ap,
                      REAL_TYPE* m_d_av,
                      REAL_TYPE* m_d_ae,
                      REAL_TYPE* m_d_dv,
                      int* m_d_bcd,
                      int* m_d_cdf,
                      double* m_mat_tbl,
                      int* d_iblk);
  
  
  /**
   * @brief sphファイルの書き出し（内部領域のみ）
   * @param [in] vf               スカラデータ
   * @param [in] sz               配列サイズ
   * @param [in] gc               ガイドセル
   * @param [in] gc_out           出力するガイドセル数
   * @param [in] org              基点
   * @param [in] ddx              ピッチ
   * @param [in] m_ModePrecision  浮動小数点の精度
   * @note 標記上，long 対応になっているが，ファイルフォーマットとの対応を確認のこと
   */
  void writeRawSPH(const REAL_TYPE *vf,
                   const int* sz,
                   const int gc,
                   const int gc_out,
                   const REAL_TYPE* org,
                   const REAL_TYPE* ddx,
                   const int m_ModePrecision);
  
  /**
   * @brief 例題のモデルをsvxフォーマットで出力する(体積率とID)
   * @param [in] vf 体積占有率
   * @param [in] id ID情報
   * @retval 出力 >> true，指定フォーマットがSVXでない場合にはfalse
   */
  bool writeSVX(REAL_TYPE *vf, int *id);
  
  
  /**
   * @brief 例題のモデルをsvxフォーマットで出力する(ID)
   * @param [in] bcd BCindex B
   * @retval 出力 >> true，指定フォーマットがSVXでない場合にはfalse
   */
  bool writeSVX(const int* bcd);
  
};

#endif // _FFV_IO_BASE_H_
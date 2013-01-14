#ifndef _FB_DFI_H_
#define _FB_DFI_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

//@file   dfi.h
//@brief  DFI class Header
//@author kero
//@note このクラスは、並列時のみにコールすること。MPI_Initialize(), Finalize()はクラス外で実行。
//      dfiファイルの書き出しは、マスターランクでのみ行う。

#include <string>
#include <stdio.h>
#include <mpi.h>
#include "FB_Define.h"

class DFI {
protected:
  int Num_Node;          ///< MPI並列数
  int my_id;             ///< 自ノードのランク番号（dfiファイルの出力はランク0のみなので自明だが）
  int Gsize[3];          ///< 計算領域全体の分割数
  int div_domain[3];     ///< 全計算領域のノード分割数
  int procGrp;           ///< プロセスグループ（デフォルト0）
  int guide;             ///< ガイドセル数
  int start_type;        ///< セッションのスタートモード（initial_start=0, restart=1, coarse_restart=2, restart_different_nproc=3）
  int* head;             ///< bboxの開始インデクス(C index) [3*Num_Node]
  int* tail;             ///< bboxの終端インデクス(C index) [3*Num_Node]
  REAL_TYPE RefLength;   ///< 代表長さ
  REAL_TYPE RefVelocity; ///< 代表速度
  std::string hostname;  ///< ラベル
  std::string Unit_L;    ///< 長さの単位 (NonDimensional, m)
  std::string Unit_V;    ///< 速度の単位 (NonDimensional, m/s)
  
public:
  DFI() {
    Num_Node   = 0;
    my_id      = 0;
    procGrp    = 0;
    guide      = 0;
    start_type = -1;
    head       = NULL;
    tail       = NULL;
    RefLength  = 0.0;
    RefVelocity= 0.0;
    
    for (int i=0; i<3; i++) {
      Gsize[i]      = 0;
      div_domain[i] = 0;
    }

  }
  ~DFI() {
    if ( head ) delete [] head;
    if ( tail ) delete [] tail;
  }
  
protected:
  
  // Endianを調べる
  bool chekcEndian();
  
  
  // 出力DFIファイル名を作成する
  std::string GenerateDFIname(const std::string prefix);
  
  
  // Indexファイルを出力する
  bool WriteIndex(const std::string dfi_name,
                  const std::string prefix,
                  const std::string dir,
                  const std::string fmt,
                  const unsigned step,
                  const double time,
                  int& dfi_mng,
                  const std::string shape,
                  const int compo,
                  const REAL_TYPE* minmax,
                  const bool avr_mode,
                  const unsigned a_step,
                  const double a_time);

  
  // Rankの情報を出力する
  bool WriteRank(FILE* fp, const unsigned tab, const int id);
  
  
  // 時系列情報を出力
  void WriteTimeSlice(FILE* fp,
                      const unsigned tab,
                      const std::string prefix,
                      const unsigned step,
                      const double time,
                      const REAL_TYPE* minmax,
                      const bool avr_mode,
                      const unsigned a_step,
                      const double a_time);

  
  // Tab(space２つ)を出力する
  void WriteTab(FILE* fp, const unsigned tab);
  
  
public:
  
  // 出力ディレクトリ名を作成する
  std::string GenerateDirName(const std::string prefix, const unsigned m_step);
  
  
  // ファイル名を作成する
  std::string GenerateFileName(const std::string prefix, const std::string fmt, const unsigned m_step, const int m_id, const char* mio);

  
  // 初期化
  bool init(const int* g_size,
            const int* m_div,
            const int gc,
            const int stype,
            const REAL_TYPE m_refL,
            const REAL_TYPE m_refV,
            const std::string m_UnitL,
            const std::string m_UnitV,
            const int* hidx,
            const int* tidx,
            const std::string m_host);
  
  
  // DFI indexファイルを生成する
  bool WriteDFIindex(const std::string prefix,
                     const std::string dir,
                     const std::string fmt,
                     const unsigned step,
                     const double time,
                     int& dfi_mng,
                     const std::string shape,
                     const int compo,
                     const REAL_TYPE* minmax,
                     const char* mio,
                     bool avr_mode=false,
                     unsigned a_step=0,
                     double a_time=0.0);
  
  // DFI procファイルを生成する
  bool WriteDFIproc(const REAL_TYPE* g_org, const REAL_TYPE* g_reg);
};

#endif // _FB_DFI_H_

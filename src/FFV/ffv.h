#ifndef _FFV_H_
#define _FFV_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################
//
// 以下のマクロはcpm_Define.hで定義されている
//   REAL_TYPE
//   X_MINUS, Y_MINUS, Z_MINUS, X_PLUS, Y_PLUS, Z_PLUS
//   X_DIR, Y_DIR, Z_DIR
//   PLUS2MINUS, MINUS2PLUS, BOTH

/** 
 * @file ffv.h
 * @brief FFV Class Header
 * @author kero
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>

#include "cpm_ParaManager.h"
#include "cpm_TextParserDomain.h"

#include "FB_Define.h"
#include "ffv_Define.h"
#include "mydebug.h"
#include "FBUtility.h"
#include "Control.h"
#include "Alloc.h"
#include "FileIO.h"
#include "ParseBC.h"
#include "ParseMat.h"
//#include "ffv_SetBC.h"
//#include "VoxInfo.h"

#include "TPControl.h"

#include "omp.h"


#include "IP_Duct.h"
#include "IP_PPLT2D.h"
#include "IP_SHC1D.h"
#include "IP_PMT.h"
#include "IP_Rect.h"
#include "IP_Step.h"
#include "IP_Cylinder.h"
#include "IP_Polygon.h"
#include "IP_Sphere.h"

// FX10 profiler
#if defined __K_FPCOLL
#include "fjcoll.h"
#elif defined __FX_FAPP
#include "fj_tool/fjcoll.h"
#endif

// Performance Monitor
#include "PerfMonitor.h"




using namespace std;
using namespace pm_lib;


class FFV {
private:
  
  int session_maxStep;     ///< セッションのステップ数
  int session_currentStep; ///< セッションの現在のステップ
  int ModeTiming;          ///< タイミング測定管理フラグ
  
  int G_size[3];           ///< 全ドメインの分割数
  REAL_TYPE G_org[3];      ///< 全ドメインの基点
  REAL_TYPE G_reg[3];      ///< 全ドメインのサイズ
  unsigned long Acell;     ///< グローバルなActive cell
  unsigned long Fcell;     ///< グローバルなFluid cell
  unsigned long Wcell;     ///< グローバルなSolid cell
  
  // Fortranへの引数
  int sz[3];        ///< ローカルの領域分割数
  int gc;           ///< ガイドセル数
  REAL_TYPE *dh;    ///< 格子幅（無次元）
  REAL_TYPE *dh0;   ///< 格子幅（有次元）
  REAL_TYPE v00[4]; ///< 参照速度
  
  
  // データ領域ポインタ
  
  // Vector3D
  REAL_TYPE *dc_v;
  REAL_TYPE *dc_vc;
  REAL_TYPE *dc_v0;
  REAL_TYPE *dc_wv;
  REAL_TYPE *dc_abf;
  REAL_TYPE *dc_vf0;
  REAL_TYPE *dc_av;
  REAL_TYPE *dc_wvex;
  REAL_TYPE *dc_qbc;
  
  // Scalar3D
  int *dc_mid;
  int *dc_bcd;
  int *dc_bcp;
  int *dc_bcv;
  int *dc_bh1;
  int *dc_bh2;
  REAL_TYPE  *dc_ws;
  REAL_TYPE  *dc_p;
  REAL_TYPE  *dc_wk2;
  REAL_TYPE  *dc_dp;
  REAL_TYPE  *dc_p0;
  REAL_TYPE  *dc_t;
  REAL_TYPE  *dc_t0;
  REAL_TYPE  *dc_vt;
  REAL_TYPE  *dc_vof;
  REAL_TYPE  *dc_ap;
  REAL_TYPE  *dc_at;
  float      *dc_cvf;
  
  // Coarse initial
  REAL_TYPE *dc_r_v;  ///< 粗格子の速度
  REAL_TYPE *dc_r_p;  ///< 粗格子の圧力
  REAL_TYPE *dc_r_t;  ///< 粗格子の温度
  
  // コンポーネントワーク配列のアドレス管理
  REAL_TYPE** component_array;
  
  
  // カット
  REAL_TYPE  *dc_cut; ///< 距離情報
  int        *dc_bid; ///< BC
  
  
  FILE *mp;    ///< 標準出力
  FILE *fp_b;  ///< 基本情報
  FILE *fp_w;  ///< 壁面情報
  FILE *fp_c;  ///< コンポーネント情報
  FILE *fp_d;  ///< 流量収支情報
  FILE *fp_i;  ///< 反復履歴情報
  FILE *fp_f;  ///< 力の履歴情報
  
  Control C;                 ///< 制御パラメータクラス
  FileIO F;                  ///< ファイル入出力クラス
  DTcntl DT;                 ///< 時間制御クラス
  ParseMat M;                ///< 媒質パラメータ管理クラス
  Intrinsic* Ex;             ///< pointer to a base class
  ItrCtl IC[ItrCtl::ic_END]; ///< 反復情報管理クラス
  ReferenceFrame RF;         ///< 参照座標系クラス
  MediumList* mat;           ///< 媒質リスト
  CompoList* cmp;            ///< コンポーネントリスト
  PerfMonitor PM;            ///< 性能モニタクラス
  
//  SetBC3D BC;                ///< BCクラス
  
  char tm_label_ptr[tm_END][TM_LABEL_MAX];  ///< プロファイラ用のラベル
  
public:
  cpm_ParaManager *paraMngr; ///< Cartesian Partition Maneger
  
public:
  /** コンストラクタ */
  FFV();
  
  /**　デストラクタ */
  ~FFV();
  
  
private:

  /**
   * @brief Adams-Bashforth法に用いる配列のアロケーション
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_AB2 (double &total);
  
  
  /**
   * @brief 平均値処理に用いる配列のアロケーション
   * @param [in/out] total  ソルバーに使用するメモリ量
   */
  void allocArray_Average (double &total);
  
  
  /**
   * @brief 粗格子読み込みに用いる配列のアロケーション
   * @param [in] r_size  粗格子の領域サイズ
   */
  void allocArray_CoarseMesh(const int* r_size);
  
  
  /**
   * @brief コンポーネント体積率の配列のアロケーション
   * @param [in/out] prep  前処理に使用するメモリ量
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_CompoVF(double &prep, double &total);
  
  
  /**
   * @brief カット情報の配列
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_Cut(double &total);
  
  
  /**
   * @brief 熱の主計算部分に用いる配列のアロケーション
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_Heat(double &total);
  
  
  /**
   * @brief 体積率の配列のアロケーション
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_Interface(double &total);
  
  
  /**
   * @brief LES計算に用いる配列のアロケーション
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_LES(double &total);
  
  
  /**
   * @brief 主計算部分に用いる配列のアロケーション
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_Main(double &total);
  
  
  /**
   * @brief 前処理に用いる配列のアロケーション
   * @param [in/out] prep  前処理に使用するメモリ量
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_Prep(double &prep, double &total);
  
  
  /**
   * @brief Runge-Kutta法に用いる配列のアロケーション
   * @param [in/out] total ソルバーに使用するメモリ量
   */
  void allocArray_RK(double &total);
  
  
public:
  
  /**
   * @brief 組み込み例題のインスタンス
   * @param [in] Cref Controlクラスのポインタ
   */
  void connectExample(Control* Cref);
  
  
  /** 計算領域情報を設定する
   * @param [in] dom_file  ドメインファイル名
   */
  void DomainInitialize(const string dom_file);
  
  
  /**
   * @brief 固定パラメータの設定
   */
  void fixed_parameters();
  
  
  /**
   * @brief 組み込み例題の設定
   * @param [in] Cref    コントロールクラス
   * @param [in] tpCntl  テキストパーサーのラッパー
   */
  void getExample(Control* Cref, TPControl* tpCntl);
  
  
  /**
   * @brief プロファイラのラベル取り出し
   * @param [in] key 格納番号
   * @return ラベル
   */
  inline const char* get_tm_label(const int key) 
  {
    return (const char*)tm_label_ptr[key];
  }
  
  
  /** 初期化 
   * 格子生成、ビットフラグ処理ほか
   * @param[in] argc  main関数の引数の個数
   * @param[in] argv  main関数の引数リスト
   */
  int Initialize(int argc, char **argv);
  
  
  /** 
   * @brief マスターノードのみ trueを返す
   * @ret true(Rank==0) / false(Rank!=0)
   */
  bool IsMaster() 
  {
    return ( paraMngr->GetMyRankID() == 0 ) ? true : false;
  }
  
  
  /** 1ステップのコアの処理
   * @param[in] m_step   現在のステップ数
   */
  int Loop(int m_step);
  
  
  /** シミュレーションの1ステップの処理
   *  Loop() + stepPost()
   */
  int MainLoop();
  
  
  /** シミュレーションの終了時の処理
   * プロファイルの統計処理ほか
   */
  bool Post();
  
  
  /**
   * @brief 読み込んだ領域情報のデバッグライト
   * @param [in] dInfo  領域情報クラス
   */
  void printDomainInfo( cpm_GlobalDomainInfo* dInfo );
  
  
  /** ParseMatクラスをセットアップし，媒質情報を入力ファイルから読み込み，媒質リストを作成する
   * @param [in] fp  ファイルポインタ
   */
  void setMediumList(FILE* fp);
  
  
  /**
   * @brief 並列化と分割の方法を保持
   * @return 並列モード
   */
  string setParallelism();
  
  
  /**
   * @brief タイミング測定区間にラベルを与えるラッパー
   * @param [in] key       キー番号
   * @param [in] label     ラベル
   * @param [in] type      測定対象タイプ(COMM or CALC)
   * @param [in] exclusive 排他測定フラグ(ディフォルトtrue)
   */
  void set_label(const int key, char* label, PerfMonitor::Type type, bool exclusive=true);
  
  
  /**
   * @brief タイミング測定区間にラベルを与える
   */
  void set_timing_label();
  
  
  /** 毎ステップ後に行う処理 */
  bool stepPost();
  
  
  /**
   * @brief タイミング測定開始
   * @param [in] key 格納番号
   */
  inline void TIMING_start(const int key) 
  {
    // Intrinsic profiler
    TIMING__ PM.start(key);
    
    // Venus FX profiler
#if defined __K_FPCOLL
    start_collection( get_tm_label(key) );
#elif defined __FX_FAPP
    fapp_start( get_tm_label(key), 0, 0);
#endif
  }
  
  
  /**
   * @brief タイミング測定終了
   * @param [in] key             格納番号
   * @param [in] flopPerTask    「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
   * @param [in] iterationCount  実行「タスク」数 (ディフォルト1)
   */
  inline void TIMING_stop(const int key, double flopPerTask=0.0, int iterationCount=1) 
  {
    // Venus FX profiler
#if defined __K_FPCOLL
    stop_collection( get_tm_label(key) );
#elif defined __FX_FAPP
    fapp_stop( get_tm_label(key), 0, 0);
#endif
  }
    
  
  /** コマンドラインヘルプ */
  virtual void Usage();
  
  

};

#endif // _FFV_H_
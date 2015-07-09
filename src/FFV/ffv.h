#ifndef _FFV_H_
#define _FFV_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################
//
// 以下のマクロはcpm_Define.hで定義されている
//   REAL_TYPE
//   X_MINUS, Y_MINUS, Z_MINUS, X_PLUS, Y_PLUS, Z_PLUS
//   X_DIR, Y_DIR, Z_DIR
//   PLUS2MINUS, MINUS2PLUS, BOTH

/** 
 * @file   ffv.h
 * @brief  FFV Class Header
 * @author aics
 */

#include "ffv_Alloc.h"
#include <math.h>
#include <float.h>

// FB
#include "omp.h"
#include "ParseBC.h"
#include "ParseMat.h"
#include "VoxInfo.h"
#include "PolyProperty.h"
#include "History.h"
#include "Monitor.h"

// FFV
#include "ffv_Version.h"
#include "ffv_Define.h"
#include "ffv_SetBC.h"
#include "ffv_Ffunc.h"
#include "ffv_LSfunc.h"
#include "ffv_TerminateCtrl.h"
#include "ffv_LS.h"

// Geometry
#include "Geometry.h"
#include "Glyph.h"
#include "CompoFraction.h"

// FileIO class
#include "ffv_sph.h"
#include "ffv_plot3d.h"

// Intrinsic class
#include "IP_Duct.h"
#include "IP_PPLT2D.h"
#include "IP_PMT.h"
#include "IP_Rect.h"
#include "IP_Step.h"
#include "IP_Cylinder.h"
#include "IP_Sphere.h"
#include "IP_Jet.h"


// FX10 profiler
#if defined __K_FPCOLL
#include "fjcoll.h"
#elif defined __FX_FAPP
#include "fj_tool/fjcoll.h"
#endif

// Performance Monitor
#include "PerfMonitor.h"

// Polylib
#include "Polylib.h"
#include "MPIPolylib.h"

// CDMlib
#include "cdm_DFI.h"


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
using namespace pm_lib;
using namespace PolylibNS;


class FFV : public FALLOC {
public:
  int EXEC_MODE;           ///< solver と filter の識別
  
private:
  int ffv_procGrp;         ///< プロセスグループ番号 => 0
  int ModeTiming;          ///< タイミング測定管理フラグ
  
  unsigned long G_Acell;   ///< グローバルなActive cell
  unsigned long G_Fcell;   ///< グローバルなFluid cell
  unsigned long G_Wcell;   ///< グローバルなSolid cell
  
  unsigned long L_Acell;   ///< ローカルなActive cell
  unsigned long L_Fcell;   ///< ローカルなFluid cell
  unsigned long L_Wcell;   ///< ローカルなSolid cell
  
  // セッション：1回のrun
  // ケース：ある一連の計算セッション
  double CurrentTime;           ///< 計算開始からの積算時刻（ケース）
  double CurrentTimeStat;       ///< 統計値操作の積算時間（ケース）
  unsigned CurrentStep;         ///< 計算開始からの積算ステップ（ケース）
  unsigned CurrentStepStat;     ///< 統計操作の積算ステップ数（ケース）
  unsigned Session_CurrentStep; ///< セッションの現在のステップ
  unsigned Session_LastStep;    ///< セッションの終了ステップ数
  
  double face_comm_size;       ///< 全ノードについて，ローカルノード1面・一層あたりの通信量の和（要素数）
  
  REAL_TYPE deltaT; ///< 時間積分幅（無次元）
  
  int communication_mode; ///< synchronous, asynchronous
  
  REAL_TYPE v00[4];      ///< 参照速度
  REAL_TYPE range_Ut[2]; ///< 
  REAL_TYPE range_Yp[2]; ///<
  
  DivConvergence DivC; ///< 発散収束判定
  
  // 定常収束モニタ
  typedef struct
  {
    double previous;  ///< 前回の反復の収束値
    double rate;      ///< 収束値の増減比
  } ConvergenceMonitor;
  
  
  // Polygon管理用
  PolygonProperty* PG;
  
  // 周期境界の方向
  int ensPeriodic[3];

  
  double *mat_tbl;    // Fortranでの多媒質対応ルーチンのため，rho, cp, lambdaの配列
  REAL_TYPE *vec_tbl; // Fortranでの速度境界条件参照用配列
  
  REAL_TYPE *cmp_force_local;      ///< コンポーネントの各ランクの力の成分
  REAL_TYPE *cmp_force_global;     ///< コンポーネントの各ランクのCompoListで計算した力の成分を積算する
  REAL_TYPE *cmp_force_avr;        ///< コンポーネントの力の成分の平均値
  REAL_TYPE *buffer_force;     ///< 力の積算用バッファ
  int *global_obstacle;        ///< 各コンポーネント毎のOBSTACLEの有無
  int num_obstacle;            ///< OBSTACLEの個数
  
  FILE *fp_b;  ///< 基本情報
  FILE *fp_w;  ///< 壁面情報
  FILE *fp_c;  ///< コンポーネント情報
  FILE *fp_d;  ///< 流量収支情報
  FILE *fp_i;  ///< 反復履歴情報
  FILE *fp_f;  ///< 力の履歴情報
  
  Control C;                 ///< 制御パラメータクラス
  DTcntl DT;                 ///< 時間制御クラス
  ParseMat M;                ///< 媒質パラメータ管理クラス
  Intrinsic* Ex;             ///< pointer to a base class
  ReferenceFrame RF;         ///< 参照座標系クラス
  MediumList* mat;           ///< 媒質リスト
  CompoList* cmp;            ///< コンポーネントリスト
  PerfMonitor PM;            ///< 性能モニタクラス
  VoxInfo V;                 ///< ボクセル前処理クラス
  ParseBC B;                 ///< 境界条件のパースクラス
  SetBC3D BC;                ///< BCクラス
  History* H;                ///< 履歴クラス
  MPIPolylib* PL;            ///< Polylibクラス
  POLYLIB_STAT poly_stat;    ///< Polylibの戻り値
  FBUtility U;               ///< ユーティリティクラス
  MonitorList MO;            ///< Monitorクラス
  IO_BASE* F;                ///< File IO class
  Geometry GM;               ///< Geometry class
  
  LinearSolver LS[ic_END];   ///< 反復解法
  
  ConvergenceMonitor CM_F;   ///< 流動の定常収束モニター
  ConvergenceMonitor CM_H;   ///< 熱の定常収束モニター
  
  
  char tm_label_ptr[PM_NUM_MAX][TM_LABEL_MAX];  ///< プロファイラ用のラベル
  int order_of_PM_key;      ///< PMlib用の登録番号カウンタ < PM_NUM_MAX
  
  string active_fname;      ///< Active subdomainのファイル名


  
public:
  /** コンストラクタ */
  FFV();
  
  /**　デストラクタ */
  ~FFV();
  
  
private:
  
  /**
   * @brief タイミング測定開始
   * @param [in] key ラベル
   */
  inline void TIMING_start(const string key)
  {
    // PMlib Intrinsic profiler
    TIMING__ PM.start(key);
    
    const char* s_label = key.c_str();
    
    // Venus FX profiler
#if defined __K_FPCOLL
    start_collection( s_label );
#elif defined __FX_FAPP
    fapp_start( s_label, 0, 0);
#endif
  }
  
  
  /**
   * @brief タイミング測定終了
   * @param [in] key             ラベル
   * @param [in] flopPerTask    「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
   * @param [in] iterationCount  実行「タスク」数 (ディフォルト1)
   */
  inline void TIMING_stop(const string key, double flopPerTask=0.0, int iterationCount=1)
  {
    // Venus FX profiler
    const char* s_label = key.c_str();
    
#if defined __K_FPCOLL
    stop_collection( s_label );
#elif defined __FX_FAPP
    fapp_stop( s_label, 0, 0);
#endif
    
    // PMlib Intrinsic profiler
    TIMING__ PM.stop(key, flopPerTask, (unsigned)iterationCount);
  }
  
  
  
  
  /** ffv_Initialize.C *******************************************************/
  
  // 主計算に用いる配列の確保
  void allocate_Main(double &total);
  
  
  // 全Voxelモデルの媒質数とKOSの整合性をチェック
  bool chkMediumConsistency();
  
  
  // mat[], cmp[]のテーブルを作成
  void createTable(FILE* fp);
  
  
  // CompoListの情報を表示する
  void displayCompoInfo(const int* cgb, FILE* fp);
  
  
  // 交点情報の表示（デバッグ）
  void displayCutInfo(const long long* cut, const int* bid);
  
  
  // メモリ使用量の表示
  void displayMemoryInfo(FILE* fp, double G_mem, double L_mem, const char* str);
  
  
  // 制御パラメータ，物理パラメータの表示
  void displayParameters(FILE* fp);
  
  
  // 境界条件のドライラン
  void dryrunBC();
  
  
  // BCIndexにビット情報をエンコードする
  void encodeBCindex(FILE* fp);
  
  
  // 固定パラメータの設定
  void fixedParameters();
  
  
  // 並列処理時の各ノードの分割数を集めてファイルに保存する
  void gatherDomainInfo();
  
  
  // Glyphを生成・出力
  void generateGlyph(const long long* cut, const int* bid, FILE* fp, int* m_st=NULL, int* m_ed=NULL);
  
  
  // Div反復のパラメータ
  bool getParaDiv(TextParser* tpCntl);
  
  
  // Intrinsic Classの同定
  void identifyExample(FILE* fp);
  
  
  // FielIO classの同定
  void identifyFIO(TextParser* tpCntl);
  
  
  // 線形ソルバを特定
  void identifyLinearSolver(TextParser* tpCntl);
  
  
  //距離情報の初期化
  void initCutInfo();
  
  
  // インターバルの初期化
  void initInterval();
  
  
  
  // 線形ソルバークラス関連の初期化
  void LS_initialize(double& TotalMemory, TextParser* tpCntl);
  
  
  // 線形ソルバを特定し，パラメータをセットする
  void LS_setParameter(TextParser* tpCntl, const int odr, const string label);
  
  
  /**
   * @brief 初期擾乱
   * @note x方向の流れが主方向で、それに垂直なdir方向のチャネルを想定
   */
  void perturbation();
  
  
  // 履歴の出力準備
  void prepHistoryOutput();
  
  
  // 収束判定条件の表示
  void printCriteria(FILE* fp);
  
  
  // 全領域情報
  void printGlobalDomain(FILE* fp);
  
  
  // 外部境界条件を読み込み，Controlクラスに保持する
  void setBCinfo();
  
  
  // 回転体をセット
  void setComponentSR();
  
  
  // HEX,FANコンポーネントなどの体積率とbboxなどをセット
  void setComponentVF();
  
  
  // コンポーネントのローカルなBbox情報からグローバルなBbox情報を求め，CompoListの情報を表示
  void dispGlobalCompoInfo(FILE* fp);
  
  
  // 線形ソルバー種別のパラメータ表示
  void printIteratoinParameter(FILE* fp, LinearSolver* IC);
  
  
  // 初期条件の設定
  void setInitialCondition();
  
  
  // ParseMatクラスをセットアップし，媒質情報を入力ファイルから読み込み，媒質リストを作成する
  void setMediumList(FILE* fp);
  
  
  // MonitorListのセットアップ
  void setMonitorList();
  
  
  // 並列化と分割の方法を保持
  string setParallelism();
  
  
  // 時間積分幅や物理パラメータの設定
  void setParameters();
  
  
  // パラメータのロードと計算領域を初期化し，並列モードを返す
  string setDomain(TextParser* tpf);
  
  
  // グローバルな領域情報を取得し、無次元の領域基本パラメータを返す
  int SD_getParameter(TextParser* tp_dom);
  
  
  // 計算領域情報を設定する
  void SD_Initialize(const int div_type, TextParser* tp_dom);
  
  
  // 各種例題のモデルをセットアップ
  void setModel(double& PrepMemory, double& TotalMemory, FILE* fp);
  
  
  // 幾何形状情報を準備し，交点計算を行う
  void SM_Polygon2Cut(double& m_prep, double& m_total, FILE* fp);
  
  

  
  
  /** ffv.C *******************************************************/
  
  // 時間平均操作を行う
  void Averaging(double& flop);
  
  
  // OBSTACLEコンポーネントの力の成分を計算し、集める
  void calcForce(double& flop);
  
  
  // 外部計算領域の各面における総流量と対流流出速度を計算する
  void DomainMonitor(BoundaryOuter* ptr, Control* R);
  
  
  // div(u)を計算する
  void NormDiv(REAL_TYPE* div, const REAL_TYPE dt);
  
  
  // タイミング測定区間にラベルを与えるラッパー
  void set_label(const string label, PerfMonitor::Type type, bool exclusive=true);
  
  
  // 毎ステップ後に行う処理
  bool stepPost();
  
  
  // タイミング測定区間にラベルを与える
  void set_timing_label();
  
  
  // 利用例の表示
  void Usage();
  
  
  // 空間平均操作と変動量の計算を行う
  // スカラ値は算術平均，ベクトル値は自乗和
  void VariationSpace(double* rms, double* avr, double& flop);
  
  
  
  
  /** ffv_Heat.C *******************************************************/
  
  
  /**
   * @brief 単媒質に対する熱伝導方程式を陰解法で解く
   * @param [in]  IC       LinearSolverクラス
   * @param [in]  rhs_nrm  Poisson定数項ベクトルの自乗和ノルム
   * @param [in]  r0       初期残差ベクトル
   */
  void ps_LS(LinearSolver* IC, const double rhs_nrm, const double r0);

  
  
  /**
   * @brief 単媒質に対する熱伝導方程式をEuler陰解法で解く
   * @retval ローカルノードの変化量の自乗和
   * @param [in,out] t    n+1時刻の温度場
   * @param [out]    b2   ソースベクトルの自乗和
   * @param [in]     dt   時間積分幅
   * @param [in]     qbc  境界条件熱流束
   * @param [in]     bh   BCindex B
   * @param [in]     ws   部分段階の温度
   * @param [in]     IC   IterationCtlクラス
   * @param [in,out] flop 浮動小数点演算数
   */
  double ps_Diff_SM_PSOR(REAL_TYPE* t, double& b2, const REAL_TYPE dt, const REAL_TYPE* qbc, const int* bh, const REAL_TYPE* ws, IterationCtl* IC, double& flop);
  
  
  
  
  
  
  /** 1ステップのコアの処理
   * @param [in] m_step   現在のステップ数
   */
  int Loop(const unsigned m_step);
  
  
  // @brief Fractional Step法でNavier-Stokes方程式を解く．バイナリ近似．
  void NS_FS_E_Binary();
  
  
  // @brief Fractional Step法でNavier-Stokes方程式を解く．距離場近似．
  void NS_FS_E_CDS();
  
  
  // 温度の移流拡散方程式をEuler陽解法/Adams-Bashforth法で解く
  void PS_Binary();

  
  
  /**
   * @brief VOF値を気体(0.0)と液体(1.0)で初期化
   */
  void setVOF();
  


  

  
public:
  
  // フィルタ処理初期化
  int FilterInitialize(int argc, char **argv);
  
  
  // フィルタ処理ループ
  int FilterLoop();
  
  int FilterLoop(unsigned int step);
  
  
  // SPHファイルのデータを取り出し、配列に格納する
  int FilterGetArrayFromSph( cdm_DFI *dfi, cdm_Rank *rank, int step, REAL_TYPE **pArray, int *nArray);
  
  
  /**
   * @brief CPMのポインタをコピーし、ランク情報を設定
   * @param [in] m_paraMngr  cpm_ParaManagerクラス
   * @return  エラーコード
   */
  bool importCPM(cpm_ParaManager* m_paraMngr)
  {
    if ( !m_paraMngr ) return false;
    paraMngr = m_paraMngr;
    
    setRankInfo(paraMngr, procGrp);
    
    return true;
  }
  
  
  // 初期化格子生成、ビットフラグ処理ほか
  int Initialize(int argc, char **argv);
  
  
  /** 
   * @brief マスターノードのみ trueを返す
   * @return true(Rank==0) / false(Rank!=0)
   */
  bool IsMaster() const
  {
    return ( paraMngr->GetMyRankID() == 0 ) ? true : false;
  }
  
  
  // シミュレーションの1ステップの処理
  int MainLoop();
  
  
  /** 
   * @brief シミュレーションの終了時の処理
   * @note プロファイルの統計処理ほか
   */
  bool Post();
  
};

#endif // _FFV_H_

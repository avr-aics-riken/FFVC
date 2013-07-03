#ifndef _FFV_H_
#define _FFV_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
//
// 以下のマクロはcpm_Define.hで定義されている
//   REAL_TYPE
//   X_MINUS, Y_MINUS, Z_MINUS, X_PLUS, Y_PLUS, Z_PLUS
//   X_DIR, Y_DIR, Z_DIR
//   PLUS2MINUS, MINUS2PLUS, BOTH

/** 
 * @file   ffv.h
 * @brief  FFV Class Header
 * @author kero
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>

#include "cpm_ParaManager.h"

#include "DomainInfo.h"

#include "FB_Define.h"
#include "ffv_Define.h"
#include "mydebug.h"
#include "FBUtility.h"
#include "Control.h"
#include "Alloc.h"
#include "FileIO.h"
#include "ParseBC.h"
#include "ParseMat.h"
#include "VoxInfo.h"
#include "TPControl.h"
#include "ffv_SetBC.h"
#include "CompoFraction.h"

/* 20130606 commentout 
#include "dfi.h"
#include "dfiinfo.h"
 */

#include "History.h"
#include "Monitor.h"
#include "ffv_Ffunc.h"
#include "ffv_TerminateCtrl.h"

#include "omp.h"

// Intrinsic class
#include "IP_Duct.h"
#include "IP_PPLT2D.h"
#include "IP_SHC1D.h"
#include "IP_PMT.h"
#include "IP_Rect.h"
#include "IP_Step.h"
#include "IP_Cylinder.h"
#include "IP_Polygon.h"
#include "IP_Sphere.h"
#include "IP_Jet.h"

/* 20130606 commentout  PLOT3D
#include "ffv_PLOT3D.h"
 */

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

// Cutlib
#include "Cutlib.h"

// CIOlib
#include "cio_DFI.h"


using namespace std;
using namespace pm_lib;
using namespace PolylibNS;
using namespace cutlib;



class FFV : public DomainInfo {
private:
  int ffv_procGrp;         ///< プロセスグループ番号 => 0
  int ModeTiming;          ///< タイミング測定管理フラグ
  
  unsigned long G_Acell;   ///< グローバルなActive cell
  unsigned long G_Fcell;   ///< グローバルなFluid cell
  unsigned long G_Wcell;   ///< グローバルなSolid cell
  
  unsigned long L_Acell;   ///< ローカルなActive cell
  unsigned long L_Fcell;   ///< ローカルなFluid cell
  unsigned long L_Wcell;   ///< ローカルなSolid cell
  
  
  double CurrentTime;           ///< 計算開始からの積算時刻（ケース）
  double CurrentTime_Avr;       ///< 平均値操作の積算時間（ケース）
  double Session_StartTime;     ///< セッションの開始時間
  double Session_CurrentTime;   ///< セッションの現在時間
  
  double step_start;            ///< 1stepのelapse time(sec)
  
  unsigned CurrentStep;         ///< 計算開始からの積算ステップ（ケース）
  unsigned CurrentStep_Avr;     ///< 平均操作の積算ステップ数（ケース）
  unsigned Session_StartStep;   ///< セッションの開始ステップ
  unsigned Session_CurrentStep; ///< セッションの現在のステップ
  unsigned Session_LastStep;    ///< セッションで計算するステップ数
  
  REAL_TYPE convergence_prev;  ///< 前回の反復の収束値
  REAL_TYPE convergence_rate;  ///< 収束値の増減比
  
  REAL_TYPE deltaT; ///< 時間積分幅（無次元）
  
  int communication_mode; ///< synchronous, asynchronous
  
  int cf_sz[3];     ///< SOR2SMAの反復の場合のバッファサイズ
  REAL_TYPE *cf_x;  ///< i方向のバッファ
  REAL_TYPE *cf_y;  ///< j方向のバッファ
  REAL_TYPE *cf_z;  ///< k方向のバッファ

  /* 20130611 commentout 
  // @attention CIOlib導入により，ポリシーの確認
  // dfi ファイル管理用 -> Kind_of_vars in FB_Define.h
  // 同じ解像度のリスタート時には、既にdfiファイルが存在する場合には、その内容を継続する
  // ラフリスタートの場合には、新規dfiファイルを生成する >> dfi.C
  //int dfi_mng[var_END];
  */
  
  REAL_TYPE v00[4];      ///< 参照速度
  REAL_TYPE range_Ut[2]; ///< 
  REAL_TYPE range_Yp[2]; ///< 
  
  
  // データ領域ポインタ
  
  // Vector3D
  REAL_TYPE *d_v;   ///< セルセンター速度
  REAL_TYPE *d_vf;  ///< セルフェイス速度
  REAL_TYPE *d_vc;  ///< セルセンター疑似速度
  REAL_TYPE *d_v0;  ///< n-stepの速度保持
  REAL_TYPE *d_wv;  ///< ワーク配列
  REAL_TYPE *d_abf; ///< Adams-bashforthワーク
  REAL_TYPE *d_av;  ///< 平均値
  REAL_TYPE *d_wo;  ///< 入出力のバッファワーク
  REAL_TYPE *d_qbc; ///< 熱BC flux保持
  
  // Scalar3D
  int *d_mid;
  int *d_bcd;
  int *d_bcp;
  int *d_bcv;
  int *d_bh1;
  int *d_bh2;
  
  REAL_TYPE *d_p;   ///< 圧力
  REAL_TYPE *d_p0;  ///< 圧力（1ステップ前）
  REAL_TYPE *d_ws;  ///< 反復中に固定のソース
  REAL_TYPE *d_sq;  ///< 反復中に変化するソース
  REAL_TYPE *d_dv;  ///< div(u)の保存
  REAL_TYPE *d_b;   ///< Ax=bの右辺ベクトル
  REAL_TYPE *d_t;   ///< 温度
  REAL_TYPE *d_t0;  ///< 温度（1ステップ前）
  REAL_TYPE *d_vt;
  REAL_TYPE *d_vof;
  REAL_TYPE *d_ap;  ///< 圧力（時間平均値）
  REAL_TYPE *d_at;  ///< 温度（時間平均値）
  float *d_cvf;     ///< 体積率
  
  // Coarse initial
  REAL_TYPE *d_r_v;  ///< 粗格子の速度
  REAL_TYPE *d_r_p;  ///< 粗格子の圧力
  REAL_TYPE *d_r_t;  ///< 粗格子の温度
  
  // GMRES
  REAL_TYPE * d_wg;   ///< テンポラリの配列 [size] 
  REAL_TYPE * d_res;  ///< 残差 = b - Ax
  REAL_TYPE * d_vm;   ///< Kryolov subspaceの直交基底 [size*FREQ_OF_RESTART]
  REAL_TYPE * d_zm;   ///< Right-hand side vector for the residual minimization problem [size*FREQ_OF_RESTART]
  
#define FREQ_OF_RESTART 15 // リスタート周期
  
  // PCG & PBiCGSTAB
	REAL_TYPE *d_pcg_r;
	REAL_TYPE *d_pcg_p;
  
	// PCG
	REAL_TYPE *d_pcg_q;
	REAL_TYPE *d_pcg_z;
  
	// PBiCGSTAB
	REAL_TYPE *d_pcg_r0;
	REAL_TYPE *d_pcg_p_;
	REAL_TYPE *d_pcg_q_;
	REAL_TYPE *d_pcg_s;
	REAL_TYPE *d_pcg_s_;
	REAL_TYPE *d_pcg_t_;
  
  REAL_TYPE** component_array; ///< コンポーネントワーク配列のアドレス管理
  
  int* compo_global_bbox; ///< グローバルなコンポーネントBbox 表示に利用
  
  // PolygonGroupの管理
  Control::Polygon_property* poly_prop;
  
  // カット
  CutPos32Array *cutPos;
  CutBid5Array  *cutBid;
  
  float  *d_cut; ///< 距離情報
  int    *d_bid; ///< BC
  
  
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
  VoxInfo V;                 ///< ボクセル前処理クラス
  ParseBC B;                 ///< 境界条件のパースクラス
  SetBC3D BC;                ///< BCクラス
  History* H;                ///< 履歴クラス
  MPIPolylib* PL;            ///< Polylibクラス
  POLYLIB_STAT poly_stat;    ///< Polylibの戻り値
  FBUtility U;               ///< ユーティリティクラス
  MonitorList MO;            ///< Monitorクラス
  
  // 20130611 ::DFI DFI;                 ///< 分散ファイルインデクス管理クラス
  // 20130611 Plot3D PLT3D;              ///< PLOT3Dクラス
  // 20130611 FileIO_PLOT3D_READ  FP3DR; ///< PLOT3D READクラス
  // 20130611 FileIO_PLOT3D_WRITE FP3DW; ///< PLOT3D WRITEクラス
  
  char tm_label_ptr[tm_END][TM_LABEL_MAX];  ///< プロファイラ用のラベル

  
  // pointers to CIO class
  
  // InFile
  cio_DFI *DFI_IN_PRS;      ///< Pressure
  cio_DFI *DFI_IN_VEL;      ///< Velocity
  cio_DFI *DFI_IN_FVEL;     ///< Face velocity
  cio_DFI *DFI_IN_TEMP;     ///< Temperature
  cio_DFI *DFI_IN_PRSA;     ///< Averaged pressure
  cio_DFI *DFI_IN_VELA;     ///< Averaged velocity
  cio_DFI *DFI_IN_TEMPA;    ///< Averaged temperature
  
  // OutFile
  cio_DFI *DFI_OUT_PRS;     ///< Pressure
  cio_DFI *DFI_OUT_VEL;     ///< Velocity
  cio_DFI *DFI_OUT_FVEL;    ///< Face velocity
  cio_DFI *DFI_OUT_TEMP;    ///< Temperature
  cio_DFI *DFI_OUT_PRSA;    ///< Averaged Pressure
  cio_DFI *DFI_OUT_VELA;    ///< Averaged velocity
  cio_DFI *DFI_OUT_TEMPA;   ///< Averaged temperature
  cio_DFI *DFI_OUT_TP;      ///< Total Pressure
  cio_DFI *DFI_OUT_VRT;     ///< Vorticity
  cio_DFI *DFI_OUT_I2VGT;   ///< 2nd Invariant of Velocity Gradient Tensor
  cio_DFI *DFI_OUT_HLT;     ///< Helicity
  cio_DFI *DFI_OUT_DIV;     ///< Divergence for debug
  
public:
  /** コンストラクタ */
  FFV();
  
  /**　デストラクタ */
  ~FFV();
  
  
private:
  
  /** ffv_Alloc.C *******************************************************/
  
  // Adams-Bashforth法に用いる配列のアロケーション
  void allocArray_AB2 (double &total);
  
  
  // 平均値処理に用いる配列のアロケーション
  void allocArray_Average (double &total);
  
  
  // 粗格子読み込みに用いる配列のアロケーション
  void allocArray_CoarseMesh(const int* r_size, double &prep);
  
  
  // コンポーネント体積率の配列のアロケーション
  void allocArray_CompoVF(double &prep, double &total);
  
  
  // カット情報の配列
  void allocArray_Cut(double &total);
  
  
  // コンポーネントのワーク用配列のアロケート
  void allocArray_Forcing(double& m_prep, double& m_total, FILE* fp);
  
  
  // 熱の主計算部分に用いる配列のアロケーション
  void allocArray_Heat(double &total);
  
  
  // 体積率の配列のアロケーション
  void allocArray_Interface(double &total);
  
  
  // Krylov-subspace法に用いる配列のアロケーション
  void allocArray_Krylov(double &total);
  
  
  // LES計算に用いる配列のアロケーション
  void allocArray_LES(double &total);
  
  
  // 主計算部分に用いる配列のアロケーション
  void allocArray_Main(double &total);
  
  
  // PCG法に用いる配列のアロケーション
  void allocArray_PCG(double &total);
  
  
  // PBiCGSTAB法に用いる配列のアロケーション
  void allocArray_PBiCGSTAB(double &total);
  
  
  // 前処理に用いる配列のアロケーション
  void allocArray_Prep(double &prep, double &total);
  
  
  // SOR2SMAのバッファ確保
  void allocate_SOR2SMA_buffer(double &total);
  
  
  
  
  
  /** ffv_Initialize.C *******************************************************/
  
  // 主計算に用いる配列の確保
  void allocate_Main(double &total);
  
  
  // ポリゴンのカット情報からVBCのboxをセット
  void Bbox_IBC();
  
  
  // 全Voxelモデルの媒質数とKOSの整合性をチェック
  bool chkMediumConsistency();
  
  
  // 時刻をRFクラスからv00[4]にコピーする
  void copyV00fromRF(double m_time);
  
  
  // コンポーネントの内容リストを表示する
  void display_Compo_Info(FILE* fp);
  
  
  // CompoListの内容とセル数の情報を表示する
  void display_CompoList(FILE* fp);
  
  
  // 制御パラメータ，物理パラメータの表示
  void display_Parameters(FILE* fp);
  
  
  // 計算領域情報を設定する
  void DomainInitialize(TPControl* tp_dom);
  
  
  //初期インデクスの情報を元に，一層拡大したインデクス値を返す
  void EnlargeIndex(int& m_st, int& m_ed, const int st_i, const int len, const int m_x, const int dir, const int m_id);
  
  
  // ポリゴンの場合のフィル操作
  void fill(FILE* fp);
  
  
  // 固定パラメータの設定
  void fixed_parameters();
  
  
  // メモリ使用量の表示
  void display_memory_info(FILE* fp, double G_mem, double L_mem, const char* str);
  
  
  // 並列処理時の各ノードの分割数を集めてファイルに保存する
  void gather_DomainInfo();
  
  
  // Binary voxelをカット情報から生成
  void generate_Solid(FILE* fp);
  
  
  // コンポーネントの面積を計算する
  void get_Compo_Area();
  
  
  // グローバルな領域情報を取得
  int get_DomainInfo(TPControl* tp_dom);
  
  
  // 出力ファイルの初期化
  void init_FileOut();
  
  
  // インターバルの初期化
  void init_Interval();
  
  
  // 距離の最小値を求める
  void min_distance(float* cut, FILE* fp);
  
  
  // 履歴の出力準備
  void prep_HistoryOutput();
  
  
  // 読み込んだ領域情報のデバッグライト
  void printDomainInfo();
  
  
  // コンポーネントリストに登録されたセル要素BCのBV情報をリサイズする
  void resizeBVface(const int* st, const int* ed, const int n, const int* bx);
  
  
  // コンポーネントリストに登録されたセル要素BCのBV情報をリサイズする
  void resizeBVcell(const int* st, const int* ed, const int n, const int* bx);
  
  
  // コンポーネントリストに登録されたBV情報をリサイズする
  void resizeCompoBV(const int kos, const bool isHeat);
  
  
  // リスタートモードの選択
  void selectRestartMode();
  
  
  // 外部境界条件を読み込み，Controlクラスに保持する
  void setBCinfo();
  
  
  // HEX,FANコンポーネントなどの体積率とbboxなどをセット
  void setComponentVF();
  
  
  // コンポーネントが存在するかを保持しておく
  void setEnsComponent();
  
  
  // コンポーネントのローカルなBbox情報からグローバルなBbox情報を求める
  void setGlobalCmpIdx();
  
  
  // 初期条件の設定
  void setInitialCondition();
  
  
  // midの情報から各BCコンポーネントのローカルなインデクスを取得する
  void setLocalCmpIdx_Binary();
  
  
  // ParseMatクラスをセットアップし，媒質情報を入力ファイルから読み込み，媒質リストを作成する
  void setMediumList(FILE* fp);
  
  
  // 各種例題のモデルをセットアップ
  void setModel(double& PrepMemory, double& TotalMemory, FILE* fp);
  
  
  // 並列化と分割の方法を保持
  string setParallelism();
  
  
  // 時間積分幅や物理パラメータの設定
  void setParameters();
  
  
  // IP用にカット領域をアロケートする
  void setup_CutInfo4IP(double& m_prep, double& m_total, FILE* fp);
  
  
  // パラメータのロードと計算領域を初期化し，並列モードを返す
  string setupDomain(TPControl* tpf, FILE* fp);
  
  
  // 幾何形状情報を準備し，交点計算を行う
  void setup_Polygon2CutInfo(double& m_prep, double& m_total, FILE* fp);
  
  
  // BCIndexにビット情報をエンコードする
  void VoxEncode();
  
  
  // ボクセルをスキャンし情報を表示する
  void VoxScan(FILE* fp);
  
  
  
  
  /** ffv.C *******************************************************/
  
  // 時間平均値のファイル出力
  void AverageOutput(double& flop);
  
  
  // 時間平均操作を行う
  void Averaging_Time(double& flop);
  
  
  
  /**
   * @brief Boussinesq浮力項の計算
   * @param [out]    v    速度
   * @param [in]     dgr  係数
   * @param [in]     t    温度
   * @param [in]     bd   BCindex ID
   * @param [in,out] flop 浮動小数点演算数
   */
  void Buoyancy(REAL_TYPE* v, const REAL_TYPE dgr, const REAL_TYPE* t, const int* bd, double& flop);
  
  
  /**
   * @brief SOR2SMAの非同期通信処理
   * @param [in,out] d_x  解ベクトル
   * @param [in]     col  オーダリングカラーの番号
   * @param [in]     ip   オーダリングカラー0の最初のインデクス
   * @param [out]    key  送信ID
   */
  void comm_SOR2SMA(REAL_TYPE* d_x, const int col, const int ip, MPI_Request* key);
  
  
  
  // 全ノードについて，ローカルノード1面・一層あたりの通信量の和を返す
  double count_comm_size(const int sz[3], const int guide);
  
  
  // 外部計算領域の各面における総流量と対流流出速度を計算する
  void DomainMonitor(BoundaryOuter* ptr, Control* R);
  
  
  //ファイル出力
  void FileOutput(double& flop, const bool crs_restart=false);
  
  
  /**
   * @brief 種類Lの線形ソルバを利用する場合，trueを返す
   * @param [in] L 線形ソルバの種類
   */
  bool hasLinearSolver(const int L);
  
  
  /**
   * @brief プロファイラのラベル取り出し
   * @param [in] key 格納番号
   * @return ラベル
   */
  inline const char* get_tm_label(const int key) 
  {
    return (const char*)tm_label_ptr[key];
  }
  
  
  /**
   * @brief FGMRES
   * @param [in]     IC      ItrCtlクラス
   * @param [in]     rhs_nrm RHS vectorのL2ノルム
   * @param [in]     r0      初期残差ベクトル
   */
  void Fgmres(ItrCtl* IC, const double rhs_nrm, const double r0);
  

  
  /**
   * @brief  FRBGS
   * @retval 反復数
   * @param [in]     IC      IterationCtlクラス
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]     rhs_nrm RHS vector
   * @param [in]     r0      初期残差ベクトル
   */
  int Frbgs(ItrCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0);
  
  /**
   * @brief  FPCG
   * @retval 反復数
   * @param [in]     IC      IterationCtlクラス
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]     rhs_nrm RHS vector
   * @param [in]     r0      初期残差ベクトル
   */
  int Fpcg(ItrCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0);
  
  /**
   * @brief  FPBiCGSTAB
   * @retval 反復数
   * @param [in]     IC      IterationCtlクラス
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]     rhs_nrm RHS vector
   * @param [in]     r0      初期残差ベクトル
   */
  int Fpbicgstab(ItrCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0);
  
  /**
   * @brief  Fcheck
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]       RHS  vector
   */
	bool Fcheck(ItrCtl* IC, REAL_TYPE res, const double rhs_nrm, const double r0);
  
  /**
   * @brief  Fpreconditioner
   * @retval 反復数
   * @param [in]     IC      IterationCtlクラス
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   */
	int Fpreconditioner(ItrCtl* IC, REAL_TYPE* x, REAL_TYPE* b);
  
  /**
   * @brief  Fsmoother
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]       RHS  vector
   */
	void Fsmoother(REAL_TYPE* x, REAL_TYPE* b, REAL_TYPE omg);
  
  /**
   * @brief  Fdot
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]       RHS  vector
   */
	void Fdot(REAL_TYPE* xy, REAL_TYPE* x, REAL_TYPE* y);
  
  
  /** 1ステップのコアの処理
   * @param [in] m_step   現在のステップ数
   */
  int Loop(const unsigned m_step);
  
  
  
  /**
   * @brief 線形ソルバーの選択実行
   * @param [in]  IC       ItrCtlクラス
   * @param [in]  rhs_nrm  Poisson定数項ベクトルの自乗和ノルム
   * @param [in]  res_init 初期残差ベクトル
   */
  void LS_Binary(ItrCtl* IC, const double rhs_nrm, const double rhs_init);
  
  

  
  
  // V-P反復のdiv(u)ノルムを計算する
  FB::Vec3i Norm_Div(ItrCtl* IC);
  
  
  /**
   * @brief Fractional Step法でNavier-Stokes方程式を解く．バイナリ近似．
   */
  void NS_FS_E_Binary();
  
  /**
   * @brief Fractional Step法でNavier-Stokes方程式を解く．距離場近似．
   */
  void NS_FS_E_CDS();
  
  
  /**
   * @brief 圧力の引き戻し操作を行う
   */
  void Pressure_Shift();


  
  /* 温度の移流拡散方程式をEuler陽解法/Adams-Bashforth法で解く
   */
  void PS_Binary();
  
  
  /**
   * @brief 移流項のEuler陽解法による時間積分
   * @param [in,out] tc      対流項の流束の和/部分段階の温度
   * @param [in]     delta_t 時間積分幅
   * @param [in]     bd      BCindex ID
   * @param [in]     t0      nステップの温度
   * @param [in,out] flop    浮動小数演算数
   * @note tc = t0 + dt/dh*sum_flux(n)
   */
  void ps_ConvectionEE(REAL_TYPE* tc, const REAL_TYPE delta_t, const int* bd, const REAL_TYPE* t0, double& flop);
  
  
  
  /**
   * @brief 単媒質に対する熱伝導方程式をEuler陽解法で解く
   * @retval 拡散項の変化量Δθの絶対値
   * @param [in,out] t    n+1時刻の温度場
   * @param [in]     dt   時間積分幅
   * @param [in]     qbc  境界条件熱流束
   * @param [in]     bh2  BCindex H2
   * @param [in]     ws   部分段階の温度
   * @param [in,out] flop 浮動小数点演算数
   */
  REAL_TYPE ps_Diff_SM_EE(REAL_TYPE* t, const REAL_TYPE dt, const REAL_TYPE* qbc, const int* bh2, const REAL_TYPE* ws, double& flop);
  
  
  /**
   * @brief 単媒質に対する熱伝導方程式をEuler陰解法で解く
   * @retval ローカルノードの変化量の自乗和
   * @param [in,out] t    n+1時刻の温度場
   * @param [out]    b2   ソースベクトルの自乗和
   * @param [in]     dt   時間積分幅
   * @param [in]     qbc  境界条件熱流束
   * @param [in]     bh2  BCindex H2
   * @param [in]     ws   部分段階の温度
   * @param [in]     IC   ItrCtlクラス
   * @param [in,out] flop 浮動小数点演算数
   */
  double ps_Diff_SM_PSOR(REAL_TYPE* t, double& b2, const REAL_TYPE dt, const REAL_TYPE* qbc, const int* bh2, const REAL_TYPE* ws, ItrCtl* IC, double& flop);
  
  
  /** SOR法
   * @retval 反復数
   * @param [in]     IC      IterationCtlクラス
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]     rhs_nrm RHS vector
   * @param [in]     r0      初期残差ベクトル
   */
  int Point_SOR(ItrCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0);
  
  /**
   * @brief 単媒質に対する熱伝導方程式を陰解法で解く
   * @param [in]  IC       IterationCtlクラス
   * @param [in]  rhs_nrm  Poisson定数項ベクトルの自乗和ノルム
   * @param [in]  r0       初期残差ベクトル
   */
  void ps_LS(ItrCtl* IC, const double rhs_nrm, const double r0);
  
  
  // リスタートプロセス
  void Restart(FILE* fp);
  
  
  // リスタート時の瞬時値ファイル読み込み
  void Restart_instantaneous(FILE* fp, double& flop);
  
  
  // リスタート時の平均値ファイル読み込み
  void Restart_avrerage(FILE* fp, double& flop);
  
  
  // リスタートの最大値と最小値の表示
  void Restart_display_minmax(FILE* fp, double& flop);

  
  
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
  
  
  /**
   * @brief VOF値を気体(0.0)と液体(1.0)で初期化
   */
  void setVOF();
  
  
  /** 毎ステップ後に行う処理 */
  bool stepPost();
  
  
  /** 2色オーダリングSORのストライドメモリアクセス版
   * @retval 反復数
   * @param [in]     IC      IterationCtlクラス
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]     rhs_nrm RHS vector
   * @param [in]     r0      初期残差ベクトル
   */
  int SOR_2_SMA(ItrCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0);
  
  
  /**
   * @brief 反復の同期処理
   * @param [in]     IC        ItrCtlクラス
   * @param [in,out] d_class   対象データ
   * @param [in]     num_layer 通信の袖数
   */
  void Sync_Scalar(ItrCtl* IC, REAL_TYPE* d_class, const int num_layer);
  
  
  /**
   * @brief タイミング測定開始
   * @param [in] key 格納番号
   */
  inline void TIMING_start(const int key) 
  {
    // Intrinsic profiler
    TIMING__ PM.start((unsigned)key);
    
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
    
    // Intrinsic profiler
    TIMING__ PM.stop((unsigned)key, flopPerTask, (unsigned)iterationCount);
  }
    
  
  /**
   * @brief コマンドラインヘルプ 
   */
  void Usage();
  
  
  /**
   * @brief 空間平均操作と変動量の計算を行う
   * @param [out]    avr  平均値
   * @param [out]    rms  変動値
   * @param [in,out] flop 浮動小数演算数
   */
  void Variation_Space(double* avr, double* rms, double& flop);
  
  
  /**
   * @brief SOR2SMAの非同期通信処理
   * @param [in,out] d_x 同期する変数
   * @param [in]     col オーダリングカラーの番号
   * @param [in]     ip  オーダリングカラー0の最初のインデクス
   * @param [in,out] key 送信ID
   */
  void wait_SOR2SMA(REAL_TYPE* d_x, const int col, const int ip, MPI_Request* key);
  

  
public:
  
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
  
  
  
  
//####################################################################
/* CIOlib導入で利用しなくなったメソッド
//####################################################################
 
  // @brief 並列分散時のファイル名の管理を行う
  void setDFI();
  
   // @brief ファイルのオープンチェック
   bool checkFile(string fname);
  
  
   // @brief 粗い格子を用いたリスタート
   void Restart_coarse(FILE* fp, double& flop);
   
   
   // @brief リスタート時の瞬時値ファイル読み込み（並列数が異なる場合）
   void Restart_different(FILE* fp, double& flop);
   
   
   // @brief ファイル読み込み＋オーバーラップを移しこみ（圧力）
   void ReadOverlap_Pressure(FILE* fp,
   double& flop,
   DifferentRestartInfo* DRI,
   DfiInfo* DI,
   REAL_TYPE* d_wk,
   int* rank_list,
   int recv_rank,
   int* assign,
   int nassign);
   
   
   // @brief ファイル読み込み＋オーバーラップを移しこみ（流速）
   void ReadOverlap_Velocity(FILE* fp,
   double& flop,
   DifferentRestartInfo* DRI,
   DfiInfo* DI,
   REAL_TYPE* d_wk,
   int* rank_list,
   int recv_rank,
   int* assign,
   int nassign);
   
   
   // @brief ファイル読み込み＋オーバーラップを移しこみ（境界流速）
   void ReadOverlap_FVelocity(FILE* fp,
   double& flop,
   DifferentRestartInfo* DRI,
   DfiInfo* DI,
   REAL_TYPE* d_wk,
   int* rank_list,
   int recv_rank,
   int* assign,
   int nassign);
   
   
   // @brief ファイル読み込み＋オーバーラップを移しこみ（温度）
   void ReadOverlap_Temperature(FILE* fp,
   double& flop,
   DifferentRestartInfo* DRI,
   DfiInfo* DI,
   REAL_TYPE* d_wk,
   int* rank_list,
   int recv_rank,
   int* assign,
   int nassign);
   
   
   // @brief ファイル読み込み＋オーバーラップを移しこみ
   void ReadOverlap(FILE* fp,
   double& flop,
   DifferentRestartInfo* DRI,
   REAL_TYPE* d_wk);
   
   // @brief オーバーラップ領域を計算
   void SetOverlap(REAL_TYPE* write_wk,
   REAL_TYPE* read_wk,
   int dim,
   int gd,
   int* h,
   int* s,
   int* overlap_h,
   int* overlap_t,
   int* head,
   int* size);
   
   // 2倍密格子の領域開始インデクス番号から、その領域が属する粗格子計算結果ファイル名と、その計算結果ファイルの開始インデクス番号を取得する
   bool getCoarseResult (int i, int j, int k,
   std::string& coarse_dfi_fname,
   std::string& coarse_prefix,
   const int m_step,
   std::string& coarse_sph_fname,
   int* c_size,
   int* coarse,
   int* block
   );
   
   
   // 2倍密格子の領域開始インデクス番号から、その領域が属する粗格子計算結果ファイル名と、その計算結果ファイルの開始インデクス番号を取得する
   bool getCoarseResult2(int i, int j, int k,
   std::string& coarse_dfi_fname,
   std::string& coarse_prefix,
   const int m_step,
   std::string& coarse_sph_fname,
   int* c_size,
   int* coarse,
   int* block
   );
   
   // @brief オーバーラップ領域を計算
   void CalOverlap(int* overlap_h, int* overlap_t, int* h, int* t, int* head, int* tail);
   
   
   // @brief 粗格子から密格子へ内挿
   void Interpolation_from_coarse_initial(const int* m_st, const int* m_bk);
   
   
  */
};

#endif // _FFV_H_

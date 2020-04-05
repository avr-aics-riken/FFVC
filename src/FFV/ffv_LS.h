#ifndef _FFV_LS_H_
#define _FFV_LS_H_
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
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
 * @file   ffv_LS.h
 * @brief  LS Class header
 * @author aics
 */

#include "cpm_ParaManager.h"

#include "FB_Define.h"
#include "DomainInfo.h"
#include "IterationControl.h"
#include "Control.h"
#include "ffv_Ffunc.h"
#include "ffv_LSfunc.h"
#include "ffv_SetBC.h"
#include "FBUtility.h"

// FX10 profiler
#if defined __K_FPCOLL
#include "fjcoll.h"
#elif defined __FX_FAPP
#include "fj_tool/fjcoll.h"
#endif

#include "PerfMonitor.h"

using namespace pm_lib;

class LinearSolver : public DomainInfo, public IterationCtl {
  
private:
  int ensPeriodic[3];
  int ModeTiming;    ///< タイミング測定管理フラグ
  double face_comm_size; ///< 通信量
  
  Control* C;        ///< Controlクラス
  SetBC3D* BC;       ///< BCクラス
  PerfMonitor* PM;   ///< PerfMonitor class
  int* bcp;          ///< BCindex P
  int* bcd;          ///< BCindex ID
  
  REAL_TYPE* pcg_p;  ///< work for BiCGstab
  REAL_TYPE* pcg_p_; ///< work for BiCGstab
  REAL_TYPE* pcg_r;  ///< work for BiCGstab
  REAL_TYPE* pcg_r0; ///< work for BiCGstab
  REAL_TYPE* pcg_q ; ///< work for BiCGstab
  REAL_TYPE* pcg_s;  ///< work for BiCGstab
  REAL_TYPE* pcg_s_; ///< work for BiCGstab
  REAL_TYPE* pcg_t ; ///< work for BiCGstab
  REAL_TYPE* pcg_t_; ///< work for BiCGstab
  
  int cf_sz[3];     ///< SOR2SMAの反復の場合のバッファサイズ
  REAL_TYPE *cf_x;  ///< i方向のバッファ
  REAL_TYPE *cf_y;  ///< j方向のバッファ
  REAL_TYPE *cf_z;  ///< k方向のバッファ
  
public:
  
  /** コンストラクタ */
  LinearSolver() : DomainInfo(), IterationCtl() {
    C   = NULL;
    BC  = NULL;
    PM  = NULL;
    bcp = NULL;
    bcd = NULL;
    pcg_p  = NULL;
    pcg_p_ = NULL;
    pcg_r  = NULL;
    pcg_r0 = NULL;
    pcg_q  = NULL;
    pcg_s  = NULL;
    pcg_s_ = NULL;
    pcg_t  = NULL;
    pcg_t_ = NULL;
    cf_x = NULL;
    cf_y = NULL;
    cf_z = NULL;
    
    ModeTiming = 0;
    face_comm_size = 0.0;
    
    for (int i=0; i<3; i++)
    {
      ensPeriodic[i] = 0;
      cf_sz[i] = 0;
    }
  }
  
  /**　デストラクタ */
  ~LinearSolver() {}
  

  
protected:
  
  /**
   * @brief  Fcheck 非Div反復
   * @retval 収束したら true
   * @param [in]  var    誤差、残差、解ベクトルのL2ノルム
   * @param [in]  b_l2   右辺ベクトルのL2ノルム
   * @param [in]  r0_l2  初期残差ベクトルのL2ノルム
   */
  bool Fcheck(double* var, const double b_l2, const double r0_l2);
  
  
  /**
   * @brief Fdot for 1 array
   * @retval  内積値
   * @param [in]   x   vector1
   */
  double Fdot1(REAL_TYPE* x);
  
  
  /**
   * @brief Fdot for 2 arrays
   * @retval  内積値
   * @param [in]   x   vector1
   * @param [in]   y   vector2
   */
  double Fdot2(REAL_TYPE* x, REAL_TYPE* y);
  
  
  /**
   * @brief Preconditioner
   * @param [in,out] x   解ベクトル
   * @param [in]     b   RHS vector
   * @param [in]     dt  時間積分幅
   */
  void Preconditioner(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt);
  
  
  
  /**
   * @brief 反復の同期処理
   * @param [in,out] d_class   対象データ
   * @param [in]     num_layer 通信の袖数
   */
  void SyncScalar(REAL_TYPE* d_class, const int num_layer);

  
  
  /**
   * @brief タイミング測定開始
   * @param [in] key ラベル
   */
  inline void TIMING_start(const string key)
  {
    // PMlib Intrinsic profiler
    TIMING__ PM->start(key);
    
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
    TIMING__ PM->stop(key, flopPerTask, (unsigned)iterationCount);
  }
  
  
public:
  
  /**
   * @brief 初期化
   * @param [in]  C      Controlクラス
   * @param [in]  BC     SetBC3Dクラス
   * @param [in]  ModeT  タイミング測定モード
   * @param [in]  f_comm 通信量
   * @param [in]  PM     PerfMonitorクラス
   * @param [in]  bcp    BCindex P
   * @param [in]  bcd    BCindex ID
   * @param [in]  pcg_p  array for BiCGstab
   * @param [in]  pcg_p_ array for BiCGstab
   * @param [in]  pcg_r  array for BiCGstab
   * @param [in]  pcg_r0 array for BiCGstab
   * @param [in]  pcg_q  array for BiCGstab
   * @param [in]  pcg_s  array for BiCGstab
   * @param [in]  pcg_s_ array for BiCGstab
   * @param [in]  pcg_t  array for BiCGstab
   * @param [in]  pcg_t_ array for BiCGstab
   * @param [in]  ensP   周期境界の存在
   * @param [in]  cf_sz  バッファサイズ
   * @param [in]  cf_x   バッファ x方向
   * @param [in]  cf_y   バッファ y方向
   * @param [in]  cf_z   バッファ z方向
   */
  void Initialize(Control* C,
                  SetBC3D* BC,
                  int ModeT,
                  double f_comm,
                  PerfMonitor* PM,
                  int* bcp,
                  int* bcd,
                  REAL_TYPE* pcg_p,
                  REAL_TYPE* pcg_p_,
                  REAL_TYPE* pcg_r,
                  REAL_TYPE* pcg_r0,
                  REAL_TYPE* pcg_q,
                  REAL_TYPE* pcg_s,
                  REAL_TYPE* pcg_s_,
                  REAL_TYPE* pcg_t,
                  REAL_TYPE* pcg_t_,
                  const int* ensP,
                  const int* cf_sz,
                  REAL_TYPE* cf_x,
                  REAL_TYPE* cf_y,
                  REAL_TYPE* cf_z);
  
  
  /** 
   * @brief SOR法
   * @retval 反復数
   * @param [in,out] x              解ベクトル
   * @param [in]     b              RHS vector
   * @param [in]     dt             時間積分幅
   * @param [in]     itrMax         反復最大値
   * @param [in]     b_l2           L2 norm of b vector
   * @param [in]     r0_l2          初期残差ベクトルのL2ノルム
   * @param [in]     converge_check 収束判定を行う(true)
   */
  int PointSOR(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt, const int itrMax, const double b_l2, const double r0_l2, bool converge_check=true);
  
  
  /**
   * @brief SSOR法
   * @retval 反復数
   * @param [in,out] x              解ベクトル
   * @param [in]     b              RHS vector
   * @param [in]     dt             時間積分幅
   * @param [in]     itrMax         反復最大値
   * @param [in]     b_l2           L2 norm of b vector
   * @param [in]     r0_l2          初期残差ベクトルのL2ノルム
   * @param [in]     converge_check 収束判定を行う(true)
   */
  int PointSSOR(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt, const int itrMax, const double b_l2, const double r0_l2, bool converge_check=true);
  
  
  /**
   * @brief 速度の発散の打ち切り誤差が4次となるSOR法
   * @retval 反復数
   * @param [in,out] x      解ベクトル
   * @param [in]     b      RHS vector
   * @param [in]     u_sum  \sum{u^*}
   * @param [in]     w1     ワーク配列
   * @param [in]     w2     ワーク配列
   * @param [in]     dt     時間積分幅
   * @param [in]     dh     格子幅
   * @param [in]     b_l2   L2 norm of b vector
   * @param [in]     r0_l2  初期残差ベクトルのL2ノルム
   */
  int PointSOR_4th(REAL_TYPE* x, REAL_TYPE* b, REAL_TYPE* u_sum, REAL_TYPE* w1, REAL_TYPE* w2, REAL_TYPE dt, REAL_TYPE dh, const double b_l2, const double r0_l2);
  
  
  
  /**
   * @brief Residual Cut SOR法
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     poi_rhs RHS vector of Poisson
   * @param [in]     b       RHS of Ax=b
   * @param [in]     bcp     BCindex P
   * @param [in]     r0_l2   初期残差ベクトルのL2ノルム
   */
  int RC_sor(REAL_TYPE* x, REAL_TYPE* pos_rhs, REAL_TYPE* b, int* bcp, const double r0_l2);
  
  
  /**
   * @brief 2色オーダリングSORのストライドメモリアクセス版
   * @retval 反復数
   * @param [in,out] x              解ベクトル
   * @param [in]     b              RHS vector
   * @param [in]     dt             時間積分幅
   * @param [in]     itrMax         反復最大値
   * @param [in]     b_l2           L2 norm of b vector
   * @param [in]     r0_l2          初期残差ベクトルのL2ノルム
   * @param [in]     converge_check 収束判定を行う(true)
   */
  int SOR2_SMA(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt, const int itrMax, const double b_l2, const double r0_l2, bool converge_check=true);
  

  /**
   * @brief 前処理つきBiCGstab
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     b       RHS vector
   * @param [in]     dt      時間積分幅
   * @param [in]     b_l2    L2 norm of b vector
   * @param [in]     r0_l2   初期残差ベクトルのL2ノルム
   */
  int PBiCGstab(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt, const double b_l2, const double r0_l2);
  
  
  /**
   * @brief  FPCG
   * @retval 反復数
   * @param [in]     IC      IterationCtlクラス
   * @param [in,out] x       解ベクトル
   * @param [in]     b  RHS  vector
   * @param [in]     rhs_nrm RHS vector
   * @param [in]     r0      初期残差ベクトル
   */
  //int Fpcg(IterationCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0);

};

#endif // _FFV_LS_H_

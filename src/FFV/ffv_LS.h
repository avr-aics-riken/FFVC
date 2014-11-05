#ifndef _FFV_LS_H_
#define _FFV_LS_H_
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

/**
 * @file   ffv_LS.h
 * @brief  LS Class header
 * @author aics
 */

#include "cpm_ParaManager.h"

#include "../FB/FB_Define.h"
#include "../FB/DomainInfo.h"
#include "../FB/IterationControl.h"
#include "../FB/Control.h"
#include "../F_CORE/ffv_Ffunc.h"
#include "../F_LS/ffv_LSfunc.h"
#include "ffv_SetBC.h"


class LinearSolver : public DomainInfo, public IterationCtl {
  
private:
  int ensPeriodic[3];
  
  Control* C;        ///< Controlクラス
  SetBC3D* BC;       ///< BCクラス
  int* bcp;          ///< BCindex P
  int* bcd;          ///< BCindex ID
  REAL_TYPE* pni;    ///< coef fot naive implementation test
  
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
    bcp = NULL;
    bcd = NULL;
    pni = NULL;
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
   * @brief Fdot1
   * @retval  内積値
   * @param [in]   x   vector1
   */
  double Fdot1(REAL_TYPE* x);
  
  
  /**
   * @brief Fdot2
   * @retval  内積値
   * @param [in]   x   vector1
   * @param [in]   y   vector2
   */
  double Fdot2(REAL_TYPE* x, REAL_TYPE* y);
  
  
  /**
   * @brief Preconditioner
   * @param [in,out] x  解ベクトル
   * @param [in]     b  RHS vector
   */
  void Preconditioner(REAL_TYPE* x, REAL_TYPE* b);
  
  
  /**
   * @brief Smoother
   * @param [in,out] x    解ベクトル
   * @param [in]     b    RHS vector
   */
  void Smoother(REAL_TYPE* x, REAL_TYPE* b);
  
  
  /**
   * @brief 反復の同期処理
   * @param [in,out] d_class   対象データ
   * @param [in]     num_layer 通信の袖数
   */
  void SyncScalar(REAL_TYPE* d_class, const int num_layer);
  
  
  int Fpreconditioner(REAL_TYPE* x, REAL_TYPE* b);
  void Fsmoother(REAL_TYPE* x, REAL_TYPE* b, REAL_TYPE omg);
  void Fdot(REAL_TYPE* xy, REAL_TYPE* x, REAL_TYPE* y);
  
  
public:
  
  /**
   * @brief  BiCGstab
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     b       RHS vector
   * @param [in]     b_l2    L2 norm of b vector
   * @param [in]     r0_l2   初期残差ベクトルのL2ノルム
   */
  int BiCGstab(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2);
  
  
  /**
   * @brief  PBiCGstab
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     b       RHS vector
   * @param [in]     b_l2    L2 norm of b vector
   * @param [in]     r0_l2   初期残差ベクトルのL2ノルム
   */
  int PBiCGstab(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2);
  
  
  /**
   * @brief 初期化
   * @param [in]  C      Controlクラス
   * @param [in]  BC     SetBC3Dクラス
   * @param [in]  bcp    BCindex P
   * @param [in]  bcd    BCindex ID
   * @param [in]  pni    array for naive implementation
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
                  int* bcp,
                  int* bcd,
                  REAL_TYPE* pni,
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
   * @param [in,out] x      解ベクトル
   * @param [in]     b      RHS vector
   * @param [in]     b_l2   L2 norm of b vector
   * @param [in]     r0_l2  初期残差ベクトルのL2ノルム
   */
  int PointSOR(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2);
  
  
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
   * @brief 2色オーダリングSORのストライドメモリアクセス版
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     b       RHS vector
   * @param [in]     itrMax  反復最大値
   * @param [in]     b_l2    L2 norm of b vector
   * @param [in]     r0_l2   初期残差ベクトルのL2ノルム
   */
  int SOR2_SMA(REAL_TYPE* x, REAL_TYPE* b, const int itrMax, const double b_l2, const double r0_l2);
  

  /**
   * @brief FGMRES
   * @param [in]     IC      IterationCtlクラス
   * @param [in]     rhs_nrm RHS vectorのL2ノルム
   * @param [in]     r0      初期残差ベクトル
   */
  //void Fgmres(IterationCtl* IC, const double rhs_nrm, const double r0);
  
  
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
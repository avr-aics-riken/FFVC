#ifndef _SKL_FB_ITRCTL_H_
#define _SKL_FB_ITRCTL_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IterationCtl.h
//@brief FlowBase ItrCtl class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "FBDefine.h"
#include "mydebug.h"

class ItrCtl {
private:
  unsigned NormType;     /// ノルムの種類
  unsigned SubType;      /// SKIP LOOP or MASK LOOP
  unsigned ItrMax;       /// 最大反復数
  unsigned LinearSolver; /// 線形ソルバーの種類
  REAL_TYPE eps;          /// 収束閾値
  REAL_TYPE omg;          /// 加速/緩和係数
  REAL_TYPE NormValue;    /// ノルムの値
  
public:
  /// 反復制御リスト
  enum itr_cntl_key {
    ic_prs_pr,
    ic_prs_cr,
    ic_vis_cn,
    ic_tdf_ei,
    ic_tdf_cn,
    ic_END
  };
  
  /// 反復法の収束基準種別
  enum norm_type { 
    v_div_max,
    v_div_max_dbg,
    v_div_l2,
    v_res_l2_a,
    v_res_l2_r,
    p_res_l2_a,
    p_res_l2_r,
    t_res_l2_a,
    t_res_l2_r
  };
  
  unsigned LoopCount;     /// 反復回数

  ItrCtl() {
    NormType = 0;
    ItrMax = LoopCount = LinearSolver = SubType = 0;
    eps = omg = 0.0;
  }
  ~ItrCtl() {}
  
public:
  // @fn unsigned get_LS(void) const
  // @brief 線形ソルバの種類を返す
  unsigned get_LS(void) const { return LinearSolver; }
  
  // @fn unsigned get_ItrMax(void) const
  // @brief 最大反復回数を返す
  unsigned get_ItrMax(void) const { return ItrMax; }
  
  // @fn unsigned get_LoopType(void) const
  // @brief ループ実装の種類を返す
  unsigned get_LoopType(void) const { return SubType; }
  
  // @fn REAL_TYPE get_omg(void) const
  // @brief 緩和/加速係数を返す
  REAL_TYPE get_omg(void) const { return omg; }
  
  // @fn REAL_TYPE get_eps(void) const
  // @brief 収束閾値を返す
  REAL_TYPE get_eps(void) const { return eps; }
  
  // @fn unsigned get_normType(void) const
  // @brief ノルムのタイプを返す
  unsigned get_normType(void) const { return NormType; }
  
  // @fn REAL_TYPE get_normValue(void) const
  // @brief keyに対応するノルムの値を返す
  REAL_TYPE get_normValue(void) const { return NormValue; }
  
  // @fn void set_LS(unsigned key)
  // @brief 線形ソルバの種類を設定する
  void set_LS(unsigned key) { LinearSolver=key; }
  
  // @fn void set_ItrMax(unsigned key)
  // @brief 最大反復回数を設定する
  void set_ItrMax(unsigned key) { ItrMax=key; }
  
  // @fn void set_LoopType(unsigned key)
  // @brief ループ実装の種類を設定する
  void set_LoopType(unsigned key) { SubType=key; }
  
  // @fn void set_omg(REAL_TYPE r)
  // @brief 緩和/加速係数を保持
  void set_omg(REAL_TYPE r) { omg = r; }
  
  // @fn void set_eps(REAL_TYPE r)
  // @brief 収束閾値を保持
  void set_eps(REAL_TYPE r) { eps = r; }
  
  // @fn void set_normType(unsigned n) 
  // @brief ノルムのタイプを保持
  void set_normType(unsigned n) { NormType = n; }
  
  // @fn void set_normValue(REAL_TYPE r)
  // @brief ノルム値を保持
  void set_normValue(REAL_TYPE r) { NormValue = r; }
};

#endif // _SKL_FB_ITRCTL_H_

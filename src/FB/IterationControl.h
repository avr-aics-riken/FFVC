#ifndef _FB_ITERATION_H_
#define _FB_ITERATION_H_

//##################################################################################
//
// Flow Base class
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

/**
 * @file   IterationControl.h
 * @brief  FlowBase IterationCtl class Header
 * @author kero
 */

#include <stdlib.h>
#include <string>
#include "FB_Define.h"


// #################################################################
class IterationCtl {
  
private:
  double NormValue;     ///< ノルムの値 （計算実行中に利用）
  double eps;           ///< 収束閾値
  double omg;           ///< 加速/緩和係数
  int NormType;         ///< ノルムの種類
  int MaxIteration;     ///< 最大反復数
  int LinearSolver;     ///< 線形ソルバーの種類
  int LoopCount;        ///< 反復回数 （計算実行中に利用）
  int valid;            ///< 有効フラグ
  int Sync;             ///< 同期モード (comm_sync, comm_async)
  std::string alias;    ///< 別名

  
public:
  
  /** コンストラクタ */
  IterationCtl() {
    NormType = 0;
    MaxIteration = 0;
    LoopCount = 0;
    LinearSolver = 0;
    eps = 0.0;
    NormValue = 0.0;
    valid = -1;
    omg = 0.0;
    Sync = -1;
  }
  
  /**　デストラクタ */
  ~IterationCtl() {}
  
  
public:
  
  // @brief 基本メンバー変数のコピー
  // @param [in] src コピー元
  void copy(IterationCtl* src);
  
  // @brief 別名を返す
  std::string getAlias() const
  {
    return alias;
  }
  
  // @brief 収束閾値を返す
  double getCriterion() const
  {
    return eps;
  }
  
  // @brief 反復カウントを返す
  int getLoopCount() const
  {
    return LoopCount;
  }
  
  // @brief 線形ソルバの種類を返す
  int getLS() const
  {
    return LinearSolver;
  }
  
  // @brief 最大反復回数を返す
  int getMaxIteration() const
  {
    return MaxIteration;
  }
  
  // @brief ノルムの文字列を返す
  std::string getNormString();
  
  // @brief ノルムのタイプを返す
  int getNormType() const
  {
    return NormType;
  }
  
  // @brief keyに対応するノルムの値を返す
  double getNormValue() const
  {
    return NormValue;
  }
  
  // @brief 有効フラグを返す
  int getValid() const
  {
    return valid;
  }
  
  // @brief 反復ループカウントを設定する
  void incLoopCount()
  {
    LoopCount++;
  }
  
  // @brief 収束しているか > ture
  bool isConverged()
  {
    return (NormValue < eps) ? true : false;
  }
  
  // @brief aliasを設定する
  void setAlias(std::string key)
  {
    alias = key;
  }
  
  // @brief 収束閾値を保持
  void setCriterion(const double r)
  {
    eps = r;
  }
  
  // @brief LoopCountを設定する
  void setLoopCount(const int key)
  {
    LoopCount = key;
  }
  
  // @brief 線形ソルバの種類を設定する
  void setLS(const int key)
  {
    LinearSolver = key;
  }
  
  // @brief 最大反復回数を設定する
  void setMaxIteration(const int key)
  {
    MaxIteration = key;
  }
  
  // @brief ノルムのタイプを保持
  void setNormType(const int n)
  {
    NormType = n;
  }
  
  // @brief ノルム値を保持
  void setNormValue(const double r)
  {
    NormValue = r;
  }
  
  // @brief validフラグを設定
  void setValid(const int n)
  {
    valid = n;
  }
  
  // @brief 緩和/加速係数を返す
  double getOmega() const
  {
    return omg;
  }
  
  // @brief 同期モードを返す
  int getSyncMode() const
  {
    return Sync;
  }
  
  // @brief 緩和/加速係数を保持
  void setOmega(const double r)
  {
    omg = r;
  }
  
  // @brief 同期モードを保持
  void setSyncMode(const Synch_Mode r)
  {
    Sync = r;
  }
  
};

#endif // _FB_ITERATION_H_
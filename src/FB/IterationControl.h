#ifndef _FB_ITERATION_H_
#define _FB_ITERATION_H_

//##################################################################################
//
// Flow Base class
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

/**
 * @file   IterationControl.h
 * @brief  FlowBase IterationCtl class Header
 * @author aics
 * @todo 先に、pressure="sor"の文字列を確認して、リストアップされたパラメータを読み込むようにした方がいいかも
 */

#include "FB_Define.h"
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <strings.h>
#include "TextParser.h"

using namespace std;

// #################################################################
class IterationCtl {
  
protected:
  
  double residual;      ///< 残差
  double error;         ///< 誤差
  double eps_res;       ///< 残差の収束閾値
  double eps_err;       ///< 誤差の収束閾値
  double omg;           ///< 加速/緩和係数
  int ErrNorm;          ///< 誤差ノルムの種類
  int ResNorm;          ///< 残差ノルムの種類
  int MaxIteration;     ///< 最大反復数
  int LinearSolver;     ///< 線形ソルバーの種類
  int smoother;         ///< 前処理法
  int LoopCount;        ///< 反復回数 （計算実行中に利用）
  int Sync;             ///< 同期モード (comm_sync, comm_async)
  int precondition;     ///< 前処理mode
  int InnerItr;         ///< 内部反復回数
  string alias;         ///< 別名
  
  
public:
  
  /** コンストラクタ */
  IterationCtl() {
    precondition = OFF;
    ErrNorm = 0;
    ResNorm = 0;
    MaxIteration = 0;
    LoopCount = 0;
    LinearSolver = 0; // zero indicates undefined state.
    eps_res = 0.0;
    eps_err = 0.0;
    residual = 0.0;
    error = 0.0;
    omg = 0.0;
    Sync = -1;
    InnerItr = 0;
    smoother = -1;
    
    eps_err = ( sizeof(REAL_TYPE) == 4 ) ? 4.0*SINGLE_EPSILON : 4.0*DOUBLE_EPSILON;
  }
  
  /**　デストラクタ */
  virtual ~IterationCtl() {}
  
  
  // @brief 基本メンバー変数のコピー
  // @param [in] src コピー元
  void copy(IterationCtl* src);
  
  
  // @brief 別名を返す
  string getAlias() const
  {
    return alias;
  }
  
  
  // @brief keyに対応する誤差を返す
  double getError() const
  {
    return error;
  }
  
  
  // @brief 誤差の収束閾値を返す
  double getErrCriterion() const
  {
    return eps_err;
  }
  
  
  // @brief 誤差ノルムの文字列を返す
  string getErrNormString();
  
  
  // @brief 誤差ノルムの種類を返す
  int getErrType() const
  {
    return ErrNorm;
  }
  
  
  /**
   * @brief 固有パラメータを取得
   * @param [in]  tpCntl   TextParser pointer
   * @param [in]  base     ラベル
   */
  bool getInherentPara(TextParser* tpCntl, const string base);
  
  
  // @brief 内部反復数を返す
  int getInnerItr() const
  {
    return InnerItr;
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
  
  
  // @brief 緩和/加速係数を返す
  double getOmega() const
  {
    return omg;
  }
  
  
  // Gmres反復固有のパラメータを指定する
  void getParaGmres(TextParser* tpCntl, const string base);
  
  
  // Jacobi反復固有のパラメータを指定する
  void getParaJacobi(TextParser* tpCntl, const string base);
  
  
  // BiCGSTAB反復固有のパラメータを指定する
  void getParaBiCGSTAB(TextParser* tpCntl, const string base);
  
  
  // RB-SOR反復固有のパラメータを指定する
  void getParaSOR2(TextParser* tpCntl, const string base);
  
  
  // Div反復固有のパラメータを指定する
  bool getParaVP(TextParser* tpCntl);
  
  
  // @brief 残差の収束閾値を返す
  double getResCriterion() const
  {
    return eps_res;
  }
  
  
  // @brief 残差ノルムの文字列を返す
  string getResNormString();
  
  
  // @brief 残差ノルムの種類を返す
  int getResType() const
  {
    return ResNorm;
  }
  
  
  // @brief keyに対応する残差を返す
  double getResidual() const
  {
    return residual;
  }
  
  
  // @brief スムーザーの種類を返す
  int getSmoother() const
  {
    return smoother;
  }
  
  
  // @brief 同期モードを返す
  int getSyncMode() const
  {
    return Sync;
  }
  
  
  // @brief 前処理の有無を返す
  // @retval true -> preconditioned
  bool isPreconditioned() const
  {
    return (precondition == ON) ? true : false;
  }
  
  
  // @brief 反復ループカウントを設定する
  void incLoopCount()
  {
    LoopCount++;
  }
  
  
  // @brief 残差が収束しているか > ture
  bool isResConverged() const
  {
    return (residual < eps_res) ? true : false;
  }
  
  
  // @brief 残差が収束しているか > ture
  bool isErrConverged() const
  {
    return (error < eps_err) ? true : false;
  }
  
  
  // @brief aliasを設定する
  void setAlias(std::string key)
  {
    alias = key;
  }
  
  
  // @brief 収束閾値を保持
  void setResCriterion(const double r)
  {
    eps_res = r;
  }
  
  
  // @brief LoopCountを設定する
  void setLoopCount(const int key)
  {
    LoopCount = key;
  }
  
  
  /**
   * @brief 線形ソルバの種類を設定する
   * @param [in] str 反復法の指定文字列
   */
  bool setLS(const string str);
  
  
  // @brief 最大反復回数を設定する
  void setMaxIteration(const int key)
  {
    MaxIteration = key;
  }
  
  
  /**
   * @brief 誤差ノルムのタイプを保持
   * @param [in]  tpCntl   TextParser pointer
   * @param [in]  label    ラベル
   */
  bool setErrType(TextParser* tpCntl, const string label);

  
  // @brief 残差ノルムのタイプを保持
  void setResType(const int n)
  {
    ResNorm = n;
  }
  
  
  // @brief 残差を保持
  void setResidual(const double r)
  {
    residual = r;
  }
  
  
  // @brief 誤差を保持
  void setError(const double r)
  {
    error = r;
  }
  
  
  // @brief 同期モードを保持
  void setSyncMode(const Synch_Mode r)
  {
    Sync = r;
  }
  
};

#endif // _FB_ITERATION_H_

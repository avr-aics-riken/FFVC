#ifndef _FB_ITERATION_H_
#define _FB_ITERATION_H_

//##################################################################################
//
// Flow Base class
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
 * @file   IterationControl.h
 * @brief  FlowBase IterationCtl class Header
 * @author aics
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
  int Naive;            ///< Naive Implementation >> on/off
  int Bit3option;       ///< bit option
  string alias;         ///< 別名
  
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
    Naive = OFF;
    Bit3option = OFF;
  }
  
  /**　デストラクタ */
  ~IterationCtl() {}
  
  
public:
  
  // @brief 基本メンバー変数のコピー
  // @param [in] src コピー元
  void copy(IterationCtl* src);
  
  
  // @brief 別名を返す
  string getAlias() const
  {
    return alias;
  }
  
  
  // @brief 収束閾値を返す
  double getCriterion() const
  {
    return eps;
  }
  
  
  /**
   * @brief 固有パラメータを取得
   * @param [in]  tpCntl   TextParser pointer
   * @param [in]  base     ラベル
   * @param [out] m_naive  Naive option
   */
  bool getInherentPara(TextParser* tpCntl, const string base, int& m_naive);
  
  
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
  
  
  // @brief NaiveOptionを返す
  int getNaive() const
  {
    return Naive;
  }
  
  // @brief Bit3Optionを返す
  int getBit3() const
  {
    return Bit3option;
  }
  
  
  // @brief ノルムの文字列を返す
  string getNormString();
  
  
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
  
private:

  // Gmres反復固有のパラメータを指定する
  void getParaGmres(TextParser* tpCntl, const string base);
  
  
  // Jacobi反復固有のパラメータを指定する
  void getParaJacobi(TextParser* tpCntl, const string base);
  
  
  // PBiCGSTAB反復固有のパラメータを指定する
  void getParaPBiCGSTAB(TextParser* tpCntl, const string base);
  
  
  // PCG反復固有のパラメータを指定する
  void getParaPCG(TextParser* tpCntl, const string base);
  
  
  // RBGS反復固有のパラメータを指定する
  void getParaRBGS(TextParser* tpCntl, const string base);
  
  
  // SOR反復固有のパラメータを指定する
  void getParaSOR(TextParser* tpCntl, const string base);
  
  
  // RB-SOR反復固有のパラメータを指定する
  void getParaSOR2(TextParser* tpCntl, const string base);
  
};

#endif // _FB_ITERATION_H_